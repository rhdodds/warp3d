c ****************************************************************
c *                                                              *
c *        domain integral for 3-d isoparametric elements        *
c *        supports finite strains-rotations                     *
c *                 body forces (including inertia)              *
c *                 crack-face tractions                         *
c *                 temperature loads                            *
c *                 kinetic energy terms                         *
c *                 anisotropic thermal expansion coefficients   *
c *                 nonhomogeneous material properties           *
c *                                                              *
c *        interaction integral for 3-d isoparametric elements   *
c *        supports linear-elastic material behavior             *
c *                 crack-face tractions                         *
c *                 nonhomogeneous material properties           *
c *                                                              *
c *         element name      type no.         description       *
c *         ------------      --------         -----------       *
c *                                                              *
c *         q3disop              1         20 node brick(*)      *
c *         l3dsiop              2          8 node brick         *
c *         ts12isiop            3         12 node brick         *
c *         ts15isiop            4         15 node brick         *
c *         ts9isiop             5          9 node brick         *
c *                                                              *
c *                                                              *
c *         (*) not fully implemented                            *
c *                                                              *
c *           strain-stress ordering in warp3d vectors:          *
c *           ----------------------------------------           *
c *                                                              *
c *              eps-x, eps-y, eps-z, gam-xy, gam-yz, gam-xz     *
c *              sig-x, sig-y, sig-z, tau-xy, tau-yz, tau-xz     *
c *                                                              *
c *                                                              *
c ****************************************************************
c
      subroutine dielem ( coord, qvals, edispl, evel, eaccel, feload,
     a                    glb_sig_gp, props, iprops, lprops, erots,
     b                    rotate, e_jresults, e_node_temps,
     c                    temperatures, enode_alpha_ij, ierr, elemno,
     d                    geonl, numrows_stress, snodes, e_W_is,
     e                    enode_swd, glb_eps_enode,
     f                    glb_displ_grad_enode,
     g                    omit_J7_J8_front_elems, fgm_e,
     h                    fgm_nu, glb_eps_tens_gp, ym_nodes,
     i                    ym_front_node,nu_nodes, nu_front_node,
     j                    e_iresults, cf_traction_flags, cf_tractions,
     k                    front_elem_flag, seg_curves_flag )
c
      use global_data, only : mxgp, mxndel, myid, numprocs, out
      use main_data, only : fgm_node_values_defined,
     a                      initial_stresses_input
      use j_data, only : nj => size_j_values, ni => size_i_values,
     a                   one_point_rule, gdebug => debug_elements,
     b                   front_nodes, num_front_nodes, front_coords,
     c                   domain_origin, domain_type, front_order,
     d                   expanded_front_nodes, crack_curvature,
     e                   face_loading, process_temperatures,
     f                   max_exp_front,  comput_i, comput_j,
     g                   j_linear_formulation, j_geonl_formulation,
     h                   process_initial_state

c
      implicit none
c
c                    parameters
c
      integer :: elemno, iout, numrows_stress
      integer :: snodes(*), iprops(*) ! props for this element only
c
      double precision :: coord(3,*), qvals(*), feload(3,*),
     a   edispl(3,*), evel(3,*), eaccel(3,*), glb_sig_gp(*),
     b   e_node_temps(*), rotate(3,3), e_jresults(nj), erots(9,*),
     c   enode_alpha_ij(6,*), enode_swd(*), glb_eps_enode(6,*),
     d   glb_displ_grad_enode(9,*), glb_eps_tens_gp(9,*),
     e   ym_nodes(*), e_iresults(ni,ni), cf_tractions(*), nu_nodes(*),
     f   e_W_is(*)
      double precision :: nu_front_node
c
      real ::    props(*)
c
      logical :: lprops(*), geonl, temperatures, fgm_e,
     &           fgm_nu, cf_traction_flags(*), omit_J7_J8_front_elems,
     &           front_elem_flag, seg_curves_flag
c
c                    locals. double are group by similar sizing
c
      integer ::  etype, ptno, enode, fnode, order, faceno, flag,
     &            nfnode, order_face, gpn, i, j, k, ii, ierr, ngpts,
     &            numipts, nnode, numrow
      logical ::  debug, linear, linear_displ, debug_i,
     &            debug_j, qp_node, fgm_alphas

c
      double precision ::
     a   crk_displ(3,mxndel), crk_vel(3,mxndel), crk_accel(3,mxndel),
     b   sf(mxndel), glb_eps_tens_enode(9,mxndel),
     c   crk_eps_tens_enode(9,mxndel), crk_displ_grad_enode(9,mxndel),
     d   detF_enodes(mxndel)
c
      double precision ::
     a   crk_sig_tens_gp(10,mxgp+1), jacob(3,3,mxgp+1),
     b   jacobi(3,3,mxgp+1), gp_det_F(mxgp+1),
     c   glb_sig_tens_gp(10,mxgp+1), detvol(mxgp+1),
     d   crk_eps_tens_gp(9,mxgp+1)
c
      double precision ::
     &   qn(6,6), elem_alpha(6), dalpha_x1(6), dstrain_x(9),
     &   dgrad_x(9), eaverage(19), ! stress, strain tensors + W
     &   F_vec(9), F_tens(3,3)
      equivalence( F_vec, F_tens )
c
      double precision ::  jterm(nj), iterm(ni,ni)
c
      double precision :: dsf(mxndel,3,mxgp+1)
c
      double precision ::
     a   aux_stress(9,8), aux_strain(9,8),
     b   daux_strain_x1(9,8), du11_aux(8), du21_aux(8),
     c   du31_aux(8), du111_aux(8), du112_aux(8),
     d   du113_aux(8), du211_aux(8), du212_aux(8),
     e   du213_aux(8), du311_aux(8), du312_aux(8),
     f   du313_aux(8), dcijkl_x1(3), sijkl(3), dsijkl_x1(3)
c
      double precision :: e, nu, nx, ny, nz, kin_energy, r, t, toler,
     a   termu_aux, termv_aux, termw_aux, gpq, point_x,
     b   xsi, eta, zeta, weight, evol, point_q, point_y, point_z,
     c   point_velocity, point_accel_x, point_accel_y, point_accel_z,
     d   point_vel_x, point_vel_y, point_vel_z, point_ym, point_nu,
     e   dux, dvx, dwx, dqx, dqy, dtheta_x, dswd_x,
     f   point_temp, dqz, dym_x, dnu_x,
     g   ym_front_node, dudotx, dvdotx, dwdotx
c
      real, parameter :: neg_99=-99.0, fgm_tol=1.0
      double precision, parameter :: zero=0.0d0, half=0.5d0, two=2.0d0,
     &                               four=4.0d0, eight=8.0d0,
     &                               one=1.0d0
      double precision, external :: didet
c
c             set up basic element properties. load all six thermal
c             expansion coefficients from material associated with
c             the element. if they are all zero, we skip processing of
c             element for temperature loading. all temperature terms
c             in domain involve derivatives of temperature at gauss points
c             and expansion coefficients at gauss points. if specified
c             coefficients for the element are zero, we skip
c             temperature processing of elements.
c
      debug    = gdebug
      debug_j  = .false.
      debug_i  = .false.
      ierr     = 0
      iout     = out
      e        = props(7)
      nu       = props(8)
      etype    = iprops(1)
      linear_displ = etype .eq. 2  .or. ! at least1 linear edge
     &               etype .eq. 3  .or.
     &               etype .eq. 4  .or.
     &               etype .eq. 5  .or.
     &               etype .eq. 8  .or.
     &               etype .eq. 11
c
c             for 'fgm' alpha values, elem_alpha(1)-elem_alpha(3) will
c             receive a value of -99.0. in the loop over integration
c             points, subroutine dielav changes elem_alpha to the value
c             of alpha at the integration point by interpolating nodal
c             alpha values. to be consistent with the solution of the
c             boundary-value problem, for linear-displacement elements, dielav
c             will assign to elem_alpha the average of nodal alpha values.
c             order: xx, yy, zz, xy, yz, xz.  Must be symmetric tensor
c
      elem_alpha(1)   = props(9)
      elem_alpha(2)   = props(13)
      elem_alpha(3:6) = props(34:37)
      fgm_alphas = abs(props(9)-neg_99) .le. fgm_tol
c
c             set integration order info for volume integrals.
c
      order      = iprops(5)
      ngpts      = iprops(6)
      numipts    = ngpts
c
      linear    = .false.
      nnode     = iprops(2)
c
      e_jresults  = zero  ! all terms
      e_iresults  = zero

c
      if( debug ) call dielem_debug_1
c
c             set number of stress values per strain point in arrays
c
      numrow = numrows_stress
c
c             set up data needed at all integration points and nodes
c             needed for the volume integrals and crack face
c             traction integrals.
c
      if( debug ) write(iout,9013)
c
c        1:   transform unrotated Cauchy stresses to (global) Cauchy
c             stresses for geonl formulation. glb_sig_gp input contains
c             unrotated stresses. first convert to Cauchy stresses
c             then to 1 st PK stresses. these are in global coordinates.
c             note: 1 PK stresses are non-symmetric so we keep the
c             3x3 tensor stored as a vector. position 10 is for the
c             work density. the work density is scaled by det F to refer
c             to t=0 configuration.
c
      if( geonl ) call dielem_stresses_geonl
c
c        2.  for small-displacement formulation, copy the input
c            (symmetric) stresses (6x1) into 3x3 tensor form.
c
      if( .not. geonl ) call dielem_stresses_linear_form
c
c        3.   copy engineering strains into 3x3 tensor form.
c             change gamma strains to tensor strains for
c             off-diagonal terms. these strains not used
c             in J computations for geonl.
c
      call dielem_copy_strains
c
c        4:   rotate nodal coordinates, displacements, velocities
c             and accelerations to crack front normal coordinate
c             system.
c
      call dielrv( coord, crk_displ, crk_vel, crk_accel,
     &             edispl, evel, eaccel, rotate, nnode, iout, debug )
c
c        5:   modify the nodal q-values to reflect linear
c             interpolations for 20, 15, 12-node elements.
c
      call dieliq( qvals, debug, iout, nnode, etype, elemno, snodes,
     &             coord, front_elem_flag, num_front_nodes, front_nodes,
     &             front_coords, expanded_front_nodes, domain_type,
     &             domain_origin, front_order, qp_node, crack_curvature,
     &             max_exp_front )
c
c        6:   compute coordinate jacobian at all gauss points
c             and the center point of element. these are now in
c             cracked coordinates for configuration at t=0
c
      call dielem_jacob_crack_coors

c        7:  rotate stresses and strains to crack-front coordinates
c            at all gauss points. get the stress work density from
c            warp3d results.
c
      do ptno = 1, ngpts
         call dielrt( rotate, glb_sig_tens_gp(1,ptno),
     &                crk_sig_tens_gp(1,ptno) )
         crk_sig_tens_gp(10,ptno) = glb_sig_tens_gp(10,ptno)
         call dielrt( rotate, glb_eps_tens_gp(1,ptno),
     &                crk_eps_tens_gp(1,ptno) )
      end do
c
c        8:  rotate: global strains, global displacement gradients
c            to crack-front coordinates at element nodes. get
c            determinant of F (or Fbar) for element nodes as needed.
c
      detF_enodes = one  ! for small displacements
      associate( F => F_tens )
c
      do enode = 1, nnode
         call dielrt( rotate, glb_eps_tens_enode(1,enode),
     &                crk_eps_tens_enode(1,enode) )
         call dielrt( rotate, glb_displ_grad_enode(1,enode),
     &                crk_displ_grad_enode(1,enode) )
         if( j_geonl_formulation ) then
           F_vec(1:9)  = crk_displ_grad_enode(1:9,enode)
           F_tens(1,1) = F_tens(1,1) + one
           F_tens(2,2) = F_tens(2,2) + one
           F_tens(3,3) = F_tens(3,3) + one
           detF_enodes(enode) = didet( F_tens )
         end if
      end do
      end associate


c
c        9:   process temperature (inital) strains if required.
c             rotate the thermal expansion coefficients from
c             from global->crack coordinates in case they are
c             anisotropic. elem_alpha(1:6) in global is replaced
c             by elem_alpha(1:6) in crack. rotate values that are
c             constant over the element (elem_alpha) and values
c             at each node of the element (enode_alpha_ij). (nodal
c             values enable computation of the x1 derivative of
c             the alpha_ij term in J.) for isotropic CTEs this
c             is some unnecessary work.
c
      if( process_temperatures ) then
         call dippie( rotate, iout, debug, elem_alpha )
         do enode = 1, nnode
            call dippie( rotate, iout, debug, enode_alpha_ij(1,enode) )
         end do
      end if
      if( debug ) then
         write(iout,9008)
         do enode = 1, nnode
            write(iout,9007) enode, snodes(enode),
     &                       (enode_alpha_ij(i,enode),i=1,6)
         end do
      end if
c
c        10:   for single point integration of bricks,
c              compute average stresses, energy density, and
c              strain at center of element (front coords).
c              set number of gauss points = one for subsequent loop.
c
      if( one_point_rule ) call dielem_set_one_pt_integration
c
c        11.  evaluate volume integrals for domain
c             ------------------------------------
c
c             set up completed. loop over all gauss points and compute
c             each term of the domain integral for j, and the interaction
c             integral for i.
c
c             jterm(1) :  work density
c             jterm(2) :  traction - displacement gradient
c             jterm(3) :  kinetic energy denisty
c             jterm(4) :  accelerations
c             jterm(5) :  crack face loading (handled separately)
c             jterm(6) :  thermal strain (2 parts. set to zero for an fgm)
c             jterm(7) :  stress times partial of strain wrt x1 (for fgms)
c             jterm(8) :  partial of stress work density wrt x1 (for fgms)
c             (jterms 7 & 8 are used to replace the explicit partial
c              derivative of stress work density wrt x1 (for fgms))
c
c             iterm(1) :  stress * derivative of aux displacement
c             iterm(2) :  aux stress * derivative of displacement
c             iterm(3) :  mixed strain energy density (aux stress * strain)
c             iterm(4) :  first term of incompatibility (stress * 2nd deriv
c                         of aux displacement
c             iterm(5) :  second term of incompatibility (stress * deriv
c                         of aux strain)
c             iterm(6) :  deriv of constitutive tensor * strain * aux strain
c             iterm(7) :  (not yet verified)
c                         aux stress * deriv of cte * relative change in temp
c                       + aux stress * cte * deriv of relative change in temp
c             iterm(8) :  crack face traction * aux disp. derivative
c
      jterm = zero ! all terms
      iterm = zero
      evol  = zero
      if( debug ) write(iout,9100)
c
      do ii = 1, numipts
        ptno = ii
        call dielem_get_J_I_ingredients_for_point
        if( comput_j ) call dielem_J_terms( crk_sig_tens_gp(1,ptno) )
        if( comput_i ) call dielem_I_terms
      end do !! over ii points
c
c
c
c        12.  perform di evaluations for traction loaded faces
c             for jterm(5) and then for iterm(8).
c
      call di_calc_surface_integrals( elemno, etype, nnode, snodes,
     &              feload, cf_traction_flags,
     &              cf_tractions, rotate, dsf, jacobi,
     &              crk_displ, qvals, coord, front_nodes,
     &              front_coords, domain_type,
     &              domain_origin, num_front_nodes,
     &              front_order, ym_front_node,
     &              nu_front_node, comput_j,
     &              comput_i, jterm, iterm,
     &              front_elem_flag, qp_node,
     &              crack_curvature, face_loading, iout,
     &              debug )
c
c        13. save di values in result vectors. last debug
c
      e_jresults = jterm    ! (1:8)
      e_iresults = iterm    ! (1:8,1:8)
c
      call dielem_debug_2
c
      return
c
 9007 format(3x,i3,2x,i3,2x,6e14.6)
 9008 format(' >>>>> node alpha(1->6) for nodes in crack x-y-z:' )
 9013 format(' >>>>> start of data set-up for edi' )
 9100 format(/,' >>>>> start of integration loop:' )
c
      contains
c     ========
c
c***************************************************************
c                                                              *
c      dielem -> contains dielem_stresses_geonl                *
c                                                              *
c        written by: rhd                                       *
c        last modified: 6/18/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_stresses_geonl
      implicit none
c
      integer :: ptno, loc, ierr, num_enodes, num_gpts, int_order
      double precision :: det_F, dg(9,mxgp) ! displ gradients
c
      num_enodes = nnode
      num_gpts   = ngpts
      int_order  = order
c
c              get displacement gradients at all element Gauss
c              points. dg = F - I.
c
c              for 8-node elements, F is replaced by F-bar
c
      call di_displ_grad_one_elem( elemno, etype, num_enodes,
     &                             num_gpts, int_order, coord,
     &                             edispl, dg )
c
      do ptno = 1, ngpts
c
       loc = ( ptno-1 ) * numrow + 1
c
c           unrotated cauchy -> cauchy. adjust energy density by
c           user-defined state as required.
c           position 10 is work density
c
       call digetr( qn, erots(1,ptno) )
       call diqmp1( qn, glb_sig_gp(loc), glb_sig_tens_gp(1,ptno) )
       glb_sig_tens_gp(10,ptno) = glb_sig_gp(loc+6) - e_W_is(ptno)
c
c           Cauchy -> 1PK (non-symmetric).
c           Scale work density to volume at t = 0.
c
       call digrad( glb_sig_tens_gp(1,ptno), dg(1,ptno), iout,
     &              debug, ptno, elemno, det_F )
c
       gp_det_F(ptno) = det_F
c
      end do ! over points
c
      return
      end subroutine dielem_stresses_geonl

c***************************************************************
c                                                              *
c      dielem -> contains dielem_stresses_linear_form          *
c                                                              *
c        written by: rhd                                       *
c        last modified: 6/18/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_stresses_linear_form
      implicit none
c
      integer :: ptno, loc
c
c           adjust energy density for user-defined initial
c           state as required.
c           position 10 is work density
c
c
!DIR$ VECTOR ALIGNED
      do ptno = 1, ngpts
          loc = ( ptno-1 )  * numrow
          glb_sig_tens_gp(1,ptno)  = glb_sig_gp(loc+1)
          glb_sig_tens_gp(2,ptno)  = glb_sig_gp(loc+4)
          glb_sig_tens_gp(3,ptno)  = glb_sig_gp(loc+6)
          glb_sig_tens_gp(4,ptno)  = glb_sig_gp(loc+4)
          glb_sig_tens_gp(5,ptno)  = glb_sig_gp(loc+2)
          glb_sig_tens_gp(6,ptno)  = glb_sig_gp(loc+5)
          glb_sig_tens_gp(7,ptno)  = glb_sig_gp(loc+6)
          glb_sig_tens_gp(8,ptno)  = glb_sig_gp(loc+5)
          glb_sig_tens_gp(9,ptno)  = glb_sig_gp(loc+3)
          glb_sig_tens_gp(10,ptno) = glb_sig_gp(loc+7) - e_W_is(ptno)
          gp_det_F(ptno) = one
      end do
c
      return
      end subroutine dielem_stresses_linear_form
c***************************************************************
c                                                              *
c      dielem -> contains dielem_copy_strains                  *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_copy_strains
      implicit none
c
      integer :: enode
      double precision :: f
c
      f = 2.0d0
c
!DIR$ VECTOR ALIGNED
      do enode = 1, nnode
         glb_eps_tens_enode(1,enode) = glb_eps_enode(1,enode)
         glb_eps_tens_enode(2,enode) = glb_eps_enode(4,enode)/f
         glb_eps_tens_enode(3,enode) = glb_eps_enode(6,enode)/f
         glb_eps_tens_enode(4,enode) = glb_eps_enode(4,enode)/f
         glb_eps_tens_enode(5,enode) = glb_eps_enode(2,enode)
         glb_eps_tens_enode(6,enode) = glb_eps_enode(5,enode)/f
         glb_eps_tens_enode(7,enode) = glb_eps_enode(6,enode)/f
         glb_eps_tens_enode(8,enode) = glb_eps_enode(5,enode)/f
         glb_eps_tens_enode(9,enode) = glb_eps_enode(3,enode)
      end do
c
      return
      end subroutine dielem_copy_strains
c***************************************************************
c                                                              *
c      dielem -> contains dielem_jacob_crack_coors             *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_jacob_crack_coors
      implicit none
c
      integer :: ptno, ierr
      double precision :: xsi, eta, zeta
c
      do ptno = 1, ngpts + 1  !  center point last
c
c             isoparametric coordinates of gauss point
c
        if( ptno .eq. ngpts+1 ) then
           xsi  = zero
           eta  = zero
           zeta = zero
        else
           call getgpts( etype, order, ptno, xsi, eta, zeta, weight )
        end if
        call derivs( etype, xsi, eta, zeta, dsf(1,1,ptno),
     &               dsf(1,2,ptno), dsf(1,3,ptno) )
        if( debug ) write(iout,9110)  ptno, xsi, eta, zeta, weight
c
c             coordinate jacobian, it's inverse, and determinant.
c
        call dielcj( dsf(1,1,ptno), dsf(1,2,ptno), dsf(1,3,ptno),
     &               coord, nnode, jacob(1,1,ptno), jacobi(1,1,ptno),
     &               detvol(ptno), ierr, iout, debug )
        if( ierr .ne. 0 ) then
           call dieler( iout,ierr,5 ); return
        end if
      end do
c
      return
 9110 format(/,'  >> ptno, xsi, eta, zeta, weight: ',i2,4(2x,f10.4))
      end subroutine dielem_jacob_crack_coors
c***************************************************************
c                                                              *
c      dielem -> contains dielem_set_one_pt_integration        *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_set_one_pt_integration
      implicit none
c
      integer :: ptno, i
c
c              for single point integration of bricks,
c              compute average stresses, energy density, and
c              strain at center of element (front coords).
c              set number of gauss points = one for subsequent loop.
      numipts  = 1
      eaverage = zero  ! all terms
      do ptno = 1, ngpts
         do i = 1, 10
            eaverage(i) = eaverage(i) + crk_sig_tens_gp(i,ptno)
         end do
         do i = 11,19
            eaverage(i) = eaverage(i) + crk_eps_tens_gp(i,ptno)
         end do
      end do
      do i = 1, 10
        crk_sig_tens_gp(i,ngpts+1) = eaverage(i) / dble(ngpts)
      end do
      do i = 11, 19
        crk_eps_tens_gp(i,ngpts+1) = eaverage(i) / dble(ngpts)
      end do
c
      return
      end subroutine dielem_set_one_pt_integration
c***************************************************************
c                                                              *
c      dielem -> contains dielem_get_j_i_ingredients_for_point *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/31/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_get_j_i_ingredients_for_point
      implicit none
c
      integer :: element_debug, ptno_debug
      logical :: ldebug
c
      element_debug = 0
      ptno_debug    = 0
      ldebug = elemno .eq. element_debug .and. ptno .eq. ptno_debug
c
c             isoparametric coordinates of gauss point and weight.
c
      if( one_point_rule ) then
         xsi    = zero
         eta    = zero
         zeta   = zero
         weight = eight
         ptno   = ngpts + 1
      else
         call getgpts( etype, order, ptno, xsi, eta, zeta, weight )
      end if
      if( ldebug ) then
        write(iout,9000) elemno, ptno
        write(iout,9110) xsi, eta, zeta, weight
      end if
c
c           evaluate shape functions of all nodes at this point.
c           used to find value of q function at integration point,
c           the total velocity at the point, values of acceleration
c           at the point, and alpha values at the point.
c
      call shapef( etype, xsi, eta, zeta, sf )
      call dielav( sf, crk_vel, crk_accel, point_velocity,
     &             point_accel_x, point_accel_y, point_accel_z,
     &             point_vel_x, point_vel_y, point_vel_z,
     &             qvals, point_q,
     &             nnode, point_temp, e_node_temps, elem_alpha,
     &             enode_alpha_ij, linear_displ, fgm_alphas, elemno,
     &             ym_nodes, nu_nodes, point_ym, point_nu, fgm_e,
     &             fgm_nu, coord, point_x, point_y, point_z,
     &             seg_curves_flag, iout )
c
      if( ldebug ) then
        write(iout,9120) point_q, point_velocity, point_accel_x,
     &                   point_accel_y, point_accel_z, point_temp
      end if
c
c           for this integration point, compute displacement
c           derivatives, q-function derivatives and temperature
c           derivative in crack coordinate system.
c
c           scale element nodal values of W by average det F over
c           element (=1.0 for small strains)
c
c           strains will be used in the order eps11, eps22, eps33,
c           eps12, eps23, eps13
c
      dux              = zero
      dvx              = zero
      dwx              = zero
      dudotx           = zero  ! partial velocity / partial x_1
      dvdotx           = zero  !  "
      dwdotx           = zero  !  "
      dqx              = zero
      dqy              = zero
      dqz              = zero
      dtheta_x         = zero
      dswd_x           = zero
      dym_x            = zero
      dnu_x            = zero
c
      dalpha_x1        = zero ! all terms
      dstrain_x        = zero !   "
      dgrad_x          = zero !   "
c
!DIR$ VECTOR ALIGNED
      do enode = 1, nnode
       nx = dsf(enode,1,ptno) * jacobi(1,1,ptno) +
     &      dsf(enode,2,ptno) * jacobi(1,2,ptno) +
     &      dsf(enode,3,ptno) * jacobi(1,3,ptno)
       ny = dsf(enode,1,ptno) * jacobi(2,1,ptno) +
     &      dsf(enode,2,ptno) * jacobi(2,2,ptno) +
     &      dsf(enode,3,ptno) * jacobi(2,3,ptno)
       nz = dsf(enode,1,ptno) * jacobi(3,1,ptno) +
     &      dsf(enode,2,ptno) * jacobi(3,2,ptno) +
     &      dsf(enode,3,ptno) * jacobi(3,3,ptno)
       dux      = dux + nx * crk_displ(1,enode)
       dvx      = dvx + nx * crk_displ(2,enode)
       dwx      = dwx + nx * crk_displ(3,enode)
       dudotx   = dudotx + nx * crk_vel(1,enode)
       dvdotx   = dvdotx + nx * crk_vel(2,enode)
       dwdotx   = dwdotx + nx * crk_vel(3,enode)
       dqx      = dqx + nx * qvals(enode)
       dqy      = dqy + ny * qvals(enode)
       dqz      = dqz + nz * qvals(enode)
       dtheta_x = dtheta_x + nx * e_node_temps(enode)
       dswd_x   = dswd_x + nx * ( enode_swd(enode) *
     &                            detF_enodes(enode) )
       dym_x    = dym_x  + nx * ym_nodes(enode)
       dnu_x    = dnu_x  + nx * nu_nodes(enode)
c
!DIR$ VECTOR ALIGNED
       dstrain_x(1:9) = dstrain_x(1:9) +
     &                   nx * crk_eps_tens_enode(1:9,enode)
!DIR$ VECTOR ALIGNED
       dgrad_x(1:9) = dgrad_x(1:9) +
     &                   nx * crk_displ_grad_enode(1:9,enode)
c
      end do ! enode
c
      if( process_temperatures ) then
        do enode = 1, nnode
           nx = dsf(enode,1,ptno) * jacobi(1,1,ptno) +
     &          dsf(enode,2,ptno) * jacobi(1,2,ptno) +
     &          dsf(enode,3,ptno) * jacobi(1,3,ptno)
           dalpha_x1(1:6) = dalpha_x1(1:6) + nx *
     &                      enode_alpha_ij(1:6,enode)
        end do
      end if
c
c             include det coord J for simplicity
c
      weight  = weight * detvol(ptno)
c
      if( .not. ldebug ) return
      write(iout,*)
     &   '...crk_eps_tens_enode, crk_displ_grad_enode. el=21, pt=5'
      do enode = 1, nnode
       write(iout,*) '     enode: ', enode
       write(iout,9200) crk_eps_tens_enode(1:9,enode)
       write(iout,9200) crk_displ_grad_enode(1:9,enode)
      end do
c
      return
c
 9000 format(/'... debug in dielem_get_j_i_ingredients_for_point.',
     &        ' elem:',i8,i3)
 9110 format(/,15x,' >> xsi, eta, zeta, weight: ',i2,4(2x,f10.4))
 9120 format(  15x,' >> q-value, velocity           ',f4.2,d14.6,
     &       /,15x,' >> accelerations:              ',3d14.6,
     &       /,15x,' >> temperature:                ',d14.6)
 9200 format(3(12x,3d14.6/))
c
      end subroutine dielem_get_j_i_ingredients_for_point
c***************************************************************
c                                                              *
c      dielem -> contains dielem_J_terms                       *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/31/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dielem_J_terms( crk_sig_tens )
c
      implicit none
c
c             parameters
c
c              crk_sig_tens -> 1PK stress tensor in front coords (geonl)
c                           -> symm stess tensor small-strain theory
c
      double precision :: crk_sig_tens(10) ! int. point
c
c             locals
c
      integer :: i, element_debug, ptno_debug, ptype
      double precision :: temp1, temp2, temp3, pk1(3,3), pk1_vec(9),
     &                    alpha(3,3), alpha_vec(9), j4_a, j4_b, swd,
     &                    product, rho, j1_incr, j2_incr, j3_incr,
     &                    j4_incr, j5_incr, j6_incr, j7_incr, j8_incr
      double precision, parameter :: zero=0.0d0, half=0.5d0
      logical :: ldebug, neglect_term_j6, include_term_j6,
     &           neglect_terms_j7_j8
      equivalence( pk1, pk1_vec ), ( alpha, alpha_vec )
c
      element_debug = 0
      ptno_debug    = 0
      ldebug = elemno .eq. element_debug .and. ptno .eq. ptno_debug
c
      if( ldebug ) then
         if( numprocs .gt. 1 ) write(iout,*) "myid = ",myid
         write(iout,890) elemno, ptno
      end if
c
      pk1_vec(1:9) = crk_sig_tens(1:9)  ! defines pk1 for geonl or
c                                       ! symm small-eps stress
c
c             form jterm 1 of di = total work density * derivative
c             of q function w.r.t. crack x. also accumulate total
c             element volume for debugging.
c
c             dux, dvx, dwx -> displ derivatives wrt crack local x
c             crk_sig_tens -> stress in crack local. 1st PK for geonl
c
c             W (work density) has already be scaled by det F for
c             geonl solutions
c
      evol     = evol + weight ! weight includes det coord Jacobian
      swd      = crk_sig_tens(10)  ! stress work density
      j1_incr  = - weight * dqx * swd
      jterm(1) = jterm(1) + j1_incr
      if( ldebug ) write(iout,892) weight, dqx, swd, j1_incr
c
c            form jterm 3 of di = total kinetic energy density *
c            derivative of q function w.r.t. crack x.
c
      rho        = props(10)
      kin_energy = half * rho * point_velocity * point_velocity
      j3_incr    = - weight * dqx * kin_energy
      jterm(3)   = jterm(3) + j3_incr
      if( ldebug ) write(iout,950) rho, point_velocity, kin_energy,
     &                             j3_incr
c
c            form jterm 4 of di = acceleration and velocity gradient
c
c             j4_a = dotdot u_i u_{i,1}
c             j4_b = dot u_i  dot u_{i,1}
c
      j4_a = point_accel_x * dux + point_accel_y * dvx +
     &       point_accel_z * dwx
      j4_b = point_vel_x * dudotx +  point_vel_y * dvdotx +
     &       point_vel_z * dwdotx
      j4_incr  = ( j4_a + j4_b ) * weight * rho * point_q
      jterm(4) = jterm(4) + j4_incr
     &
      if( ldebug ) write(iout,902) dudotx, dvdotx, dwdotx, j4_a, j4_b,
     &                             j4_incr
c
c            form jterm 2 of di = stress * displacement derivatives
c            use full tensor since we can have 1st PK stresses.
c
c              J_2 = P_{ij} u_{i,1} q_{,j}
c
      temp1 = ( pk1(1,1)*dqx + pk1(1,2)*dqy + pk1(1,3)*dqz ) * dux
      temp2 = ( pk1(2,1)*dqx + pk1(2,2)*dqy + pk1(2,3)*dqz ) * dvx
      temp3 = ( pk1(3,1)*dqx + pk1(3,2)*dqy + pk1(3,3)*dqz ) * dwx
      j2_incr  =  weight*( temp1 + temp2 + temp3 )
      jterm(2) = jterm(2) + j2_incr
      if( ldebug ) write(iout,900) dux, dvx, dwx, temp1, temp2,
     &                             temp3, j2_incr
c
c            form jterm 6 of di caused by temperature loads.
c
c            there are two parts both dotted into stresses:
c             (a) expansion coefficients * derivative of temperature
c                 w.r.t. crack x
c             (b) temperature * derivative of expansion coefficients
c                 w.r.t. crack x
d
c            stresses can be non-symmetric. thermal expansion coeffs.
c            have been rotated to crack coordinates and are symmetric.
c
c            (a) we use the specified alpha expansion coefficients for
c            the element directly.
c
c            (b) we compute derivatives from node values of alpha for
c            element (computed by averaging
c            element values for elements incident on the node). off-
c            diagonal alpha terms are multiplied by 'half' to obtain the
c            tensorial component corresponding to tensorial strain. (see
c            subroutine dippie.)
c
c            note: jterm6 is not used when e and nu
c            are entered using the 'fgm' option. this is because jterm6
c            is implicitly included in jterm 7-8 = -W,1.
c
c            temperature dependent stress-strain curves. E, nu, ...
c            can vary spatially from temperature dependence.
c            jterm6 not used for this condition (seg_curves_flag)
c
      alpha(1,1) = elem_alpha(1)         ! xx, alpha equived -> alpha_vec
      alpha(2,1) = elem_alpha(4) * half  ! yx
      alpha(3,1) = elem_alpha(6) * half  ! zx
      alpha(1,2) = elem_alpha(4) * half  ! xy
      alpha(2,2) = elem_alpha(2)         ! yy
      alpha(3,2) = elem_alpha(5) * half  ! zy
      alpha(1,3) = elem_alpha(6) * half  ! xz
      alpha(2,3) = elem_alpha(5) * half  ! yz
      alpha(3,3) = elem_alpha(3)         ! zz
      temp1 = dot_product( alpha_vec, pk1_vec ) * dtheta_x * weight
     &                                          * point_q
      temp2 = pk1(1,1)*dalpha_x1(1) +
     &        pk1(2,2)*dalpha_x1(2) +
     &        pk1(3,3)*dalpha_x1(3) +
     &        dalpha_x1(4) * half * ( pk1(2,1)+pk1(1,2) ) +
     &        dalpha_x1(5) * half * ( pk1(3,2)+pk1(2,3) ) +
     &        dalpha_x1(6) * half * ( pk1(3,1)+pk1(1,3) )
      temp2 = temp2 * weight * point_q * point_temp
c
      j6_incr = temp1 + temp2
      jterm(6) = jterm(6) + temp1 + temp2
      neglect_term_j6 = fgm_node_values_defined .or. seg_curves_flag
     &                  .or. process_initial_state 
     &                  .or. initial_stresses_input
      include_term_j6 = .not. neglect_term_j6
      if( neglect_term_j6 ) jterm(6) = zero
      if( ldebug )
     &  write(iout,955) process_temperatures, include_term_j6,
     &                  dtheta_x, point_temp, dalpha_x1(1),
     &                  point_q, temp1, temp2, j6_incr
c
c              form jterm 7 & jterm 8 of J for explicit partial
c              derivative terms
c
c                   jterm7 = sigma_{ij} * eps_{ij},1  small strains
c                          = 1PK_{ij} * u_{i,j},1  geonl
c                   jterm8 = -W,1
c                   strains, gradients are in tensorial form.
c
c                   swd  scaled by det F when F computed
c
      if( geonl ) then
        product = dot_product( pk1_vec, dgrad_x )
        ptype = 1
      else
        ptype = 2    !  pk1_vec has small-strain, symm. stresses
        product = dot_product( pk1_vec, dstrain_x )
      end if
      j7_incr  = product * weight * point_q
      j8_incr  = - dswd_x * weight * point_q
      if( omit_J7_J8_front_elems ) then
         j7_incr = zero
         j8_incr = zero
      end if
      jterm(7) = jterm(7) + j7_incr
      jterm(8) = jterm(8) + j8_incr
      if( ldebug )
     &  write(iout,960) omit_J7_J8_front_elems, ptype, product,
     &                  j7_incr, j8_incr
c
      neglect_terms_j7_j8 = include_term_j6
      if( neglect_terms_j7_j8 ) jterm(7:8) = zero
c
      if( ldebug ) then
         write(iout,970)
         if( numprocs .gt. 1 ) write(iout,*) "   myid = ", myid
         write(iout,920) pk1_vec(1:9)
         write(iout,921) dstrain_x(1:9)
         write(iout,922) dgrad_x(1:9)
         write(iout,906) elem_alpha(1:6)
         write(iout,907) dalpha_x1(1:6)
       end if
c
c              done with j-integral calculations for this integration point.
c              calculation of jterm(5) for the surface integral occurs later.
c
      return
c
c
 890  format(///,"J-integral terms: element",2x,i7,2x,"point",2x,i2)
 892  format(' >>> weight, dqx, swd, j1_incre: ',4(1x,d13.6))
 900  format(' >>> dux, dvx, dwx: ',3d14.6,/,
     &       '     temp1, 2, 3  : ',3d14.6,/,
     &       '     j2_incr      : ',d14.6 )
 902  format(' >>> dudotx, dvdotx, dwdotx: ',3d14.6,/,
     &       '     j4_a, j_4b:             ',2d14.6,/,
     &       '     j4_incr:                ',d14.6 )
 905  format(//,10x,'stress tensor:',/,3(10x,3(e11.4,2x),/))
 906  format(10x,'elem_alpha(1:6): ',6d14.6)
 907  format(10x,'dalpha_x1(1:6) : ',6d14.6)
!  908  format(/,10x,'dtheta_x       : ',e13.6,
!      &       /,10x,'weight         : ',e13.6,
!      &       /,10x,'point_q        : ',e13.6,
!      &       /,10x,'point_temp     : ',e13.6 )
 910  format(' >>> jterm6 contributions: ',2(1x,e13.6),/,
     &       '     jterm(6)            : ',1x,e13.6 )
 920  format(10x,'stress tensor (1PK, or small-strain):',
     &        /,3(15x,3d14.6/))
 921  format(10x,'dstrain_x tensor:'/,3(15x,3d14.6/))
 922  format(10x,'drad_x tensor:'/,3(15x,3d14.6/))
 950  format(' >>> rho, pt. velocity, KE: ',3d14.6,/,
     &       '     j3_incr:               ',d14.6 )
 955  format(' >>> proc_temps, include_j6:         ',2l2,/,
     &       '     dtheta_x, point_temp, dalpha_x1(1):',3d14.6,/,
     &       '     point_q, temp1, temp2:',3d14.6,/,
     &       '     j6_incr:                       :',d14.6 )
 960  format(' >>> omit_J7_J8_front_elems, ptype, product:      ',
     &  l2,i2,d14.6,/,
     &       '     j7_incr, j8_incr:                  ',2d14.6)
 970  format(/,' Details of quantities to makeup J terms:')
c
      end subroutine dielem_J_terms
c
c***************************************************************
c                                                              *
c      dielem -> contains dielem_I_terms                       *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c

      subroutine dielem_I_terms
      implicit none
c
c             for straight element edges:
c
c             compute distance from integration point to the closest
c             line that connects two adjacent front nodes.
c             compute angle between integration point, line connecting
c             front nodes, and projection of integration point onto
c             crack plane.
c
c             for curved element edges:
c
c             compute distance from integration point to the curve
c             fitted through adjacent front nodes.
c             compute angle between integration point, curve, and
c             projection of integration point onto crack plane.
c
      call di_calc_r_theta( 2, front_nodes, num_front_nodes,
     &                   front_coords, domain_type, domain_origin,
     &                   front_order, point_x, point_y, point_z,
     &                   elemno, ptno, r, t, crack_curvature,
     &                   debug, iout )
c
c        compute constitutive tensor components
c
      call di_calc_constitutive( dcijkl_x1, sijkl, dsijkl_x1,
     &                        point_ym, point_nu, dym_x, dnu_x,
     &                        elemno, debug, iout )
c
c        compute auxiliary fields for stress intensity factors
c
      call di_calc_aux_fields_k( elemno, ptno, r, t, ym_front_node,
     &                        nu_front_node, dcijkl_x1, sijkl,
     &                        dsijkl_x1, aux_stress,
     &                        aux_strain, daux_strain_x1,
     &                        du11_aux,  du21_aux,  du31_aux,
     &                        du111_aux, du112_aux, du113_aux,
     &                        du211_aux, du212_aux, du213_aux,
     &                        du311_aux, du312_aux, du313_aux,
     &                        iout )
c
c        compute auxiliary fields for t-stresses
c
      call di_calc_aux_fields_t( elemno, ptno, r, t, ym_front_node,
     &                        nu_front_node, dcijkl_x1, sijkl,
     &                        dsijkl_x1, aux_stress,
     &                        aux_strain, daux_strain_x1,
     &                        du11_aux,  du21_aux,  du31_aux,
     &                        du111_aux, du112_aux, du113_aux,
     &                        du211_aux, du212_aux, du213_aux,
     &                        du311_aux, du312_aux, du313_aux,
     &                        iout )
c
c        compute interaction integral terms for stress intensity factors
c
      call di_calc_i_terms( ptno, dqx, dqy, dqz, dux, dvx, dwx,
     &                   dtheta_x, crk_sig_tens_gp, aux_stress,
     &                   crk_eps_tens_gp,
     &                   aux_strain, dstrain_x,
     &                   daux_strain_x1, dcijkl_x1,
     &                   du11_aux,  du21_aux,  du31_aux,
     &                   du111_aux, du211_aux, du311_aux,
     &                   du112_aux, du212_aux, du312_aux,
     &                   du113_aux, du213_aux, du313_aux,
     &                   process_temperatures, elem_alpha,
     &                   dalpha_x1, point_temp, point_q, weight,
     &                   elemno, fgm_e, fgm_nu, iterm,
     &                   iout, debug)
c
      return
      end subroutine dielem_I_terms
c***************************************************************
c                                                              *
c      dielem -> contains dielem_debug_2                       *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dielem_debug_2
      implicit none
c
      if( debug_j ) write (iout,9003) elemno,   jterm(1), jterm(2),
     &                                jterm(3), jterm(4), jterm(5),
     &                                jterm(6), jterm(7), jterm(8),
     &                                evol
      if( debug_i ) then
         write(iout,9009) elemno
         do j = 1, 8
            if( j.eq.1 ) write(iout,9010) "KI,   plane stress     :"
            if( j.eq.2 ) write(iout,9010) "KI,   plane strain     :"
            if( j.eq.3 ) write(iout,9010) "KII,  plane stress     :"
            if( j.eq.4 ) write(iout,9010) "KII,  plane strain     :"
            if( j.eq.5 ) write(iout,9010) "KIII, anti-plane shear :"
            if( j.eq.6 ) write(iout,9010) "T11,  plane stress     :"
            if( j.eq.7 ) write(iout,9010) "T11,  plane strain     :"
            if( j.eq.8 ) write(iout,9010) "T13,  anti-plane shear :"
            write(iout,9011) (iterm(k,j),k=1,8)
         end do
         write(iout,9012) evol
      end if
c
      return
c
 9003 format(//,' >>>>> final j-integral contributions for element ',
     & i6,':',
     & /,15x,'jterm1: ',e13.6,'   jterm2: ',e13.6,
     & /,15x,'jterm3: ',e13.6,'   jterm4: ',e13.6,
     & /,15x,'jterm5: ',e13.6,'   jterm6: ',e13.6,
     & /,15x,'jterm7: ',e13.6,'   jterm8: ',e13.6,
     & /,15x,'element volume: ',e13.6 )
 9009 format(//,' >>>>> final i-integral contributions for element ',i7)
 9010 format(/,a)
 9011 format('iterm1: ',e13.6,' iterm2: ',e13.6,' iterm3: ',e13.6,
     &       ' iterm4: ',e13.6,' iterm5: ',e13.6,' iterm6: ',e13.6,
     &       ' iterm7: ',e13.6,' iterm8: ',e13.6 )
 9012 format(//,15x,'element volume: ',e13.6 )
c
      end subroutine dielem_debug_2
c
c***************************************************************
c                                                              *
c      dielem -> contains dielem_debug_1                       *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dielem_debug_1
      implicit none
c
      write(iout,9002) elemno, e, nu, etype, order, ngpts, nnode,
     &                   linear, (rotate(1,i),i=1,3),
     &                   (rotate(2,j),j=1,3), (rotate(3,k), k=1,3)
      write(iout,9001) (qvals(i),i=1,nnode)
      write(iout,9005) process_temperatures
      if( temperatures ) then
          write(iout,9004) (e_node_temps(i),i=1,nnode)
      end if
      if( process_temperatures ) then
        write(iout,9006)
        do enode = 1, nnode
         write(iout,9007) enode, (enode_alpha_ij(i,enode),i=1,6)
        end do
      end if
c
      return
c
 9001 format(' >>>>> supplied values of s (q) function at element',
     &       ' nodes', /,10(/,5x,3f10.3) )
 9002 format(//,
     &       ' >>>>> domain computations for element: ',i7,
     & //,   '       e, nu, etype, order, ngpts, nnode, linear,: ',
     & /,15x, f10.1, f5.2, i5, i5, i5, i5, l5,
     & /,    ' >>>>> global->crack rotation matrix:',
     & 3(/,15x,3f10.6) )
 9004 format(' >>>>> nodal temperatures:',
     &  /,15x,6f10.3,/,15x,6f10.3 )
 9005 format(' >>>>> process temperature loading: ',l1 )
 9006 format(' >>>>> node alpha(1->6) for element in global x-y-z:' )
 9007 format(3x,i3,2x,i3,2x,6e14.6)
c
      end subroutine dielem_debug_1
c
      end subroutine dielem


c *******************************************************************
c *                                                                 *
c *                                                                 *
c *                                                                 *
c *            standalone routines outside dielem                   *
c *                                                                 *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
c
c
c
c
c *******************************************************************
c *                                                                 *
c *    error messages for 3-d edi evaluation                        *
c *                                                                 *
c *******************************************************************
c
       subroutine dieler( iout, ierr, erno )
       implicit none
c
       integer :: iout, ierr, erno
c
c          write an error message
c
c          return error flag
c                      ierr = 0 non-fatal
c                      ierr = 2 terminate element computations
c
      select case( erno )
c
      case( 1 )
      write(iout,1001)
      ierr = 2
c
      case( 2 )
      ierr = 2
c
      case( 3 )
      ierr = 2
c
      case( 4 )
      ierr = 2
c
      case( 5 )
      write (iout,1005)
      ierr = 2
c
      case( 6 )
      ierr = 0
c
      case( 7 )
      write(iout,1007)
      ierr = 0
c
      case( 8 )
      ierr = 0
c
      end select
      return
c
 1001 format(' >>>>> integration order invalid for brick element' )
 1005 format(' >>>>> the determinant of the coordinate jacobian is',
     &/,     '       not positive' )
 1007 format(' >>>>> equivalent loads detected on more than one',
     &/,     '       element face. only the lowest numbered face with',
     &/,     '       loads will be processed in edi computations')
      end
c *******************************************************************
c *                                                                 *
c *   get q, velocity, accelerations, temperature and alpha at      *
c *   integration point all in crack front coords                   *
c *                                                                 *
c *                    last modified: 3/29/2018 rhd                 *
c *                                                                 *
c *******************************************************************
c
      subroutine dielav( sf, evel, eaccel, velocity, x_accel, y_accel,
     &                   z_accel, x_vel, y_vel, z_vel,
     &                   eqvalues, q, nnode, point_temp,
     &                   enode_temps, elem_alpha, enode_alpha_ij,
     &                   linear_displ, fgm_alphas, elemno, ym_nodes,
     &                   nu_nodes, point_ym, point_nu, fgm_e, fgm_nu,
     &                   coord, point_x, point_y, point_z,
     &                   seg_curves_flag, out )
      implicit none
c
c          parameters -- already in crack front coords
c
      integer :: nnode, elemno, mxndel, out
      logical :: linear_displ, fgm_e, fgm_nu, seg_curves_flag,
     &           fgm_alphas
      double precision :: sf(*), evel(3,*), eaccel(3,*), velocity,
     &   x_accel, y_accel, z_accel, eqvalues(*), q, point_temp,
     &   enode_temps(*), elem_alpha(*), enode_alpha_ij(6,*),
     &   ym_nodes(*), nu_nodes(*), point_ym, point_nu, coord(3,*),
     &   point_x, point_y, point_z, x_vel, y_vel, z_vel
c
c          locals
c
      integer :: enode, i, j, k
      logical :: debug
      double precision :: gp_temp, avg_temp,
     &                    gp_alpha, avg_alpha, compl_tens(3,2),
     &                    avg_ym, avg_nu
      double precision, parameter :: zero=0.0d0
c
      debug = .false.
c
      if( debug ) write(out,100) (elem_alpha(i),i=1,6)
c
      x_vel                   = zero
      y_vel                   = zero
      z_vel                   = zero
      x_accel                 = zero
      y_accel                 = zero
      z_accel                 = zero
      q                       = zero
      gp_temp                 = zero
      avg_temp                = zero
      gp_alpha                = zero
      avg_alpha               = zero
      point_ym                = zero
      point_nu                = zero
      avg_ym                  = zero
      avg_nu                  = zero
      point_x                 = zero
      point_y                 = zero
      point_z                 = zero
c
!DIR$ VECTOR ALIGNED
      do enode = 1, nnode
         q         = q         + sf(enode) * eqvalues(enode)
         x_vel     = x_vel     + sf(enode) * evel(1,enode)
         y_vel     = y_vel     + sf(enode) * evel(2,enode)
         z_vel     = z_vel     + sf(enode) * evel(3,enode)
         x_accel   = x_accel   + sf(enode) * eaccel(1,enode)
         y_accel   = y_accel   + sf(enode) * eaccel(2,enode)
         z_accel   = z_accel   + sf(enode) * eaccel(3,enode)
         gp_temp   = gp_temp   + sf(enode) * enode_temps(enode)
         gp_alpha  = gp_alpha  + sf(enode) * enode_alpha_ij(1,enode)
         avg_temp  = avg_temp  + enode_temps(enode)
         avg_alpha = avg_alpha + enode_alpha_ij(1,enode)
         point_ym  = point_ym  + sf(enode) * ym_nodes(enode)
         point_nu  = point_nu  + sf(enode) * nu_nodes(enode)
         avg_ym    = avg_ym    + ym_nodes(enode)
         avg_nu    = avg_nu    + nu_nodes(enode)
         point_x   = point_x   + sf(enode) * coord(1,enode)
         point_y   = point_y   + sf(enode) * coord(2,enode)
         point_z   = point_z   + sf(enode) * coord(3,enode)
      end do
c
c               for linear elements, stiffnesses and stresses were
c               computed using temperatures, alphas and constitutive
c               properties that were constant over the element,
c               being the average of nodal values. thus when J is
c               computed for linear elements, temperature and alpha
c               values will also be constant over the element, being
c               the average of nodal values.
c
      point_temp = gp_temp
      if( linear_displ ) point_temp = avg_temp / dble( nnode )
c
c             for element constant alphas, elem_alpha already
c             contains the right alpha values. only adjust
c             values for fgm or temperature-dependent alphas.
c
      if( fgm_alphas .or. seg_curves_flag ) then
         if( linear_displ ) then
            elem_alpha(1) = avg_alpha / dble( nnode )
            elem_alpha(2) = elem_alpha(1)
            elem_alpha(3) = elem_alpha(1)
            elem_alpha(4) = zero
            elem_alpha(5) = zero
            elem_alpha(6) = zero
         else
            elem_alpha(1) = gp_alpha
            elem_alpha(2) = gp_alpha
            elem_alpha(3) = gp_alpha
            elem_alpha(4) = zero
            elem_alpha(5) = zero
            elem_alpha(6) = zero
         end if
      end if
c
      if( fgm_e .and. linear_displ )
     &    point_ym = avg_ym / dble( nnode )
c
      if( fgm_nu .and. linear_displ )
     &   point_nu = avg_nu / dble( nnode )
c
      if( seg_curves_flag .and. linear_displ ) then
        point_ym = avg_ym / dble( nnode )
        point_nu = avg_nu / dble( nnode )
      end if
c
c               return the scalar velocity of point
c
      velocity = sqrt( x_vel**2 + y_vel**2 + z_vel**2 )
c
      if( debug ) then
         if( linear_displ ) write(out,200)
         if( fgm_alphas ) write(out,300)
         if( fgm_e ) write(out,400)
         if( fgm_nu ) write(out,500)
         if( seg_curves_flag ) write(out,600)
         if( debug ) write(out,700) (elem_alpha(i),i=1,6)
         write(out,800) gp_temp, avg_temp, gp_alpha, avg_alpha,
     &                  point_ym, avg_ym, point_nu, avg_nu
      end if
c
      return
c
 100  format(/,10x,'elem_alpha(1:6) before: ',6(e13.6,2x))
 200  format(/,10x,'linear_displ    = true')
 300  format(/,10x,'fgm_alphas      = true')
 400  format(/,10x,'fgm_e           = true')
 500  format(/,10x,'fgm_nu          = true')
 600  format(/,10x,'seg_curves_flag = true')
 700  format(/,10x,'elem_alpha(1:6) after : ',6(e13.6,2x))
 800  format(/,10x,'gp_temp        : ',e13.6,
     &       /,10x,'avg_temp       : ',e13.6,
     &       /,10x,'gp_alpha       : ',e13.6,
     &       /,10x,'avg_alpha      : ',e13.6,
     &       /,10x,'point_ym       : ',e13.6,
     &       /,10x,'avg_ym         : ',e13.6,
     &       /,10x,'point_nu       : ',e13.6,
     &       /,10x,'avg_nu         : ',e13.6 )
c
      end
c *******************************************************************
c *                                                                 *
c *   rotate nodal vector values to crack x-y-z                     *
c *                                                                 *
c *   last modified: 3/29/2018 rhd                                  *
c *                                                                 *
c *******************************************************************
c
      subroutine dielrv( coord, cdispl, cvel, caccel, edispl, evel,
     &                   eaccel, rotate, nnode, iout, debug )
      implicit none
c
c
c             rotate coordinates, displacements, velocities and
c             accelerations of element nodes from the global
c             system to the crack reference frame.
c
c
      integer :: iout, nnode
      double precision :: coord(3,*), cdispl(3,*), cvel(3,*),
     &                    caccel(3,*), edispl(3,*), evel(3,*),
     &                    eaccel(3,*), rotate(3,3)
      logical :: debug
c
      integer :: inode, i, j
      double precision :: tvec(3)
c
      do inode = 1, nnode
        tvec = coord(1:3,inode)
        call dielrv_a( coord(1,inode), tvec, rotate )
        call dielrv_a( cdispl(1,inode), edispl(1,inode), rotate )
        call dielrv_a( cvel(1,inode),   evel(1,inode),   rotate )
        call dielrv_a( caccel(1,inode), eaccel(1,inode), rotate )
      end do
c
      return
c
      contains
c     ========
c
      subroutine dielrv_a( crack, global, rot )
      implicit none
c
      double precision :: crack(3), global(3), rot(3,3)
c
!DIR$ VECTOR ALIGNED
      crack = matmul( rot, global )
      return
      end subroutine dielrv_a
      end subroutine dielrv
c *******************************************************************
c *                                                                 *
c *  rotate a tensor from global -> crack x-y-z                     *
c *                                                                 *
c *   last modified: 3/30/2018 rhd                                  *
c *                                                                 *
c *******************************************************************
c
c
      subroutine dielrt ( rotate, global, crack )
      implicit none
c
      double precision :: rotate(3,3), global(3,3), crack(3,3),
     &                    space(3,3)
c
      double precision :: trotate(3,3)
c
c            rotation is performed in tensor notation
c
c            crack = rotate * global * rotate(transpose)
c
c            for stresses:
c            global are element global stresses - 1st PK for
c            geonl nonlinear formulation. see digrad for component
c            ordering of 1PK
c
c            for strains:
c            global are element global strains
c
!DIR$ VECTOR ALIGNED
      trotate = transpose( rotate )
!DIR$ VECTOR ALIGNED
      space   = matmul( rotate, global )
!DIR$ VECTOR ALIGNED
      crack   = matmul( space, trotate )
c
      return
c
      end
c
c *******************************************************************
c *                                                                 *
c *    compute coordinate jacobian, its determinate, and inverse    *
c *                                                                 *
c *******************************************************************
c
c
      subroutine dielcj( nxi, neta, nzeta, coord, nnode, cj, cjinv,
     &                   det, ierr, iout, debug )
      implicit none
c
c              compute the 3 x 3 jacobian, its determinate and
c              inverse for the 3-d isoparametrics. nxi, neta,
c              nzeta are the 3 cols of dsf( )
c
c              ierr = 0  ok, = 1 failed
c
      integer :: nnode, ierr, iout
      double precision :: nxi(*), neta(*), nzeta(*), coord(3,*),
     &                    cj(3,3), cjinv(3,3)
c
      integer :: i, j
      double precision :: det, deti
      double precision, external :: dieldp, didet
      double precision, parameter :: zero=0.0d0, one=1.0d0
      logical :: debug
c
c              compute jacobian at the point. use a dot product
c              support function.
c
      cj(1,1) = dieldp( nxi, coord(1,1), nnode, 1, 3 )
      cj(1,2) = dieldp( nxi, coord(2,1), nnode, 1, 3 )
      cj(1,3) = dieldp( nxi, coord(3,1), nnode, 1, 3 )
      cj(2,1) = dieldp( neta, coord(1,1), nnode, 1, 3 )
      cj(2,2) = dieldp( neta, coord(2,1), nnode, 1, 3 )
      cj(2,3) = dieldp( neta, coord(3,1), nnode, 1, 3 )
      cj(3,1) = dieldp( nzeta, coord(1,1), nnode, 1, 3 )
      cj(3,2) = dieldp( nzeta, coord(2,1), nnode, 1, 3 )
      cj(3,3) = dieldp( nzeta, coord(3,1), nnode, 1, 3 )
      if( debug ) write(iout,9001) (( cj(j,i), i = 1,3 ), j=1,3)
c
      det = didet( cj )
      if( det .le. zero ) then
         ierr = 1;  return
      end if
c
      deti = one / det
      cjinv(1,1) =  (cj(2,2) * cj(3,3) - cj(3,2) * cj(2,3)) * deti
      cjinv(2,2) =  (cj(1,1) * cj(3,3) - cj(3,1) * cj(1,3)) * deti
      cjinv(3,3) =  (cj(1,1) * cj(2,2) - cj(2,1) * cj(1,2)) * deti
      cjinv(2,1) = -(cj(2,1) * cj(3,3) - cj(3,1) * cj(2,3)) * deti
      cjinv(3,1) =  (cj(2,1) * cj(3,2) - cj(3,1) * cj(2,2)) * deti
      cjinv(1,2) = -(cj(1,2) * cj(3,3) - cj(3,2) * cj(1,3)) * deti
      cjinv(3,2) = -(cj(1,1) * cj(3,2) - cj(3,1) * cj(1,2)) * deti
      cjinv(1,3) =  (cj(1,2) * cj(2,3) - cj(2,2) * cj(1,3)) * deti
      cjinv(2,3) = -(cj(1,1) * cj(2,3) - cj(2,1) * cj(1,3)) * deti
      if( debug ) write(iout,9002) det,((cjinv(j,i),i=1,3 ), j=1,3)
c
      ierr = 0  ! good result
c
      return
c
 9001 format(/,15x,' jacobian at point', /,3(/,7x,3f15.5) )
 9002 format(15x,' determinant',f15.9, /,
     &       15x,' jacobian inverse', /,3(/,7x,3f15.5) )
c
      end
c
c ******************************************************************
c *                                                                *
c *         interpolate q-values at side nodes                     *
c *                                                                *
c *         if any 1/4-point node is detected, special integration *
c *         for crack-face traction on crack-front element faces   *
c *         is not used.                                           *
c *                                                                *
c *                last modified: 11/20/03                         *
c *                     by: mcw                                    *
c *                                                                *
c ******************************************************************
c
      subroutine dieliq( qvals, debug, out, nnode, etype, elemno,
     &                   snodes, coord, front_elem_flag,
     &                   num_front_nodes, front_nodes, front_coords,
     &                   expanded_front_nodes, domain_type,
     &                   domain_origin, front_order, qp_node,
     &                   crack_curvature, max_exp_front )
      implicit none
c
      double precision :: qvals(*), coord(3,*), front_coords(3,*),
     &                    crack_curvature(*)
      logical :: debug, front_elem_flag, qp_node
      integer :: out, nnode, etype, elemno, snodes(*), num_front_nodes,
     &           max_exp_front, front_nodes(*),
     &           expanded_front_nodes(0:max_exp_front,*),
     &           domain_type, domain_origin, front_order
c
c             local variables
c
      integer :: i, j, k, snode, num_qp_nodes
      logical :: found
      double precision :: r_max, r, t, rs(nnode), node_x,
     &                    node_y, node_z
      double precision, parameter :: zero=0.0d0, toler=1.0d-5,
     &                               half=0.5d0, p4=0.4d0, p75=0.75d0
c
      qp_node      = .false.
      num_qp_nodes = 0
c
c                set q-value of interior nodes by
c                interpolating linearly from corresponding
c                corner nodes for each edge of the element.
c
      select case( etype )
c
c             20 node elements
c
      case( 1 )
        qvals(17) = ( qvals(1) + qvals(5) ) * half
        qvals(18) = ( qvals(2) + qvals(6) ) * half
        qvals(9)  = ( qvals(1) + qvals(2) ) * half
        qvals(10) = ( qvals(2) + qvals(3) ) * half
        qvals(20) = ( qvals(4) + qvals(8) ) * half
        qvals(19) = ( qvals(3) + qvals(7) ) * half
        qvals(12) = ( qvals(4) + qvals(1) ) * half
        qvals(11) = ( qvals(3) + qvals(4) ) * half
        qvals(13) = ( qvals(5) + qvals(6) ) * half
        qvals(14) = ( qvals(6) + qvals(7) ) * half
        qvals(16) = ( qvals(5) + qvals(8) ) * half
        qvals(15) = ( qvals(8) + qvals(7) ) * half
c
c             8-nodes
c
      case( 2 )
        return
c
c             12 node elements
c
      case( 3  )
        qvals(9)  = ( qvals(1) + qvals(2) ) * half
        qvals(10) = ( qvals(2) + qvals(3) ) * half
        qvals(11) = ( qvals(3) + qvals(4) ) * half
        qvals(12) = ( qvals(4) + qvals(1) ) * half
c
c             15 node elements
c
      case( 4 )
        qvals(9)  = ( qvals(1) + qvals(2) ) * half
        qvals(10) = ( qvals(2) + qvals(3) ) * half
        qvals(11) = ( qvals(3) + qvals(4) ) * half
        qvals(12) = ( qvals(4) + qvals(1) ) * half
        qvals(13) = ( qvals(5) + qvals(6) ) * half
        qvals(14) = ( qvals(1) + qvals(5) ) * half
        qvals(15) = ( qvals(2) + qvals(6) ) * half
c
c             9 node elements
c
      case( 5 )
        qvals(9)  = ( qvals(1) + qvals(2) ) * half
c
      case default
        write(out,9000)
        call die_abort
      end select
c
c             set the q-value to 0.75 for 1/4-point nodes.
c             if any 1/4-point node is detected, special integration
c             for crack-face traction on crack-front element faces
c             is not used.
c             skip process when element not on crack front.
c
      if( .not. front_elem_flag ) then
          call dieliq_debug
          return
      end if
c
c             compute distance from each element node to crack front.
c             for curved elements, di_calc_r_theta currently gives an
c             approximation: the distance is measured to a line
c             connecting adjacent crack front nodes.
c
      do i = 1, nnode
         node_x = coord(1,i)
         node_y = coord(2,i)
         node_z = coord(3,i)
         r = zero
         t = zero
         call di_calc_r_theta( 1, front_nodes, num_front_nodes,
     &                      front_coords, domain_type, domain_origin,
     &                      front_order, node_x, node_y, node_z,
     &                      elemno, i, r, t, crack_curvature,
     &                      debug, out )
         rs(i) = r
      end do
c
c             find node farthest from front
c
      r_max = zero
      do i = 1, nnode
         r = rs(i)
         if( r.gt.r_max ) r_max = r
      end do
c
c             identify 1/4-point elements as follows:
c                1. node must not be on crack front.
c                2. q-value must currently be 0.5
c
      do i = 1, nnode
         snode = snodes(i)
c
c             skip node if it's on the crack front
c
        found = .false.
        do j = 1, num_front_nodes
            if( found ) exit
            do k = 1, max_exp_front
               if( expanded_front_nodes(k,j) .eq. snode ) then
                  found = .true.
                  exit
               end if
            end do
        end do
        if( found ) cycle
c
        if( abs(qvals(i) - half) .lt. toler ) then
c
c              the current node is a mid-side node on an element edge
c              'normal' to the crack front.
c
c              if the distance from current node to the crack front
c              is less than 0.4 times the maximum distance of any
c              element node to the crack front (approximately measured
c              using straight lines), consider the node as a 1/4-point
c              node. this identification will work when element
c              distortion or curvature is not too large.
c
            if( rs(i) .lt. r_max * p4 ) then
               qvals(i)     = p75
               qp_node      = .true.
               num_qp_nodes = num_qp_nodes + 1
            end if
c
         end if
      end do ! over i
c
      call dieliq_debug
c
      return
c
 9000 format('>> FATAL ERROR: bad code dieliq',//)
c
      contains
c     ========
c
      subroutine dieliq_debug
      implicit none
c
      if( debug ) then
         if( front_elem_flag ) write(out,9012) r_max, num_qp_nodes
         write(out,9014) elemno
         do i=1,nnode
            write(out,9016) i, snodes(i), qvals(i)
         end do
      end if
      return
c
 9012 format(/,' >>>>> max. elem. node dist. from front:     ',e13.6,
     &       /,' >>>>> number of 1/4-point nodes on element: ',i7)
 9014 format(' >>>>> final nodal q-function values for element',
     &       2x,i7,/,5x,'node',5x,'snode',8x,'q')
 9016 format(5x,i2,5x,i7,2x,f10.3)
c
      end subroutine dieliq_debug
c
      end subroutine dieliq
c
c *******************************************************************
c *                                                                 *
c *   build 6x6 rotation matrix to convert unrotated cauchy         *
c *   stresses (6x1) to cauchy stresses (6x1)                       *
c *                                                                 *
c *   last updated: 3/12/2018 rhd                                   *
c *   code extracted from getrm1 in gtmat1.f                        *
c *                                                                 *
c *******************************************************************
c
      subroutine digetr( q, r  )
      implicit none
c
c           parameters
c
      double precision :: q(6,6), r(3,3)
c
c           locals
c
      double precision, parameter :: two=2.0d0
c
c       cauchy stress {T} = [q] * (unrotated) cauchy stress {t}.
c       in tensor form:
c
c                 [T] = [r] [t] trans([r])
c
c       both [T] and [t] are symmetric, [r] is orthogonal rotation.
c       vector ordering is {x,y,z,xy,yz,xz}.
c
      q(1,1) = r(1,1)**2
      q(1,2) = r(1,2)**2
      q(1,3) = r(1,3)**2
      q(1,4) = two*r(1,1)*r(1,2)
      q(1,5) = two*r(1,3)*r(1,2)
      q(1,6) = two*r(1,1)*r(1,3)
      q(2,1) = r(2,1)**2
      q(2,2) = r(2,2)**2
      q(2,3) = r(2,3)**2
      q(2,4) = two*r(2,1)*r(2,2)
      q(2,5) = two*r(2,3)*r(2,2)
      q(2,6) = two*r(2,1)*r(2,3)
      q(3,1) = r(3,1)**2
      q(3,2) = r(3,2)**2
      q(3,3) = r(3,3)**2
      q(3,4) = two*r(3,1)*r(3,2)
      q(3,5) = two*r(3,3)*r(3,2)
      q(3,6) = two*r(3,1)*r(3,3)
      q(4,1) = r(1,1)*r(2,1)
      q(4,2) = r(1,2)*r(2,2)
      q(4,3) = r(1,3)*r(2,3)
      q(4,4) = r(1,1)*r(2,2)+r(2,1)*r(1,2)
      q(4,5) = r(1,2)*r(2,3)+r(1,3)*r(2,2)
      q(4,6) = r(1,1)*r(2,3)+r(1,3)*r(2,1)
      q(5,1) = r(2,1)*r(3,1)
      q(5,2) = r(3,2)*r(2,2)
      q(5,3) = r(2,3)*r(3,3)
      q(5,4) = r(2,1)*r(3,2)+r(2,2)*r(3,1)
      q(5,5) = r(2,2)*r(3,3)+r(3,2)*r(2,3)
      q(5,6) = r(2,1)*r(3,3)+r(2,3)*r(3,1)
      q(6,1) = r(1,1)*r(3,1)
      q(6,2) = r(1,2)*r(3,2)
      q(6,3) = r(1,3)*r(3,3)
      q(6,4) = r(1,1)*r(3,2)+r(1,2)*r(3,1)
      q(6,5) = r(1,2)*r(3,3)+r(1,3)*r(3,2)
      q(6,6) = r(1,1)*r(3,3)+r(3,1)*r(1,3)
c
      return
      end
c *******************************************************************
c *                                                                 *
c *   rotate unrotated cauchy stresses (6x1) to global              *
c *   cauchy stresses in tensor form                                *
c *                                                                 *
c *******************************************************************
c
      subroutine diqmp1( q, unrotated_stress, cauchy_tensor )
      implicit none
c
c                    parameter declarations
c
      double precision :: q(6,6), unrotated_stress(6),
     &                    cauchy_tensor(3,3)
c
c                    locals
c
      double precision :: cauchy_vector(6)
c
c            {cstress} = [q] * {stress}; 6x1 vectors and 6x6 q
c             (Cauchy)         (u.r. Cauchy)
c
!DIR$ VECTOR ALIGNED
      cauchy_vector = matmul( q, unrotated_stress )
c
c            store cauchy stresses in 3x3 tensor form.
c
      cauchy_tensor(1,1) = cauchy_vector(1) ! xx
      cauchy_tensor(2,1) = cauchy_vector(4) ! xy, yx
      cauchy_tensor(3,1) = cauchy_vector(6) ! xz, zx
      cauchy_tensor(1,2) = cauchy_vector(4) ! xy, yx
      cauchy_tensor(2,2) = cauchy_vector(2) ! yy
      cauchy_tensor(3,2) = cauchy_vector(5) ! yz, zy
      cauchy_tensor(1,3) = cauchy_vector(6) ! xz, zx
      cauchy_tensor(2,3) = cauchy_vector(5) ! yz, zy
      cauchy_tensor(3,3) = cauchy_vector(3) ! zz
c
      return
      end
c *******************************************************************
c *                                                                 *
c *   form deformation gradient, determinant, inverse transpose     *
c *   at a gauss point. convert Cauchy stress to 1st PK stress      *
c *   scale energy density to volume density at t=0                 *
c *                                                                 *
c *   last modified: 4/22/2018  rhd                                 *
c *                                                                 *
c *******************************************************************
c
      subroutine digrad( sig, displ_grad, iout, gdebug, ptno, elemno,
     &                   gp_det_F )
c
      implicit none
c
c                    parameters
c                     sig on input is Cauchy stress tensor in global
c                     coordinates + work density.
c                     output: 1PK stress tensor + W for t=0 volume dens
c
      integer :: iout, ptno, elemno
      double precision ::  sig(10), gp_det_F, displ_grad(3,3)
      logical :: gdebug
c
c                    locals
c
      logical :: debug
      double precision :: det_F, F(3,3), pk1(3,3),
     &                    pk1vec(9), deti, sigtens(3,3), sigvec(9),
     &                    Fit(3,3)
      double precision, parameter :: zero=0.0d0, one=1.0d0
      double precision, external :: didet
      equivalence( sigtens, sigvec ), ( pk1, pk1vec )
c
c             X -> coordinates at t = 0
c
      debug = gdebug
c
!DIR$ VECTOR ALIGNED
      F = displ_grad
      F(1,1) = F(1,1) + one
      F(2,2) = F(2,2) + one
      F(3,3) = F(3,3) + one
c
c             compute determinant and inverse *transpose*
c             of 3 x 3 matrix.
c             symmetry of input array not required.
c             use method of cofactors.
c
      det_F = didet( F )
      gp_det_F = det_F
      if( det_F .le. zero ) then
         write(iout,9200) elemno, ptno, det_F
         call die_abort
      end if
c
c             omit / by det_F since the Cauchy -> 1PK has * det_F
c
      Fit(1,1) =  (F(2,2) * F(3,3) - F(3,2) * F(2,3))
      Fit(2,2) =  (F(1,1) * F(3,3) - F(3,1) * F(1,3))
      Fit(3,3) =  (F(1,1) * F(2,2) - F(2,1) * F(1,2))
      Fit(1,2) = -(F(2,1) * F(3,3) - F(3,1) * F(2,3))
      Fit(1,3) =  (F(2,1) * F(3,2) - F(3,1) * F(2,2))
      Fit(2,1) = -(F(1,2) * F(3,3) - F(3,2) * F(1,3))
      Fit(2,3) = -(F(1,1) * F(3,2) - F(3,1) * F(1,2))
      Fit(3,1) =  (F(1,2) * F(2,3) - F(2,2) * F(1,3))
      Fit(3,2) = -(F(1,1) * F(2,3) - F(2,1) * F(1,3))
c
c              1st PK = det * Cauchy * trans(fi)
c              ordering of non-symmetric 1PK_{ij} -> stress in i
c                   direction on +face j
c               xx xy xz
c               yx yy yz
c               zx zy zz
c
      sigvec   = sig(1:9)     ! see equivalence
      pk1      = matmul( sigtens, Fit )
      sig(1:9) = pk1vec(1:9) ! see equivalence
c
c              scale work density to volume at t=0
c
      sig(10) = sig(10) * det_F
c
      if( .not. debug ) return
      write(iout,9110) elemno, ptno
      write(iout,*) '>> inside digrad: '
      write(iout,*) '   > deformation gradient:'
      write(iout,9000) f(1,1), f(1,2), f(1,3),
     &                 f(2,1), f(2,2), f(2,3),
     &                 f(3,1), f(3,2), f(3,3)
      write(iout,*) '   > transpose( inv(F):'
      write(iout,9000) fit(1,1), fit(1,2), fit(1,3),
     &                 fit(2,1), fit(2,2), fit(2,3),
     &                 fit(3,1), fit(3,2), fit(3,3)
      write(iout,*) '   > Cauchy stresses:'
      write(iout,9000) sig(1:9)
      write(iout,*) '   > 1 st PK stresses:'
      write(iout,9000) pk1(1,1), pk1(1,2), pk1(1,3),
     &                 pk1(2,1), pk1(2,2), pk1(2,3),
     &                 pk1(3,1), pk1(3,2), pk1(3,3)
      write(iout,9100) det_F, sig(10), sig(10)/det_F
c
      return
c
 9000 format(3(5x,3f20.6,/))
 9100 format(5x,'det_F: ',f20.6,/,5x,'W(0):',f20.6,
     &      /,5x,'W(n):',f20.6 )
 9110 format(5x,'   .. get 1PK stresses, digrad. elem, pt: ',i3,i10)
 9200 format('>> FATAL ERROR: digrad. bad det[F]. ptno, elem: ',i3,i8,
     &  /,   '                det: ',d14.6,//)
c
      end

c *******************************************************************
c *                                                                 *
c *  set up to process temperature effects for domain integral      *
c *                                                                 *
c *   last modified: 3/12/2018 rhd                                  *
c *                                                                 *
c *******************************************************************
c
c
      subroutine dippie( rotate, iout, debug, elem_alpha )
      implicit none
c
c                    parameters
c
      integer :: iout
      double precision :: rotate(3,3),  elem_alpha(*)
      logical ::  debug
c
c                    locals
c
      double precision :: work(3,3), tensor(3,3), trot(3,3)
      double precision, parameter :: zero=0.0d0, half=0.5d0
c
c                    the 6 alpha values for thermal expansion define
c                    a symmetric 3x3 tensor to define initial
c                    strains due to a temperature change. rotate the
c                    alphas to crack coordinates and put back into
c                    the 6x1 elem_alpha vector. The tensor alphas
c                    transform the same as strains.
c
c                    alpha order: xx, yy, zz, xy, yz, xz
c
      tensor(1,1) = elem_alpha(1)
      tensor(2,2) = elem_alpha(2)
      tensor(3,3) = elem_alpha(3)
      tensor(1,2) = elem_alpha(4) * half
      tensor(2,1) = tensor(1,2)
      tensor(2,3) = elem_alpha(5) * half
      tensor(3,2) = tensor(2,3)
      tensor(1,3) = elem_alpha(6) * half
      tensor(3,1) = tensor(1,3)
c
c            rotation is performed in tensor notation
c              crack_tensor = rotate * tensor * rotate(transpose)
c
      trot   = transpose( rotate )
      work   = matmul( rotate, tensor )
      tensor = matmul( work, trot )
c
      elem_alpha(1) =  tensor(1,1)
      elem_alpha(2) =  tensor(2,2)
      elem_alpha(3) =  tensor(3,3)
      elem_alpha(4) =  tensor(1,2) / half
      elem_alpha(5) =  tensor(2,3) / half
      elem_alpha(6) =  tensor(1,3) / half
c
      return
      end
c
c   *******************************************************************
c   *                                                                 *
c   * support routine:  dot product: non-unit stride                  *
c   *                                                                 *
c   *******************************************************************
c
      function dieldp( veca, vecb, nterms, stepa, stepb ) result( res )
c
      implicit none
c
      double precision :: veca(*), vecb(*), res
      integer :: stepa, stepb, nterms
c
      integer :: indexa, indexb, i
      double precision, parameter :: zero=0.0d0
c
      indexa = 1
      indexb = 1
      res = zero
c
      if( stepa .eq. 1 .and. stepb .eq. 1 ) then
!DIR$ VECTOR ALIGNED
       do i = 1, nterms
          res = res + veca(i)*vecb(i)
        end do
      else
        indexa = 1
        indexb = 1
!DIR$ VECTOR ALIGNED
        do i = 1, nterms
          res = res + veca(indexa)*vecb(indexb)
          indexa = indexa + stepa
          indexb = indexb + stepb
        end do
      end if
c
      return
      end
c

c   *******************************************************************
c   *                                                                 *
c   * support routine:  determinaant of 3x3 matrix                    *
c   *                                                                 *
c   *******************************************************************
c
      function didet( a ) result( det )
c
      implicit none
c
      double precision :: a(3,3), det
c
      det =   a(1,1) * a(2,2) * a(3,3) + a(2,1) * a(3,2) * a(1,3)
     &      + a(3,1) * a(1,2) * a(2,3) - a(1,1) * a(3,2) * a(2,3)
     &      - a(2,1) * a(1,2) * a(3,3) - a(3,1) * a(2,2) * a(1,3)
c
      return
      end

