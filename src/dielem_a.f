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
c *                                                              *              
c * dielem calling tree:                                         *              
c *                                                              *              
c *    (dielem_a.f)                                              *              
c *       -dielem                                                *              
c *           -digetr                                            *              
c *           -diqmp1                                            *              
c *           -getgpts                                           *              
c *           -derivs                                            *              
c *           -dielcj                                            *              
c *           -dieler                                            *              
c *           -digrad                                            *              
c *           -dielrv                                            *              
c *           -dieliq                                            *              
c *           -dielrt                                            *              
c *           -dippie                                            *              
c *           -shapef                                            *              
c *           -dielav                                            *              
c *           -di_calc_j_terms                                   *              
c *                                                              *              
c *    (dielem_b.f)                                              *              
c *           -di_calc_r_theta                                   *              
c *           -di_calc_constitutive                              *              
c *           -di_calc_aux_fields_k                              *              
c *           -di_calu_aux_fields_t                              *              
c *           -di_calc_i_terms                                   *              
c *                                                              *              
c *    (dielem_c.f)                                              *              
c *           -di_calc_surface_integrals                         *              
c *               -dielwf                                        *              
c *               -dielrl                                        *              
c *               -di_calc_surface_j                             *              
c *               -di_calc_surface_i                             *              
c *                   -di_reorder_nodes                          *              
c *                                                              *              
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
      subroutine dielem ( coord, qvals, edispl, evel, eaccel,                   
     &                    feload, estress, history,                             
     &                    props, iprops, lprops, erots, nrow_hist,              
     &                    rotate, iout, e_jresults, e_node_temps,               
     &                    temperatures, enode_alpha_ij, ierr, elemno,           
     &                    gdebug, one_point, geonl, numrows_stress,             
     &                    snodes, elem_nod_swd, elem_nod_strains,               
     &                    omit_ct_elem, fgm_e, fgm_nu, e_strain,                
     &                    ym_nodes, ym_front_node, nu_nodes,                    
     &                    nu_front_node, front_nodes, num_front_nodes,          
     &                    front_coords, domain_origin, domain_type,             
     &                    front_order, e_iresults, comput_i, comput_j,          
     &                    cf_traction_flags, cf_tractions,                      
     &                    front_elem_flag, expanded_front_nodes,                
     &                    myid, numprocs, crack_curvature, face_loading,        
     &                    seg_curves_flag, process_temperatures,                
     &                    max_exp_front )                                       
c                                                                               
      implicit double precision ( a-h, o-z )                                    
c                                                                               
c                    parameter declarations                                     
c                                                                               
      dimension  coord(3,*), qvals(*), feload(3,*), edispl(3,*),                
     &           evel(3,*), eaccel(3,*), estress(*), e_node_temps(*),           
     &           history(nrow_hist,*), rotate(3,3), e_jresults(*),              
     &           erots(9,*), enode_alpha_ij(6,*), elem_nod_swd(*),              
     &           elem_nod_strains(6,*), iprops(*), e_strain(9,*),               
     &           ym_nodes(*), front_coords(3,*), e_iresults(8,8),               
     &           cf_tractions(*), crack_curvature(*)                            
      double precision                                                          
     &  nu_nodes(*), nu_front_node                                              
      real       props(*)                                                       
      logical    lprops(*), one_point, geonl, gdebug, temperatures,             
     &           omit_ct_elem, fgm_e, fgm_nu, comput_i, comput_j,               
     &           cf_traction_flags(*), front_elem_flag, seg_curves_flag,        
     &           process_temperatures, face_loading                             
      integer    elemno, snodes(*), front_nodes(*), num_front_nodes,            
     &           domain_origin, domain_type, front_order, iout,                 
     &           expanded_front_nodes(0:max_exp_front,*), myid, numprocs        
c                                                                               
c                    local declarations                                         
c                                                                               
      dimension                                                                 
     &   csig(10,27), cdispl(3,20), cvel(3,20), caccel(3,20),                   
     &   jacob(3,3,28), jacobi(3,3,28), dsf(20,3,28), sf(20),                   
     &   eaverage(19),                                                          
     &   qn(6,6), stresses(10,28), detvol(28), elem_alpha(6),                   
     &   dalpha_x1(6), dstrain_x1(9), strains(9,20),                            
     &   aux_stress(9,8), aux_strain(9,8), daux_strain_x1(9,8),                 
     &   ceps_n(9,20), ceps_gp(9,27), du11_aux(8), du21_aux(8),                 
     &   du31_aux(8), du111_aux(8), du112_aux(8),                               
     &   du113_aux(8), du211_aux(8), du212_aux(8),                              
     &   du213_aux(8), du311_aux(8), du312_aux(8),                              
     &   du313_aux(8), dcijkl_x1(3), sijkl(3), dsijkl_x1(3)                     
      double precision                                                          
     &  e, nu, jacob, jacobi, nx, ny, nz, kin_energy, lg(28), four,             
     &  jterm(8), iterm(8,8), r, t, toler, termu_aux, termv_aux,                
     &  termw_aux, gpq                                                          
      real neg_99, fgm_tol                                                      
      integer  etype, ptno, enode, fnode, order, faceno, flag, nfnode,          
     &         order_face, gpn, i, j, k                                         
      logical  debug, linear, linear_displ, debug_i,                            
     &         debug_j, qp_node, fgm_alphas                                     
      data  neg_99, fgm_tol / -99.0, 1.0 /                                      
      data zero, half, two, four, eight                                         
     & / 0.d0, 0.5d0, 2.d0, 4.0d0, 8.0d0 /                                      
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
      e        = props(7)                                                       
      nu       = props(8)                                                       
      etype    = iprops(1)                                                      
      linear_displ = etype .eq. 2  .or.                                         
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
c                                                                               
      elem_alpha(1) = props(9)                                                  
      elem_alpha(2) = props(13)                                                 
      elem_alpha(3) = props(34)                                                 
      elem_alpha(4) = props(35)                                                 
      elem_alpha(5) = props(36)                                                 
      elem_alpha(6) = props(37)                                                 
c                                                                               
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
      rho       = props(10)                                                     
c                                                                               
      e_jresults(1:8)     = zero                                                
      e_iresults(1:8,1:8) = zero                                                
c                                                                               
      if ( debug ) then                                                         
        write(iout,9002) elemno, e, nu, etype, order, ngpts, nnode,             
     &                   linear, (rotate(1,i),i=1,3),                           
     &                   (rotate(2,j),j=1,3), (rotate(3,k), k=1,3)              
        write(iout,9001) (qvals(i),i=1,nnode)                                   
        write(iout,9005) process_temperatures                                   
        if ( temperatures ) then                                                
          write(iout,9004) (e_node_temps(i),i=1,nnode)                          
        end if                                                                  
        if ( process_temperatures ) then                                        
          write(iout,9006)                                                      
          do enode = 1, nnode                                                   
           write(iout,9007) enode, (enode_alpha_ij(i,enode),i=1,6)              
          end do                                                                
        end if                                                                  
      end if                                                                    
c                                                                               
c             set number of stress values per strain point in arrays            
c             passed by warp                                                    
c                                                                               
      numrow = numrows_stress                                                   
c                                                                               
c             set up data needed at all integration points and nodes            
c             needed for the volume integrals and crack face                    
c             traction integrals.                                               
c                                                                               
      if ( debug ) write(iout,9013)                                             
c                                                                               
c        1:   transform unrotated Cauchy stresses to Cauchy                     
c             stresses for geonl formulation. estress on input contains         
c             unrotated stresses. first convert to Cauchy stresses              
c             then to 1 st PK stresses. these are in global coordinates.        
c             note: 1 PK stresses are non-symmetric so we keep the              
c             3x3 tensor stored as a vector. position 10 is for the             
c             work density. the work density is scaled by det F to refer        
c             to t=0 configuration.                                             
c                                                                               
      if ( geonl ) then                                                         
         do ptno = 1, ngpts                                                     
          loc = ( ptno-1 ) * numrow + 1                                         
          call digetr( qn, erots(1,ptno) )                                      
          call diqmp1( qn, estress(loc), stresses(1,ptno) )                     
          stresses(10,ptno) = estress(loc+6)                                    
          call getgpts( etype, order, ptno, xsi, eta, zeta, weight )            
          call derivs( etype, xsi, eta, zeta, dsf(1,1,ptno),                    
     &                 dsf(1,2,ptno), dsf(1,3,ptno) )                           
          call dielcj( dsf(1,1,ptno), coord, nnode, jacob(1,1,ptno),            
     &                 jacobi(1,1,ptno), detvol(ptno), ierr, iout,              
     &                 debug )                                                  
          if ( ierr .ne. 0 ) then                                               
           call dieler( iout,ierr,5 )                                           
           return                                                               
          end if                                                                
          call digrad( nnode, dsf(1,1,ptno), jacobi(1,1,ptno),                  
     &                 edispl, stresses(1,ptno), iout, debug )                  
         end do                                                                 
      end if                                                                    
c                                                                               
c        2.  for small-displacement formulation, copy the input                 
c            (symmetric) stresses (6x1) into 3x3 tensor form.                   
c                                                                               
      if ( .not. geonl ) then                                                   
        do ptno = 1, ngpts                                                      
          loc = ( ptno-1 )  * numrow                                            
          stresses(1,ptno)  = estress(loc+1)                                    
          stresses(2,ptno)  = estress(loc+4)                                    
          stresses(3,ptno)  = estress(loc+6)                                    
          stresses(4,ptno)  = estress(loc+4)                                    
          stresses(5,ptno)  = estress(loc+2)                                    
          stresses(6,ptno)  = estress(loc+5)                                    
          stresses(7,ptno)  = estress(loc+6)                                    
          stresses(8,ptno)  = estress(loc+5)                                    
          stresses(9,ptno)  = estress(loc+3)                                    
          stresses(10,ptno) = estress(loc+7)                                    
        end do                                                                  
      end if                                                                    
c                                                                               
c        3.   copy engineering strains into 3x3 tensor form.                    
c             change gamma strains to tensor strains for                        
c             off-diagonal terms. this formulation is not                       
c             correct for large strains.                                        
c                                                                               
      do enode = 1, nnode                                                       
         strains(1,enode) = elem_nod_strains(1,enode)                           
         strains(2,enode) = elem_nod_strains(4,enode)/2.0                       
         strains(3,enode) = elem_nod_strains(6,enode)/2.0                       
         strains(4,enode) = elem_nod_strains(4,enode)/2.0                       
         strains(5,enode) = elem_nod_strains(2,enode)                           
         strains(6,enode) = elem_nod_strains(5,enode)/2.0                       
         strains(7,enode) = elem_nod_strains(6,enode)/2.0                       
         strains(8,enode) = elem_nod_strains(5,enode)/2.0                       
         strains(9,enode) = elem_nod_strains(3,enode)                           
      end do                                                                    
c                                                                               
c        4:   rotate nodal coordinates, displacements, velocities               
c             and accelerations to crack front normal coordinate                
c             system.                                                           
c                                                                               
      call dielrv( coord, cdispl, cvel, caccel, edispl, evel, eaccel,           
     &             rotate, nnode, iout, debug )                                 
c                                                                               
c        5:   modify the nodal q-values to reflect linear                       
c             interpolations for 20, 15, 12-node elements.                      
c                                                                               
      call dieliq( qvals, debug, iout, nnode, etype, elemno, snodes,            
     &             coord, front_elem_flag, num_front_nodes, front_nodes,        
     &             front_coords, expanded_front_nodes, domain_type,             
     &             domain_origin, front_order, qp_node, crack_curvature,        
     &             max_exp_front)                                               
c                                                                               
c        6:   compute coordinate jacobian at all gauss points                   
c             and the center point of element. these are now in                 
c             cracked coordinates.                                              
c                                                                               
      do ptno = 1, ngpts + 1                                                    
c                                                                               
c             6a: isoparametric coordinates of gauss point                      
c                                                                               
        if ( ptno .eq. ngpts+1 ) then                                           
           xsi  = zero                                                          
           eta  = zero                                                          
           zeta = zero                                                          
        else                                                                    
           call getgpts( etype, order, ptno, xsi, eta, zeta, weight )           
        end if                                                                  
        call derivs( etype, xsi, eta, zeta, dsf(1,1,ptno),                      
     &               dsf(1,2,ptno), dsf(1,3,ptno) )                             
        if ( debug ) write(iout,9110)  ptno, xsi, eta, zeta, weight             
c                                                                               
c             6b: coordinate jacobian, it's inverse, and determinant.           
c                                                                               
        call dielcj( dsf(1,1,ptno), coord, nnode, jacob(1,1,ptno),              
     &               jacobi(1,1,ptno), detvol(ptno), ierr, iout, debug )        
        if ( ierr .ne. 0 ) then                                                 
           call dieler( iout,ierr,5 )                                           
           return                                                               
        end if                                                                  
      end do                                                                    
c                                                                               
c        7:  rotate stresses and strains to crack-front coordinates             
c             at all gauss points. get the stress work density from             
c             warp results.                                                     
c                                                                               
      do ptno = 1, ngpts                                                        
         call dielrt( ptno, rotate, stresses(1,ptno), csig(1,ptno),             
     &               1, iout, debug )                                           
         csig(10,ptno) = stresses(10,ptno)                                      
         call dielrt( ptno, rotate, e_strain(1,ptno), ceps_gp(1,ptno),          
     &               2, iout, debug )                                           
      end do                                                                    
c                                                                               
c        8:  rotate strains to crack-front coordinates at all nodes.            
c                                                                               
      do enode = 1, nnode                                                       
         call dielrt( enode, rotate, strains(1,enode), ceps_n(1,enode),         
     &        3, iout, debug )                                                  
      end do                                                                    
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
      if ( process_temperatures ) then                                          
         call dippie( rotate, iout, debug, elem_alpha )                         
         do enode = 1, nnode                                                    
            call dippie( rotate, iout, debug, enode_alpha_ij(1,enode))          
         end do                                                                 
      end if                                                                    
      if ( debug ) then                                                         
         write(iout,9008)                                                       
         do enode = 1, nnode                                                    
            write(iout,9007) enode, snodes(enode),                              
     &                       (enode_alpha_ij(i,enode),i=1,6)                    
         end do                                                                 
      end if                                                                    
c                                                                               
c        10:   for single point integration of bricks,                          
c              compute average stresses, energy density, and                    
c              strain at center of element. set number of gauss                 
c              points to one to control subsequent loop.                        
c                                                                               
      if ( one_point ) then                                                     
        numipts         = 1                                                     
        eaverage(1:19)  = zero                                                  
        do ptno = 1, ngpts                                                      
           do i = 1, 10                                                         
              eaverage(i) = eaverage(i) + csig(i,ptno)                          
           end do                                                               
           do i = 11,19                                                         
              eaverage(i) = eaverage(i) + ceps_gp(i,ptno)                       
           end do                                                               
        end do                                                                  
        do i = 1, 10                                                            
       csig(i,ngpts+1) = eaverage(i) / dble(ngpts)                              
        end do                                                                  
        do i = 11, 19                                                           
       ceps_gp(i,ngpts+1) = eaverage(i) / dble(ngpts)                           
        end do                                                                  
      end if                                                                    
c                                                                               
c             evaluate volume integrals for domain                              
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
      jterm(1:8)     = zero                                                     
      iterm(1:8,1:8) = zero                                                     
      evol           = zero                                                     
      if ( debug ) write(iout,9100)                                             
c                                                                               
      do ii = 1, numipts                                                        
c                                                                               
c             isoparametric coordinates of gauss point and weight.              
c                                                                               
        ptno = ii                                                               
        if ( one_point ) then                                                   
           xsi    = zero                                                        
           eta    = zero                                                        
           zeta   = zero                                                        
           weight = eight                                                       
           ptno   = ngpts + 1                                                   
        else                                                                    
           call getgpts( etype, order, ptno, xsi, eta, zeta, weight )           
        end if                                                                  
        if ( debug ) write(iout,9110)  ptno, xsi, eta, zeta, weight             
c                                                                               
c             evaluate shape functions of all nodes at this point.              
c             used to find value of q function at integration point,            
c             the total velocity at the point, values of acceleration           
c             at the point, and alpha values at the point.                      
c                                                                               
        call shapef( etype, xsi, eta, zeta, sf )                                
        call dielav( sf, evel, eaccel, point_velocity, point_accel_x,           
     &               point_accel_y, point_accel_z, qvals, point_q,              
     &               nnode, point_temp, e_node_temps, elem_alpha,               
     &               enode_alpha_ij, linear_displ, fgm_alphas, elemno,          
     &               ym_nodes, nu_nodes, point_ym, point_nu, fgm_e,             
     &               fgm_nu, coord, point_x, point_y, point_z,                  
     &               seg_curves_flag, iout )                                    
c                                                                               
        kin_energy = half * rho * point_velocity * point_velocity               
        if ( debug ) then                                                       
          write(iout,9120) point_q, point_velocity, kin_energy,                 
     &                     point_accel_x, point_accel_y, point_accel_z,         
     &                     point_temp                                           
        endif                                                                   
c                                                                               
c             for this integration point, compute displacement                  
c             derivatives, q-function derivatives and temperature               
c             derivative in crack coordinate system.                            
c             strains will be used in the order eps11, eps22, eps33,            
c             eps12, eps23, eps13                                               
c                                                                               
        dux              = zero                                                 
        dvx              = zero                                                 
        dwx              = zero                                                 
        dqx              = zero                                                 
        dqy              = zero                                                 
        dqz              = zero                                                 
        dtx              = zero                                                 
        dalpha_x1(1:6)   = zero                                                 
        dswd_x1          = zero                                                 
        dstrain_x1(1:9)  = zero                                                 
        csig_dstrain_x1  = zero                                                 
        de_x1            = zero                                                 
        dnu_x1           = zero                                                 
c                                                                               
        do enode = 1, nnode                                                     
           nx = dsf(enode,1,ptno) * jacobi(1,1,ptno) +                          
     &          dsf(enode,2,ptno) * jacobi(1,2,ptno) +                          
     &          dsf(enode,3,ptno) * jacobi(1,3,ptno)                            
           ny = dsf(enode,1,ptno) * jacobi(2,1,ptno) +                          
     &          dsf(enode,2,ptno) * jacobi(2,2,ptno) +                          
     &          dsf(enode,3,ptno) * jacobi(2,3,ptno)                            
           nz = dsf(enode,1,ptno) * jacobi(3,1,ptno) +                          
     &          dsf(enode,2,ptno) * jacobi(3,2,ptno) +                          
     &          dsf(enode,3,ptno) * jacobi(3,3,ptno)                            
           dux = dux + nx * cdispl(1,enode)                                     
           dvx = dvx + nx * cdispl(2,enode)                                     
           dwx = dwx + nx * cdispl(3,enode)                                     
           dqx = dqx + nx * qvals(enode)                                        
           dqy = dqy + ny * qvals(enode)                                        
           dqz = dqz + nz * qvals(enode)                                        
           dtx = dtx + nx * e_node_temps(enode)                                 
           dswd_x1 = dswd_x1 + nx * elem_nod_swd(enode)                         
           de_x1   = de_x1  + nx * ym_nodes(enode)                              
           dnu_x1  = dnu_x1 + nx * nu_nodes(enode)                              
           dstrain_x1(1) = dstrain_x1(1) + nx * ceps_n(1,enode)                 
           dstrain_x1(2) = dstrain_x1(2) + nx * ceps_n(2,enode)                 
           dstrain_x1(3) = dstrain_x1(3) + nx * ceps_n(3,enode)                 
           dstrain_x1(4) = dstrain_x1(4) + nx * ceps_n(4,enode)                 
           dstrain_x1(5) = dstrain_x1(5) + nx * ceps_n(5,enode)                 
           dstrain_x1(6) = dstrain_x1(6) + nx * ceps_n(6,enode)                 
           dstrain_x1(7) = dstrain_x1(7) + nx * ceps_n(7,enode)                 
           dstrain_x1(8) = dstrain_x1(8) + nx * ceps_n(8,enode)                 
           dstrain_x1(9) = dstrain_x1(9) + nx * ceps_n(9,enode)                 
c                                                                               
           if ( process_temperatures ) then                                     
              dalpha_x1(1) = dalpha_x1(1) + nx * enode_alpha_ij(1,enode)        
              dalpha_x1(2) = dalpha_x1(2) + nx * enode_alpha_ij(2,enode)        
              dalpha_x1(3) = dalpha_x1(3) + nx * enode_alpha_ij(3,enode)        
              dalpha_x1(4) = dalpha_x1(4) + nx * enode_alpha_ij(4,enode)        
              dalpha_x1(5) = dalpha_x1(5) + nx * enode_alpha_ij(5,enode)        
              dalpha_x1(6) = dalpha_x1(6) + nx * enode_alpha_ij(6,enode)        
           end if                                                               
        end do                                                                  
c                                                                               
c             combine terms for simplicity                                      
c                                                                               
        weight = weight * detvol(ptno)                                          
c                                                                               
c                                                                               
c ----------- call routines to compute domain-integral terms ------------       
c ------------         for current integration point.           -----------     
c                                                                               
c                                                                               
        if( comput_j ) then                                                     
           call di_calc_j_terms(weight, evol, jterm, ptno, dqx,                 
     &                       dqy, dqz, dux, dvx, dwx, dtx, csig,                
     &                       kin_energy, rho, point_q, point_accel_x,           
     &                       point_accel_y, point_accel_z, point_temp,          
     &                       process_temperatures, elem_alpha,                  
     &                       dalpha_x1, csig_dstrain_x1, dstrain_x1,            
     &                       dswd_x1, omit_ct_elem, fgm_e, fgm_nu,              
     &                       elemno, myid, numprocs, seg_curves_flag,           
     &                       debug, iout)                                       
        end if                                                                  
c                                                                               
        if( comput_i ) then                                                     
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
     &                        front_coords, domain_type, domain_origin,         
     &                        front_order, point_x, point_y, point_z,           
     &                        elemno, ptno, r, t, crack_curvature,              
     &                        debug, iout )                                     
c                                                                               
c             compute constitutive tensor components                            
c                                                                               
           call di_calc_constitutive( dcijkl_x1, sijkl, dsijkl_x1,              
     &                             point_ym, point_nu, de_x1, dnu_x1,           
     &                             elemno, debug, iout )                        
c                                                                               
c             compute auxiliary fields for stress intensity factors             
c                                                                               
           call di_calc_aux_fields_k( elemno, ptno, r, t, ym_front_node,        
     &                             nu_front_node, dcijkl_x1, sijkl,             
     &                             dsijkl_x1, aux_stress,                       
     &                             aux_strain, daux_strain_x1,                  
     &                             du11_aux,  du21_aux,  du31_aux,              
     &                             du111_aux, du112_aux, du113_aux,             
     &                             du211_aux, du212_aux, du213_aux,             
     &                             du311_aux, du312_aux, du313_aux,             
     &                             iout )                                       
c                                                                               
c             compute auxiliary fields for t-stresses                           
c                                                                               
           call di_calc_aux_fields_t( elemno, ptno, r, t, ym_front_node,        
     &                             nu_front_node, dcijkl_x1, sijkl,             
     &                             dsijkl_x1, aux_stress,                       
     &                             aux_strain, daux_strain_x1,                  
     &                             du11_aux,  du21_aux,  du31_aux,              
     &                             du111_aux, du112_aux, du113_aux,             
     &                             du211_aux, du212_aux, du213_aux,             
     &                             du311_aux, du312_aux, du313_aux,             
     &                             iout )                                       
c                                                                               
c             compute interaction integral terms for stress intensity factors   
c                                                                               
           call di_calc_i_terms( ptno, dqx, dqy, dqz, dux, dvx, dwx,            
     &                        dtx, csig, aux_stress, ceps_gp,                   
     &                        aux_strain, dstrain_x1,                           
     &                        daux_strain_x1, dcijkl_x1,                        
     &                        du11_aux,  du21_aux,  du31_aux,                   
     &                        du111_aux, du211_aux, du311_aux,                  
     &                        du112_aux, du212_aux, du312_aux,                  
     &                        du113_aux, du213_aux, du313_aux,                  
     &                        process_temperatures, elem_alpha,                 
     &                        dalpha_x1, point_temp, point_q, weight,           
     &                        elemno, fgm_e, fgm_nu, iterm,                     
     &                        iout, debug)                                      
c                                                                               
c                                                                               
        end if                                                                  
c                                                                               
c     -----------------------------------------------------------------         
c     ---------------end of loop over integration points---------------         
c     -----------------------------------------------------------------         
c                                                                               
      end do                                                                    
c                                                                               
c             perform edi evaluations for traction loaded faces                 
c             for jterm(5) and then for iterm(8).                               
c                                                                               
      call di_calc_surface_integrals( elemno, etype, nnode, snodes,             
     &                             feload, cf_traction_flags,                   
     &                             cf_tractions, rotate, dsf, jacobi,           
     &                             cdispl, qvals, coord, front_nodes,           
     &                             front_coords, domain_type,                   
     &                             domain_origin, num_front_nodes,              
     &                             front_order, ym_front_node,                  
     &                             nu_front_node, comput_j,                     
     &                             comput_i, jterm, iterm,                      
     &                             front_elem_flag, qp_node,                    
     &                             crack_curvature, face_loading, iout,         
     &                             debug)                                       
c                                                                               
c            save edi values in result vectors                                  
c                                                                               
      e_jresults(1) = jterm(1)                                                  
      e_jresults(2) = jterm(2)                                                  
      e_jresults(3) = jterm(3)                                                  
      e_jresults(4) = jterm(4)                                                  
      e_jresults(5) = jterm(5)                                                  
      e_jresults(6) = jterm(6)                                                  
      e_jresults(7) = jterm(7)                                                  
      e_jresults(8) = jterm(8)                                                  
c                                                                               
      do j=1,8                                                                  
         e_iresults(1,j) = iterm(1,j)                                           
         e_iresults(2,j) = iterm(2,j)                                           
         e_iresults(3,j) = iterm(3,j)                                           
         e_iresults(4,j) = iterm(4,j)                                           
         e_iresults(5,j) = iterm(5,j)                                           
         e_iresults(6,j) = iterm(6,j)                                           
         e_iresults(7,j) = iterm(7,j)                                           
         e_iresults(8,j) = iterm(8,j)                                           
      end do                                                                    
c                                                                               
      if ( debug_j ) write (iout,9003) elemno,   jterm(1), jterm(2),            
     &                                 jterm(3), jterm(4), jterm(5),            
     &                                 jterm(6), jterm(7), jterm(8),            
     &                                 evol                                     
      if ( debug_i ) then                                                       
         write(iout,9009) elemno                                                
         do j=1,8                                                               
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
 9001 format(' >>>>> supplied values of s (q) function at element',             
     &       ' nodes', /,10(/,5x,3f10.3) )                                      
 9002 format(//,                                                                
     &       ' >>>>> domain computations for element: ',i7,                     
     & //,   '       e, nu, etype, order, ngpts, nnode, linear,: ',             
     & /,15x, f10.1, f5.2, i5, i5, i5, i5, l5,                                  
     & /,    ' >>>>> global->crack rotation matrix:',                           
     & 3(/,15x,3f10.6) )                                                        
 9003 format(//,' >>>>> final j-integral contributions for element ',           
     & i6,':',                                                                  
     & /,15x,'jterm1: ',e13.6,'   jterm2: ',e13.6,                              
     & /,15x,'jterm3: ',e13.6,'   jterm4: ',e13.6,                              
     & /,15x,'jterm5: ',e13.6,'   jterm6: ',e13.6,                              
     & /,15x,'jterm7: ',e13.6,'   jterm8: ',e13.6,                              
     & /,15x,'element volume: ',e13.6 )                                         
 9004 format(' >>>>> nodal temperatures:',                                      
     &  /,15x,6f10.3,/,15x,6f10.3 )                                             
 9005 format(' >>>>> process temperature loading: ',l1 )                        
 9006 format(' >>>>> node alpha(1->6) for element in global x-y-z:' )           
 9007 format(3x,i3,2x,i3,2x,6e14.6)                                             
 9008 format(' >>>>> node alpha(1->6) for nodes in crack x-y-z:' )              
 9009 format(//,' >>>>> final i-integral contributions for element ',i7)        
 9010 format(/,a)                                                               
 9011 format('iterm1: ',e13.6,' iterm2: ',e13.6,' iterm3: ',e13.6,              
     &       ' iterm4: ',e13.6,' iterm5: ',e13.6,' iterm6: ',e13.6,             
     &       ' iterm7: ',e13.6,' iterm8: ',e13.6 )                              
 9012 format(//,15x,'element volume: ',e13.6 )                                  
 9013 format(' >>>>> start of data set-up for edi' )                            
 9100 format(/,' >>>>> start of integration loop:' )                            
 9110 format(/,'     >> ptno, xsi, eta, zeta, weight: ',i2,4(2x,f10.4))         
 9120 format(/,'     >> q-value, velocity, k. ener: ',3e14.6,                   
     &       /,'        accelerations:              ',3e14.6,                   
     &       /,'     >> int. pt. temperature:       ',e14.6)                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *    error messages for 3-d edi evaluation                        *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
       subroutine dieler( iout,ierr,erno )                                      
       implicit integer (a-z)                                                   
c                                                                               
c                                                                               
c          write an error message                                               
c                                                                               
c          return error flag                                                    
c                      ierr = 0 non-fatal                                       
c                      ierr = 2 terminate element computations                  
c                                                                               
c                                                                               
       go to (10,20,30,40,50,60,70,80),  erno                                   
c                                                                               
   10 write(iout,1001)                                                          
      ierr = 2                                                                  
      return                                                                    
c                                                                               
   20 continue                                                                  
      ierr = 2                                                                  
      return                                                                    
c                                                                               
   30 continue                                                                  
      ierr = 2                                                                  
      return                                                                    
c                                                                               
   40 continue                                                                  
      ierr = 2                                                                  
      return                                                                    
c                                                                               
   50 write (iout,1005)                                                         
      ierr = 2                                                                  
      return                                                                    
c                                                                               
   60 continue                                                                  
      ierr = 0                                                                  
      return                                                                    
c                                                                               
   70 write(iout,1007)                                                          
      ierr = 0                                                                  
      return                                                                    
c                                                                               
   80 continue                                                                  
      ierr = 0                                                                  
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
c *   gauss point.                                                  *           
c *                                                                 *           
c *                    last modified: 9/5/01    by: mcw             *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine dielav( sf, evel, eaccel, velocity, x_accel, y_accel,          
     &                   z_accel, eqvalues, q, nnode, point_temp,               
     &                   enode_temps, elem_alpha, enode_alpha_ij,               
     &                   linear_displ, fgm_alphas, elemno, ym_nodes,            
     &                   nu_nodes, point_ym, point_nu, fgm_e, fgm_nu,           
     &                   coord, point_x, point_y, point_z,                      
     &                   seg_curves_flag, out )                                 
      implicit none                                                             
c                                                                               
c          dummy declarations                                                   
c                                                                               
      double precision                                                          
     &         sf(*), evel(3,*), eaccel(3,*), velocity, x_accel,                
     &         y_accel, z_accel, eqvalues(*), q, point_temp,                    
     &         enode_temps(*), elem_alpha(*), enode_alpha_ij(6,*),              
     &         ym_nodes(*), nu_nodes(*), point_ym, point_nu,                    
     &         coord(3,*), point_x, point_y, point_z                            
      integer nnode, elemno, mxndel, out                                        
      logical linear_displ, fgm_e, fgm_nu, seg_curves_flag, fgm_alphas          
c                                                                               
c          local declarations                                                   
c                                                                               
      double precision                                                          
     &         x_vel, y_vel, z_vel, gp_temp, avg_temp, gp_alpha,                
     &         avg_alpha, zero, compl_tens(3,2), avg_ym, avg_nu                 
      integer enode, i, j, k                                                    
      logical debug                                                             
      data zero / 0.0d0 /                                                       
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
      if( linear_displ ) then                                                   
        point_temp = avg_temp / dble( nnode )                                   
      else                                                                      
            point_temp = gp_temp                                                
      end if                                                                    
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
      if( fgm_e .and. linear_displ ) then                                       
        point_ym = avg_ym / dble( nnode )                                       
      end if                                                                    
c                                                                               
      if( fgm_nu .and. linear_displ ) then                                      
        point_nu = avg_nu / dble( nnode )                                       
      end if                                                                    
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
c *******************************************************************           
c                                                                               
      subroutine dielrv( coord, cdispl, cvel, caccel, edispl, evel,             
     &                   eaccel, rotate, nnode, iout, debug )                   
      implicit double precision (a-h,o-z)                                       
c                                                                               
c                                                                               
c             rotate coordinates, displacements, velocities and                 
c             accelerations of element nodes from the global                    
c             system to the crack reference frame.                              
c                                                                               
c                                                                               
      integer i                                                                 
      dimension coord(3,*), cdispl(3,*), cvel(3,*), caccel(3,*),                
     &          edispl(3,*), evel(3,*), eaccel(3,*), rotate(3,3)                
      logical debug                                                             
c                                                                               
      if ( debug ) write(iout,9001) (inode,(coord(i,inode),i=1,3),              
     &                               inode=1,nnode)                             
c                                                                               
      r11 = rotate(1,1)                                                         
      r21 = rotate(2,1)                                                         
      r31 = rotate(3,1)                                                         
      r12 = rotate(1,2)                                                         
      r22 = rotate(2,2)                                                         
      r32 = rotate(3,2)                                                         
      r13 = rotate(1,3)                                                         
      r23 = rotate(2,3)                                                         
      r33 = rotate(3,3)                                                         
c                                                                               
      do inode = 1, nnode                                                       
        x =   r11 * coord(1,inode) + r12 * coord(2,inode)                       
     &      + r13 * coord(3,inode)                                              
        y =   r21 * coord(1,inode) + r22 * coord(2,inode)                       
     &      + r23 * coord(3,inode)                                              
        z =   r31 * coord(1,inode) + r32 * coord(2,inode)                       
     &      + r33 * coord(3,inode)                                              
        coord(1,inode) = x                                                      
        coord(2,inode) = y                                                      
        coord(3,inode) = z                                                      
        cdispl(1,inode) = r11*edispl(1,inode) + r12*edispl(2,inode) +           
     &                    r13*edispl(3,inode)                                   
        cdispl(2,inode) = r21*edispl(1,inode) + r22*edispl(2,inode) +           
     &                    r23*edispl(3,inode)                                   
        cdispl(3,inode) = r31*edispl(1,inode) + r32*edispl(2,inode) +           
     &                    r33*edispl(3,inode)                                   
        cvel(1,inode)   = r11*evel(1,inode)   + r12*evel(2,inode) +             
     &                    r13*evel(3,inode)                                     
        cvel(2,inode)   = r21*evel(1,inode)   + r22*evel(2,inode) +             
     &                    r23*evel(3,inode)                                     
        cvel(3,inode)   = r31*evel(1,inode)   + r32*evel(2,inode) +             
     &                    r33*evel(3,inode)                                     
        caccel(1,inode) = r11*eaccel(1,inode) + r12*eaccel(2,inode) +           
     &                    r13*eaccel(3,inode)                                   
        caccel(2,inode) = r21*eaccel(1,inode) + r22*eaccel(2,inode) +           
     &                    r23*eaccel(3,inode)                                   
        caccel(3,inode) = r31*eaccel(1,inode) + r32*eaccel(2,inode) +           
     &                    r33*eaccel(3,inode)                                   
      end do                                                                    
c                                                                               
      if ( debug ) then                                                         
         write(iout,9002) (inode,(coord(i,inode),i=1,3),inode=1,nnode)          
         write(iout,9003) (inode,(edispl(i,inode),                              
     &                     i=1,3),(cdispl(j,inode),j=1,3),inode=1,nnode)        
         write(iout,9004) (inode,(evel(i,inode),                                
     &                     i=1,3),(cvel(j,inode),j=1,3),inode=1,nnode)          
         write(iout,9005) (inode,(eaccel(i,inode),                              
     &                     i=1,3),(caccel(j,inode),j=1,3),inode=1,nnode)        
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
c                                                                               
 9001 format (' >>>>>  element global nodal coordinates',//,                    
     & 5x,'node',15x,'x',20x,'y',20x,'z',27(/,i5,3(7x,e13.6)))                  
 9002 format (' >>>>>  nodal coordinates in crack reference frame',//,          
     & 5x,'node',15x,'x',19x,'y',19x,'z',27(/,i5,3(7x,e13.6)))                  
 9003 format (' >>>>> nodal displacements',//,                                  
     & 3x,'node',10x,'global u',13x,'global v',13x,'global w',13x,              
     & 'crack u',13x,'crack v',13x,'crack w',/,32(/,i5,6(7x,e13.6)))            
 9004 format (' >>>>> nodal velocities',//,                                     
     & 3x,'node',10x,'global u',13x,'global v',13x,'global w',13x,              
     & 'crack u',13x,'crack v',13x,'crack w',/,32(/,i5,6(7x,e13.6)))            
 9005 format (' >>>>> nodal accelerations',//,                                  
     & 3x,'node',10x,'global u',13x,'global v',13x,'global w',13x,              
     & 'crack u',13x,'crack v',13x,'crack w',/,32(/,i5,6(7x,e13.6)))            
c                                                                               
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *  rotate stresses or strains from global->crack x-y-z            *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine dielrt ( ptno, rotate, global, crack, id, iout,                
     &                    debug )                                               
      implicit double precision (a-h,o-z)                                       
      integer ptno, id                                                          
      dimension rotate(3,3), global(3,3), crack(3,3), space(3,3)                
      logical   debug                                                           
c                                                                               
c            rotation is performed in tensor notation                           
c                                                                               
c            crack = rotate * global * rotate(transpose)                        
c                                                                               
c            for stresses:                                                      
c            global are element global stresses - 1st PK for                    
c            geonl nonlinear formulation                                        
c                                                                               
c            for strains:                                                       
c            global are element global strains                                  
c                                                                               
c            space = rotate * global                                            
c                                                                               
      space(1,1) = rotate(1,1)*global(1,1) + rotate(1,2)*global(2,1) +          
     &             rotate(1,3)*global(3,1)                                      
      space(1,2) = rotate(1,1)*global(1,2) + rotate(1,2)*global(2,2) +          
     &             rotate(1,3)*global(3,2)                                      
      space(1,3) = rotate(1,1)*global(1,3) + rotate(1,2)*global(2,3) +          
     &             rotate(1,3)*global(3,3)                                      
      space(2,1) = rotate(2,1)*global(1,1) + rotate(2,2)*global(2,1) +          
     &             rotate(2,3)*global(3,1)                                      
      space(2,2) = rotate(2,1)*global(1,2) + rotate(2,2)*global(2,2) +          
     &             rotate(2,3)*global(3,2)                                      
      space(2,3) = rotate(2,1)*global(1,3) + rotate(2,2)*global(2,3) +          
     &             rotate(2,3)*global(3,3)                                      
      space(3,1) = rotate(3,1)*global(1,1) + rotate(3,2)*global(2,1) +          
     &             rotate(3,3)*global(3,1)                                      
      space(3,2) = rotate(3,1)*global(1,2) + rotate(3,2)*global(2,2) +          
     &             rotate(3,3)*global(3,2)                                      
      space(3,3) = rotate(3,1)*global(1,3) + rotate(3,2)*global(2,3) +          
     &             rotate(3,3)*global(3,3)                                      
c                                                                               
c            crack = space * trans(rotate)                                      
c                                                                               
      crack(1,1) = space(1,1)*rotate(1,1) + space(1,2)*rotate(1,2) +            
     &             space(1,3)*rotate(1,3)                                       
      crack(1,2) = space(1,1)*rotate(2,1) + space(1,2)*rotate(2,2) +            
     &             space(1,3)*rotate(2,3)                                       
      crack(1,3) = space(1,1)*rotate(3,1) + space(1,2)*rotate(3,2) +            
     &             space(1,3)*rotate(3,3)                                       
      crack(2,1) = space(2,1)*rotate(1,1) + space(2,2)*rotate(1,2) +            
     &             space(2,3)*rotate(1,3)                                       
      crack(2,2) = space(2,1)*rotate(2,1) + space(2,2)*rotate(2,2) +            
     &             space(2,3)*rotate(2,3)                                       
      crack(2,3) = space(2,1)*rotate(3,1) + space(2,2)*rotate(3,2) +            
     &             space(2,3)*rotate(3,3)                                       
      crack(3,1) = space(3,1)*rotate(1,1) + space(3,2)*rotate(1,2) +            
     &             space(3,3)*rotate(1,3)                                       
      crack(3,2) = space(3,1)*rotate(2,1) + space(3,2)*rotate(2,2) +            
     &             space(3,3)*rotate(2,3)                                       
      crack(3,3) = space(3,1)*rotate(3,1) + space(3,2)*rotate(3,2) +            
     &             space(3,3)*rotate(3,3)                                       
c                                                                               
      if ( .not. debug ) return                                                 
      if( id .eq. 1 ) write ( iout,9001 ) ptno                                  
      if( id .eq. 2 ) write ( iout,9002 ) ptno                                  
      if( id .eq. 3 ) write ( iout,9003 ) ptno                                  
      write ( iout,9100 ) global(1,1), global(1,2), global(1,3),                
     &                    global(2,1), global(2,2), global(2,3),                
     &                    global(3,1), global(3,2), global(3,3)                 
      if( id .eq. 1 ) write ( iout,9004 ) ptno                                  
      if( id .eq. 2 ) write ( iout,9005 ) ptno                                  
      if( id .eq. 3 ) write ( iout,9006 ) ptno                                  
      write ( iout,9100 ) crack(1,1), crack(1,2), crack(1,3),                   
     &                    crack(2,1), crack(2,2), crack(2,3),                   
     &                    crack(3,1), crack(3,2), crack(3,3)                    
c                                                                               
      return                                                                    
c                                                                               
 9001 format (' >>>>> global element stresses. point:',i3)                      
 9002 format (' >>>>> global element strains. point:',i3)                       
 9003 format (' >>>>> global element strains. node:',i3)                        
 9004 format (' >>>>> crack reference frame stresses. point:',i3)               
 9005 format (' >>>>> crack reference frame strains. point:',i3)                
 9006 format (' >>>>> crack reference frame strains. node:',i3)                 
 9100 format ('  point number:',3(/10x,3e15.6) )                                
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *    compute coordinate jacobian, its determinate, and inverse    *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine dielcj( dsf, coord, nnode, cjacob, cjinv, det, ierr,           
     &                   iout, debug )                                          
      implicit double precision (a-h,o-z)                                       
c                                                                               
c                                                                               
c              compute the 3 x 3 jacobian, its determinate and                  
c              inverse for the 3-d isoparametrics.                              
c                                                                               
c                                                                               
      dimension  dsf(20,*), coord(3,1), cjacob(3,3), cjinv(3,3)                 
      logical  debug                                                            
      data zero, one / 0.0d0, 1.0d0 /                                           
c                                                                               
c              compute jacobian at the point. use a dot product                 
c              support function.                                                
c                                                                               
      cjacob(1,1) = dieldp( dsf(1,1), coord(1,1), nnode, 1, 3 )                 
      cjacob(1,2) = dieldp( dsf(1,1), coord(2,1), nnode, 1, 3 )                 
      cjacob(1,3) = dieldp( dsf(1,1), coord(3,1), nnode, 1, 3 )                 
      cjacob(2,1) = dieldp( dsf(1,2), coord(1,1), nnode, 1, 3 )                 
      cjacob(2,2) = dieldp( dsf(1,2), coord(2,1), nnode, 1, 3 )                 
      cjacob(2,3) = dieldp( dsf(1,2), coord(3,1), nnode, 1, 3 )                 
      cjacob(3,1) = dieldp( dsf(1,3), coord(1,1), nnode, 1, 3 )                 
      cjacob(3,2) = dieldp( dsf(1,3), coord(2,1), nnode, 1, 3 )                 
      cjacob(3,3) = dieldp( dsf(1,3), coord(3,1), nnode, 1, 3 )                 
      if ( debug ) write(iout,9001) (( cjacob(j,i), i = 1,3 ), j=1,3)           
c                                                                               
      det = cjacob(1,1) * cjacob(2,2) * cjacob(3,3)                             
     &    + cjacob(2,1) * cjacob(3,2) * cjacob(1,3)                             
     &    + cjacob(3,1) * cjacob(1,2) * cjacob(2,3)                             
     &    - cjacob(1,1) * cjacob(3,2) * cjacob(2,3)                             
     &    - cjacob(2,1) * cjacob(1,2) * cjacob(3,3)                             
     &    - cjacob(3,1) * cjacob(2,2) * cjacob(1,3)                             
      if ( det .le. zero ) then                                                 
         ierr = 1                                                               
         return                                                                 
      end if                                                                    
c                                                                               
      deti = one / det                                                          
      cjinv(1,1) =  (cjacob(2,2) * cjacob(3,3)                                  
     &              - cjacob(3,2) * cjacob(2,3)) * deti                         
      cjinv(2,2) =  (cjacob(1,1) * cjacob(3,3)                                  
     &              - cjacob(3,1) * cjacob(1,3)) * deti                         
      cjinv(3,3) =  (cjacob(1,1) * cjacob(2,2)                                  
     &              - cjacob(2,1) * cjacob(1,2)) * deti                         
      cjinv(2,1) = -(cjacob(2,1) * cjacob(3,3)                                  
     &              - cjacob(3,1) * cjacob(2,3)) * deti                         
      cjinv(3,1) =  (cjacob(2,1) * cjacob(3,2)                                  
     &              - cjacob(3,1) * cjacob(2,2)) * deti                         
      cjinv(1,2) = -(cjacob(1,2) * cjacob(3,3)                                  
     &              - cjacob(3,2) * cjacob(1,3)) * deti                         
      cjinv(3,2) = -(cjacob(1,1) * cjacob(3,2)                                  
     &             -  cjacob(3,1) * cjacob(1,2)) * deti                         
      cjinv(1,3) =  (cjacob(1,2) * cjacob(2,3)                                  
     &              - cjacob(2,2) * cjacob(1,3)) * deti                         
      cjinv(2,3) = -(cjacob(1,1) * cjacob(2,3)                                  
     &              - cjacob(2,1) * cjacob(1,3)) * deti                         
      if ( debug ) write(iout,9002) det,((cjinv(j,i),i=1,3 ), j=1,3)            
c                                                                               
      return                                                                    
c                                                                               
 9001 format(/,15x,' jacobian at point', /,                                     
     &                 3(/,7x,3f15.5) )                                         
 9002 format(15x,' determinant',f15.9, /,                                       
     &       15x,' jacobian inverse', /,3(/,7x,3f15.5) )                        
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
c *                                 last modified: 11/20/03        *            
c *                                            by: mcw             *            
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
      double precision                                                          
     &         qvals(*), coord(3,*), front_coords(3,*),                         
     &         crack_curvature(*)                                               
      logical debug, front_elem_flag, qp_node                                   
      integer out, nnode, etype, elemno, snodes(*), num_front_nodes,            
     &        max_exp_front,                                                    
     &        front_nodes(*), expanded_front_nodes(0:max_exp_front,*),          
     &        domain_type, domain_origin, front_order                           
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, j, k, snode, num_qp_nodes                                      
      double precision                                                          
     &         r_max, r, t, rs(nnode), zero, toler, half, p75, node_x,          
     &         node_y, node_z, p4                                               
      data zero, toler, half, p4, p75 / 0.0d0, 1.0d-5, 0.5d0,                   
     &                                  0.4d0, 0.75d0 /                         
c                                                                               
      qp_node      = .false.                                                    
      num_qp_nodes = 0                                                          
c                                                                               
c                set q-value of interior nodes by                               
c                interpolating linearly from corresponding                      
c                corner nodes for each edge of the element.                     
c                                                                               
      go to ( 100, 200, 300, 400, 500 ), etype                                  
c                                                                               
c             20 node elements                                                  
c                                                                               
 100  continue                                                                  
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
      goto 600                                                                  
c                                                                               
c             8-nodes                                                           
c                                                                               
 200  continue                                                                  
      return                                                                    
c                                                                               
c             12 node elements                                                  
c                                                                               
 300  continue                                                                  
      qvals(9)  = ( qvals(1) + qvals(2) ) * half                                
      qvals(10) = ( qvals(2) + qvals(3) ) * half                                
      qvals(11) = ( qvals(3) + qvals(4) ) * half                                
      qvals(12) = ( qvals(4) + qvals(1) ) * half                                
      goto 600                                                                  
c                                                                               
c             15 node elements                                                  
c                                                                               
 400  continue                                                                  
      qvals(9)  = ( qvals(1) + qvals(2) ) * half                                
      qvals(10) = ( qvals(2) + qvals(3) ) * half                                
      qvals(11) = ( qvals(3) + qvals(4) ) * half                                
      qvals(12) = ( qvals(4) + qvals(1) ) * half                                
      qvals(13) = ( qvals(5) + qvals(6) ) * half                                
      qvals(14) = ( qvals(1) + qvals(5) ) * half                                
      qvals(15) = ( qvals(2) + qvals(6) ) * half                                
      goto 600                                                                  
c                                                                               
c             9 node elements                                                   
c                                                                               
 500  continue                                                                  
      qvals(9)  = ( qvals(1) + qvals(2) ) * half                                
c                                                                               
 600  continue                                                                  
c                                                                               
c             set the q-value to 0.75 for 1/4-point nodes.                      
c             if any 1/4-point node is detected, special integration            
c             for crack-face traction on crack-front element faces              
c             is not used.                                                      
c             skip process when element not on crack front.                     
c                                                                               
      if( .not. front_elem_flag ) goto 700                                      
c                                                                               
c             compute distance from each element node to crack front.           
c             for curved elements, di_calc_r_theta currently gives an           
c             approximation: the distance is measured to a line                 
c             connecting adjacent crack front nodes.                            
c                                                                               
      do i=1,nnode                                                              
c                                                                               
         node_x = coord(1,i)                                                    
         node_y = coord(2,i)                                                    
         node_z = coord(3,i)                                                    
c                                                                               
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
      do i=1,nnode                                                              
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
        do j = 1, num_front_nodes                                               
            do k = 1, max_exp_front                                             
               if( expanded_front_nodes(k,j) .eq. snode ) goto 650              
            end do                                                              
         end do                                                                 
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
 650     continue                                                               
      end do                                                                    
c                                                                               
 700  continue                                                                  
c                                                                               
      if( debug ) then                                                          
         if( front_elem_flag ) write(out,9012) r_max, num_qp_nodes              
         write(out,9014) elemno                                                 
         do i=1,nnode                                                           
            write(out,9016) i, snodes(i), qvals(i)                              
         end do                                                                 
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9012 format(/,' >>>>> max. elem. node dist. from front:     ',e13.6,           
     &       /,' >>>>> number of 1/4-point nodes on element: ',i7)              
 9014 format(' >>>>> final nodal q-function values for element',                
     &       2x,i7,/,5x,'node',5x,'snode',8x,'q')                               
 9016 format(5x,i2,5x,i7,2x,f10.3)                                              
c                                                                               
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *   build 6x6 rotation matrix to convert unrotated cauchy         *           
c *   stresses to cauchy stresses                                   *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine digetr( q, r  )                                                
      implicit integer (a-z)                                                    
c                                                                               
c           parameter declarations                                              
c                                                                               
      double precision                                                          
     &  q(6,6), r(3,3)                                                          
c                                                                               
c           locals                                                              
c                                                                               
      double precision                                                          
     & two                                                                      
      data two / 2.0d0 /                                                        
c                                                                               
c       cauchy stress {T} = [q] * (rotated) cauchy stress {t}.                  
c       in tensor form:                                                         
c                                                                               
c                 [T] = [r] [t] trans([r])                                      
c                                                                               
c       both [T] and [t] are symmetric, [r] is orthogonal rotation.             
c       vector ordering is {x,y,z,xy,yz,xz}.                                    
c                                                                               
      q(1,1)= r(1,1)**2                                                         
      q(1,2)= r(1,2)**2                                                         
      q(1,3)= r(1,3)**2                                                         
      q(1,4)= two*r(1,1)*r(1,2)                                                 
      q(1,5)= two*r(1,3)*r(1,2)                                                 
      q(1,6)= two*r(1,1)*r(1,3)                                                 
      q(2,1)= r(2,1)**2                                                         
      q(2,2)= r(2,2)**2                                                         
      q(2,3)= r(2,3)**2                                                         
      q(2,4)= two*r(2,1)*r(2,2)                                                 
      q(2,5)= two*r(2,3)*r(2,2)                                                 
      q(2,6)= two*r(2,1)*r(2,3)                                                 
      q(3,1)= r(3,1)**2                                                         
      q(3,2)= r(3,2)**2                                                         
      q(3,3)= r(3,3)**2                                                         
      q(3,4)= two*r(3,1)*r(3,2)                                                 
      q(3,5)= two*r(3,3)*r(3,2)                                                 
      q(3,6)= two*r(3,1)*r(3,3)                                                 
      q(4,1)= r(1,1)*r(2,1)                                                     
      q(4,2)= r(1,2)*r(2,2)                                                     
      q(4,3)= r(1,3)*r(2,3)                                                     
      q(4,4)= r(1,1)*r(2,2)+r(2,1)*r(1,2)                                       
      q(4,5)= r(1,2)*r(2,3)+r(1,3)*r(2,2)                                       
      q(4,6)= r(1,1)*r(2,3)+r(1,3)*r(2,1)                                       
      q(5,1)= r(2,1)*r(3,1)                                                     
      q(5,2)= r(3,2)*r(2,2)                                                     
      q(5,3)= r(2,3)*r(3,3)                                                     
      q(5,4)= r(2,1)*r(3,2)+r(2,2)*r(3,1)                                       
      q(5,5)= r(2,2)*r(3,3)+r(3,2)*r(2,3)                                       
      q(5,6)= r(2,1)*r(3,3)+r(2,3)*r(3,1)                                       
      q(6,1)= r(1,1)*r(3,1)                                                     
      q(6,2)= r(1,2)*r(3,2)                                                     
      q(6,3)= r(1,3)*r(3,3)                                                     
      q(6,4)= r(1,1)*r(3,2)+r(1,2)*r(3,1)                                       
      q(6,5)= r(1,2)*r(3,3)+r(1,3)*r(3,2)                                       
      q(6,6)= r(1,1)*r(3,3)+r(3,1)*r(1,3)                                       
c                                                                               
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *   rotate unrotated cauchy stresses to cauchy stresses           *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine diqmp1( q, stress, causig )                                    
      implicit integer (a-z)                                                    
c                                                                               
c                    parameter declarations                                     
c                                                                               
      double precision                                                          
     & q(6,6), stress(*), cstress(6), causig(3,3)                               
c                                                                               
c            {cstress} = [q] * {stress}; 6x1 vectors and 6x6 q                  
c             (Cauchy)         (u.r. Cauchy)                                    
c                                                                               
      cstress(1) = q(1,1)*stress(1)+q(1,2)*stress(2)+                           
     &             q(1,3)*stress(3)+q(1,4)*stress(4)+                           
     &             q(1,5)*stress(5)+q(1,6)*stress(6)                            
      cstress(2) = q(2,1)*stress(1)+q(2,2)*stress(2)+                           
     &             q(2,3)*stress(3)+q(2,4)*stress(4)+                           
     &             q(2,5)*stress(5)+q(2,6)*stress(6)                            
      cstress(3) = q(3,1)*stress(1)+q(3,2)*stress(2)+                           
     &             q(3,3)*stress(3)+q(3,4)*stress(4)+                           
     &             q(3,5)*stress(5)+q(3,6)*stress(6)                            
      cstress(4) = q(4,1)*stress(1)+q(4,2)*stress(2)+                           
     &             q(4,3)*stress(3)+q(4,4)*stress(4)+                           
     &             q(4,5)*stress(5)+q(4,6)*stress(6)                            
      cstress(5) = q(5,1)*stress(1)+q(5,2)*stress(2)+                           
     &             q(5,3)*stress(3)+q(5,4)*stress(4)+                           
     &             q(5,5)*stress(5)+q(5,6)*stress(6)                            
      cstress(6) = q(6,1)*stress(1)+q(6,2)*stress(2)+                           
     &             q(6,3)*stress(3)+q(6,4)*stress(4)+                           
     &             q(6,5)*stress(5)+q(6,6)*stress(6)                            
c                                                                               
c            store cauchy stresses in 3x3 tensor form.                          
c                                                                               
      causig(1,1) = cstress(1)                                                  
      causig(2,1) = cstress(4)                                                  
      causig(3,1) = cstress(6)                                                  
      causig(1,2) = cstress(4)                                                  
      causig(2,2) = cstress(2)                                                  
      causig(3,2) = cstress(5)                                                  
      causig(1,3) = cstress(6)                                                  
      causig(2,3) = cstress(5)                                                  
      causig(3,3) = cstress(3)                                                  
c                                                                               
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *   form deformation gradient, determinant, inverse transpose     *           
c *   at a gauss point. convert cauchy stress to 1st PK stress      *           
c *   scale energy density to volume density at t=0                 *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine digrad( nnode, dsf, jacobi, edispl, stress, iout,              
     &                   gdebug )                                               
      implicit integer (a-z)                                                    
c                                                                               
c                    parameter declarations                                     
c                                                                               
      double precision                                                          
     & dsf(20,3), jacobi(3,3), edispl(3,*), stress(*)                           
       logical gdebug, debug                                                    
c                                                                               
c                    local declarations                                         
c                                                                               
      double precision                                                          
     &  dux, dvx, dwx, duy, dvy, dwy, duz, dvz, dwz, zero, nx, ny, nz,          
     &  one, det, f(3,3), fi(3,3), pk1(3,3), deti                               
      data zero, one / 0.0d0, 1.0d0 /                                           
c                                                                               
c             for this integration point, compute displacement                  
c             derivatives in global coordinate system                           
c                                                                               
      debug = gdebug                                                            
      dux   = zero                                                              
      dvx   = zero                                                              
      dwx   = zero                                                              
      duy   = zero                                                              
      dvy   = zero                                                              
      dwy   = zero                                                              
      duz   = zero                                                              
      dvz   = zero                                                              
      dwz   = zero                                                              
      do enode = 1, nnode                                                       
         nx = dsf(enode,1) * jacobi(1,1) + dsf(enode,2) * jacobi(1,2) +         
     &        dsf(enode,3) * jacobi(1,3)                                        
         ny = dsf(enode,1) * jacobi(2,1) + dsf(enode,2) * jacobi(2,2) +         
     &        dsf(enode,3) * jacobi(2,3)                                        
         nz = dsf(enode,1) * jacobi(3,1) + dsf(enode,2) * jacobi(3,2) +         
     &        dsf(enode,3) * jacobi(3,3)                                        
         dux = dux + nx * edispl(1,enode)                                       
         dvx = dvx + nx * edispl(2,enode)                                       
         dwx = dwx + nx * edispl(3,enode)                                       
         duy = duy + ny * edispl(1,enode)                                       
         dvy = dvy + ny * edispl(2,enode)                                       
         dwy = dwy + ny * edispl(3,enode)                                       
         duz = duz + nz * edispl(1,enode)                                       
         dvz = dvz + nz * edispl(2,enode)                                       
         dwz = dwz + nz * edispl(3,enode)                                       
      end do                                                                    
c                                                                               
      f(1,1) = one + dux                                                        
      f(2,1) = dvx                                                              
      f(3,1) = dwx                                                              
      f(1,2) = duy                                                              
      f(2,2) = one + dvy                                                        
      f(3,2) = dwy                                                              
      f(1,3) = duz                                                              
      f(2,3) = dvz                                                              
      f(3,3) = one + dwz                                                        
c                                                                               
c             compute determinant and inverse of 3 x 3 matrix.                  
c             symmetry of input array not required.                             
c             use method of cofactors.                                          
c                                                                               
      det =   f(1,1) * f(2,2) * f(3,3)                                          
     &      + f(2,1) * f(3,2) * f(1,3)                                          
     &      + f(3,1) * f(1,2) * f(2,3)                                          
     &      - f(1,1) * f(3,2) * f(2,3)                                          
     &      - f(2,1) * f(1,2) * f(3,3)                                          
     &      - f(3,1) * f(2,2) * f(1,3)                                          
      if ( det .le. zero ) then                                                 
         write(*,*) '>> stop 1. det < 0 in digrad'                              
         call die_abort                                                         
         stop                                                                   
      end if                                                                    
c                                                                               
      deti    = one / det                                                       
      fi(1,1) =  (f(2,2) * f(3,3) - f(3,2) * f(2,3)) * deti                     
      fi(2,2) =  (f(1,1) * f(3,3) - f(3,1) * f(1,3)) * deti                     
      fi(3,3) =  (f(1,1) * f(2,2) - f(2,1) * f(1,2)) * deti                     
      fi(2,1) = -(f(2,1) * f(3,3) - f(3,1) * f(2,3)) * deti                     
      fi(3,1) =  (f(2,1) * f(3,2) - f(3,1) * f(2,2)) * deti                     
      fi(1,2) = -(f(1,2) * f(3,3) - f(3,2) * f(1,3)) * deti                     
      fi(3,2) = -(f(1,1) * f(3,2) - f(3,1) * f(1,2)) * deti                     
      fi(1,3) =  (f(1,2) * f(2,3) - f(2,2) * f(1,3)) * deti                     
      fi(2,3) = -(f(1,1) * f(2,3) - f(2,1) * f(1,3)) * deti                     
c                                                                               
      if ( debug ) then                                                         
        write(iout,*) '>> inside digrad: '                                      
        write(iout,*) '   > deformation gradient:'                              
        write(iout,9000) f(1,1), f(1,2), f(1,3),                                
     &                   f(2,1), f(2,2), f(2,3),                                
     &                   f(3,1), f(3,2), f(3,3)                                 
        write(iout,*) '   > inverse of deformation gradient:'                   
        write(iout,9000) fi(1,1), fi(1,2), fi(1,3),                             
     &                   fi(2,1), fi(2,2), fi(2,3),                             
     &                   fi(3,1), fi(3,2), fi(3,3)                              
      end if                                                                    
c                                                                               
c              1st PK = det * Cauchy * trans(fi)                                
c                                                                               
      pk1(1,1) = det * ( stress(1)*fi(1,1) + stress(4)*fi(1,2) +                
     &                   stress(7)*fi(1,3) )                                    
      pk1(1,2) = det * ( stress(1)*fi(2,1) + stress(4)*fi(2,2) +                
     &                   stress(7)*fi(2,3) )                                    
      pk1(1,3) = det * ( stress(1)*fi(3,1) + stress(4)*fi(3,2) +                
     &                   stress(7)*fi(3,3) )                                    
      pk1(2,1) = det * ( stress(2)*fi(1,1) + stress(5)*fi(1,2) +                
     &                   stress(8)*fi(1,3) )                                    
      pk1(2,2) = det * ( stress(2)*fi(2,1) + stress(5)*fi(2,2) +                
     &                   stress(8)*fi(2,3) )                                    
      pk1(2,3) = det * ( stress(2)*fi(3,1) + stress(5)*fi(3,2) +                
     &                   stress(8)*fi(3,3) )                                    
      pk1(3,1) = det * ( stress(3)*fi(1,1) + stress(6)*fi(1,2) +                
     &                   stress(9)*fi(1,3) )                                    
      pk1(3,2) = det * ( stress(3)*fi(2,1) + stress(6)*fi(2,2) +                
     &                   stress(9)*fi(2,3) )                                    
      pk1(3,3) = det * ( stress(3)*fi(3,1) + stress(6)*fi(3,2) +                
     &                   stress(9)*fi(3,3) )                                    
c                                                                               
      if ( debug ) then                                                         
        write(iout,*) '   > Cauchy stresses:'                                   
        write(iout,9000) stress(1), stress(4), stress(7),                       
     &                   stress(2), stress(5), stress(8),                       
     &                   stress(3), stress(6), stress(9)                        
        write(iout,*) '   > 1 st PK stresses:'                                  
        write(iout,9000) pk1(1,1), pk1(1,2), pk1(1,3),                          
     &                   pk1(2,1), pk1(2,2), pk1(2,3),                          
     &                   pk1(3,1), pk1(3,2), pk1(3,3)                           
      end if                                                                    
c                                                                               
      stress(1) = pk1(1,1)                                                      
      stress(2) = pk1(2,1)                                                      
      stress(3) = pk1(3,1)                                                      
      stress(4) = pk1(1,2)                                                      
      stress(5) = pk1(2,2)                                                      
      stress(6) = pk1(3,2)                                                      
      stress(7) = pk1(1,3)                                                      
      stress(8) = pk1(2,3)                                                      
      stress(9) = pk1(3,3)                                                      
c                                                                               
c              scale work density to volume at t=0                              
c                                                                               
      stress(10) = stress(10) * det                                             
c                                                                               
      if ( debug ) then                                                         
        write(iout,9100) det, stress(10), stress(10)/det                        
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(3(5x,3f20.6,/))                                                    
 9100 format(5x,'det: ',f20.6,/,5x,'w(0):',f20.6,                               
     &      /,5x,'w(n):',f20.6 )                                                
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *  set up to process temperature effects for domain integral      *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine dippie( rotate, iout, debug, elem_alpha )                      
      implicit double precision (a-h,o-z)                                       
c                                                                               
c                    parameter declarations                                     
c                                                                               
      dimension rotate(3,3),  elem_alpha(*)                                     
      logical   debug                                                           
c                                                                               
c                    local declarations                                         
c                                                                               
      dimension  work(3,3), tensor(3,3)                                         
      data zero, half / 0.0d0, 0.5d0 /                                          
c                                                                               
c                    the 6 alpha values for thermal expansion define            
c                    a symmetric 3x3 tensor to define initial                   
c                    strains due to a temperature change. rotate the            
c                    alphas to crack coordinates and put back into              
c                    the 6x1 elem_alpha vector. The tensor alphas               
c                    transform the same as strains.                             
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
      work(1,1) = rotate(1,1)*tensor(1,1) + rotate(1,2)*tensor(2,1) +           
     &            rotate(1,3)*tensor(3,1)                                       
      work(1,2) = rotate(1,1)*tensor(1,2) + rotate(1,2)*tensor(2,2) +           
     &            rotate(1,3)*tensor(3,2)                                       
      work(1,3) = rotate(1,1)*tensor(1,3) + rotate(1,2)*tensor(2,3) +           
     &            rotate(1,3)*tensor(3,3)                                       
      work(2,1) = rotate(2,1)*tensor(1,1) + rotate(2,2)*tensor(2,1) +           
     &            rotate(2,3)*tensor(3,1)                                       
      work(2,2) = rotate(2,1)*tensor(1,2) + rotate(2,2)*tensor(2,2) +           
     &            rotate(2,3)*tensor(3,2)                                       
      work(2,3) = rotate(2,1)*tensor(1,3) + rotate(2,2)*tensor(2,3) +           
     &            rotate(2,3)*tensor(3,3)                                       
      work(3,1) = rotate(3,1)*tensor(1,1) + rotate(3,2)*tensor(2,1) +           
     &            rotate(3,3)*tensor(3,1)                                       
      work(3,2) = rotate(3,1)*tensor(1,2) + rotate(3,2)*tensor(2,2) +           
     &            rotate(3,3)*tensor(3,2)                                       
      work(3,3) = rotate(3,1)*tensor(1,3) + rotate(3,2)*tensor(2,3) +           
     &            rotate(3,3)*tensor(3,3)                                       
c                                                                               
c            alpha_crack = work * trans(rotate)                                 
c                                                                               
      tensor(1,1) = work(1,1)*rotate(1,1) + work(1,2)*rotate(1,2) +             
     &              work(1,3)*rotate(1,3)                                       
      tensor(1,2) = work(1,1)*rotate(2,1) + work(1,2)*rotate(2,2) +             
     &              work(1,3)*rotate(2,3)                                       
      tensor(1,3) = work(1,1)*rotate(3,1) + work(1,2)*rotate(3,2) +             
     &              work(1,3)*rotate(3,3)                                       
      tensor(2,1) = work(2,1)*rotate(1,1) + work(2,2)*rotate(1,2) +             
     &              work(2,3)*rotate(1,3)                                       
      tensor(2,2) = work(2,1)*rotate(2,1) + work(2,2)*rotate(2,2) +             
     &              work(2,3)*rotate(2,3)                                       
      tensor(2,3) = work(2,1)*rotate(3,1) + work(2,2)*rotate(3,2) +             
     &              work(2,3)*rotate(3,3)                                       
      tensor(3,1) = work(3,1)*rotate(1,1) + work(3,2)*rotate(1,2) +             
     &              work(3,3)*rotate(1,3)                                       
      tensor(3,2) = work(3,1)*rotate(2,1) + work(3,2)*rotate(2,2) +             
     &              work(3,3)*rotate(2,3)                                       
      tensor(3,3) = work(3,1)*rotate(3,1) + work(3,2)*rotate(3,2) +             
     &              work(3,3)*rotate(3,3)                                       
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
c                                                                               
c   *******************************************************************         
c   *                                                                 *         
c   * support routine:  dot product: non-unit stride                  *         
c   *                                                                 *         
c   *******************************************************************         
c                                                                               
      double precision function dieldp( veca, vecb, nterms, stepa,              
     &                                  stepb )                                 
      implicit double precision (a-h,o-z)                                       
      dimension veca(1), vecb(1)                                                
      integer stepa, stepb                                                      
      data zero / 0.0d0 /                                                       
c                                                                               
      indexa = 1                                                                
      indexb = 1                                                                
      dieldp = zero                                                             
c                                                                               
      if ( stepa .eq. 1 .and. stepb .eq. 1 ) then                               
        do i = 1, nterms                                                        
          dieldp = dieldp + veca(i)*vecb(i)                                     
        end do                                                                  
      else                                                                      
        indexa = 1                                                              
        indexb = 1                                                              
        do i = 1, nterms                                                        
          dieldp = dieldp + veca(indexa)*vecb(indexb)                           
          indexa = indexa + stepa                                               
          indexb = indexb + stepb                                               
        end do                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c                                                                               
c***************************************************************                
c                                                              *                
c    subroutine to compute terms of volume integrals for J     *                
c                                                              *                
c                                 written by: rhd              *                
c                                   modified: mcw              *                
c                              last modified: 9/18/03          *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_calc_j_terms(weight, evol, jterm, ptno, dqx,                
     &                        dqy, dqz, dux, dvx, dwx, dtx, csig,               
     &                        kin_energy, rho, point_q, point_accel_x,          
     &                        point_accel_y, point_accel_z, point_temp,         
     &                        process_temperatures, elem_alpha,                 
     &                        dalpha_x1, csig_dstrain_x1, dstrain_x1,           
     &                        dswd_x1, omit_ct_elem, fgm_e, fgm_nu,             
     &                        elemno, myid, numprocs, seg_curves_flag,          
     &                        debug, out)                                       
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer ptno, elemno, myid, numprocs, out                                 
      double precision                                                          
     &     weight, evol, jterm(*), dqx, dqy, dqz, dux, dvx,                     
     &     dwx, dtx, csig(10,*), kin_energy, rho, point_q, dswd_x1,             
     &     point_accel_x, point_accel_y, point_accel_z, elem_alpha(*),          
     &     dalpha_x1(*), csig_dstrain_x1, dstrain_x1(*), point_temp             
      logical process_temperatures, omit_ct_elem, fgm_e, fgm_nu,                
     &        seg_curves_flag, debug                                            
c                                                                               
c             local variables                                                   
c                                                                               
      integer i                                                                 
      double precision                                                          
     &  temp1, temp2, temp3, zero, half                                         
c                                                                               
      data  zero, half                                                          
     &     / 0.0d0, 0.5d0 /                                                     
c                                                                               
c                                                                               
c      debug = .true.                                                           
c                                                                               
      if( debug ) then                                                          
         if( numprocs .gt. 1 ) write(out,*) "myid = ",myid                      
         write(out,890) elemno, ptno                                            
      end if                                                                    
c                                                                               
c             form jterm 1 of edi = total work density * derivative             
c             of q function w.r.t. crack x. also accumulate total               
c             element volume for debugging.                                     
c                                                                               
      evol     = evol + weight                                                  
      jterm(1) = jterm(1) - weight * dqx * csig(10,ptno)                        
      if( debug ) write(out,892) weight, dqx, csig(10,ptno), jterm(1)           
c                                                                               
c            form jterm 3 of edi = total kinetic energy density *               
c            derivative  of q function w.r.t. crack x.                          
c                                                                               
      jterm(3) = jterm(3) - weight * dqx * kin_energy                           
c                                                                               
c            form jterm 4 of edi = forces due to acceleration.                  
c                                                                               
      jterm(4) = jterm(4) + weight * rho * point_q * ( point_accel_x            
     &         * dux + point_accel_y * dvx + point_accel_z * dwx )              
c                                                                               
c            form jterm 2 of edi = stress * displacement derivatives            
c            in crack x * derivatives of q- function in crack x-y-z.            
c            use full tensor since we can have 1st PK stresses.                 
c                                                                               
      temp1 = csig(1,ptno)*dqx + csig(2,ptno)*dqy + csig(3,ptno)*dqz            
      temp2 = csig(4,ptno)*dqx + csig(5,ptno)*dqy + csig(6,ptno)*dqz            
      temp3 = csig(7,ptno)*dqx + csig(8,ptno)*dqy + csig(9,ptno)*dqz            
      jterm(2) = jterm(2) + weight*( dux*temp1 + dvx*temp2 + dwx*temp3 )        
      if ( debug )write(out,900) dux,dvx,dwx,temp1,temp2,temp3,jterm(2)         
c                                                                               
c            form jterm 6 of edi caused by temperature loads.                   
c            there are two parts both dotted into stresses:                     
c             (a) expansion coefficients * derivative of temperature            
c                 w.r.t. crack x                                                
c             (b) temperature * derivative of expansion coefficients            
c                 w.r.t. crack x                                                
c            stresses can be non-symmetric. thermal expansion coeffs.           
c            have been rotated to crack coordinates and are symmetric.          
c                                                                               
c            for (a) we use the specified alpha expansion coefficients for      
c            the element directly. for (b) we compute derivatives               
c            from node values of alpha for element (computed by averaging       
c            element values for elements incident on the node). off-            
c            diagonal alpha terms are multiplied by 'half' to obtain the        
c            tensorial component corresponding to tensorial strain. (see        
c            subroutine dippie.) note: jterm6 is not used when e and nu         
c            are entered using the 'fgm' option. this is because jterm6         
c            is implicitly included in jterm8 = -W,1. for the same reason,      
c            jterm6 is also not used when temperature-dependent curves          
c            are used to input material properties.                             
c                                                                               
      temp1 = zero                                                              
      temp2 = zero                                                              
      if ( process_temperatures ) then                                          
         temp1 = csig(1,ptno)*elem_alpha(1) +                                   
     &           csig(5,ptno)*elem_alpha(2) +                                   
     &           csig(9,ptno)*elem_alpha(3) +                                   
     &           elem_alpha(4) * half * ( csig(2,ptno)+csig(4,ptno) ) +         
     &           elem_alpha(5) * half * ( csig(6,ptno)+csig(8,ptno) ) +         
     &           elem_alpha(6) * half * ( csig(3,ptno)+csig(7,ptno) )           
         temp1 = temp1 * dtx * weight * point_q                                 
         temp2 = csig(1,ptno)*dalpha_x1(1) +                                    
     &           csig(5,ptno)*dalpha_x1(2) +                                    
     &           csig(9,ptno)*dalpha_x1(3) +                                    
     &           dalpha_x1(4) * half * ( csig(2,ptno)+csig(4,ptno) ) +          
     &           dalpha_x1(5) * half * ( csig(6,ptno)+csig(8,ptno) ) +          
     &           dalpha_x1(6) * half * ( csig(3,ptno)+csig(7,ptno) )            
         temp2 = temp2 * weight * point_q * point_temp                          
      end if                                                                    
      jterm(6) = jterm(6) + temp1 + temp2                                       
      if( fgm_e .or. fgm_nu ) jterm(6) = zero                                   
      if( seg_curves_flag )   jterm(6) = zero                                   
c                                                                               
      if( debug ) then                                                          
         if( numprocs .gt. 1 ) write(out,*) "myid = ",myid                      
         write(out,905) (csig(i,ptno),i=1,9)                                    
         write(out,906) (elem_alpha(i),i=1,6)                                   
         write(out,907) (dalpha_x1(i),i=1,6)                                    
         write(out,908) dtx, weight, point_q, point_temp, half                  
         write(out,910)  temp1, temp2, jterm(6)                                 
      end if                                                                    
c                                                                               
c              form jterm 7 & jterm 8 of J to account for gradients in          
c              temperature and thermal expansion coefficient when               
c              material is nonhomogeneous.                                      
c                                                                               
c                   jterm7 = sigma_ij * eps_ij,1                                
c                   jterm8 = -W,1                                               
c                   strains are in tensorial form.                              
c                                                                               
      csig_dstrain_x1 = csig(1,ptno) * dstrain_x1(1)                            
     &                + csig(2,ptno) * dstrain_x1(2)                            
     &                + csig(3,ptno) * dstrain_x1(3)                            
     &                + csig(4,ptno) * dstrain_x1(4)                            
     &                + csig(5,ptno) * dstrain_x1(5)                            
     &                + csig(6,ptno) * dstrain_x1(6)                            
     &                + csig(7,ptno) * dstrain_x1(7)                            
     &                + csig(8,ptno) * dstrain_x1(8)                            
     &                + csig(9,ptno) * dstrain_x1(9)                            
c                                                                               
      jterm(7) = jterm(7) + csig_dstrain_x1 * weight * point_q                  
      jterm(8) = jterm(8) - dswd_x1 * weight * point_q                          
c                                                                               
      if( omit_ct_elem ) jterm(7:8) = zero                                      
      if( .not.fgm_e .and. .not.fgm_nu .and. .not.seg_curves_flag )             
     &    jterm(7:8) = zero                                                     
c                                                                               
      if( debug ) then                                                          
         if( numprocs .gt. 1 ) write(out,*) "myid = ",myid                      
         write(out,920) (csig(i,ptno),i=1,9), (dstrain_x1(i),i=1,9)             
         write(out,930) csig_dstrain_x1, dswd_x1, weight,                       
     &                  point_q, jterm(7), jterm(8)                             
      end if                                                                    
c                                                                               
c              done with j-integral calculations for this integration point.    
c              calculation of jterm(5) for the surface integral occurs later.   
c                                                                               
      return                                                                    
c                                                                               
c                                                                               
 890  format(///,"J-integral terms: element",2x,i7,2x,"point",2x,i2)            
 892  format(' >>> weight, dqx, csig(10,ptno),jterm(1): ',4(1x,e13.6))          
 900  format(' >>> dux, dvx, dwx: ',3(1x,e13.6),/,                              
     &       '     temp1, 2, 3  : ',3(1x,e13.6 ),/,                             
     &       '     jterm(2)     : ',1x,e13.6 )                                  
 905  format(//,10x,'stress tensor:',/,3(10x,3(e11.4,2x),/))                    
 906  format(/,10x,'elem_alpha(1:6): ',6(e13.6,2x))                             
 907  format(/,10x,'dalpha_x1(1:6) : ',6(e13.6,2x))                             
 908  format(/,10x,'dtx            : ',e13.6,                                   
     &       /,10x,'weight         : ',e13.6,                                   
     &       /,10x,'point_q        : ',e13.6,                                   
     &       /,10x,'point_temp     : ',e13.6,                                   
     &       /,10x,'half           : ',e13.6)                                   
 910  format(' >>> jterm6 contributions: ',2(1x,e13.6),/,                       
     &       '     jterm(6)            : ',1x,e13.6 )                           
 920  format(//,10x,'stress tensor:',/,3(10x,3(e11.4,2x),/),                    
     &       /,10x,'dstrain_x1 tensor:',/,3(10x,3(e11.4,2x),/) )                
 930  format(' >>> jterm7, jterm8 contributions: ',                             
     &       /,10x,'csig_dstrain_x1: ',e13.6,                                   
     &       /,10x,'dswd_x1:         ',e13.6,                                   
     &       /,10x,'weight:          ',e13.6,                                   
     &       /,10x,'point_q:         ',e13.6,                                   
     &       /,10x,'jterm(7):        ',e13.6,                                   
     &       /,10x,'jterm(8):        ',e13.6  )                                 
c                                                                               
      end                                                                       
c                                                                               
