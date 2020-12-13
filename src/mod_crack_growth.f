c     ****************************************************************          
c     *                                                              *          
c     *                      f-90 module crack_growth_data           *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 11/20/2020 rhd             *          
c     *                                                              *          
c     *     define the data structures for crack growth              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c           module elem_extinct_data contains allocated arrays for              
c           the element extinction algorithms                                   
c                                                                               
c                                                                               
      module elem_extinct_data                                                  
c                                                                               
c                                                                               
      integer, allocatable, dimension(:)  :: dam_node_elecnt                    
      integer, allocatable, dimension(:)  :: dam_state                          
      double precision, allocatable       :: dam_ifv(:,:)                       
      double precision, allocatable       :: dam_dbar_elems(:,:)                
      logical, allocatable, dimension(:)  :: dam_blk_killed                     
      integer, allocatable, dimension(:)  :: dam_print_list                     
      integer, allocatable, dimension(:)  :: kill_order_list                    
      integer, allocatable                :: dam_face_nodes(:,:)                
      double precision, allocatable       :: old_porosity(:)                    
      double precision, allocatable       :: old_deff(:)                   
      double precision, allocatable       :: old_plast_strain(:)                
      double precision, allocatable       :: old_mises(:)                       
      double precision, allocatable       :: old_mean(:)    
      double precision, allocatable       :: smcs_weighted_T(:),
     &                                       smcs_old_epsplas(:),
     &                                       smcs_weighted_zeta(:),
     &                                       smcs_weighted_bar_theta(:)
c                                                                               
      end module  elem_extinct_data                                                            
c                                                                               
c                                                                               
c           module node_release_data contains allocated arrays for              
c           the node release algorithms                                         
c                                                                               
c                                                                               
      module node_release_data                                                  
c                                                                               
c                                                                               
      integer, allocatable, dimension(:)  :: crack_plane_nodes                  
      integer, allocatable, dimension(:)  :: inv_crkpln_nodes                   
      integer, allocatable, dimension(:)  :: num_neighbors                      
      integer, allocatable                :: neighbor_nodes(:,:)                
      integer, allocatable                :: crack_front_nodes(:,:)             
      allocatable                         :: crkpln_nodes_react(:)              
      integer, allocatable, dimension(:)  :: crkpln_nodes_state                 
      allocatable                         :: node_release_frac(:)               
      allocatable                         :: old_angles_at_front(:,:)           
      integer, allocatable, dimension(:)  :: master_nodes                       
      integer, allocatable, dimension(:,:):: crack_front_list                   
      integer, allocatable, dimension(:,:):: master_lines                       
c                                                                               
      double precision                                                          
     &  crkpln_nodes_react, node_release_frac, old_angles_at_front              
c                                                                               
c                                                                               
      end module node_release_data                                                               
                 
      module dam_param_code       

      contains                   
c                                          
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_param                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/5/20 rhd                *          
c     *                                                              *          
c     *     for a killable element not yet killed, determine if the  *          
c     *     element should be killed now. isolating decision here    *          
c     *     makes possible different kinds of killing criteria       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_param( elem, kill_now, debug, porosity,                    
     &                      eps_plas, eps_crit, sig_mean, sig_mises,            
     &                      ext_gurson, ext_shape, ext_spacing,
     &                      avg_triaxiality, avg_zeta, avg_bar_theta )
      use global_data, only : iprops, out  ! old common.main
c                                                                               
c        elem      -- (input)   element number to be checked for                
c                               killing                                         
c        kill_now  -- (output)  set .true. if element should be                 
c                               killed now                                      
c        debug     -- (input)   .true. if the calling routine wants             
c                               debugging info output                           
c        porosity  -- (output)  for Gurson-type material models,                
c                               this is the average element porosity.           
c                               Just used for output messages                   
c        eps_plas  -- (output)  for SMCS type models. this the                  
c                               average plastic strain over element.            
c                               Just used for output                            
c        eps_crit  -- (output)  for SMCS type models. this the                  
c                               critical plastic strain value.                  
c                               Just used for output                            
c        sig_mean  -- (output)  all models. average mean stress                 
c                               over element. just used for output              
c                               messages                                        
c        sig_mises -- (ouput)   all models. average mises stress                
c                               over element. just used for output.             
c        ext_gurson -- (output) logical = .true. if element death using         
c                               the extended gurson model                       
c        ext_shape -- (output)  W (shape) parameter for extended gurson model   
c        ext_spacing -- (output) X (spacing) parameter for extended             
c                                gurson model    
c        avg_triaxiality -- (output) plastic strain weighted average of
c                                triaxiality over loading history                               
c        avg_zeta  --  (output) plastic strain weighted average of
c                                Lode parameter (-1<=mean_zeta<=1)
c                                over loading history                               
c        avg_bar_theta --  (output) plastic strain weighted average of
c                                MIT MMC Lode parameter
c                                (-1<=mean_bar_theta<=1)
c                                over loading history                               
c                                                                               
c                                                                               
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
      use damage_data, only     : crack_growth_type 
c                                                         
      implicit none                                                    
c                                                                               
c              parameter declarations                                           
c                             
      integer, intent(in) :: elem      !  not optional                                             
      logical, optional,intent(out) :: kill_now, ext_gurson 
      logical, optional,intent(in) ::  debug 
      double precision, optional, intent(out) ::
     &      porosity, eps_plas, eps_crit, sig_mean, sig_mises, 
     &      ext_shape, ext_spacing, avg_triaxiality, avg_zeta,
     &      avg_bar_theta 
c                                                                               
c              local declarations                                               
c       
      integer :: mat_model
      logical :: kill_now_local
      double precision :: mean_zeta, mean_omega, triaxiality,
     &                    mean_bar_theta, sig_mean_local,
     &                    sig_mises_local, eps_plas_local,
     &                    eps_crit_local
      double precision, parameter :: zero = 0.d0, one = 1.0d0    
c                     
      double precision, dimension(:), pointer :: history, urcs_n, eps_n                         
c   
      if( present( ext_gurson ) ) ext_gurson = .false.                                                     
      if( present( ext_shape ) )  ext_shape = one                                                         
      if( present( ext_shape ) )  ext_shape= zero                                                        
c                                                                               
      select case( crack_growth_type )                         
c                                                                               
c                                                                               
c           1 gurson model- the material model 3 is the standard         
c                           Gurson model, 6 is the extended                
c                           Gurson model. call model dependent routine     
c                           to assess element status for killing and       
c                           to compute parameters.                         
c                                                                               
      case( 1) 
       mat_model = iprops(25,elem)                                               
       if( mat_model .eq. 3 ) then                                               
         call dam_param_gt( elem, kill_now_local, debug, porosity,                    
     &                      sig_mean_local, sig_mises_local )   
         if( present( kill_now ) ) kill_now = kill_now_local
         if( present( sig_mean ) ) sig_mean = sig_mean_local
         if( present( sig_mises ) ) sig_mises = sig_mises_local                            
       else if( mat_model .eq. 6 ) then                                          
         ext_gurson = .true.                                                    
         call dam_param_agt( elem, kill_now_local, debug, porosity,                   
     &                       sig_mean_local, sig_mises_local,
     &                       ext_spacing, ext_shape )                                        
         if( present( kill_now ) ) kill_now = kill_now_local
         if( present( sig_mean ) ) sig_mean = sig_mean_local
         if( present( sig_mises ) ) sig_mises = sig_mises_local                            
       else                                                                      
         write(out,9100) 1; write(out,9120)                                                        
         call die_gracefully                                                       
       end if                                                                    
c                                                                               
c           2  ctoa - not processed here                  
c                                                                               
      case(2) 
       write(out,9100) 2; write(out,9120)                                                      
       call die_gracefully                                                       
c                                                                               
c           3 smcs - failure based on a stress modified critical strain     
c
      case(3)   
       call dam_param_3_get_values( 
     &    elem, debug, eps_plas_local, eps_crit_local, sig_mean_local,
     &    sig_mises_local, triaxiality, mean_zeta, mean_omega,
     &    mean_bar_theta, 1, kill_now_local )   
       if( present( sig_mean ) ) sig_mean = sig_mean_local
       if( present( sig_mises ) ) sig_mises = sig_mises_local                            
       if( present( eps_plas ) ) eps_plas = eps_plas_local
       if( present( eps_crit ) ) eps_crit = eps_crit_local
       if( present( avg_triaxiality ) ) avg_triaxiality = triaxiality      
       if( present( avg_zeta ) )  avg_zeta = mean_zeta
       if( present( avg_bar_theta ) ) avg_bar_theta = mean_bar_theta
       if( present( kill_now ) ) kill_now = kill_now_local
c                                                                               
c           4 cohesive - not processed here               
c                                                                                                                                                             
      case(4)
        write(out,9100) 4; write(out,9120)                                                          
        call die_gracefully                                                       
c
      end select  
c
      return                                                              
c                                                                               
 9100 format('>>>> FATAL ERROR: invalid model no. in dam_param @:',i1 )                
 9120 format('                  job terminated....')                
c    
      end subroutine dam_param
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine dam_param_3_get_values               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 11/21/20 rhd               *          
c     *                                                              *          
c     *             SMCS compute values for this element             *
c     *                                                              *          
c     ****************************************************************          
c
                                                                     
      subroutine dam_param_3_get_values( elem, debug,                     
     &                      eps_plas, eps_crit, sig_mean, sig_mises,            
     &                      triaxiality, mean_zeta, mean_omega,
     &                      mean_bar_theta, dowhat, kill_now )  
c              
      use global_data, only : iprops, nstr, nstrs, out ! old common.main
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
      use damage_data, only     : dam_ptr, smcs_gamma, smcs_alpha, 
     &                            smcs_beta,  smcs_type, 
     &                            smcs_beta_1,smcs_beta_2, 
     &                            smcs_a_plus, smcs_a_minus,
     &                            smcs_kappa, smcs_alpha_1, 
     &                            smcs_alpha_2, max_eps_critical,
     &                            smcs_cutoff_triaxiality,
     &                            smcs_type_4_A, smcs_type_4_n,
     &                            smcs_type_4_c1, smcs_type_4_c2,
     &                            smcs_type_4_c3                                                          
      use elem_extinct_data, only : smcs_old_epsplas, smcs_weighted_T,
     &                              smcs_weighted_zeta,
     &                              smcs_weighted_bar_theta
c                                                                               
c        elem      -- (input)   element number to be checked for                
c                               killing                                         
c        debug     -- (input)   .true. if the calling routine wants             
c                               debugging info output                           
c        eps_plas  -- (output)  for SMCS type models. this the current                 
c                               average plastic strain over element.            
c                               Just used for output                            
c        eps_crit  -- (output)  for SMCS type models. this the                  
c                               current, predicted critical 
c                               plastic strain value.                  
c                               Just used for output                            
c        sig_mean  -- (output)  all models. average, current mean stress                 
c                               over element. just used for output              
c                               messages                                        
c        sig_mises -- (ouput)   all models. average, current mises 
c                               stress over element. for output.             
c        triaxiality -- (output) plastic strain weighted average of
c                                triaxiality over loading history 
c        mean_zeta -- (output)  Lode angle parameter average weighted
c                               by plastic strain
c        mean_omega -- (output)  Nahshon-Hutchinson Lode angle parameter
c                                average weighted by plastic strain
c        mean_bar_theta(output)  bar_theta angle in the MMC critical
c                                strain function from MIT group
c                                average weighted by plastic strain
c        dowhat    --  (input)   = 1 return damage variables - no update.
c                                    set kill now flag
c                                = 2 just return values - no update 
c                                = 3 return instantaneous values for 
c                                  non-killable elem
c                                = 4 compute/update damage variables  
c        kill_now  --  (output)  if dowhat = 1, set flag if critical
c                                condition reached                           
c                                                                                                                                                              
      implicit none                                                    
c                                                                               
c               parameter declarations                                           
c                             
      integer, intent(in) :: elem, dowhat                                                  
      logical,intent(in) :: debug
      logical, intent(out),optional :: kill_now                                       
      double precision,intent(out) :: eps_crit, sig_mean, sig_mises,
     &                                triaxiality, eps_plas, mean_zeta,
     &                                mean_omega, mean_bar_theta
c                                                                               
c               local declarations                                               
c
      double precision, dimension(:), pointer :: urcs_n, eps_n                         
      integer :: gp, blk, epsoffset, mat_model, elem_type,
     &           ngp, rel_elem, sigoffset, j, elem_ptr
      double precision ::  sig(6), eps(6), fpngp,
     &                     mean_triaxiality,
     &                     epspls_vec(27),s11, s22, s33, s12, s13, s23,
     &                     j2, j3, zeta, omega, t1, t2, t3, t4, 
     &                     p_omega, omega2, bar_theta
      logical :: update, local_debug, threed_solid_elem, l1, l2, l3,
     &           l4, l5, l6, l7, l8, l9, l10
     &          
      double precision, parameter :: zero = 0.d0,
     &     third = 1.d0/3.0d0, iroot2 = 1.0d0/dsqrt(2.0d0), six = 6.0d0,
     &     one = 1.0d0, two = 2.0d0, onept5 = 1.5d0, three = 3.d0,
     &     root3 = dsqrt(3.0d0), half = 0.5d0, thirteenpt5 = 13.5d0,
     &     small_plastic_strain = 1.0d-07, small_mises = 1.0d-07,
     &     type_2_tol = 1.0d-06, pi = 3.14159265359d0,
     &     type_4_toler = 0.0001d0
c                              
c               stress/strain vectors are (x,y,z,xy,yz,xz).                
c
      update      = dowhat .eq. 4
      local_debug = .false.  !   elem == 16
      if( local_debug )then
         write(out,*) '... dam_param_3_get_values, dowhat: ...', dowhat
      end if
      if( present(kill_now) ) kill_now = .false.
c
      elem_type   = iprops(1,elem)
      elem_ptr    = dam_ptr(elem)
      ngp         = iprops(6,elem)                                           
      blk         = elems_to_blocks(elem,1)                                  
      rel_elem    = elems_to_blocks(elem,2)                                  
      epsoffset   = (rel_elem-1)*nstr*ngp                                    
      sigoffset   = (rel_elem-1)*nstrs*ngp                                   
      urcs_n      => urcs_n_blocks(blk)%ptr                                  
      eps_n       => eps_n_blocks(blk)%ptr    
c
      eps_plas    = zero
      eps_crit    = max_eps_critical
      sig_mean    = zero
      sig_mises   = zero
      triaxiality = zero
      mean_zeta   = zero
      mean_omega  = zero  
      mean_bar_theta = zero                             
      sig         = zero
      eps         = zero  !   strains not used current code       
      call set_element_type( elem_type, threed_solid_elem, l1, l2, l3,
     &                       l4, l5, l6, l7, l8, l9, l10 )  
      if( .not. threed_solid_elem ) return                                 
c
      call mm_return_values( "plastic_strain", elem, epspls_vec, ngp )
c
c               average stresses over all integration points
c
      if( ngp == 8 ) then  ! for optimizer
       do gp = 1, 8    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
         eps(j) = eps(j) + eps_n(epsoffset+j)
        end do                                      
        eps_plas  = eps_plas + epspls_vec(gp)
        epsoffset = epsoffset + nstr                                         
        sigoffset = sigoffset + nstrs                                        
       end do   
      else  ! general case
       do gp = 1, ngp    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
         eps(j) = eps(j) + eps_n(epsoffset+j)
        end do                                      
        eps_plas  = eps_plas + epspls_vec(gp)
        epsoffset = epsoffset + nstr                                         
        sigoffset = sigoffset + nstrs                                        
       end do 
      end if  
c
c               average mean stress, mises stress, plastic strain,
c               Lode parameters
c
      fpngp     = one / dble( ngp )
      eps_plas  = fpngp * eps_plas     !  averaged at element center
      sig = sig * fpngp                !    ""
      sig_mean  = (sig(1) + sig(2) + sig(3)) * third                    
      s11 = sig(1) - sig_mean
      s22 = sig(2) - sig_mean
      s33 = sig(3) - sig_mean
      s12 = sig(4)
      s13 = sig(5)
      s23 = sig(6)
      j2 = half * ( s11**2 + s22**2 + s33**2 + 
     &              two*(s12**2 + s23**2 + s13**2))
      j3 = s11*s22*s33 + two*s12*s23*s13 - s11*s23*s23 -
     &     s22*s13*s13 - s33*s12*s12
      sig_mises = sqrt( three * j2 )
      zeta      = zero  !  called \xi in manual
      omega     = zero
      bar_theta = zero
      if( sig_mises > small_mises ) then
        zeta  = onept5*root3 * j3 / j2**onept5
        omega = one - zeta*zeta 
        t1    = thirteenpt5 * j3 / sig_mises**3
        bar_theta = one - ( acos( t1 ) * two / pi )
        if( abs(bar_theta) > one+type_4_toler ) then
          write(out,9120) elem, bar_theta  ! improve message
          call die_abort
        endif
      end if
      if( local_debug ) write(out,9100) sig_mean, sig_mises, eps_plas
c
      if( dowhat == 3 ) then ! compute values for non-killable elem
         eps_crit       = max_eps_critical
         mean_zeta      = zeta
         mean_omega     = omega
         mean_bar_theta = bar_theta
         triaxiality = 100.0d0   ! in cases mises -> 0                                             
         if( sig_mises .gt. small_mises )
     &         triaxiality = sig_mean / sig_mises 
         return
       end if
c
c               update plastic strain weighted triaxiality, zeta,
c               omega, bar_theta. update saved plastic strain.
c               only update if instantaneous triaxiality exceeds
c               user-specified cutoff value
c              
      if( update ) call dam_param_3_get_values_a
c
c               if no plastic strain yet, set defaults and return.
c               cannot compute weighed averages over history (divide
c               by zero plastic strain).
c
      if( eps_plas < small_plastic_strain ) then 
        triaxiality = zero
        eps_crit    = max_eps_critical
        mean_zeta   = zero
        mean_omega  = zero
        mean_bar_theta = zero
        return
      end if
c
c               use weighted triaxiality, zeta, omega, bar_theta
c               values to set critical plastic strain
c
      mean_triaxiality =  smcs_weighted_T(elem_ptr) / eps_plas
      triaxiality    = mean_triaxiality  ! to return
      mean_zeta      = smcs_weighted_zeta(elem_ptr) / eps_plas
      mean_omega     = one - mean_zeta**2
      mean_bar_theta = smcs_weighted_bar_theta(elem_ptr) / eps_plas
      eps_crit       = max_eps_critical
      if( triaxiality <= smcs_cutoff_triaxiality ) return
c
      if( smcs_type .eq. 1 ) then
c
          t2 = exp(-smcs_kappa*abs(mean_zeta))
          t1 = smcs_alpha * exp(-smcs_beta * mean_triaxiality)
          eps_crit = smcs_gamma + t1*t2
          if( local_debug ) write(out,9110) mean_triaxiality,eps_crit
c
      elseif( smcs_type .eq. 2 ) then
c
          t3 = exp(-smcs_kappa*abs(mean_zeta))
          t2 = smcs_beta_2 * exp(-smcs_a_minus * mean_triaxiality)
          t1 = smcs_beta_1 * exp(smcs_a_plus * mean_triaxiality)
          t4 = t1 - t2
          if( abs(t4) < type_2_tol .or. t4 < zero ) then
            eps_crit = max_eps_critical
          else
            eps_crit = smcs_gamma + t3 / t4
          end if
c
      elseif( smcs_type .eq. 3 ) then
c
          p_omega = (one - mean_omega)**2
          t1 = smcs_alpha_1*exp(-smcs_beta_1 * mean_triaxiality )
          t2 = smcs_alpha_2*exp(-smcs_beta_2 * mean_triaxiality )
          eps_crit = (one-p_omega)*t1 + p_omega*t2
c
      elseif( smcs_type .eq. 4 ) then
c
        block
          double precision :: A, n, c1, c2, c3, eta, t1, t2, ta, tb,
     &                        t1c, t1s, t1sec
          A     = smcs_type_4_A
          n     = smcs_type_4_n
          c1    = smcs_type_4_c1
          c2    = smcs_type_4_c2
          c3    = smcs_type_4_c3
          eta   = mean_triaxiality
          t1    = mean_bar_theta * pi / six
          t1c   = cos( t1 )
          t1s   = sin( t1 )
          t1sec = one / t1c
          t2 = root3 / ( two - root3 )
          ta = (A/c2) * ( c3 + t2*(one-c3)*(t1sec-one) )
          tb = sqrt( third*(one+c1*c1) ) * t1c + 
     &         c1*(eta+third*t1s)
          eps_crit = ( ta * tb )**(-one/n)
        end block
c
      else
          write(out,9130) elem, smcs_type
          call die_abort
      end if
c      
      if( present( kill_now ) ) kill_now = 
     &             ( eps_plas - eps_crit ) >= zero .or. 
     &             eps_plas >= max_eps_critical
c
      return
c
 9100 format(10x,'mean, mises, eps_pls:',2f10.3,f12.6)
 9110 format(10x,'T_mean, eps crit: ',f10.3,f12.6) 
 9120 format(/1x,'>>>>> error: inconsistency in bar_theta computation',
     & ' for smcs type 4',/,
     & 14x,'element: ',i8,' bar_theta value: ',f6.2,/,
     & 14x,'value must lie in range (-1.0,1.0)',/,
     & 14x,'job aborted by routine dam_param_3_get_values',/)
 9130 format(/1x,'>>>>> error: inconsistency in smcs_type',/,
     & 14x,'element: ',i8,' smcs_type: ',i5,/,
     & 14x,'job aborted by routine dam_param_3_get_values',/)
c
      contains
c     ========
c
      subroutine dam_param_3_get_values_a
      implicit none
c
      double precision :: deps_plas, now_triaxiality
c
      deps_plas = eps_plas - smcs_old_epsplas(elem_ptr)
      smcs_old_epsplas(elem_ptr) = eps_plas 
c
      if( sig_mises <= small_mises ) return
c
      now_triaxiality = sig_mean / sig_mises 
      if( now_triaxiality < smcs_cutoff_triaxiality ) return
c
      smcs_weighted_T(elem_ptr) = 
     &      smcs_weighted_T(elem_ptr) + deps_plas * now_triaxiality
      smcs_weighted_zeta(elem_ptr) = 
     &      smcs_weighted_zeta(elem_ptr) + deps_plas * zeta
      smcs_weighted_bar_theta(elem_ptr) = 
     &      smcs_weighted_bar_theta(elem_ptr) + deps_plas * bar_theta
c
      return
      end subroutine dam_param_3_get_values_a

      end subroutine dam_param_3_get_values

c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine growth_set_dbar                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/6/2020 rhd              *          
c     *                                                              *          
c     *  supports traction-separation law release option. compute    *          
c     *  elongation of element normal to crack plane when first      *          
c     *  killed (option=0) or at some load step later (option=1)     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine growth_set_dbar( elem, elem_ptr, debug, action,                
     &                            dbar_now )  
c                                  
      use global_data, only : out, u, du, dstmap
      use elem_extinct_data, only : dam_face_nodes, dam_dbar_elems              
      use damage_data, only : gurson_cell_size, crk_pln_normal_idx
c                                                                               
      implicit none
c                
      integer, intent(in) :: elem, elem_ptr, action                                    
      logical, intent(in) :: debug
      double precision, optional, intent(out) :: dbar_now          
c
      integer :: i, snode
      integer, parameter :: num_face_nodes = 4
      logical, parameter :: local_debug = .false.                                                
      double precision :: sum, node_displ(3), dbar_now_local        
      double precision, parameter :: zero = 0.d0                           
c                                                                               
      if( debug ) write (out,*) '>>>> in growth_set_dbar'                      
c                                                                               
      sum = zero       
c                                                         
      if( action .lt. 0 ) then                                                 
        do i = 1, num_face_nodes                                                
          snode         = dam_face_nodes(i,elem_ptr)                            
          node_displ(1) = u(dstmap(snode)+0)                                    
          node_displ(2) = u(dstmap(snode)+1)                                    
          node_displ(3) = u(dstmap(snode)+2)                                    
          sum           = sum + node_displ(crk_pln_normal_idx)                  
        end do                                                                  
      end if   
c                                                                 
      if( action .gt. 0 ) then                                                 
        do i = 1, num_face_nodes                                                
          snode         = dam_face_nodes(i,elem_ptr)                            
          node_displ(1) = u(dstmap(snode)+0) + du(dstmap(snode)+0)              
          node_displ(2) = u(dstmap(snode)+1) + du(dstmap(snode)+1)              
          node_displ(3) = u(dstmap(snode)+2) + du(dstmap(snode)+2)              
          sum           = sum + node_displ(crk_pln_normal_idx)                  
        end do                                                                  
      end if                                                                    
c                                                                               
      sum = sum / dble(num_face_nodes)                                          
      dbar_now_local = gurson_cell_size + sum     
      if( present( dbar_now ) ) dbar_now = dbar_now_local
c                                  
      if( abs(action) .eq. 1 ) then                                            
       dam_dbar_elems(1,elem_ptr) = dbar_now_local                                    
       dam_dbar_elems(2,elem_ptr) = zero                                        
      end if                                                                    
c                                                                               
      if( local_debug ) write(out,9000)  action, elem, sum,
     &                                   dbar_now_local                            
c                                                                               
      return                                                                    
c                                                                               
 9000 format(' > element release type 2. Cell height computation.',             
     & /,5x,'action: ',i1,10x,'element: ',i7,' avg. face displ: ',f10.6,        
     & ' D-bar now: ',f10.6 )                                                   
c                                                                               
      end subroutine growth_set_dbar        
                                                               
     
                                                               
      end module   dam_param_code                                                            
