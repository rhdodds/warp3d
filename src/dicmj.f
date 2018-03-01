c ***************************************************************               
c *                                                             *               
c * domain_compute - drive execution of element routine to      *               
c *                  compute j and i-integrals for a single     *               
c *                  domain                                     *               
c *                                                             *               
c ***************************************************************               
c                                                                               
c                                                                               
      subroutine dicmj                                                          
      use global_data ! old common.main
      use main_data                                                             
      use elem_block_data, only : history_blocks, rot_n1_blocks,                
     &                            urcs_n_blocks, cdest_blocks,                  
     &                            edest_blocks, eps_n_blocks,                   
     &                            history_blk_list                              
      use j_data                                                                
      implicit integer (a-z)                                                    
c                                                                               
c                 parameter declarations (none)                                 
c                                                                               
c                                                                               
c                 local declarations. includes pointers to simplify             
c                 addressing into blocks.                                       
c                                                                               
      double precision                                                          
     & e_coord(3,mxndel), e_displ(3,mxndel), e_vel(3,mxndel),                   
     & e_accel(3,mxndel), q_element(mxndel), e_jresults(8),                     
     & e_stress(nstrs*mxgp), e_rots(9*mxgp),                                    
     & e_force(3,mxndel), diterms(8),                                           
     & zero, jtotal, itotal(8), one, two, symm_factor, di_value_j,              
     & sum_force, sum_accel, sum_velocities, toler, static_di,                  
     & static_total, toler_static_j, e_node_temps(mxndel),                      
     & elem_uniform_temp, e_alpha_ij(6,mxndel), eq_load_modifier,               
     & dummy, elem_nod_strains(6,mxndel), elem_nod_swd(mxndel),                 
     & iiterms(8,8), e_strain(9,mxgp), di_value_i(8), sign,                     
     & e_iresults(8,8), Ki, T11, T33, T13, e(mxndel), nu(mxndel),               
     & cf_load(3)                                                               
c                                                                               
       real e_props(mxelpr)                                                     
       logical dibmck, q_same, skip_element, geonl, chk_killed,                 
     &         elem_temps, omit_ct_elem, face_loads_possible,                   
     &         front_elem_flag, seg_curves_flag                                 
       integer eps_offset, gpn, snodes(mxndel), node, node_id,                  
     &         orig_node                                                        
       integer, dimension (:,:), pointer :: edest, cdest                        
       double precision,                                                        
     & dimension(:), pointer :: urcs_n, history, rot_n1, eps_n                  
       double precision,                                                        
     & dimension(:), allocatable :: elem_hist                                   
c                                                                               
       data zero, one, two, toler_static_j / 0.0d0, 1.0d0, 2.0d0,               
     &                                           1.0d-5 /                       
c                                                                               
c          if running parallel MPI:                                             
c             alert all slaves to enter this routine, then send all             
c             data required for J calculations that is only on root.            
c          if running serial:                                                   
c             this routines are dummy routines which return immerdiately        
c                                                                               
c          dummy allocate element history vector to simplify logic              
c          in element loop                                                      
c                                                                               
      call wmpi_alert_slaves( 27 )                                              
      call wmpi_send_jint                                                       
c                                                                               
      iout = out                                                                
      if ( debug_driver ) write(iout,*) ' >>> entered element driver'           
c                                                                               
c             keep track of ring_count for storing and retrieving output        
c                                                                               
      ring_count = ring_count + 1                                               
c                                                                               
c          initialize totals for domain and set symmetry factor.                
c          print necessary column headers.                                      
c                                                                               
      diterms(1:8)     = zero                                                   
      iiterms(1:8,1:8) = zero                                                   
      hist_size_old    = 0                                                      
      skipped_killed   = 0                                                      
      symm_factor      = one                                                    
      if( symmetric_domain ) symm_factor = two                                  
      if ( myid .eq. 0 ) then                                                   
         if( print_elem_values ) write(iout,9100)                               
         if( print_elem_values ) write(iout,9110)                               
      end if                                                                    
      allocate( elem_hist( 1 ) )                                                
c                                                                               
c          loop over all structure elements. if element is involved             
c          in integral computations, call element dependent routine.            
c                                                                               
      do elemno = 1, noelem                                                     
c                                                                               
c          handle only elements whose data is owned by this processor.          
c                                                                               
         blk = elems_to_blocks(elemno,1)                                        
         if( elblks(2,blk) .ne. myid ) cycle                                    
         if( .not. dibmck( elemno, q_element_maps, bits ) ) go to 10            
         if( chk_killed( elemno ) ) then                                        
            skipped_killed = skipped_killed + 1                                 
            go to 10                                                            
         end if                                                                 
         num_enodes = iprops(2,elemno)                                          
         num_gpts   = iprops(6,elemno)                                          
         geonl      = lprops(18,elemno)                                         
         incptr     = incmap(elemno)                                            
         k1         = num_enodes                                                
         k2         = 2*num_enodes                                              
         blk        = elems_to_blocks(elemno,1)                                 
         rel_elem   = elems_to_blocks(elemno,2)                                 
         cdest      => cdest_blocks(blk)%ptr                                    
         edest      => edest_blocks(blk)%ptr                                    
         hist_size  = history_blk_list(blk)                                     
c                                                                               
c           build copy of element coordinates, displacements,                   
c           velocities, accelerations, q-values from global data.               
c           zero load terms that correspond to constrained dof.                 
c           if we have temperature loadings, get temps for element              
c           nodes (node values + element uniform value). rotate                 
c           displacements, velocities and accelerations from                    
c           constraint-compatible to global coordinates. obtain                 
c           element nodal values of stress work density and strain              
c           from structure node-value arrays.                                   
c                                                                               
         seg_curves_flag = .false.                                              
         do enode = 1, num_enodes                                               
            e_coord(1,enode) = c(cdest(enode,rel_elem))                         
            e_coord(2,enode) = c(cdest(enode+k1,rel_elem))                      
            e_coord(3,enode) = c(cdest(enode+k2,rel_elem))                      
            snode            = incid(incptr+enode-1)                            
            snodes(enode)    = snode                                            
            q_element(enode) = q_values(snode)                                  
            i1               = edest(enode,   rel_elem)                         
            i2               = edest(k1+enode,rel_elem)                         
            i3               = edest(k2+enode,rel_elem)                         
            e_displ(1,enode) = u(i1)                                            
            e_displ(2,enode) = u(i2)                                            
            e_displ(3,enode) = u(i3)                                            
            e_vel(1,enode)   = v(i1)                                            
            e_vel(2,enode)   = v(i2)                                            
            e_vel(3,enode)   = v(i3)                                            
            e_accel(1,enode) = a(i1)                                            
            e_accel(2,enode) = a(i2)                                            
            e_accel(3,enode) = a(i3)                                            
            e(enode)         = props(7,elemno)                                  
            nu(enode)        = props(8,elemno)                                  
            if( fgm_e ) then                                                    
               e(enode)  = fgm_node_values(snode,1)                             
            end if                                                              
            if( fgm_nu ) then                                                   
               nu(enode) = fgm_node_values(snode,2)                             
            end if                                                              
            if( temperatures .or. temperatures_ref ) then                       
               if( block_seg_curves(blk) ) then                                 
                  seg_curves_flag = .true.                                      
                  e(enode)        = seg_snode_e(snode)                          
                  nu(enode)       = seg_snode_nu(snode)                         
               end if                                                           
            end if                                                              
c                                                                               
c            if necessary, transform nodal values from                          
c            constraint-compatible coordinates to global                        
c            coordinates.                                                       
c                                                                               
            if( trn(snode) ) then                                               
               call di_trans_nodvals( e_displ(1,enode),                         
     &                                trnmat(snode)%mat(1,1))                   
               call di_trans_nodvals( e_vel(1,enode),                           
     &                                trnmat(snode)%mat(1,1))                   
               call di_trans_nodvals( e_accel(1,enode),                         
     &                                trnmat(snode)%mat(1,1))                   
            end if                                                              
c                                                                               
            e_force(1,enode)          = zero                                    
            e_force(2,enode)          = zero                                    
            e_force(3,enode)          = zero                                    
            e_node_temps(enode)       = zero                                    
            e_alpha_ij(1,enode)       = zero                                    
            e_alpha_ij(2,enode)       = zero                                    
            e_alpha_ij(3,enode)       = zero                                    
            e_alpha_ij(4,enode)       = zero                                    
            e_alpha_ij(5,enode)       = zero                                    
            e_alpha_ij(6,enode)       = zero                                    
            elem_nod_swd(enode)       = swd_at_nodes(snode)                     
            elem_nod_strains(1,enode) = strain_at_nodes(1,snode)                
            elem_nod_strains(2,enode) = strain_at_nodes(2,snode)                
            elem_nod_strains(3,enode) = strain_at_nodes(3,snode)                
            elem_nod_strains(4,enode) = strain_at_nodes(4,snode)                
            elem_nod_strains(5,enode) = strain_at_nodes(5,snode)                
            elem_nod_strains(6,enode) = strain_at_nodes(6,snode)                
         end do                                                                 
c                                                                               
c           for problems with temperature loads, pull out                       
c           temperatures of element nodes and set the thermal                   
c           expansion coefficients at element nodes. remove reference           
c           temperatures of model nodes.                                        
c                                                                               
         elem_temps = .false.                                                   
         if( temperatures .or. temperatures_ref ) then                          
            elem_uniform_temp = temper_elems(elemno)                            
            do enode = 1, num_enodes                                            
               snode = incid(incptr + enode - 1)                                
               e_node_temps(enode) = temper_nodes(snode)                        
     &                             + elem_uniform_temp                          
     &                             - temper_nodes_ref(snode)                    
               if( abs(e_node_temps(enode)).gt.zero ) elem_temps=.true.         
            end do                                                              
            do enode = 1, num_enodes                                            
               snode = incid( incptr + enode - 1 )                              
               e_alpha_ij(1,enode) = snode_alpha_ij(snode,1)                    
               e_alpha_ij(2,enode) = snode_alpha_ij(snode,2)                    
               e_alpha_ij(3,enode) = snode_alpha_ij(snode,3)                    
               e_alpha_ij(4,enode) = snode_alpha_ij(snode,4)                    
               e_alpha_ij(5,enode) = snode_alpha_ij(snode,5)                    
               e_alpha_ij(6,enode) = snode_alpha_ij(snode,6)                    
            end do                                                              
         end if                                                                 
c                                                                               
c                 crack-face tractions. the contribution to J arising           
c                 from crack-face tractions is calculated element-by-element    
c                 using the equivalent nodal loads on each element face         
c                 computed during setup of load step for applied pressures,     
c                 surface tractions and body forces. here we just get           
c                 the element equiv. nodal forces. they are passed to           
c                 element J routine which sorts out the loaded face.            
c                                                                               
c                 a similar procedure is used to compute the contribution       
c                 to I.                                                         
c                                                                               
         face_loads_possible = .not. ignore_face_loads .and.                    
     &                               eq_node_force_len .gt. 0                   
         if( face_loads_possible ) then                                         
            if( eq_node_force_indexes(elemno) .ne. 0 ) then                     
               call vec_ops( e_force,                                           
     &              eq_node_forces(eq_node_force_indexes(elemno)),              
     &              dummy, 3*num_enodes, 5 )                                    
            end if                                                              
         end if                                                                 
c                                                                               
c                  determine if we really need to process element               
c                  based on q-values, nodal velocities, accelerations           
c                  and specified temperature loadings.                          
c                                                                               
         q_same = .true.                                                        
         do enode = 1, num_enodes                                               
            if( q_element(enode) .ne. one ) q_same = .false.                    
         end do                                                                 
         sum_velocities = zero                                                  
         sum_accel      = zero                                                  
         sum_force      = zero                                                  
         do enode = 1, num_enodes                                               
            sum_velocities = sum_velocities + abs(e_vel(1,enode)) +             
     &                       abs(e_vel(2,enode)) + abs(e_vel(3,enode))          
            sum_accel      = sum_accel + abs(e_accel(1,enode)) +                
     &                       abs(e_accel(2,enode))+abs(e_accel(3,enode))        
            sum_force      = sum_force + abs(e_force(1,enode)) +                
     &                       abs(e_force(2,enode))+abs(e_force(3,enode))        
         end do                                                                 
         toler = 1.0e-04                                                        
         skip_element = q_same .and. sum_velocities .lt. toler                  
     &                  .and. sum_accel .lt. toler                              
     &                  .and. sum_force .lt. toler                              
     &                  .and. (.not. elem_temps)                                
     &                  .and. (.not. fgm_e )                                    
     &                  .and. (.not. fgm_nu )                                   
         omit_ct_elem = .false.                                                 
         if( fgm_e .and. omit_crack_front_elems) then                           
            omit_ct_elem = crack_front_elem(elemno)                             
         end if                                                                 
         if( fgm_nu .and. omit_crack_front_elems) then                          
            omit_ct_elem = crack_front_elem(elemno)                             
         end if                                                                 
         if ( skip_element ) go to 10                                           
c                                                                               
c                  gather element stresses                                      
c                                                                               
         offset = (rel_elem-1) * nstrs * num_gpts                               
         urcs_n =>  urcs_n_blocks(blk)%ptr                                      
         do i = 1, nstrs * num_gpts                                             
            e_stress(i) = urcs_n(offset+i)                                      
         end do                                                                 
c                                                                               
c                  gather element histories                                     
c                                                                               
        offset = (rel_elem-1) * hist_size * num_gpts                            
        history => history_blocks(blk)%ptr                                      
        hist_size_new = hist_size * num_gpts                                    
        if( hist_size_new .ne. hist_size_old ) then                             
            deallocate( elem_hist )                                             
            allocate( elem_hist(hist_size_new) )                                
        end if                                                                  
        hist_size_old = hist_size_new                                           
        do i = 1, hist_size_new                                                 
          elem_hist(i) = history(offset+i)                                      
        end do                                                                  
c                                                                               
c                  gather element properties                                    
c                                                                               
        do i = 1, mxelpr                                                        
           e_props(i) = props(i,elemno)                                         
        end do                                                                  
c                                                                               
c                  gather 3x3 rotation matrices from polar decompositions       
c                  at each gauss point of element for geonl                     
c                                                                               
         if ( geonl ) then                                                      
            rot_n1 => rot_n1_blocks(blk)%ptr                                    
            offset =  (rel_elem-1)*9*num_gpts                                   
            do i = 1, 9 * num_gpts                                              
               e_rots(i) = rot_n1(offset+i)                                     
            end do                                                              
         end if                                                                 
c                                                                               
         e_strain(1:9,1:mxgp) = zero                                            
         if( comput_i ) then                                                    
c                                                                               
c                  gather element strains at integration points.                
c                  and store in tensor form. convert engineering                
c                  (total) strains of off-diagonal terms to tensor              
c                  strains.                                                     
c                                                                               
            eps_n => eps_n_blocks(blk)%ptr                                      
            eps_offset = (rel_elem - 1) * nstr * num_gpts                       
c                                                                               
            do gpn = 1, num_gpts                                                
               e_strain(1,gpn) = eps_n(eps_offset + 1)                          
               e_strain(2,gpn) = eps_n(eps_offset + 4)/two                      
               e_strain(3,gpn) = eps_n(eps_offset + 6)/two                      
               e_strain(4,gpn) = eps_n(eps_offset + 4)/two                      
               e_strain(5,gpn) = eps_n(eps_offset + 2)                          
               e_strain(6,gpn) = eps_n(eps_offset + 5)/two                      
               e_strain(7,gpn) = eps_n(eps_offset + 6)/two                      
               e_strain(8,gpn) = eps_n(eps_offset + 5)/two                      
               e_strain(9,gpn) = eps_n(eps_offset + 3)                          
               eps_offset      = eps_offset + nstr                              
            end do                                                              
            if( debug_driver ) then                                             
               do gpn = 1, num_gpts                                             
                  write(iout,9130) elemno, gpn, (e_strain(i,gpn),i=1,9)         
               end do                                                           
            end if                                                              
c                                                                               
c                  for crack-face loading, make a copy of the tractions         
c                  input by the user in the domain definition. if none          
c                  was input, dielem_c.f automatically uses the equivalent      
c                  nodal loads used to solve the boundary value problem.        
c                                                                               
            cf_load(1) = cf_tractions(1)                                        
            cf_load(2) = cf_tractions(2)                                        
            cf_load(3) = cf_tractions(3)                                        
c                                                                               
         end if                                                                 
c                                                                               
c                  set flag for crack-front element.                            
c                                                                               
         front_elem_flag = .false.                                              
         if( crack_front_elem(elemno) ) front_elem_flag = .true.                
c                                                                               
c                  call the element routine to compute contribution             
c                  to the j-integral and or i-integral for domain               
c                                                                               
         element    = elemno                                                    
         numrow_sig = nstrs                                                     
         call dielem( e_coord, q_element, e_displ, e_vel, e_accel,              
     &                e_force, e_stress, elem_hist, e_props, e_props,           
     &                e_props, e_rots, hist_size, domain_rot, iout,             
     &                e_jresults, e_node_temps, elem_temps, e_alpha_ij,         
     &                ierr, element, debug_elements, one_point_rule,            
     &                geonl, numrow_sig, snodes, elem_nod_swd,                  
     &                elem_nod_strains, omit_ct_elem, fgm_e, fgm_nu,            
     &                e_strain, e, e_front, nu, nu_front, front_nodes,          
     &                num_front_nodes, front_coords, domain_origin,             
     &                domain_type, front_order, e_iresults, comput_i,           
     &                comput_j, cf_traction_flags, cf_load,                     
     &                front_elem_flag, expanded_front_nodes,                    
     &                myid, numprocs, crack_curvature, face_loading,            
     &                seg_curves_flag, process_temperatures,                    
     &                max_exp_front )                                           
c                                                                               
         if( comput_j ) then                                                    
            do i=1,8                                                            
               e_jresults(i) = e_jresults(i) * symm_factor                      
               diterms(i)    = diterms(i) + e_jresults(i)                       
            end do                                                              
            if ( print_elem_values .and. comput_j ) then                        
               write(iout,9140) "J     ", elemno, (e_jresults(i),i=1,8)         
            end if                                                              
         end if                                                                 
c                                                                               
         if( comput_i ) then                                                    
            do i=1,8                                                            
               do j=1,8                                                         
                  e_iresults(i,j) = e_iresults(i,j) * symm_factor               
                  iiterms(i,j)    = iiterms(i,j) + e_iresults(i,j)              
               end do                                                           
            end do                                                              
         end if                                                                 
c                                                                               
         if( comput_i .and. symmetric_domain ) then                             
            iiterms(1:8,3:5) = zero                                             
            iiterms(1:8,8)   = zero                                             
         end if                                                                 
c                                                                               
         if ( print_elem_values .and. comput_i ) then                           
            write(iout,9140) "I_KI  ", elemno, (e_iresults(i,1),i=1,8)          
            write(iout,9140) "I_KI  ", elemno, (e_iresults(i,2),i=1,8)          
            write(iout,9140) "I_KII ", elemno, (e_iresults(i,3),i=1,8)          
            write(iout,9140) "I_KII ", elemno, (e_iresults(i,4),i=1,8)          
            write(iout,9140) "I_KIII", elemno, (e_iresults(i,5),i=1,8)          
            write(iout,9140) "I_T11 ", elemno, (e_iresults(i,6),i=1,8)          
            write(iout,9140) "I_T11 ", elemno, (e_iresults(i,7),i=1,8)          
            write(iout,9140) "I_T13 ", elemno, (e_iresults(i,8),i=1,8)          
            write(iout,9150)                                                    
         end if                                                                 
 10      continue                                                               
      end do                                                                    
c                                                                               
c             release local block for element history data                      
c                                                                               
      deallocate( elem_hist )                                                   
c                                                                               
c             for parallel computation, reduce all processor contributions      
c             to the J-value and I-value terms for this domain.                 
c             slaves are finished after this procedure.                         
c                                                                               
c      do proc = 0, numprocs - 1                                                
c        if ( proc .eq. myid) then                                              
c           write (*,*) myid,':  diterms:',(diterms(i),i=1,8)                   
c        endif                                                                  
c      enddo                                                                    
c                                                                               
      call wmpi_reduce_vec( diterms, 8 )                                        
      call wmpi_reduce_vec( iiterms(1,1), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,2), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,3), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,4), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,5), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,6), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,7), 8 )                                   
      call wmpi_reduce_vec( iiterms(1,8), 8 )                                   
      call wmpi_redlog( face_loading )                                          
c                                                                               
      if ( slave_processor ) then                                               
         deallocate( q_values, q_element_maps, crack_front_elem,                
     &               expanded_front_nodes )                                     
          return                                                                
      end if                                                                    
c                                                                               
c      if ( 0 .eq. myid) then                                                   
c         write (*,*) myid,':  collected diterms:',(diterms(i),i=1,8)           
c      endif                                                                    
c                                                                               
c             done with this domain computation. print values                   
c             as required based on user specified flags.                        
c                                                                               
c             J-integral results:                                               
c                                                                               
      if( comput_j ) then                                                       
         jtotal = diterms(1) + diterms(2) + diterms(3) + diterms(4) +           
     &            diterms(5) + diterms(6) + diterms(7) + diterms(8)             
         if ( abs(diterms(3)+diterms(4)).lt.toler_static_j*abs(jtotal))         
     &        static_j = .true.                                                 
         static_total = diterms(1) + diterms(2) + diterms(5)                    
     &                + diterms(6) + diterms(7) + diterms(8)                    
         di_value_j   = jtotal / front_q_area                                   
         static_di    = static_total / front_q_area                             
         domain_min_j = min( domain_min_j, di_value_j )                         
         domain_max_j = max( domain_max_j, di_value_j )                         
         domain_avg_j = domain_avg_j + di_value_j                               
         static_min   = min( static_min, static_di  )                           
         static_max   = max( static_max, static_di  )                           
         static_avg   = static_avg + static_di                                  
c                                                                               
c             calculate stress intensity factor K, from J, for:                 
c                                                                               
c               1. plane stress where loading is pure mode I or mode II.        
c               2. plane strain where loading is pure mode I or mode II.        
c               3. anti-plane shear where loading is pure mode III.             
c                                                                               
c             if J is negative, calculate K from |J|, and then change sign.     
c                                                                               
         sign = di_value_j / abs(di_value_j)                                    
c                                                                               
         ks_from_j(ring_count,1) = sqrt(abs(di_value_j) * e_front)*sign         
c                                                                               
         ks_from_j(ring_count,2) =                                              
     &        sqrt(abs(di_value_j) * e_front/(one - nu_front**2))*sign          
c                                                                               
         ks_from_j(ring_count,3) =                                              
     &        sqrt(abs(di_value_j) * e_front/(one + nu_front)) * sign           
c                                                                               
c             store J-integral values for output in didriv.f                    
c                                                                               
      j_storage(ring_count,1)   = dble( nowring )                               
          j_storage(ring_count,2:9) = diterms(1:8)                              
          j_storage(ring_count,10)  = di_value_j                                
      j_storage(ring_count,11)  = dble(skipped_killed)                          
      end if                                                                    
c                                                                               
c             interaction-integral results:                                     
c             calculate K from relationship with interaction integral           
c                                                                               
c             the 8 terms for I-integral results are stored                     
c             in array e_iresults(8,7) as follows:                              
c                                                                               
c                               value    auxiliary field                        
c                                                                               
c                  iiterms(i,1): KI       plane stress                          
c                  iiterms(i,2): KI       plane stress                          
c                  iiterms(i,3): KII      plane stress                          
c                  iiterms(i,4): KII      plane stress                          
c                  iiterms(i,5): KIII     anti-plane shear                      
c                  iiterms(i,6): T11,T33  plane stress                          
c                  iiterms(i,7): T11,T33  plane strain                          
c                  iiterms(i,8): T13      anti-plane strain                     
c                                                                               
      if( comput_i ) then                                                       
         do j=1,8                                                               
            itotal(j) = iiterms(1,j) + iiterms(2,j) + iiterms(3,j)              
     &                + iiterms(4,j) + iiterms(5,j) + iiterms(6,j)              
     &                + iiterms(7,j) + iiterms(8,j)                             
            di_value_i(j) = itotal(j) / front_q_area                            
c                                                                               
c             calculate stress intensity factors and T-stress                   
c                                                                               
            Ki  = zero                                                          
            T11 = zero                                                          
            T33 = zero                                                          
            T13 = zero                                                          
c                                                                               
c             KI (plane stress auxiliary fields)                                
c                                                                               
            if( j .eq. 1 ) Ki = di_value_i(j) * e_front / two                   
c                                                                               
c             KI (plane strain auxiliary fields)                                
c                                                                               
            if( j .eq. 2 ) Ki = di_value_i(j) * e_front                         
     &                        / ( two * (one - nu_front**2) )                   
c                                                                               
c             KII (plane stress auxiliary fields)                               
c                                                                               
            if( j .eq. 3 ) Ki = di_value_i(j) * e_front / two                   
c                                                                               
c             KII (plane strain auxiliary fields)                               
c                                                                               
            if( j .eq. 4 ) Ki = di_value_i(j) * e_front                         
     &                        / ( two * (one - nu_front**2) )                   
c                                                                               
c             KIII (anti-plane shear auxiliary fields)                          
c                                                                               
            if( j .eq. 5 ) Ki = di_value_i(j) * e_front                         
     &                        / (two * (one + nu_front))                        
c                                                                               
c             T11, T33 (plane stress auxiliary fields)                          
c                                                                               
            if( j .eq. 6 ) then                                                 
               T11 = e_front * di_value_i(j)                                    
               Ki  = T11                                                        
               T33 = zero                                                       
               domain_min_i(9) = T33                                            
               domain_max_i(9) = T33                                            
               domain_avg_i(9) = T33                                            
            end if                                                              
c                                                                               
c             T11, T33 (plane strain auxiliary fields)                          
c                                                                               
            if( j .eq. 7 ) then                                                 
               T11 = e_front / (one - nu_front**2 )                             
     &             * ( di_value_i(j) + nu_front * e33_front )                   
               T33 = e33_front * e_front + nu_front * T11                       
               Ki  = T11                                                        
               domain_min_i(10) = min( domain_min_i(10), T33 )                  
               domain_max_i(10) = max( domain_max_i(10), T33 )                  
               domain_avg_i(10) = domain_avg_i(10) + T33                        
            end if                                                              
c                                                                               
c             T13 (anti-plane shear auxiliary fields)                           
c                                                                               
            if( j .eq. 8 ) then                                                 
               T13 = di_value_i(j) * e_front / (two * (one + nu_front ))        
               Ki  = T13                                                        
            end if                                                              
c                                                                               
c             determine maximum and minimum values                              
c                                                                               
            domain_min_i(j) = min( domain_min_i(j), Ki )                        
            domain_max_i(j) = max( domain_max_i(j), Ki )                        
            domain_avg_i(j) = domain_avg_i(j) + Ki                              
c                                                                               
c             store i-integral values for output in didriv.f                    
c                                                                               
        i_storage(ring_count,1,j)   = dble( nowring )                           
            i_storage(ring_count,2:9,j) = iiterms(1:8,j)                        
            i_storage(ring_count,10,j)  = di_value_i(j)                         
            i_storage(ring_count,11,j)  = Ki                                    
            if( j.eq.7 ) i_storage(ring_count,12,j) = T33                       
        i_storage(ring_count,13,j)  = dble(skipped_killed)                      
         end do                                                                 
      end if                                                                    
c                                                                               
c             compute J-values from K-values...                                 
c                                                                               
c             including KI plane stress, KII plane stress, and KIII             
c                                                                               
      j_from_ks(ring_count,1) = i_storage(ring_count,11,1)**2 / e_front         
     &                        + i_storage(ring_count,11,3)**2 / e_front         
     &                        + i_storage(ring_count,11,5)**2                   
     &                        * (one + nu_front ) / e_front                     
c                                                                               
c             including KI plane strain, KII plane strain, and KIII             
c                                                                               
      j_from_ks(ring_count,2) = i_storage(ring_count,11,2)**2                   
     &                          * (one - nu_front**2 ) / e_front                
     &                        + i_storage(ring_count,11,4)**2                   
     &                          * (one - nu_front**2 ) / e_front                
     &                        + i_storage(ring_count,11,5)**2                   
     &                        * (one + nu_front ) / e_front                     
c                                                                               
c             print sum of values from all elements in domain                   
c                                                                               
      if ( print_elem_values ) then                                             
         write(iout,9160)                                                       
         if( comput_j ) write(iout,9170) "J     ", (diterms(i),i=1,8),          
     &        jtotal                                                            
         if( comput_i ) then                                                    
            write(iout,9170) "I_KI  ", (iiterms(i,1),i=1,8), itotal(1)          
            write(iout,9170) "I_KI  ", (iiterms(i,2),i=1,8), itotal(2)          
            write(iout,9170) "I_KII ", (iiterms(i,3),i=1,8), itotal(3)          
            write(iout,9170) "I_KII ", (iiterms(i,4),i=1,8), itotal(4)          
            write(iout,9170) "I_KIII", (iiterms(i,5),i=1,8), itotal(5)          
            write(iout,9170) "I_T11 ", (iiterms(i,6),i=1,8), itotal(6)          
            write(iout,9170) "I_T11 ", (iiterms(i,7),i=1,8), itotal(7)          
            write(iout,9170) "I_T13 ", (iiterms(i,8),i=1,8), itotal(8)          
         end if                                                                 
      end if                                                                    
c                                                                               
c             print sum of values from all elements in domain,                  
c             and J, I for domain                                               
c                                                                               
      if( print_totals ) then                                                   
         write(iout,9180)                                                       
         if( nowring.eq.0 ) then                                                
            if( comput_j .and. comput_i ) write(iout,9185)                      
            if( comput_j .and. .not. comput_i ) write(iout,9186)                
            if( comput_i .and. .not. comput_j ) write(iout,9187)                
            if( comput_j ) then                                                 
               write(iout,9220) "J     ", domain_id(1:8),                       
     &              (diterms(j),j=1,8), di_value_j, skipped_killed              
            end if                                                              
            if( comput_i ) then                                                 
               write(iout,9220) "I_KI  ", domain_id(1:8),                       
     &              (iiterms(i,1),i=1,8), di_value_i(1), skipped_killed         
               write(iout,9220) "I_KI  ", domain_id(1:8),                       
     &              (iiterms(i,2),i=1,8), di_value_i(2), skipped_killed         
               write(iout,9220) "I_KII ", domain_id(1:8),                       
     &              (iiterms(i,3),i=1,8), di_value_i(3), skipped_killed         
               write(iout,9220) "I_KII ", domain_id(1:8),                       
     &              (iiterms(i,4),i=1,8), di_value_i(4), skipped_killed         
               write(iout,9220) "I_KIII", domain_id(1:8),                       
     &              (iiterms(i,5),i=1,8), di_value_i(5), skipped_killed         
               write(iout,9220) "I_T11 ", domain_id(1:8),                       
     &              (iiterms(i,6),i=1,8), di_value_i(6), skipped_killed         
               write(iout,9220) "I_T11 ", domain_id(1:8),                       
     &              (iiterms(i,7),i=1,8), di_value_i(7), skipped_killed         
               write(iout,9220) "I_T13 ", domain_id(1:8),                       
     &              (iiterms(i,8),i=1,8), di_value_i(8), skipped_killed         
            end if                                                              
         end if                                                                 
c                                                                               
         if( nowring.gt.0 ) then                                                
            if( comput_j .and.       comput_i ) write(iout,9250)                
            if( comput_j .and. .not. comput_i ) write(iout,9270)                
            if( comput_i .and. .not. comput_j ) write(iout,9290)                
            if( comput_j ) then                                                 
               write(iout,9300) "J     ", nowring, (diterms(i),i=1,8),          
     &              di_value_j, skipped_killed                                  
            end if                                                              
            if( comput_i ) then                                                 
               write(iout,9300) "I_KI  ", nowring,                              
     &              (iiterms(i,1),i=1,8), di_value_i(1), skipped_killed         
               write(iout,9300) "I_KI  ", nowring,                              
     &              (iiterms(i,2),i=1,8), di_value_i(2), skipped_killed         
               write(iout,9300) "I_KII ", nowring,                              
     &              (iiterms(i,3),i=1,8), di_value_i(3), skipped_killed         
               write(iout,9300) "I_KII ", nowring,                              
     &              (iiterms(i,4),i=1,8), di_value_i(4), skipped_killed         
               write(iout,9300) "I_KIII", nowring,                              
     &              (iiterms(i,5),i=1,8), di_value_i(5), skipped_killed         
               write(iout,9300) "I_T11 ", nowring,                              
     &              (iiterms(i,6),i=1,8), di_value_i(6), skipped_killed         
               write(iout,9300) "I_T11 ", nowring,                              
     &              (iiterms(i,7),i=1,8), di_value_i(7), skipped_killed         
               write(iout,9300) "I_T13 ", nowring,                              
     &              (iiterms(i,8),i=1,8), di_value_i(8), skipped_killed         
            end if                                                              
         end if                                                                 
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9100 format(/,1x,25x,'domain integral components',                             
     &  /,      1x,25x,'--------------------------')                            
 9110 format(7x,'element',5x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                  
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8' )                                 
 9130 format(/,5x,' gauss point strains for element ',i7,', gpt ',i2,           
     &       ':',/, 3(10x,3(e11.4,2x),/))                                       
 9140 format(1x,a,1x,i7,8(1x,e11.4))                                            
 9150 format('')                                                                
 9160 format(/,1x,'element totals:',3x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',        
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total')                       
 9170 format(1x,a,7x,9(1x,e11.4))                                               
 9180 format(/,1x,25x,'domain integral components',                             
     &  /,      1x,25x,'--------------------------')                            
 9185 format(8x,'domain',8x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                   
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',5x,'total J,I',                   
     &    2x,'killed ele' )                                                     
 9186 format(8x,'domain',8x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                   
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total J',                     
     &    2x,'killed ele' )                                                     
 9187 format(8x,'domain',8x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                   
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total I',                     
     &    2x,'killed ele' )                                                     
 9220 format(1x,a,1x,a8,9(1x,e11.4),2x,'(',i3,')')                              
 9250 format(8x,'domain',6x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                   
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',5x,'total J,I',                   
     &    2x,'killed ele' )                                                     
 9270 format(8x,'domain',6x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                   
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total J',                     
     &    2x,'killed ele' )                                                     
 9290 format(8x,'domain',6x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',                   
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total I',                     
     &    2x,'killed ele' )                                                     
 9300 format(1x,a,1x,i7,9(1x,e11.4),2x,'(',i3,')')                              
 9400 format(/,1x,'area under q-function along crack front:  ',e11.4,           
     &       /,1x,'length along crack front for this domain: ',e11.4,/)         
c                                                                               
      end                                                                       
c                                                                               
c***************************************************************                
c                                                              *                
c subroutine to transform a 3x1 vector of constraint-          *                
c compatible-coordinate-system values to global-coordinate-    *                
c system values. the routine premultiplies the incoming        *                
c constraint-compatible values by the transpose of the stored  *                
c global-to-constraint coordinate system rotation matrix.      *                
c                                                              *                
c                                 written by: mcw              *                
c                              last modified: 02/01            *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_trans_nodvals(vec, transmat)                                
      implicit none                                                             
c                                                                               
c     dummy variables                                                           
c                                                                               
      double precision                                                          
     & vec(3), transmat(3,3)                                                    
c                                                                               
c     local variables                                                           
c                                                                               
      double precision                                                          
     & tempvec(3)                                                               
c                                                                               
c                                                                               
      tempvec(1) = transmat(1,1) * vec(1)                                       
     &           + transmat(2,1) * vec(2)                                       
     &           + transmat(3,1) * vec(3)                                       
      tempvec(2) = transmat(1,2) * vec(1)                                       
     &           + transmat(2,2) * vec(2)                                       
     &           + transmat(3,2) * vec(3)                                       
      tempvec(3) = transmat(1,3) * vec(1)                                       
     &           + transmat(2,3) * vec(2)                                       
     &           + transmat(3,3) * vec(3)                                       
c                                                                               
      vec(1) = tempvec(1)                                                       
      vec(2) = tempvec(2)                                                       
      vec(3) = tempvec(3)                                                       
c                                                                               
      return                                                                    
      end                                                                       
