c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine initdm                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/7/10 by rhd              *          
c     *                                   add tangent vector         *          
c     *                                                              *          
c     *     this subroutine initializes various variables and arrays *          
c     *     that define a domain for j-integral computations         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine initdm                                                         
      use j_data                                                                
      implicit integer (a-z)                                                    
      double precision                                                          
     &   zero, one                                                              
      real rword, rzero                                                         
      equivalence ( iword, rword )                                              
      data zero, rzero, one / 0.d0, 0.0, 1.d0 /                                 
c                                                                               
      if ( allocated( compr_q_list ) ) deallocate( compr_q_list )               
      if ( allocated( q_element_maps ) ) deallocate( q_element_maps )           
      if ( allocated( node_set ) ) deallocate( node_set )                       
c                                                                               
      crack_plane_normal(1)  = zero                                             
      crack_plane_normal(2)  = zero                                             
      crack_plane_normal(3)  = zero                                             
      crack_front_tangent(1)  = zero                                            
      crack_front_tangent(2)  = zero                                            
      crack_front_tangent(3)  = zero                                            
      tangent_vector_defined = .false.                                          
      symmetric_domain       = .false.                                          
      one_point_rule         = .false.                                          
      num_front_nodes        = 0                                                
      front_nodes(1:30)      = 0                                                
      verify_front           = .false.                                          
      front_order            = 0                                                
      print_elem_values      = .false.                                          
      print_totals           = .false.                                          
      qvals_given            = .false.                                          
      rings_given            = .false.                                          
      q_vals_linear          = .true.                                           
      domain_type            = 0                                                
      debug_driver           = .false.                                          
      debug_elements         = .false.                                          
      domain_id(1:24)        = ' '                                              
      last_compr             = 0                                                
      ring_list(1:300)       = 0                                                
      ignore_face_loads      = .false.                                          
      omit_crack_front_elems = .false.                                          
      max_node_set_id        = 0                                                
      num_auto_rings         = 0                                                
      output_packet_j        = .false.                                          
      comput_j               = .false.                                          
      comput_i               = .false.                                          
      ring_count             = 0                                                
      box_tol_relat      = 0.0003d0                                             
      out_pstress            = .true.                                           
      out_pstrain            = .true.                                           
      output_packet_i        = .false.                                          
      cf_traction_flags(1:3) = .false.                                          
      cf_tractions(1:3)      = zero                                             
      face_loading           = .false.                                          
      e33_front              = zero                                             
      j_to_k                 = .false.                                          
      crack_curvature(1:7)   = zero                                             
c                                                                               
      return                                                                    
      end                                                                       
