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
c***********************************************************************        
c                                                                      *        
c    subroutine to drive computation of surface-traction               *        
c    integrals                                                         *        
c                                                                      *        
c                                 written by: rhd/mcw                  *        
c                                   modified: mcw                      *        
c                              last modified: 9/18/03                  *        
c                                                                      *        
c***********************************************************************        
c                                                                               
      subroutine di_calc_surface_integrals( elemno, etype, nnode,               
     &     snodes, feload, cf_traction_flags, cf_tractions, rotate, dsf,        
     &     jacobi, cdispl, qvals, coord, front_nodes, front_coords,             
     &     domain_type, domain_origin, num_front_nodes, front_order,            
     &     ym_front_node, nu_front_node, comput_j, comput_i, jterm,             
     &     iterm, front_elem_flag, qp_node, crack_curvature,                    
     &     face_loading, out, debug )                                           
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer elemno, etype, nnode, snodes(*), front_nodes(*),                  
     &        domain_type, domain_origin, num_front_nodes, front_order,         
     &        out                                                               
      double precision                                                          
     &     feload(3,*), cf_tractions(*), rotate(3,3), dsf(20,3,28),             
     &     jacobi(3,3,28), cdispl(3,20), qvals(*), jterm(9),                    
     &     iterm(8,7), coord(3,*), front_coords(3,*), ym_front_node,            
     &     nu_front_node, crack_curvature(*)                                    
      logical cf_traction_flags(*), comput_j, comput_i, front_elem_flag,        
     &     qp_node, face_loading, debug                                         
c                                                                               
c             local variables                                                   
c                                                                               
      integer j, enode, faceno, nfnode, fnodes(20), sfnodes(20), flag           
      double precision                                                          
     &     sum_load, eloads(3,20), zero                                         
c                                                                               
      data zero                                                                 
     & / 0.d0 /                                                                 
c                                                                               
      if( debug ) write(out,100) elemno                                         
c                                                                               
      sum_load = zero                                                           
      do enode = 1, nnode                                                       
        sum_load = sum_load + abs( feload(1,enode) ) +                          
     &                        abs( feload(2,enode) ) +                          
     &                        abs( feload(3,enode) )                            
      end do                                                                    
      if ( sum_load .le. 1.0e-8 ) return                                        
c                                                                               
      if( debug ) write(out,110) elemno                                         
      do enode = 1, nnode                                                       
        if( debug ) then                                                        
           write(out,120) snodes(enode), feload(1,enode),                       
     &                    feload(2,enode), feload(3,enode)                      
        end if                                                                  
      end do                                                                    
      if( debug ) write(out, 130) sum_load                                      
c                                                                               
c             determine which face is loaded. only one loaded face              
c             is permitted. rotate equivalent loads from element                
c             local coordinates to crack coordinates. if all nodes              
c             of element have non-zero equiv loads, then loading                
c             is a body force or a thermal loading. thermals are                
c             handled above. We ignore body force loadings. also                
c             rotate crack-face tractions input in domain definition            
c             to local crack coordinates                                        
c                                                                               
      if ( debug ) write(out,140)                                               
      call dielwf( feload, etype, nnode, faceno, out, debug,                    
     &             flag, snodes, nfnode, fnodes, sfnodes )                      
      if ( flag .eq. 1 .or. faceno .eq. 0 ) return                              
      face_loading = .true.                                                     
      if ( debug ) write(out,150) faceno, (fnodes(j),j=1,nfnode)                
      call dielrl( feload, eloads, cf_tractions, rotate, nnode, debug,          
     &             out )                                                        
c                                                                               
c             set-up is complete. calculate integrals.                          
c                                                                               
      if( comput_j ) then                                                       
         call di_calc_surface_j( etype, nnode, nfnode, fnodes,                  
     &                        dsf, jacobi, cdispl, eloads, qvals,               
     &                        jterm, elemno, out, debug )                       
      end if                                                                    
c                                                                               
      if( comput_i ) then                                                       
         call di_calc_surface_i( elemno, etype, nnode, nfnode, fnodes,          
     &                        sfnodes, coord, front_nodes, front_coords,        
     &                        domain_type, domain_origin,                       
     &                        num_front_nodes, front_order,                     
     &                        cf_traction_flags, cf_tractions,                  
     &                        cdispl, eloads, qvals, iterm,                     
     &                        ym_front_node, nu_front_node,                     
     &                        front_elem_flag, qp_node, crack_curvature,        
     &                        out, debug )                                      
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 100  format(//,'---begin evaluation of surface-traction integral:',            
     &       ' element ',i7)                                                    
 110  format(//,'---element ',i7,' has nodal loads',                            
     &       /,1x,'snode',9x,'fx',13x,'fy',13x,'fz')                            
 120  format(i8,3(2x,e13.6))                                                    
 130  format('sum of nodal loads',2x,e13.6)                                     
 140  format(//,' >>> begin face-traction integration: ' )                      
 150  format(/,' >>>>> face: ',i2,' face nodes: ',8i3)                          
c                                                                               
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *    determine which face of an element is loaded. if face is     *           
c *    collapsed, ignore it.                                        *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine dielwf( eqload, etype, nnode, faceno, iout, debug,             
     &                   flag, snodes, nfnode, fnodes, sfnodes )                
      implicit double precision (a-h,o-z)                                       
c                                                                               
c             parameter declarations                                            
c                                                                               
      integer etype, faceno, flag, snodes(*), nfnode, fnodes(*),                
     &        sfnodes(*)                                                        
      dimension eqload(3,*)                                                     
      logical debug                                                             
c                                                                               
c             local declarations                                                
c                                                                               
      integer fnode, enode, face, counter                                       
      logical ldface(6)                                                         
      data toler / 1.0d-10 /                                                    
c                                                                               
c             by examining the non-zero entries in the equivalent               
c             loads for the element, determine which face has an                
c             applied surface traction. only one loaded face will               
c             be processed. routine eqelfn returns the nodes on a               
c             face for any element.                                             
c                                                                               
      if ( debug ) write(iout,9000)                                             
c                                                                               
c             check first for loads on all element nodes indicating             
c             a body force type loading.  skip face loads                       
c             processing if we have this case.                                  
c                                                                               
      flag = 0                                                                  
      do enode = 1, nnode                                                       
          if ( (abs(eqload(1,enode)) + abs(eqload(2,enode)) +                   
     &          abs(eqload(3,enode))) .le. toler ) go to 80                     
      end do                                                                    
      flag = 1                                                                  
      return                                                                    
c                                                                               
c             at least one node has no loading applied.                         
c             for a face to be considered loaded with a surface                 
c             traction, every node on the face must have a non-zero             
c             component of equivalent nodal load.                               
c                                                                               
 80   continue                                                                  
      numfac = 6                                                                
      do 100 face = 1, numfac                                                   
        ldface(face) = .false.                                                  
        call eqelfn( fnodes, etype, face, nfnode )                              
c                                                                               
c             obtain the structure node numbers for the face being              
c             considered.                                                       
c                                                                               
        do i= 1, nfnode                                                         
           fnode      = fnodes(i)                                               
           sfnodes(i) = snodes(fnode)                                           
        end do                                                                  
c                                                                               
        if( debug ) write(iout, 9500) (sfnodes(i),i=1,nfnode)                   
c                                                                               
c             if any nodes on the face are repeated, consider it a              
c             collapsed face and cycle the "do 100" loop.                       
c                                                                               
        do i = 1, nfnode                                                        
           counter = i                                                          
           do 90 j = 1, nfnode                                                  
              if( j .eq. counter) go to 90                                      
              if( sfnodes(i) .eq. sfnodes(j)) go to 100                         
 90        continue                                                             
        end do                                                                  
c                                                                               
c             if any node has a load smaller than toler, consider the           
c             face traction-free and cycle do loop.                             
c                                                                               
        do i = 1, nfnode                                                        
          fnode = fnodes(i)                                                     
          if ( (abs(eqload(1,fnode)) + abs(eqload(2,fnode)) +                   
     &          abs(eqload(3,fnode))) .le. toler ) go to 100                    
        end do                                                                  
        ldface(face) = .true.                                                   
 100  continue                                                                  
c                                                                               
      if ( debug ) write(iout,9010) (ldface(i),i=1,numfac)                      
c                                                                               
c             if more than one face is loaded, send                             
c             message to user. we only process the first loaded                 
c             face.                                                             
c                                                                               
      faceno = 0                                                                
      do face = 1, numfac                                                       
         if ( ldface(face) ) then                                               
            if ( faceno .gt. 0 ) then                                           
               call dieler( iout, ierr, 7 )                                     
               return                                                           
            else                                                                
               faceno = face                                                    
               call eqelfn( fnodes, etype, faceno, nfnode )                     
c                                                                               
c             obtain the structure node numbers for the loaded face.            
c                                                                               
               do i= 1, nfnode                                                  
                  fnode      = fnodes(i)                                        
                  sfnodes(i) = snodes(fnode)                                    
               end do                                                           
            end if                                                              
         end if                                                                 
      end do                                                                    
c                                                                               
      if ( debug ) then                                                         
         write(iout,9030) faceno                                                
         do i=1,nfnode                                                          
            write(iout,9040) fnodes(i), sfnodes(i)                              
         end do                                                                 
      end if                                                                    
                                                                                
      return                                                                    
c                                                                               
 9000 format( ' >>> entered find loaded face' )                                 
 9010 format( //,'  >> loaded face flags: ',6l1 )                               
 9020 format( ' >>>>> more than one element face has applied surface',          
     & /,     '       traction. lowest numbered face processed' )               
 9030 format( '  >> loaded face number: ',i4,// )                               
 9040 format(' node, structure node: ',i4,2x,i8 )                               
 9500 format(//,'snodes for face:',8(2x,i5))                                    
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *   rotate equivalent nodal loads to crack coordinate system      *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine dielrl( feload, eloads, cf_tractions, rotate, nnode,           
     &                   debug, iout )                                          
      implicit double precision (a-h,o-z)                                       
c                                                                               
c             rotate equivalent nodal loads from element local                  
c             to the crack coordinate system.                                   
c                                                                               
      dimension feload(3,*), eloads(3,*), cf_tractions(3), rotate(3,3)          
      logical debug                                                             
c                                                                               
c             local declarations                                                
c                                                                               
c                                                                               
      do inode = 1,nnode                                                        
        eloads(1,inode) = rotate(1,1) * feload(1,inode) +                       
     &                    rotate(1,2) * feload(2,inode) +                       
     &                    rotate(1,3) * feload(3,inode)                         
        eloads(2,inode) = rotate(2,1) * feload(1,inode) +                       
     &                    rotate(2,2) * feload(2,inode) +                       
     &                    rotate(2,3) * feload(3,inode)                         
        eloads(3,inode) = rotate(3,1) * feload(1,inode) +                       
     &                    rotate(3,2) * feload(2,inode) +                       
     &                    rotate(3,3) * feload(3,inode)                         
      end do                                                                    
      tx =   rotate(1,1) * cf_tractions(1)                                      
     &     + rotate(1,2) * cf_tractions(2)                                      
     &     + rotate(1,3) * cf_tractions(3)                                      
      ty =   rotate(2,1) * cf_tractions(1)                                      
     &     + rotate(2,2) * cf_tractions(2)                                      
     &     + rotate(2,3) * cf_tractions(3)                                      
      tz =   rotate(3,1) * cf_tractions(1)                                      
     &     + rotate(3,2) * cf_tractions(2)                                      
     &     + rotate(3,3) * cf_tractions(3)                                      
      cf_tractions(1) = tx                                                      
      cf_tractions(2) = ty                                                      
      cf_tractions(3) = tz                                                      
      if ( debug ) write(iout,9002) (inode,(eloads(i,inode),                    
     &             i=1,3),inode=1,nnode)                                        
      if( debug ) write(iout,9004) (cf_tractions(i),i=1,3)                      
c                                                                               
      return                                                                    
c                                                                               
 9002 format (//,'>>>>>  nodal loads crack reference frame',//,                 
     & 1x,'node',9x,'x',15x,'y',15x,'z',27(/,i5,3(2x,e13.6)))                   
 9004 format ('>>>>>  crack-face tractions in crack reference frame',//,        
     & 1x,'tx, ty, tz: ',3(2x,e13.6))                                           
      end                                                                       
c                                                                               
c                                                                               
c***************************************************************                
c                                                              *                
c    subroutine to calculate surface-traction integral for     *                
c    j-integral term5                                          *                
c                                                              *                
c                                 written by: rhd              *                
c                                   modified: mcw 9/18/03      *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_calc_surface_j( etype, nnode, nfnode, fnodes,               
     &                           dsf, jacobi, cdispl, eloads, qvals,            
     &                           jterm, elemno, out, debug )                    
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer etype, nnode, nfnode, fnodes(*), elemno, out                      
      double precision                                                          
     &     dsf(20,3,28), jacobi(3,3,28), cdispl(3,20), eloads(3,20),            
     &     qvals(*), jterm(9)                                                   
      logical debug                                                             
c                                                                               
c             local variables                                                   
c                                                                               
      integer ptno, enode, k1, idumvec(1), order_face, ngpts_face               
      double precision                                                          
     &     xsi, eta, zeta, weight, dux, dvx, dwx, nx, ny, nz, dieldp,           
     &     cderiv(27,3), termu, termv, termw, dum_vec(1), lg(28), zero          
c                                                                               
      data zero                                                                 
     & / 0.d0 /                                                                 
c                                                                               
c             set integration order for face integrals for                      
c             tractions applied on crack face. we always use                    
c             a 2x2x2 order for the face loadings.                              
c                                                                               
      order_face = 1                                                            
      ngpts_face = 8                                                            
      if ( etype .eq. 1 ) order_face = 8                                        
c                                                                               
c             loop over all gauss points and compute                            
c             cartesian displacement derivatives: dux, dvx, dwx                 
c                                                                               
c             calculate ui,1 at gauss points on the interior, not the           
c             loaded face, of the brick element, based on displacements         
c             from all element nodes.                                           
c                                                                               
      if ( debug ) write(out,9105)                                              
      if ( debug ) write(out,9106) etype, order_face                            
c                                                                               
      do ptno = 1, ngpts_face                                                   
c                                                                               
c             isoparametric coordinates of gauss point and weight.              
c                                                                               
         call getgpts( etype, order_face, ptno, xsi, eta, zeta, weight )        
         if ( debug ) write(out,9110) ptno, xsi, eta, zeta                      
c                                                                               
         dux     = zero                                                         
         dvx     = zero                                                         
         dwx     = zero                                                         
         do enode = 1, nnode                                                    
            nx = dsf(enode,1,ptno) * jacobi(1,1,ptno) +                         
     &           dsf(enode,2,ptno) * jacobi(1,2,ptno) +                         
     &           dsf(enode,3,ptno) * jacobi(1,3,ptno)                           
            ny = dsf(enode,1,ptno) * jacobi(2,1,ptno) +                         
     &           dsf(enode,2,ptno) * jacobi(2,2,ptno) +                         
     &           dsf(enode,3,ptno) * jacobi(2,3,ptno)                           
            nz = dsf(enode,1,ptno) * jacobi(3,1,ptno) +                         
     &           dsf(enode,2,ptno) * jacobi(3,2,ptno) +                         
     &           dsf(enode,3,ptno) * jacobi(3,3,ptno)                           
            dux = dux + nx * cdispl(1,enode)                                    
            dvx = dvx + nx * cdispl(2,enode)                                    
            dwx = dwx + nx * cdispl(3,enode)                                    
         end do                                                                 
         cderiv(ptno,1) = dux                                                   
         cderiv(ptno,2) = dvx                                                   
         cderiv(ptno,3) = dwx                                                   
      end do                                                                    
      if( debug ) write(out,9210)                                               
      do ptno = 1, ngpts_face                                                   
         if ( debug ) write(out,9211) ptno, cderiv(ptno,1),                     
     &        cderiv(ptno,2), cderiv(ptno,3)                                    
      end do                                                                    
c                                                                               
c                 1)  loop over all nodes on the face.                          
c                 2)  get parametric coordinates of the node.                   
c                 3)  compute lagrange interpolation functions                  
c                     from gauss points to the node                             
c                 2)  extrapolate cartesian derivatives to node                 
c                 3)  add each nodes' contribution for j to the total.          
c                                                                               
      termu = zero                                                              
      termv = zero                                                              
      termw = zero                                                              
c                                                                               
      do k1 = 1, nfnode                                                         
         enode = fnodes(k1)                                                     
         call ndpts1( idumvec, 0, dum_vec, etype, enode, xsi, eta, zeta)        
         if ( etype .eq. 1 ) then                                               
            call oulgr1( xsi, eta, zeta, lg, order_face )                       
         else                                                                   
            call oulgr2( xsi, eta, zeta, lg, order_face )                       
         end if                                                                 
         if ( debug ) write(out,9112) enode, xsi, eta, zeta                     
c                                                                               
         dux = dieldp( lg, cderiv(1,1), ngpts_face, 1, 1 )                      
         dvx = dieldp( lg, cderiv(1,2), ngpts_face, 1, 1 )                      
         dwx = dieldp( lg, cderiv(1,3), ngpts_face, 1, 1 )                      
c                                                                               
c             surface traction term is given by (-) product of nodal            
c             q function values, the equivalent nodal loads,                    
c             and the cartesian crack face derivatives.                         
c                                                                               
         termu = termu - eloads(1,enode) * qvals(enode) * dux                   
         termv = termv - eloads(2,enode) * qvals(enode) * dvx                   
         termw = termw - eloads(3,enode) * qvals(enode) * dwx                   
         if ( debug ) write(out,9020) enode, dux, dvx, dwx,                     
     &                                qvals(enode)                              
      end do                                                                    
      if ( debug ) write(out,9022) termu, termv, termw                          
      jterm(5) = termu + termv + termw                                          
c                                                                               
      return                                                                    
c                                                                               
 9105 format(/,' >>>>> start of integration loop for dux,dvx,dwx',              
     &       ' of jterm(5):' )                                                  
 9106 format(/,'     >> etype',2x,i2,2x,'order_face',2x,i2)                     
 9110 format(/,'     >> ptno, xsi, eta, zeta: ',i2,3f10.4)                      
 9210 format(/,'     >> dux, dvx, dwx at gauss points:')                        
 9211 format(10x,i4,3f10.6)                                                     
 9112 format(/,'     >> enode, xsi, eta, zeta: ',i2,3f10.4)                     
 9020 format('   > at face node: ', i3,                                         
     & /,    '       dux, dvx, dwx, q: ',4(1x,e13.6) )                          
 9022 format('  >> end face integrations. x, y , z terms: ',                    
     & /,    '     ',e13.6,1x,e13.6,1x,e13.6 )                                  
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c***********************************************************************        
c                                                                      *        
c    subroutine to calculate surface-traction integral for             *        
c    i-integral, iterm(8,j).                                           *        
c                                                                      *        
c                                 written by: mcw                      *        
c                              last modified: 9/18/03                  *        
c                                                                      *        
c***********************************************************************        
c                                                                               
      subroutine di_calc_surface_i( elemno, etype, nnode, nfnode,               
     &                           fnodes, sfnodes, coord, front_nodes,           
     &                           front_coords, domain_type,                     
     &                           domain_origin, num_front_nodes,                
     &                           front_order, cf_traction_flags,                
     &                           cf_tractions, cdispl, eloads, qvals,           
     &                           iterm, e_front, nu_front,                      
     &                           front_elem_flag, qp_node,                      
     &                           crack_curvature, out, debug )                  
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer elemno, etype, nnode, nfnode, fnodes(20), sfnodes(20),            
     &        domain_type, front_nodes(*), domain_origin,                       
     &        num_front_nodes, front_order, out                                 
      double precision                                                          
     &     coord(3,*), front_coords(3,*), cf_tractions(*), cdispl(3,20),        
     &     eloads(3,20), qvals(*), iterm(8,7), e_front, nu_front,               
     &     crack_curvature(*)                                                   
      logical cf_traction_flags(*), front_elem_flag, qp_node, debug             
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, j, k, enode, order_face, ngpts_face, ptno, fnode,              
     &        idumvec(1), ierr                                                  
      double precision                                                          
     &     fnode_coords(3,20), xsi, eta, zeta,                                  
     &     weight, dux, dvx, dwx, nx, ny, nz, cderiv_aux(27,3,7),               
     &     termu, termv, termw, dum_vec(1), lg(28), sf(20),                     
     &     dsf(20,3,28), jacob(2,2), det_j, point_x, point_y,                   
     &     point_z, gpq, r, t, du11_aux(7), du21_aux(7),                        
     &     du31_aux(7), ave_y_coord, dieldp, mu_front, kappa,                   
     &     zero, half, one, two, three, pi, four, eta_old, eta_shifted,         
     &     sum, ratio                                                           
      logical approx, shift                                                     
c                                                                               
      data zero, half, one, two, three, pi, four                                
     & / 0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0, 3.14159265359d0, 4.0d0 /            
c                                                                               
c                                                                               
c             method 1: when crack-face tractions are known only                
c                       by equivalent nodal loads, we use the same              
c                       approximate integration scheme used for crack-face      
c                       contributions to the j-integral. this method            
c                       employs nodal loads instead of tractions.               
c                                                                               
c             method 2: when the user inputs uniform crack-face tractions       
c                       in the domain definition, we use standard               
c                       gauss quadrature to evaluate the integral over the      
c                       crack face area.                                        
c                                                                               
c             evaluate iterm8 using method 1. if user has not defined           
c             crack-face tractions in the domain integral definition,           
c             the result will be zero. then, use method 2. this involves        
c             a little extra work.                                              
c                                                                               
c             define element types, number of gauss points, and face order.     
c             we consider planar quadrilaterals with either 2x2 or 3x3 integrati
c                                                                               
c             order_face  2 = 2x2                                               
c             order_face  3 = 3x3 etc.                                          
c             etype  1      = 20-noded brick element                            
c             etype  2      =  8-noded brick element                            
c             etype  3      = 12-noded brick element                            
c             etype  4      = 15-noded brick element                            
c             etype  5      =  9-noded brick element                            
c             etype  7      = 15-noded wedge element                            
c             etype  9      =  8-node planar quadrilateral "quad8" element      
c             etype 16      =  4-node planar quadrilateral "quad4" element      
c                                                                               
c                                                                               
c                                                                               
c             if user has input tractions in the domain definition, use the     
c             exact integration scheme. if not, employ approximate scheme       
c             used in calc_surface_j.                                           
c                                                                               
      approx = .true.                                                           
      if( cf_traction_flags(1) ) approx = .false.                               
      if( cf_traction_flags(2) ) approx = .false.                               
      if( cf_traction_flags(3) ) approx = .false.                               
c                                                                               
c             set element type--either 4 or 8 nodes on face                     
c                                                                               
      if( nfnode .eq. 4 ) etype = 16                                            
      if( nfnode .eq. 8 ) etype =  9                                            
c                                                                               
c             set integration order.                                            
c             for approximate integration, use 2x2 quadrature on face.          
c             for exact scheme, 2x2 quadrature is adequate for elements         
c             with straight edges, but we use 4x4 to habdle elements            
c             with curved edges.                                                
c                                                                               
      if( approx ) order_face = 2                                               
      if( .not. approx ) order_face = 4                                         
      if( order_face .eq.  2 ) ngpts_face =    4                                
      if( order_face .eq.  3 ) ngpts_face =    9                                
      if( order_face .eq.  4 ) ngpts_face =   16                                
c                                                                               
c             obtain local coordinates of nodes on loaded face                  
c                                                                               
      do i=1, nfnode                                                            
         enode = fnodes(i)                                                      
         fnode_coords(1,i) = coord(1,enode)                                     
         fnode_coords(2,i) = coord(2,enode)                                     
         fnode_coords(3,i) = coord(3,enode)                                     
      end do                                                                    
c                                                                               
      if( debug ) then                                                          
         write( out,9900 )                                                      
         do i=1, nfnode                                                         
            write( out,9905 ) fnodes(i), sfnodes(i),                            
     &         fnode_coords(1,i), fnode_coords(2,i), fnode_coords(3,i)          
         end do                                                                 
      end if                                                                    
c                                                                               
c             calculate the local y-coordinate of the element centroid          
c             to determine if loaded element face lies on top or                
c             bottom crack face. cf_tractions() reflect loading on the          
c             crack face corresponding to local y >= 0. for faces               
c             corresponding to the bottom surface, apply tractions in           
c             opposite sense to those applied on top surface.                   
c                                                                               
      ave_y_coord = zero                                                        
      do enode = 1, nnode                                                       
         ave_y_coord = ave_y_coord                                              
     &               + (coord(2,enode) - front_coords(2,domain_origin))         
      end do                                                                    
         ave_y_coord = ave_y_coord / dble(nnode)                                
c                                                                               
      if( ave_y_coord .lt. zero ) then                                          
         cf_tractions(1) = - cf_tractions(1)                                    
         cf_tractions(2) = - cf_tractions(2)                                    
         cf_tractions(3) = - cf_tractions(3)                                    
      end if                                                                    
      if( debug ) write(out,9960) (cf_tractions(i),i=1,3)                       
c                                                                               
c             reorder nodes on face such that positive eta = - X1. this is      
c             necessary for correct calculation of jacobian, and for            
c             treatment of inverse square-root singularity which shifts eta.    
c                                                                               
      call di_reorder_nodes( nfnode, fnodes, sfnodes, fnode_coords,             
     &                    ave_y_coord, elemno, front_nodes,                     
     &                    num_front_nodes, front_coords, domain_type,           
     &                    domain_origin, front_order, crack_curvature,          
     &                    out )                                                 
c                                                                               
c             loop over gauss points of loaded face to evaluate integral        
c                                                                               
c             1. obtain parametric coordinates of gauss point (gp).             
c                value of zeta should be zero.                                  
c                                                                               
      do ptno = 1, ngpts_face                                                   
         call getgpts( etype, order_face, ptno, xsi, eta, zeta, weight )        
         if ( debug ) then                                                      
            write(out,9909) ptno                                                
            write(out,9910) xsi, eta, zeta, weight                              
         end if                                                                 
c                                                                               
c             2. if current element borders the crack front, to integrate       
c                the inverse square-root singularity, we use a substitution     
c                of the form                                                    
c                   int[a,b]f(x)dx = int[0,(b-a)^0.5]2t*f(t^2 + a)dt            
c                this requires us to shift eta in order to conform to           
c                the new limits of integration. we then evaluate f(x) at        
c                x = t^2 + a. this requires shape functions and derivatives     
c                to use the shifted value of eta.                               
c                                                                               
         shift = .false.                                                        
         if( front_elem_flag.and. .not.approx .and. .not.qp_node ) then         
            shift = .true.                                                      
         end if                                                                 
         if( debug .and. shift ) write(out,9907)                                
         if( shift ) then                                                       
            eta_old     = eta                                                   
            eta_shifted = (eta + one) / two**half                               
            eta         = eta_shifted**two - one                                
         end if                                                                 
         if ( debug ) write(out,9908) eta_old, eta_shifted, eta                 
c                                                                               
c             3. evaluate element shape functions, shape function               
c                derivatives.                                                   
c                                                                               
         call shapef( etype, xsi, eta, zeta, sf )                               
         call derivs( etype, xsi, eta, zeta, dsf(1,1,1),                        
     &                dsf(1,2,1), dsf(1,3,1) )                                  
         if( debug ) then                                                       
            write(out,9911)                                                     
            do j=1, nfnode                                                      
               write(out,9912) j, sf(j), dsf(j,1,1), dsf(j,2,1)                 
            end do                                                              
         end if                                                                 
c                                                                               
c             4. evaluate coordinate jacobian at current gauss point.           
c                                                                               
         jacob(1,1) = dieldp(dsf(1,1,1),fnode_coords(3,1),nfnode,1,3)           
         jacob(1,2) = dieldp(dsf(1,1,1),fnode_coords(1,1),nfnode,1,3)           
         jacob(2,1) = dieldp(dsf(1,2,1),fnode_coords(3,1),nfnode,1,3)           
         jacob(2,2) = dieldp(dsf(1,2,1),fnode_coords(1,1),nfnode,1,3)           
c                                                                               
c             5. calculate determinant of coordinate jacobian at                
c                current gauss point. for element on top crack face,            
c                parent coordinate xsi corresponds to +local z, and             
c                for elements on bottom crack face, parent xsi                  
c                corresponds to -local z. after subroutine di_reorder_nodes,    
c                parent eta always corresponds to -local x. adjust the          
c                values of jacobian terms to give a positive determinant.       
c                                                                               
         if( ave_y_coord .gt. zero ) then                                       
            jacob(1,2) = -jacob(1,2)                                            
            jacob(2,2) = -jacob(2,2)                                            
         end if                                                                 
         if( ave_y_coord .lt. zero ) then                                       
            jacob(1,1) = -jacob(1,1)                                            
            jacob(1,2) = -jacob(1,2)                                            
            jacob(2,1) = -jacob(2,1)                                            
            jacob(2,2) = -jacob(2,2)                                            
         end if                                                                 
c                                                                               
         det_j = jacob(1,1) * jacob(2,2)                                        
     &         - jacob(1,2) * jacob(2,1)                                        
c                                                                               
         if( debug ) write(out,9913) ((jacob(j,i),i=1,2),j=1,2), det_j          
c                                                                               
         if ( det_j.le.zero ) then                                              
            call dieler( out,ierr,5 )                                           
         end if                                                                 
c                                                                               
c             6. calculate local coordinates of gp, and q-value at gp.          
c                for non-uniform crack-face traction, change this section.      
c                                                                               
         point_x = zero                                                         
         point_y = zero                                                         
         point_z = zero                                                         
         gpq     = zero                                                         
         do enode = 1, nfnode                                                   
            if( debug ) write(out,9965) enode, fnodes(enode),                   
     &                                  sfnodes(enode),                         
     &                                  qvals(fnodes(enode))                    
            point_x = point_x + sf(enode) * fnode_coords(1,enode)               
            point_y = point_y + sf(enode) * fnode_coords(2,enode)               
            point_z = point_z + sf(enode) * fnode_coords(3,enode)               
            gpq     = gpq     + sf(enode) * qvals(fnodes(enode))                
         end do                                                                 
c                                                                               
c             7. calculate r and theta at gp.                                   
c                                                                               
         r = zero                                                               
         t = zero                                                               
         call di_calc_r_theta( 2, front_nodes, num_front_nodes,                 
     &                      front_coords, domain_type, domain_origin,           
     &                      front_order, point_x, point_y, point_z,             
     &                      elemno, ptno, r, t, crack_curvature,                
     &                      debug, out )                                        
c                                                                               
c             8. adjust theta when point_y is close to zero in order to         
c                obtain a value of theta with the correct sign.                 
c                this assumes a planar crack in the local x1-x3 plane.          
c                                                                               
         if( ave_y_coord .gt. zero .and. t .lt. zero ) then                     
            t = - t                                                             
         end if                                                                 
         if( ave_y_coord .lt. zero .and. t .gt. zero ) then                     
            t = - t                                                             
         end if                                                                 
c                                                                               
         if ( debug ) write(out,9915) ptno, point_x, point_y,                   
     &                                point_z, gpq, ave_y_coord                 
c                                                                               
c             9. calculate ui,1_aux at gp for all loading modes, p-stress, p-str
c                these are the same expressions for ui,1 that are in subroutine 
c                di_calc_aux_fields_k, and di_calc_aux_fields_t, but are simplif
c                for theta = +/- pi.                                            
c                                                                               
         du11_aux(1:7) = zero                                                   
         du21_aux(1:7) = zero                                                   
         du31_aux(1:7) = zero                                                   
         mu_front          = e_front / (two * ( one + nu_front ))               
         if ( debug ) write(out,9916) e_front, nu_front, mu_front               
c                                                                               
c             plane stress auxiliary fields for KI, KII, T11                    
c                                                                               
         kappa = (three - nu_front) / (one + nu_front)                          
         if ( debug ) write(out,9917) kappa                                     
c                                                                               
c                  KI                                                           
c                                                                               
         du11_aux(1) = zero                                                     
         du21_aux(1) = one / four / mu_front / sqrt( two * pi * r )             
     &               * (kappa + one) * cos(t) * sin(t/two)                      
         du31_aux(1) = zero                                                     
c                                                                               
c                  KII                                                          
c                                                                               
         du11_aux(3) = one / four / mu_front / sqrt( two * pi * r )             
     &               * (kappa + one) * cos(t) * sin(t/two)                      
         du21_aux(3) = zero                                                     
         du31_aux(3) = zero                                                     
c                                                                               
c                  T11                                                          
c                                                                               
         du11_aux(6) = one / (pi * e_front * r)                                 
         du21_aux(6) = zero                                                     
         du31_aux(6) = zero                                                     
c                                                                               
c             plane strain auxiliary fields for KI, KII, T11                    
c                                                                               
         kappa = three - four * nu_front                                        
         if ( debug ) write(out,9917) kappa                                     
c                                                                               
c                  KI                                                           
c                                                                               
         du11_aux(2) = zero                                                     
         du21_aux(2) = one / four / mu_front / sqrt( two * pi * r )             
     &                 * (kappa + one) * cos(t) * sin(t/two)                    
         du31_aux(2) = zero                                                     
c                                                                               
c                  KII                                                          
c                                                                               
         du11_aux(4) = one / four / mu_front / sqrt( two * pi * r )             
     &               * (kappa + one) * cos(t) * sin(t/two)                      
         du21_aux(4) = zero                                                     
         du31_aux(4) = zero                                                     
c                                                                               
c                  T11                                                          
c                                                                               
         du11_aux(7) = ( one - nu_front**2 ) / ( pi * e_front * r )             
         du21_aux(7) = zero                                                     
         du31_aux(7) = zero                                                     
c                                                                               
c             anti-plane shear auxiliary fields for KIII                        
c                                                                               
c                  KIII                                                         
c                                                                               
         du11_aux(5) = zero                                                     
         du21_aux(5) = zero                                                     
         du31_aux(5) = one / mu_front / sqrt( two * pi * r )                    
     &               * cos(t) * sin(t/two)                                      
c                                                                               
c             10. store values of ui,1 for use in method 2 if necessary.        
c                                                                               
         do j = 1,7                                                             
            cderiv_aux(ptno,1,j) = du11_aux(j)                                  
            cderiv_aux(ptno,2,j) = du21_aux(j)                                  
            cderiv_aux(ptno,3,j) = du31_aux(j)                                  
            if( debug ) then                                                    
               write( out,9920 ) du11_aux(j), du21_aux(j),                      
     &                           du31_aux(j)                                    
            end if                                                              
         end do                                                                 
c                                                                               
c             11. calculate surface traction integral and add to                
c                 running total                                                 
c                                                                               
         do j = 1,7                                                             
            termu = - cf_tractions(1) * du11_aux(j)                             
     &            *   gpq * weight * det_j                                      
            termv = - cf_tractions(2) * du21_aux(j)                             
     &            *   gpq * weight * det_j                                      
            termw = - cf_tractions(3) * du31_aux(j)                             
     &            *   gpq * weight * det_j                                      
c                                                                               
c             for crack-front elements, adjust final values:                    
c             multiply by the jacobian of the interval                          
c             shift, i.e. 1/2**0.5, and by the jacobian of the                  
c             substitution, i.e. 2*eta_shifted.                                 
c                                                                               
            if( shift ) then                                                    
               termu = termu * eta_shifted * two**half                          
               termv = termv * eta_shifted * two**half                          
               termw = termw * eta_shifted * two**half                          
            end if                                                              
            iterm(8,j) = iterm(8,j) + termu + termv + termw                     
            if ( debug ) then                                                   
               write(out,9945) j                                                
               write(out,9955) termu, termv, termw, iterm(8,j)                  
            end if                                                              
         end do                                                                 
c                                                                               
      end do                                                                    
c                                                                               
      if( .not. approx ) return                                                 
      if( debug ) write(out,*) "approximate surface integration..."             
c                                                                               
c ----------------------------------------------------------------------        
c ----------------------------------------------------------------------        
c                                                                               
c             method 2: for case where crack-face tractions are not             
c                       specified in the domain definition:                     
c                                                                               
c             loop over face nodes to extrapolate ui,1 from                     
c             gauss points to nodes. calculate contribution                     
c             from each node to the integral (ti * ui,1 * q)                    
c                                                                               
      iterm(8,1:7) = zero                                                       
      if( debug ) write(out,9925)                                               
      do fnode = 1, nfnode                                                      
         enode = fnodes(fnode)                                                  
         if( debug ) then                                                       
            write(out,9930) fnode, enode, sfnodes(fnode)                        
            write(out,9935) eloads(1,enode), eloads(2,enode),                   
     &                      eloads(3,enode), qvals(enode)                       
         end if                                                                 
c                                                                               
c             12. obtain parametric coordinates of current face node.           
c                                                                               
         call ndpts1( idumvec, 0, dum_vec, etype, fnode, xsi, eta, zeta)        
c                                                                               
c             13. evaluate "gauss point shape functions" with parametric        
c                coordinates of current face node.                              
c                                                                               
         call oulgr3( xsi, eta, lg, order_face )                                
         if( debug ) then                                                       
            write(out,9940) xsi, eta                                            
            do j=1, ngpts_face                                                  
               write(out,9942) j, lg(j)                                         
            end do                                                              
         end if                                                                 
c                                                                               
c             14. calculate ui,1 at current node by extrapolating gp value:     
c                 (ui,1)node = sum [ (N)gp * (ui,1)gp ]                         
c                                                                               
         do j = 1,7                                                             
            dux = zero                                                          
            dvx = zero                                                          
            dwx = zero                                                          
            dux = dieldp( lg, cderiv_aux(1,1,j), ngpts_face, 1, 1 )             
            dvx = dieldp( lg, cderiv_aux(1,2,j), ngpts_face, 1, 1 )             
            dwx = dieldp( lg, cderiv_aux(1,3,j), ngpts_face, 1, 1 )             
c                                                                               
c             surface traction term is given by (-) product of nodal            
c             q-values, the equivalent nodal loads, and auxiliary               
c             displacement derivatives.                                         
c                                                                               
            termu = - eloads(1,enode) * qvals(enode) * dux                      
            termv = - eloads(2,enode) * qvals(enode) * dvx                      
            termw = - eloads(3,enode) * qvals(enode) * dwx                      
            iterm(8,j) = iterm(8,j) + termu + termv + termw                     
            if ( debug ) then                                                   
               write(out,9945) j                                                
               write(out,9950) dux, dvx, dwx                                    
               write(out,9955) termu, termv,                                    
     &              termw, iterm(8,j)                                           
               end if                                                           
         end do                                                                 
      end do                                                                    
c                                                                               
c                                                                               
      return                                                                    
c                                                                               
 9900 format(///,'coordinates of nodes on loaded face:',                        
     &       /,3x,'node',4x,'snode',8x,'x',14x,'y',14x,'z')                     
 9905 format(2x,i4,2x,i8,3(2x,e13.6))                                           
 9907 format(/,'>>>shifting gauss-point eta value for cf element...')           
 9908 format(/,'     >> eta, eta_shifted, eta_t: ',3(2x,e13.6))                 
 9909 format(//,'---gauss point ',i2,'---',/)                                   
 9910 format(/,'     >> xsi, eta, zeta, weight: ',4(2x,e13.6))                  
 9911 format(/,'     >> shape functions and derivatives:')                      
 9912 format(10x,'node ',i2,3(2x,e13.6))                                        
 9913 format('     >> jacobian: ',2(2x,e13.6),/,18x,2(2x,e13.6),/               
     &       '     >> det|j|:',2x,e13.6)                                        
 9915 format(/,'     >> ptno, x, y, z, gpq, ave_y_coord: ',i2,                  
     &       5(2x,e13.6),/)                                                     
 9916 format(/,'     >> e_front, nu_front, mu_front: ',3(2x,e13.6),/)           
 9917 format(/,'     >> kappa: ',2x,e13.6,/)                                    
 9920 format('     >> du11_aux, du21_aux, du31_aux: ',3(2x,e13.6))              
 9925 format(//,'---loop over face nodes for interaction integral:---')         
 9930 format(//'       >> face node, hex node, snode: ',i2,2x,i2,2x,i8)         
 9935 format('       >> P1, P2, P3, qval: ',4(2x,e13.6))                        
 9940 format('       >> xsi, eta:',2(2x,e13.6),/,                               
     &       '       >> shape functions:')                                      
 9942 format(10x,'gpt ',i2,2(2x,e13.6))                                         
 9945 format(/,'           >> j: ',2x,i2)                                       
 9950 format('           >> dux, dvx, dwx:',46x,3(2x,e13.6))                    
 9955 format('           >> termu_aux, termv_aux,',                             
     &       ' termw_aux, iterm(8,j)',4(2x,e13.6),/)                            
 9960 format(/,'     >>  crack-face tractions on element face:',/,              
     &        10x,'tx, ty, tz: ',3(2x,e13.6))                                   
 9965 format(/,5x,'node, fnode, qval: ',i2,2x,i2,2x,i7,2x,e13.6)                
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c***********************************************************************        
c                                                                      *        
c    subroutine to reorder the numbering of face nodes.                *        
c    this is necessary to ensure that xsi corresponds to x1, and       *        
c    eta corresponds to x3 for jacobian calculations.                  *        
c                                                                      *        
c                                 written by: mcw                      *        
c                              last modified: 11/10/03                 *        
c                                                                      *        
c***********************************************************************        
c                                                                               
      subroutine di_reorder_nodes( nfnode, fnodes, sfnodes,                     
     &                          fnode_coords, ave_y_coord, elemno,              
     &                          front_nodes, num_front_nodes,                   
     &                          front_coords, domain_type,                      
     &                          domain_origin, front_order,                     
     &                          crack_curvature, out )                          
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer nfnode, fnodes(*), sfnodes(*), elemno, front_nodes(*),            
     &        num_front_nodes, domain_type, domain_origin, front_order,         
     &        out                                                               
      double precision                                                          
     &     fnode_coords(3,*), ave_y_coord, front_coords(3,*),                   
     &     crack_curvature(*)                                                   
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, fnode, node1, node2, sfnode1, sfnode2, pair,                   
     &        fnode_temp(20), sfnode_temp(20), index                            
      double precision                                                          
     &     x1_1, x1_2, measure, zero, half, one, two, huge, d1, d1_old,         
     &     coord_temp(3,20), node_x, node_y, node_z, r, t, rs(4)                
      logical debug                                                             
c                                                                               
      data zero, half, one, two, huge                                           
     & / 0.0d0, 0.5d0, 1.0d0, 2.0d0, 1.0d10 /                                   
c                                                                               
      debug = .false.                                                           
c                                                                               
      if( debug ) write(out,1101) elemno                                        
c                                                                               
c             determine which edge of the element face is closest to            
c             the crack front. we simply determine the edge whose               
c             corner nodes have the smallest cumulative value of                
c             distance r to the crack front.                                    
c                                                                               
      do i=1,4                                                                  
         node_x = fnode_coords(1,i)                                             
         node_y = fnode_coords(2,i)                                             
         node_z = fnode_coords(3,i)                                             
         r      = zero                                                          
         t      = zero                                                          
         call di_calc_r_theta( 1, front_nodes, num_front_nodes,                 
     &                      front_coords, domain_type, domain_origin,           
     &                      front_order, node_x, node_y, node_z,                
     &                      elemno, sfnodes(i), r, t, crack_curvature,          
     &                      debug, out )                                        
         rs(i) = r                                                              
      end do                                                                    
c                                                                               
      d1     = huge                                                             
      d1_old = huge                                                             
      do i = 1,4                                                                
         node1   = fnodes(i)                                                    
         sfnode1 = sfnodes(i)                                                   
         x1_1    = rs(i)                                                        
         if( i.lt.4 ) then                                                      
            node2   = fnodes(i+1)                                               
            sfnode2 = sfnodes(i+1)                                              
            x1_2    = rs(i+1)                                                   
         else                                                                   
            node2   = fnodes(1)                                                 
            sfnode2 = sfnodes(1)                                                
            x1_2    = rs(1)                                                     
         end if                                                                 
         d1 = x1_1 + x1_2                                                       
         if( d1.le.d1_old ) then                                                
            d1_old = d1                                                         
            pair = i                                                            
         end if                                                                 
         if( debug ) write(out,1111) node1, node2, sfnode1, sfnode2,            
     &                               x1_1, x1_2, d1_old                         
      end do                                                                    
      if( debug ) write(out,1121) d1, pair                                      
c                                                                               
c             if closest pair is 1, order of face nodes is ok. return.          
c             if closest pair is other, restructure nodal arrays.               
c             'pair' represents the starting node for the new numbering.        
c                                                                               
      if( pair.eq.1 ) goto 7000                                                 
c                                                                               
      fnode_temp(1:nfnode)     = fnodes(1:nfnode)                               
      sfnode_temp(1:nfnode)    = sfnodes(1:nfnode)                              
      coord_temp(1:3,1:nfnode) = fnode_coords(1:3,1:nfnode)                     
c                                                                               
c             reorder corner nodes                                              
c                                                                               
      do i = 1,4                                                                
         index = pair + (i - 1)                                                 
         if( index.gt.4 ) index = index - 4                                     
         fnodes(i)           = fnode_temp(index)                                
         sfnodes(i)          = sfnode_temp(index)                               
         fnode_coords(1:3,i) = coord_temp(1:3,index)                            
      end do                                                                    
      if( nfnode.eq.4 ) goto 7000                                               
c                                                                               
c             reorder mid-side nodes.                                           
c                                                                               
      do i = 5,8                                                                
         index = pair + (i - 1)                                                 
         if( index.gt.8 ) index = index - 4                                     
         fnodes(i)           = fnode_temp(index)                                
         sfnodes(i)          = sfnode_temp(index)                               
         fnode_coords(1:3,i) = coord_temp(1:3,index)                            
      end do                                                                    
c                                                                               
 7000 continue                                                                  
c                                                                               
      if( debug ) then                                                          
         write(out,1131)                                                        
         do i=1, nfnode                                                         
            write(out,1141) fnodes(i), sfnodes(i),                              
     &         fnode_coords(1,i), fnode_coords(2,i), fnode_coords(3,i)          
         end do                                                                 
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 1101 format(/,5x,' >>>> reordering face nodes on element ',i8)                 
 1111 format(/,10x,'node1, node2:',2(2x,i7),/,10x,'sfnode1, sfnode2:',          
     &       2(2x,i7),/,10x,'x1_1, x1_2, d1:',3(2x,e13.6))                      
 1121 format(/,10x,'distance, closest pair: ',2x,e13.6,2x,i8 )                  
 1131 format(/,10x,'numbering of reordered nodes:',                             
     &       /,13x,'node',4x,'snode',8x,'x',14x,'y',14x,'z')                    
 1141 format(12x,i4,2x,i8,3(2x,e13.6))                                          
c                                                                               
      end                                                                       
c                                                                               
