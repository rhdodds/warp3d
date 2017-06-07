c ****************************************************************              
c *                                                              *              
c *      equivalent nodal loads -- constant body force           *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
       subroutine body_load( element, etype, nnode, body_dir,                   
     &                       body_intens, ecoord, equiv_loads_bdf )             
      implicit double precision (a-h,o-z)                                       
c                                                                               
      dimension  equiv_loads_bdf(*), ecoord(3,*)                                
      integer element, etype, body_dir                                          
c                                                                               
c                local arrays                                                   
c                                                                               
      double precision                                                          
     &    jacob, jacobi, sixth                                                  
c                                                                               
      dimension qvec(32), dsf(32,3), jacob(3,3), jacobi(3,3),                   
     &          sf(32)                                                          
      logical local_debug, tet_elem                                             
      integer enode                                                             
      data zero, local_debug / 0.0, .false. /                                   
c                                                                               
c                                                                               
c                if element is a tetrahedron, we must use a                     
c                1/6 multiple for integration.                                  
c                                                                               
c                                                                               
      tet_elem = etype .eq. 6 .or. etype .eq. 13                                
c                                                                               
c                                                                               
      sixth = 1.0D0/6.0D0                                                       
c                                                                               
c                                                                               
c                compute the {q} vector by integrating the                      
c                shape functions over the element volume.                       
c                                                                               
c                compute element volume. We know the distribution               
c                of the total body force to each node once                      
c                we have the element volume.                                    
c                                                                               
c                 1)  loop for all integration points over the                  
c                     element;                                                  
c                 2)  compute derivatives of shape functions at                 
c                     the integration point;                                    
c                 3)  compute the jacobian and its determinant at               
c                     the integration point;                                    
c                 4)  evaluate shape functions at the integration               
c                     point;                                                    
c                 5)  form partial contribution to each {q}                     
c                     term and add to existing terms.                           
c                                                                               
c                                                                               
      if ( local_debug ) then                                                   
         write(*,*) '>> inside body load'                                       
         write(*,*) '     element coordinates'                                  
         do i = 1, nnode                                                        
           write(*,9900) i,ecoord(1,i),ecoord(2,i),                             
     &                   ecoord(3,i)                                            
         end do                                                                 
         write(*,*) '   nnode,etype: ',nnode,etype                              
      end if                                                                    
c                                                                               
      vol = zero                                                                
      do i = 1, nnode                                                           
        qvec(i) = zero                                                          
      end do                                                                    
c                                                                               
c                use full integration for elements with                         
c                quadratic faces. for tet elements, max                         
c                integration is 15 point rule.                                  
c                                                                               
c                                                                               
      if ( tet_elem ) then                                                      
           ngp    = 15                                                          
           iorder = 15                                                          
      else                                                                      
         if ( nnode .gt. 8 ) then                                               
             ngp    = 27                                                        
             iorder = 1                                                         
         else                                                                   
             ngp    = 8                                                         
             iorder = 1                                                         
         end if                                                                 
      end if                                                                    
c                                                                               
c                                                                               
c                get all relevant quantities at each                            
c                integrationt point. if tet element,                            
c                multiply scalar by 1/6.                                        
c                                                                               
c                                                                               
      do igp = 1, ngp                                                           
        call getgpts( etype, iorder, igp, xi, eta, zeta, weight )               
        call derivs( etype, xi, eta, zeta, dsf(1,1), dsf(1,2),                  
     &               dsf(1,3) )                                                 
        call eqldjb( dsf, ecoord, nnode, jacob, jacobi, det, ierr )             
        call shapef( etype, xi, eta, zeta, sf )                                 
        scale = weight * det                                                    
c                                                                               
        if (tet_elem) scale = scale * sixth                                     
c                                                                               
        do enode = 1, nnode                                                     
         qvec(enode) = qvec(enode) + sf(enode)*scale                            
        end do                                                                  
        vol = vol + scale                                                       
      end do                                                                    
c                                                                               
c              equivalent loads are distribution factors, {q},                  
c              x the intensity of the body force                                
c                                                                               
      totlod = zero                                                             
      do enode = 1, nnode                                                       
        equiv_loads_bdf(enode) = qvec(enode) * body_intens                      
        totlod = totlod + equiv_loads_bdf(enode)                                
      end do                                                                    
c                                                                               
      if ( local_debug ) then                                                   
        write(*,3040) ( equiv_loads_bdf(i),i=1,nnode )                          
        write(*,3050) totlod, vol                                               
      end if                                                                    
      return                                                                    
c                                                                               
c                                                                               
 3040 format( 1h0, 5x,32hequivalent loads for body forces,                      
     &    / 20(/,10x,f20.6 ) )                                                  
 3010 format(1x,10f10.3)                                                        
 3050 format(1x,"total load, vol: ",2f10.4 )                                    
 9900 format(i4,3f10.3)                                                         
c                                                                               
      end                                                                       
c ****************************************************************              
c *                                                              *              
c *      equivalent nodal loads -- constant face force           *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
       subroutine face_load( element, etype, nnode, face,                       
     &                       face_intens, ecoord, equiv_loads_face,             
     &                       ldtype, constant_intens )                          
c                                                                               
      implicit double precision (a-h,o-z)                                       
c                                                                               
      dimension  equiv_loads_face(*),                                           
     &           ecoord(3,*),                                                   
     &           face_intens(*)                                                 
c                                                                               
      integer element, etype, face                                              
c                                                                               
c                                                                               
c                local arrays                                                   
c                                                                               
c                                                                               
      dimension qmat(10,10), fnodes(10)                                         
c                                                                               
      logical local_debug, tet_elem                                             
      integer ptno, fnodes                                                      
c                                                                               
      data zero, one, local_debug                                               
     & / 0.d0, 1.d0, .false. /                                                  
c                                                                               
c                                                                               
c                check if we've got a tet element. if so,                       
c                call a special subroutine later. if face                       
c                greater than 4, print an error and quit.                       
c                                                                               
c                                                                               
      tet_elem = etype .eq. 6 .or. etype .eq. 13                                
      if (tet_elem .and. face .gt. 4) then                                      
         write(*,9920)                                                          
         call die_abort                                                         
      end if                                                                    
c                                                                               
c                                                                               
c                compute the {q} vector by integrating the                      
c                shape functions over the element volume.                       
c                                                                               
c                compute element volume. We know the distribution               
c                of the total body force to each node once                      
c                we have the element volume.                                    
c                                                                               
c                 1)  loop for all integration points over the                  
c                     element;                                                  
c                 2)  compute derivatives of shape functions at                 
c                     the integration point;                                    
c                 3)  compute the jacobian and its determinant at               
c                     the integration point;                                    
c                 4)  evaluate shape functions at the integration               
c                     point;                                                    
c                 5)  form partial contribution to each {q}                     
c                     term and add to existing terms.                           
c                                                                               
c                                                                               
      if ( local_debug ) then                                                   
         write(*,*) '>> inside face load'                                       
         write(*,*) '     element coordinates'                                  
         do i = 1, nnode                                                        
           write(*,9900) i,ecoord(1,i),ecoord(2,i),                             
     &                   ecoord(3,i)                                            
         end do                                                                 
         write(*,*) '   nnode,etype: ',nnode,etype                              
      end if                                                                    
c                                                                               
c               ldtype = 1 means constant intensity over the                    
c               loaded face. fill vector of values. ldtype = 2                  
c               means variable intensity at face nodes. if                      
c               all values are zero, leave (could occur in                      
c               support of pressure loads).                                     
c                                                                               
      if ( ldtype .eq. 1 ) then                                                 
         do i = 1, nnode                                                        
           face_intens(i) = constant_intens                                     
         end do                                                                 
      else if ( ldtype .eq. 2 ) then                                            
         sum = zero                                                             
         do i = 1, nnode                                                        
           sum = sum + abs(face_intens(i))                                      
         end do                                                                 
         if ( sum .eq. zero ) return                                            
      else                                                                      
         write(*,9910)                                                          
       call die_gracefully                                                      
         stop                                                                   
      end if                                                                    
c                                                                               
c                                                                               
c               tet and hex elements diverge here to compute [q].               
c               2-D gauss integration schemes are used on the                   
c               loaded face.                                                    
c                                                                               
c                                                                               
      if (tet_elem) then                                                        
         call tet_compute_q (element, etype, nnode, face,                       
     &                       face_intens, ecoord, equiv_loads_face,             
     &                       ldtype, constant_intens, qmat, fnodes,             
     &                       nfnode, area)                                      
      else                                                                      
         call hex_compute_q (element, etype, nnode, face,                       
     &                       face_intens, ecoord, equiv_loads_face,             
     &                       ldtype, constant_intens, qmat, fnodes,             
     &                       nfnode, area)                                      
      end if                                                                    
c                                                                               
c             equivalent loads for face nodes are given by                      
c             [q] * intensities.                                                
c                                                                               
      totlod = zero                                                             
      do i = 1, nfnode                                                          
       sum = zero                                                               
       do j = 1, nfnode                                                         
         sum = sum + qmat(i,j) * face_intens(fnodes(j))                         
       end do                                                                   
       equiv_loads_face(fnodes(i)) = sum                                        
       totlod = totlod + sum                                                    
      end do                                                                    
c                                                                               
      if ( .not. local_debug ) return                                           
      write(*,3000)                                                             
      do i = 1, nfnode                                                          
        write(*,3010) (qmat(i,j),j=1,nfnode)                                    
      end do                                                                    
      write(*,3020) ( i,equiv_loads_face(i),i=1,nnode )                         
      write(*,3030) totlod, area                                                
      return                                                                    
c                                                                               
c                                                                               
c                                                                               
 3020 format( 1h0,5x,38hequivalent nodal loads for edge forces,                 
     &   / 20(/,10x,i4,f20.6) )                                                 
 3000 format(/,2x,"[q] matrix:" )                                               
 3010 format(1x,10f10.3)                                                        
 3030 format(1x,"total load, area: ",2f10.4 )                                   
 9900 format(i4,3f10.3)                                                         
 9910 format('>> FATAL ERROR: routine face_load. job aborted')                  
 9920 format('>> FATAL ERROR: tet face > 4 not valid. job aborted')             
      end                                                                       
c                                                                               
c                                                                               
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *    subroutine hex_compute_q                                     *           
c *                                                                 *           
c *               compute the [q] matrix for hex elements.          *           
c *               this routine contains 2D gauss rules for          *           
c *               the quad faces.                                   *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine hex_compute_q (element, etype, nnode, face,                    
     &                       face_intens, ecoord, equiv_loads_face,             
     &                       ldtype, constant_intens, qmat, fnodes,             
     &                       nfnode, area)                                      
c                                                                               
c                                                                               
      implicit double precision (a-h,o-z)                                       
c                                                                               
      dimension  equiv_loads_face(*),                                           
     &           ecoord(3,*),                                                   
     &           face_intens(*),                                                
     &           qmat(10,*),                                                    
     &           fnodes(*)                                                      
c                                                                               
      integer element, etype, face                                              
c                                                                               
c                                                                               
c                local arrays                                                   
c                                                                               
c                                                                               
      double precision                                                          
     &    jacob, jacobi                                                         
c                                                                               
      dimension dsf(32,3), jacob(3,3), jacobi(3,3),                             
     &          sf(32), fcoor(3,9), fweigt(3),                                  
     &          gauss(3)                                                        
c                                                                               
      logical local_debug                                                       
      integer ptno, fnodes                                                      
c                                                                               
      data zero, one, local_debug                                               
     & / 0.d0, 1.d0, .false. /                                                  
c                                                                               
c                                                                               
c               gauss points and weights for 2x2 and                            
c               3x3 integration over the loaded face.                           
c                                                                               
c                                                                               
      data gp1, gp2, w1, w2                                                     
     &  / 0.57735026918962576450d0,                                             
     &    0.77459666924148337703d0,                                             
     &    0.55555555555555555555d0,                                             
     &    0.88888888888888888888d0 /                                            
c                                                                               
c               get the element nodes on the loaded face and                    
c               zero the [q] matrix.                                            
c                                                                               
      call eqelfn( fnodes, etype, face, nfnode )                                
c                                                                               
      do i = 1, nfnode                                                          
       do j = 1, nfnode                                                         
         qmat(j,i) = zero                                                       
       end do                                                                   
      end do                                                                    
c                                                                               
c               integrate the product of shape functions over                   
c               the loaded face to obtain the [q] matrix.                       
c                                                                               
c                 1)  loop over all integration points on the face;             
c                 2)  compute derivatives of shape functions at                 
c                     the integration point;                                    
c                 3)  compute jacobian matrix at integration point;             
c                 4)  evaluate shape functions at the point;                    
c                 5)  compute differential area ( ratio of area                 
c                     in parent and real element )                              
c                 6)  add the integration point contribution                    
c                     to each [q] term on the face.                             
c                                                                               
      if ( nfnode .eq. 4 ) then                                                 
         iorder = 2                                                             
         gauss(1) = -gp1                                                        
         gauss(2) =  gp1                                                        
         fweigt(1) = one                                                        
         fweigt(2) = one                                                        
      else                                                                      
         iorder    = 3                                                          
         gauss(1)  = -gp2                                                       
         gauss(2)  = zero                                                       
         gauss(3)  = gp2                                                        
         fweigt(1) = w1                                                         
         fweigt(2) = w2                                                         
         fweigt(3) = w1                                                         
      end if                                                                    
c                                                                               
      call eqfnic( fcoor, iorder, gauss, face )                                 
      ptno = 0                                                                  
      area = zero                                                               
      do i = 1, iorder                                                          
        do j = 1, iorder                                                        
          ptno = ptno + 1                                                       
          call derivs( etype, fcoor(1,ptno), fcoor(2,ptno),                     
     &                 fcoor(3,ptno), dsf(1,1), dsf(1,2), dsf(1,3) )            
          call eqldjb( dsf, ecoord, nnode, jacob, jacobi, det, ierr )           
          call shapef( etype, fcoor(1,ptno), fcoor(2,ptno),                     
     &                 fcoor(3,ptno), sf )                                      
          call eqfcda( jacob, face, darea )                                     
          weight = fweigt(i) * fweigt(j) * darea                                
          area   = area + weight                                                
          do irow = 1, nfnode                                                   
            a = sf(fnodes(irow)) * weight                                       
            do jcol = 1, nfnode                                                 
              qmat(irow,jcol) = qmat(irow,jcol) +  a * sf(fnodes(jcol))         
            end do                                                              
          end do                                                                
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *    subroutine tet_compute_q                                     *           
c *                                                                 *           
c *               compute the [q] matrix for tetrahedral            *           
c *               elements.                                         *           
c *                                                                 *           
c *******************************************************************           
      subroutine tet_compute_q (element, etype, nnode, face,                    
     &                       face_intens, ecoord, equiv_loads_face,             
     &                       ldtype, constant_intens, qmat, fnodes,             
     &                       nfnode, area)                                      
c                                                                               
c                                                                               
      implicit double precision (a-h,o-z)                                       
c                                                                               
      dimension  equiv_loads_face(*),                                           
     &           ecoord(3,*),                                                   
     &           face_intens(*),                                                
     &           qmat(10,*),                                                    
     &           fnodes(*)                                                      
c                                                                               
      integer element, etype, face, index, gp, irow, jcol, i, j                 
c                                                                               
c                                                                               
c                local arrays                                                   
c                                                                               
c                                                                               
      double precision                                                          
     &    jacob, jacobi, temp_dsf, temp_sf, vec1, vec2                          
c                                                                               
      dimension dsf(32,3), jacob(3,3), jacobi(3,3),                             
     &          sf(32), fcoor(3,9), fweight(7),                                 
     &          gauss(3), temp_dsf(32,3), temp_sf(32),                          
     &          vec1(3), vec2(3)                                                
c                                                                               
      logical local_debug                                                       
      integer ptno, fnodes                                                      
c                                                                               
      data zero, one, local_debug                                               
     & / 0.d0, 1.d0, .false. /                                                  
c                                                                               
      ngpts = 0                                                                 
c                                                                               
c                                                                               
c               get the element nodes on the loaded face and                    
c               zero the [q] matrix.                                            
c                                                                               
c                                                                               
      call tet_get_nodes ( fnodes, etype, face, nfnode )                        
c                                                                               
c                                                                               
c               get the isoparametric coordinates of the                        
c               integration points for 3 and 6 noded tri's.                     
c               each column of 'fcoor' contains the coords                      
c               for a gauss point.                                              
c                                                                               
c                                                                               
      call tet_get_gpts ( etype, fcoor, fweight, ngpts, index )                 
c                                                                               
c                                                                               
c               compute the [q] matrix for the element. note                    
c               that the shape functions and derivatives are                    
c               calculated using triangular interface element                   
c               routines. the 'index' argument determines                       
c               whether the linear of quadratic tri is used.                    
c                                                                               
c               this was necessary to simplify the computation                  
c               of 'darea' for the integragtion. problems                       
c               arose when derivatives with respect to the                      
c               dependent natural coordinate were needed.                       
c               treating the case in 2D solved these problems.                  
c                                                                               
c                                                                               
      do i = 1, nfnode                                                          
       do j = 1, nfnode                                                         
         qmat(j,i) = zero                                                       
       end do                                                                   
      end do                                                                    
c                                                                               
      area = 0.0                                                                
      do gp = 1, ngpts                                                          
c                                                                               
         call derivs( index, fcoor(1,gp), fcoor(2,gp), fcoor(3,gp),             
     &                temp_dsf(1,1), temp_dsf(1,2), temp_dsf(1,3) )             
         call shapef( index, fcoor(1,gp), fcoor(2,gp), fcoor(3,gp),             
     &                temp_sf )                                                 
         call map_tet_tri( etype, nfnode, fnodes, sf, temp_sf,                  
     &                     dsf, temp_dsf )                                      
         call tet_darea( dsf, nfnode, fnodes, ecoord, darea,                    
     &                   vec1, vec2 )                                           
c                                                                               
         area = area + fweight(gp)*darea                                        
c                                                                               
         do irow = 1, nfnode                                                    
           a = sf(fnodes(irow)) * fweight(gp) * darea                           
           do jcol = 1, nfnode                                                  
              qmat(irow,jcol) = qmat(irow,jcol) +  a * sf(fnodes(jcol))         
           end do                                                               
         end do                                                                 
c                                                                               
      end do                                                                    
c                                                                               
c                                                                               
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *         compute darea for numerical integration of tets         *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine tet_darea( dsf, nfnode, fnodes, ecoord, darea,                 
     &                      vec1, vec2 )                                        
      implicit none                                                             
c                                                                               
      double precision                                                          
     &    dsf(32,*), ecoord(3,*), darea, vec1(*), vec2(*)                       
c                                                                               
      integer fnodes(*), nfnode                                                 
c                                                                               
c                 local variables                                               
c                                                                               
      double precision                                                          
     &    a1, a2, a3                                                            
c                                                                               
      integer i, loc                                                            
c                                                                               
c                                                                               
c                 first build the vectors using the derivatives                 
c                 of shape functions and nodal coordinates on                   
c                 the loaded face.                                              
c                                                                               
c                                                                               
      do i=1,3                                                                  
      vec1(i) = 0.d0                                                            
      vec2(i) = 0.d0                                                            
      end do                                                                    
c                                                                               
      do i=1,nfnode                                                             
         loc = fnodes(i)                                                        
         vec1(1) = vec1(1) + dsf(loc,1)*ecoord(1,loc)                           
         vec1(2) = vec1(2) + dsf(loc,1)*ecoord(2,loc)                           
         vec1(3) = vec1(3) + dsf(loc,1)*ecoord(3,loc)                           
c                                                                               
         vec2(1) = vec2(1) + dsf(loc,2)*ecoord(1,loc)                           
         vec2(2) = vec2(2) + dsf(loc,2)*ecoord(2,loc)                           
         vec2(3) = vec2(3) + dsf(loc,2)*ecoord(3,loc)                           
      end do                                                                    
c                                                                               
c                                                                               
c                  darea is simply the magnitude of the                         
c                  cross product of the two vectors.                            
c                                                                               
c                                                                               
      a1 = vec1(2)*vec2(3) - vec2(2)*vec1(3)                                    
      a2 = vec2(1)*vec1(3) - vec1(1)*vec2(3)                                    
      a3 = vec1(1)*vec2(2) - vec2(1)*vec1(2)                                    
      darea = sqrt( a1*a1 + a2*a2 + a3*a3 )                                     
c                                                                               
c                                                                               
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *         map the 2D triangular nodes to 3D tet nodes.            *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine map_tet_tri( etype, nfnode, fnodes, sf, temp_sf,               
     &                     dsf, temp_dsf )                                      
      implicit none                                                             
c                                                                               
      double precision                                                          
     &    sf(*), dsf(32,*), temp_dsf(32,*), temp_sf(*)                          
c                                                                               
      integer etype, nfnode, fnodes(*)                                          
c                                                                               
c                local variables                                                
c                                                                               
      integer node, map_3node(3), map_6node(6)                                  
c                                                                               
      logical etype_tet4, etype_tet10                                           
c                                                                               
c                                                                               
c                use the map vectors to match the shape                         
c                functions and derivatives with the                             
c                proper nodes.                                                  
c                                                                               
c                example:                                                       
c                                                                               
c        tri node scheme                  tet node scheme                       
c        ---------------                  ---------------                       
c        3                                c                                     
c         o                               o                                     
c         |  \                            |  \                                  
c         |    \                          |    \                                
c         |      \  5                     |      \  e                           
c       6 o        o                    f o        o                            
c         |          \                    |          \                          
c         |            \                  |            \                        
c         |              \                |              \                      
c         o - - - o - - - o               o - - - o - - - o                     
c        1        4        2             a        d         b                   
c                                                                               
c                                                                               
c                the 'fnodes' array for this example would be:                  
c                                                                               
c                fnodes    = { a,d,b,e,c,f }                                    
c                map_6node = { 1,4,2,5,3,6 }                                    
c                                                                               
c                the node given by tri node 'map(i)' corresponds                
c                to the tet node 'fnodes(i)'.                                   
c                                                                               
c                                                                               
      data   map_3node   / 1,2,3 /                                              
      data   map_6node   / 1,4,2,5,3,6 /                                        
c                                                                               
c                                                                               
c                determine if linear or quadratic tet.                          
c                                                                               
c                                                                               
      etype_tet10 = etype .eq. 6                                                
      etype_tet4  = etype .eq. 13                                               
c                                                                               
c                                                                               
c                scatter the 2D information into the full                       
c                3D arrays.                                                     
c                                                                               
c                                                                               
      if (etype_tet4) then                                                      
         do node = 1,nfnode                                                     
            sf(fnodes(node))    = temp_sf(map_3node(node))                      
            dsf(fnodes(node),1) = temp_dsf(map_3node(node),1)                   
            dsf(fnodes(node),2) = temp_dsf(map_3node(node),2)                   
         end do                                                                 
      else if (etype_tet10) then                                                
         do node = 1,nfnode                                                     
            sf(fnodes(node))    = temp_sf(map_6node(node))                      
            dsf(fnodes(node),1) = temp_dsf(map_6node(node),1)                   
            dsf(fnodes(node),2) = temp_dsf(map_6node(node),2)                   
         end do                                                                 
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *         provide isoparametric coordinates for gauss             *           
c *         points on triangular faces.                             *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine tet_get_gpts (etype, fcoor, fweight, ngpts, index)             
      implicit none                                                             
c                                                                               
      integer etype, ngpts, order, gp, index                                    
c                                                                               
      double precision                                                          
     &    fcoor(3,*),                                                           
     &    fweight(*),                                                           
     &    s1, s2, s3,                                                           
     &    weight, zero, one                                                     
c                                                                               
      logical etype_tet4, etype_tet10                                           
c                                                                               
      data zero, one                                                            
     & / 0.d0, 1.d0 /                                                           
c                                                                               
c                                                                               
c               determine whether linear or quadratic tet.                      
c               then set the number of gauss points on the                      
c               element face.                                                   
c                                                                               
c               the 'index' variable is used when the shape                     
c               functions or derivatives are calculated. the                    
c               subroutines for triangular interface elements                   
c               are used to integrate over the triangular                       
c               face. 'index' determines whether the linear                     
c               or quadratic tri routines are called.                           
c                                                                               
c                                                                               
      etype_tet10 = etype .eq. 6                                                
      etype_tet4  = etype .eq. 13                                               
c                                                                               
      if (etype_tet4) then                                                      
        ngpts = 4                                                               
        order = 4                                                               
        index = 14                                                              
      else if (etype_tet10) then                                                
        ngpts = 7                                                               
        order = 7                                                               
        index = 15                                                              
      end if                                                                    
c                                                                               
c                                                                               
c               loop over the number of gauss points to                         
c               build the 'fcoor' and 'fweight' arrays.                         
c               use the subroutines in 'getgpts.f' created                      
c               for triangular interface elements.                              
c                                                                               
c               NOTE: the trint12_gp routine has multiplied                     
c                     the weight by 0.5 in preparation for                      
c                     integration over a triangular face.                       
c                     therefore, this scaling will not be                       
c                     necessary during the integration in                       
c                     tet_compute_q.                                            
c                                                                               
c                     also, the trint12 routine doesn't return                  
c                     's3' for reasons related to interface                     
c                     elements. 's3' will be calculated in                      
c                     the shapef and derivs routines so it                      
c                     won't be calculated here.                                 
c                                                                               
c                                                                               
      do gp = 1,ngpts                                                           
         call trint12_gp( order, gp, s1, s2, s3, weight )                       
         fcoor (1,gp) = s1                                                      
         fcoor (2,gp) = s2                                                      
         fweight (gp) = weight                                                  
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *         provide face nodes for 3d isoparametrics for            *           
c *         tetrahedral elements.                                   *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine tet_get_nodes(fnodes, etype, face, nfnode)                     
      implicit integer (a-z)                                                    
c                                                                               
c              put the node numbers for face "face" , element                   
c              type "etype" into vector fnodes and set number of                
c              nodes on the face.                                               
c                                                                               
      dimension  fnodes(*), tet_10node(6,4), tet_4node(3,4)                     
c                                                                               
      logical etype_tet10, etype_tet4                                           
c                                                                               
c                                                                               
      data   tet_10node   / 1,7,3,6,2,5,                                        
     &                      1,5,2,9,4,8,                                        
     &                      2,6,3,10,4,9,                                       
     &                      1,8,4,10,3,7 /                                      
c                                                                               
      data   tet_4node    / 1,3,2,                                              
     &                      1,2,4,                                              
     &                      2,3,4,                                              
     &                      1,4,3 /                                             
c                                                                               
c                                                                               
c               determine whether linear or quadratic tet                       
c                                                                               
c                                                                               
      etype_tet10 = etype .eq. 6                                                
      etype_tet4  = etype .eq. 13                                               
c                                                                               
c                                                                               
c               store face nodes in 'fnodes'                                    
c                                                                               
c                                                                               
      if (etype_tet10) then                                                     
         nfnode = 6                                                             
         do i = 1, nfnode                                                       
            fnodes(i) = tet_10node(i,face)                                      
         end do                                                                 
         return                                                                 
      end if                                                                    
c                                                                               
      if (etype_tet4) then                                                      
         nfnode = 3                                                             
         do i = 1, nfnode                                                       
            fnodes(i) = tet_4node(i,face)                                       
         end do                                                                 
         return                                                                 
      end if                                                                    
c                                                                               
c                                                                               
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *    compute jacobian, its determinate, and inverse               *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine  eqldjb( dsf, coord, nnode, jacob, jacobi, det, ierr )         
      implicit double precision (a-h,o-z)                                       
c                                                                               
c                                                                               
c              compute the 3 x 3 jacobian, its determinate and                  
c              inverse for the 3-d isoparametrics.                              
c                                                                               
c                                                                               
      dimension  dsf(32,*), coord(3,*), jacob(3,3), jacobi(3,3)                 
      double precision                                                          
     &  jacob, jacobi                                                           
      logical  debug                                                            
      data zero, debug / 0.0, .false. /                                         
c                                                                               
c              compute jacobian at the point.                                   
c                                                                               
      jacob(1,1) = zero                                                         
      jacob(2,1) = zero                                                         
      jacob(3,1) = zero                                                         
      jacob(1,2) = zero                                                         
      jacob(2,2) = zero                                                         
      jacob(3,2) = zero                                                         
      jacob(1,3) = zero                                                         
      jacob(2,3) = zero                                                         
      jacob(3,3) = zero                                                         
      do k = 1, nnode                                                           
       jacob(1,1) = jacob(1,1) + dsf(k,1) * coord(1,k)                          
       jacob(2,1) = jacob(2,1) + dsf(k,2) * coord(1,k)                          
       jacob(3,1) = jacob(3,1) + dsf(k,3) * coord(1,k)                          
       jacob(1,2) = jacob(1,2) + dsf(k,1) * coord(2,k)                          
       jacob(2,2) = jacob(2,2) + dsf(k,2) * coord(2,k)                          
       jacob(3,2) = jacob(3,2) + dsf(k,3) * coord(2,k)                          
       jacob(1,3) = jacob(1,3) + dsf(k,1) * coord(3,k)                          
       jacob(2,3) = jacob(2,3) + dsf(k,2) * coord(3,k)                          
       jacob(3,3) = jacob(3,3) + dsf(k,3) * coord(3,k)                          
      end do                                                                    
      if ( debug ) write(*,100) ((jacob(i,j),j=1,3),i=1,3)                      
c                                                                               
c              inverse and determinant of the jacobian.                         
c                                                                               
      det = jacob(1,1) * jacob(2,2) * jacob(3,3)                                
     &    + jacob(2,1) * jacob(3,2) * jacob(1,3)                                
     &    + jacob(3,1) * jacob(1,2) * jacob(2,3)                                
     &    - jacob(1,1) * jacob(3,2) * jacob(2,3)                                
     &    - jacob(2,1) * jacob(1,2) * jacob(3,3)                                
     &    - jacob(3,1) * jacob(2,2) * jacob(1,3)                                
      if ( det .le. zero ) then                                                 
        ierr = 1                                                                
        return                                                                  
      end if                                                                    
c                                                                               
      jacobi(1,1) =  (jacob(2,2) * jacob(3,3)                                   
     &              - jacob(3,2) * jacob(2,3)) / det                            
      jacobi(2,2) =  (jacob(1,1) * jacob(3,3)                                   
     &              - jacob(3,1) * jacob(1,3)) / det                            
      jacobi(3,3) =  (jacob(1,1) * jacob(2,2)                                   
     &              - jacob(2,1) * jacob(1,2)) / det                            
      jacobi(2,1) = -(jacob(2,1) * jacob(3,3)                                   
     &              - jacob(3,1) * jacob(2,3)) / det                            
      jacobi(3,1) =  (jacob(2,1) * jacob(3,2)                                   
     &              - jacob(3,1) * jacob(2,2)) / det                            
      jacobi(1,2) = -(jacob(1,2) * jacob(3,3)                                   
     &              - jacob(3,2) * jacob(1,3)) / det                            
      jacobi(3,2) = -(jacob(1,1) * jacob(3,2)                                   
     &             -  jacob(3,1) * jacob(1,2)) / det                            
      jacobi(1,3) =  (jacob(1,2) * jacob(2,3)                                   
     &              - jacob(2,2) * jacob(1,3)) / det                            
      jacobi(2,3) = -(jacob(1,1) * jacob(2,3)                                   
     &              - jacob(2,1) * jacob(1,3)) / det                            
c                                                                               
      if ( debug ) write(*,110) det, ((jacobi(i,j),j=1,3), i = 1,3)             
c                                                                               
      return                                                                    
 100  format(1h0,5x,18hjacobian at point  ,/,                                   
     &                 2(/,7x,3f15.5) )                                         
 110  format(1h0,5x,12hdeterminant ,f15.5,                                      
     &       /,  5x,16hjacobian inverse  /,2(/,7x,3f15.5) )                     
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *         provide face nodes for 3d isoparametrics                *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine eqelfn( fnodes, etype, face, nfnode )                          
      implicit integer (a-z)                                                    
c                                                                               
c              put the node numbers for face "face" , element                   
c              type "etype" into vector fnodes and set number of                
c              nodes on the face.                                               
c                                                                               
      dimension  fnodes(*), f8nod(4,6), f20nod(8,6), f12nod(8,6),               
     &           f15nod(8,6), nfnodes(6,5), f9nod(5,6)                          
c                                                                               
      data   f8nod   / 4,3,2,1,                                                 
     &                 5,6,7,8,                                                 
     &                 1,2,6,5,                                                 
     &                 3,4,8,7,                                                 
     &                 2,3,7,6,                                                 
     &                 4,1,5,8 /                                                
c                                                                               
c              ordering of face nodes in f20nod coincides with the              
c              counter-clockwise numbering of 2D elements in                    
c              shapef.f subroutine "shape9"                                     
c                                                                               
      data   f20nod  /  4,3,2,1,11,10,9,12,                                     
     &                  5,6,7,8,13,14,15,16,                                    
     &                  1,2,6,5,9,18,13,17,                                     
     &                  3,4,8,7,11,20,15,19,                                    
     &                  2,3,7,6,10,19,14,18,                                    
     &                  4,1,5,8,12,17,16,20 /                                   
      data   f9nod  /  4,3,2,1,9,                                               
     &                 5,6,7,8,0,                                               
     &                 1,2,6,5,9,                                               
     &                 3,4,8,7,0,                                               
     &                 2,3,7,6,0,                                               
     &                 4,1,5,8,0 /                                              
      data   f12nod  / 4,3,2,1,12,9,10,11,                                      
     &                 5,6,7,8,0,0,0,0,                                         
     &                 1,2,6,5,9,0,0,0,                                         
     &                 3,4,8,7,11,0,0,0,                                        
     &                 2,3,7,6,10,0,0,0,                                        
     &                 4,1,5,8,12,0,0,0 /                                       
      data   f15nod / 4,3,2,1,12,9,10,11,                                       
     &                5,6,7,8,13,0,0,0,                                         
     &                1,2,6,5,9,15,13,14,                                       
     &                3,4,8,7,11,0,0,0,                                         
     &                2,3,7,6,10,15,0,0,                                        
     &                4,1,5,8,12,14,0,0 /                                       
      data nfnodes / 8,8,8,8,8,8,                                               
     &               4,4,4,4,4,4,                                               
     &               8,4,5,5,5,5,                                               
     &               8,5,8,5,6,6,                                               
     &               5,4,5,4,4,4 /                                              
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600, 1300 ), etype                       
c                                                                               
c             20 node element.                                                  
c                                                                               
 100  continue                                                                  
      nfnode = 8                                                                
      do i = 1, 8                                                               
        fnodes(i) = f20nod(i,face)                                              
      end do                                                                    
      return                                                                    
c                                                                               
c              8 node element.                                                  
c                                                                               
 200  continue                                                                  
      nfnode = 4                                                                
      do i = 1, 4                                                               
        fnodes(i) = f8nod(i,face)                                               
      end do                                                                    
      return                                                                    
c                                                                               
c             12 node element                                                   
c                                                                               
 300  continue                                                                  
      nfnode = nfnodes(face,3)                                                  
      do i = 1, nfnode                                                          
       fnodes(i) = f12nod(i,face)                                               
      end do                                                                    
      return                                                                    
c                                                                               
c             15 node element                                                   
c                                                                               
 400  continue                                                                  
      nfnode = nfnodes(face,4)                                                  
      do i = 1, nfnode                                                          
       fnodes(i) = f15nod(i,face)                                               
      end do                                                                    
      return                                                                    
c                                                                               
c             9 node element                                                    
c                                                                               
 500  continue                                                                  
      nfnode = nfnodes(face,5)                                                  
      do i = 1, nfnode                                                          
       fnodes(i) = f9nod(i,face)                                                
      end do                                                                    
      return                                                                    
c                                                                               
c            10 node tet element                                                
c                                                                               
 600  continue                                                                  
      call tet_get_nodes( fnodes, etype, face, nfnode )                         
      return                                                                    
c                                                                               
c             4 node tet element                                                
c                                                                               
 1300 continue                                                                  
      call tet_get_nodes( fnodes, etype, face, nfnode )                         
      return                                                                    
c                                                                               
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *   isoparametric coordinates of face integration points          *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine eqfnic( fcoor, iorder, gauss, face )                           
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
c              return the isoparametric coordinates of integration              
c              points on an element face.  order of points                      
c              is of no importance to user or integrators.                      
c                                                                               
c                                                                               
      double precision                                                          
     &   fcoor(3,*), gauss(*)                                                   
c                                                                               
c              local variables                                                  
c                                                                               
      double precision                                                          
     &  const(6), value                                                         
      integer row1(6),  row2(6), row3(6)                                        
c                                                                               
      data    const  /   -1.,1.,  -1.,1.,  -1.,1.  /                            
      data    row1   /  1,1,  2,2,  3,3  /                                      
      data    row2   /  2,2,  1,1,  1,1  /                                      
      data    row3   /  3,3,  3,3,  2,2  /                                      
c                                                                               
c                                                                               
      value = const(face)                                                       
      irow1 = row1(face)                                                        
      irow2 = row2(face)                                                        
      irow3 = row3(face)                                                        
      npts  = iorder * iorder                                                   
c                                                                               
c                    assign the isoparametric coordinate of                     
c                    integration points that is constant over the               
c                    face.                                                      
c                                                                               
      do ptno = 1, npts                                                         
        fcoor(irow1,ptno) = value                                               
      end do                                                                    
c                                                                               
c                    assign isoparametric coordinates of                        
c                    integration points that vary over the face.                
c                                                                               
      ptno = 0                                                                  
      do i = 1, iorder                                                          
       do j = 1, iorder                                                         
        ptno = ptno + 1                                                         
        fcoor(irow2,ptno) = gauss(i)                                            
        fcoor(irow3,ptno) = gauss(j)                                            
       end do                                                                   
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c ***************************************************************               
c *                                                             *               
c *      differential surface area computation                  *               
c *                                                             *               
c ***************************************************************               
c                                                                               
c                                                                               
      subroutine eqfcda( jacob, face, darea )                                   
c                                                                               
c                                                                               
c                 compute the differential area at some point                   
c                 on an element face given the jacobian evaluated               
c                 at that point.  take length of vector normal to               
c                 surface at point as the required area.  use                   
c                 cross product of two tangential vectors on the face.          
c                                                                               
c                                                                               
      double precision                                                          
     &  jacob(3,*), veca(3), vecb(3), darea, a1, a2, a3                         
      integer   face                                                            
c                                                                               
c                                                                               
      go to ( 100,100,200,200,300,300 ), face                                   
c                                                                               
c                face 1 and 2 --  cross jbar and kbar                           
c                                                                               
 100  continue                                                                  
      do  i = 1, 3                                                              
       veca(i) = jacob(2,i)                                                     
       vecb(i) = jacob(3,i)                                                     
      end do                                                                    
      go to 500                                                                 
c                                                                               
c                face 3 and 4 -- cross kbar and ibar                            
c                                                                               
 200  continue                                                                  
      do i = 1, 3                                                               
       veca(i) = jacob(3,i)                                                     
       vecb(i) = jacob(1,i)                                                     
      end do                                                                    
      go to 500                                                                 
c                                                                               
c                face 5 and 6 -- cross ibar and jbar                            
c                                                                               
 300  continue                                                                  
      do i = 1, 3                                                               
        veca(i) = jacob(1,i)                                                    
        vecb(i) = jacob(2,i)                                                    
      end do                                                                    
      go to 500                                                                 
c                                                                               
c                                                                               
c                find length of vector defined by cross product                 
c                of veca x vecb.                                                
c                                                                               
c                                                                               
 500  continue                                                                  
      a1 = veca(2)*vecb(3) - vecb(2)*veca(3)                                    
      a2 = vecb(1)*veca(3) - veca(1)*vecb(3)                                    
      a3 = veca(1)*vecb(2) - vecb(1)*veca(2)                                    
      darea = sqrt( a1*a1 + a2*a2 + a3*a3 )                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ****************************************************************              
c *                                                              *              
c *     face intensity on hex element from piston theory         *              
c *                                                              *              
c ****************************************************************              
c                                                                               
      subroutine piston_face_intens(                                            
     &   face_intens, element, etype, nnode, face,                              
     &   p3, u3, m3, gam, fdirc, tet_elem, elem_nodes, ecoord )                 
      use global_data ! old common.main
c                                                                               
      use main_data                                                             
c                                                                               
      implicit none                                                             
c                                                                               
c     parameter declarations                                                    
c     ----------------------                                                    
c                                                                               
      integer::  element, etype, nnode, face                                    
      integer :: elem_nodes(*)                                                  
      logical :: tet_elem                                                       
c                                                                               
      double precision                                                          
     &     face_intens, ecoord(3,*), p3, gam, m3, u3, fdirc(3)                  
                                                                                
c                                                                               
c     declare local variables                                                   
c     -----------------------                                                   
c                                                                               
      double precision ::                                                       
     &     jacob(3,3), jacobi(3,3), nrmvec(3), sf(32), dsf(32,3),               
     &     xi, eta, zeta, oecoord(3,mxndel), lvec(3),                           
     &     dwdt, dv1, dv2, dv3, dwidt, du1, du2, du3, dwi,                      
     &     dsidx, dsidy, dsidz, dsidf, dwdf, fpv, g1, g2, g3, nv1,              
     &     zero, half, one, four, twlv, x, y, z, det                            
      logical ::  debug                                                         
      integer :: enode, fnodes(10), node, nfnode, fnd, ierr                     
      data zero, half, one, four, twlv, debug                                   
     &     / 0.0d0, 0.5d0, 1.0d0, 4.0d0, 12.0d0, .false. /                      
c                                                                               
      if( debug ) write(*,*) '>> In piston_face_intens'                         
c                                                                               
      face_intens = zero                                                        
c                                                                               
c     determine the center of an element face                                   
c                                                                               
      call ccfe( xi, eta, zeta, face, tet_elem )                                
      if( debug ) write(*,9010) xi, eta, zeta                                   
c                                                                               
c     determine original face normal vector (outwards)                          
c                                                                               
      do node = 1, nnode                                                        
         oecoord(1,node) = c(crdmap(elem_nodes(node)))                          
         oecoord(2,node) = c(crdmap(elem_nodes(node))+1)                        
         oecoord(3,node) = c(crdmap(elem_nodes(node))+2)                        
      end do                                                                    
c                                                                               
      if( tet_elem ) then                                                       
         call tet_face_nvec( xi, eta, zeta, oecoord, lvec,                      
     &        element, etype, nnode, face )                                     
      else                                                                      
         call hex_face_nvec( xi, eta, zeta, oecoord, lvec,                      
     &        element, etype, nnode, face )                                     
      end if                                                                    
c                                                                               
      if(debug) write(*,9020) lvec(1), lvec(2), lvec(3)                         
c                                                                               
c     determine change in face velocity                                         
c                                                                               
      call eqelfn( fnodes, etype, face, nfnode )                                
      call shapef( etype, xi, eta, zeta, sf )                                   
      dwdt = zero                                                               
      if( debug ) then                                                          
         x = zero                                                               
         y = zero                                                               
         z = zero                                                               
      end if                                                                    
c                                                                               
c        get the velocities for each node on face and then                      
c        take the dot product with org. face vector                             
c                                                                               
      do node = 1, nfnode                                                       
         fnd = fnodes(node)                                                     
         dv1 = v(3*elem_nodes(fnd)-2)                                           
         dv2 = v(3*elem_nodes(fnd)-1)                                           
         dv3 = v(3*elem_nodes(fnd)  )                                           
         dwidt = dv1*lvec(1) + dv2*lvec(2) + dv3*lvec(3)                        
         if( debug ) write(*,9030) elem_nodes(fnd),                             
     &        dv1,dv2,dv3,dwidt                                                 
         if( debug ) then                                                       
            x = x + oecoord(1,fnd)*sf(fnd)                                      
            y = y + oecoord(2,fnd)*sf(fnd)                                      
            z = z + oecoord(3,fnd)*sf(fnd)                                      
         end if                                                                 
         dwdt = dwdt + dwidt*sf(fnd)                                            
      end do                                                                    
      if( debug ) write(*,9040) dwdt                                            
c                                                                               
c     determine change in deformation wrt flow direction                        
c                                                                               
      call derivs( etype, xi, eta, zeta, dsf(1,1), dsf(1,2), dsf(1,3) )         
      call eqldjb( dsf, ecoord, nnode, jacob, jacobi, det, ierr )               
      dwdf = zero                                                               
c                                                                               
c        get the displacements for each node on face and then                   
c        take the dot product with org. face vector                             
c                                                                               
      do node = 1,nfnode                                                        
         fnd = fnodes(node)                                                     
         du1 = u(3*elem_nodes(fnd)-2)                                           
         du2 = u(3*elem_nodes(fnd)-1)                                           
         du3 = u(3*elem_nodes(fnd)  )                                           
         dwi = du1*lvec(1) + du2*lvec(2) + du3*lvec(3)                          
         if( debug ) write(*,9050) elem_nodes(fnd),                             
     &        du1,du2,du3,dwi                                                   
c          get the derivative of shape functions wrt real coordinates           
         dsidx = dsf(fnd,1)*jacobi(1,1) + dsf(fnd,2)*jacobi(1,2)                
     &        +  dsf(fnd,3)*jacobi(1,3)                                         
         dsidy = dsf(fnd,1)*jacobi(2,1) + dsf(fnd,2)*jacobi(2,2)                
     &        +  dsf(fnd,3)*jacobi(2,3)                                         
         dsidz = dsf(fnd,1)*jacobi(3,1) + dsf(fnd,2)*jacobi(3,2)                
     &        +  dsf(fnd,3)*jacobi(3,3)                                         
c          get the derivatives shape function wrt flow direction                
c          and sum                                                              
         dsidf = dsidx*fdirc(1) + dsidy*fdirc(2) + dsidz*fdirc(3)               
         if( debug ) write(*,9060) dsidx,dsidy,dsidz,dsidf                      
         dwdf  = dwdf + dwi*dsidf                                               
      end do                                                                    
      if( debug ) write(*,9070) dwdf                                            
c                                                                               
c     compute piston theory flow pressure                                       
c                                                                               
      fpv = one/u3*dwdt + dwdf                                                  
      g1  = gam*m3                                                              
      g2  = (gam+one)/four*m3                                                   
      g3  = (gam+one)/twlv*m3*m3                                                
      nv1 = g1*fpv*( one + g2*fpv + g3*fpv*fpv )                                
c                                                                               
      face_intens = p3*(one+nv1)                                                
      if( debug ) write(*,9080) face_intens                                     
      if( debug ) write(*,9090) x, y, z                                         
c                                                                               
      if( debug ) write(*,*) '>> Leaving piston_face_intens'                    
c                                                                               
 9010 format(4x,'(xi,eta,zeta)            ',5x,3f10.6 )                         
 9020 format(4x,'lvec                     ',5x,3f10.6 )                         
 9030 format(4x,'node, vel (dx,dy,dz,dwdt)',i5,4f10.6)                          
 9040 format(4x,'value of dwdt            ',5x,f10.6)                           
 9050 format(4x,'node, dis (dx,dy,dz,dw)  ',i5,4f10.6)                          
 9060 format(4x,'dsf wrt (x,y,z,f)        ',5x,4f10.6)                          
 9070 format(4x,'value of dwdf            ',5x,f10.6)                           
 9080 format(4x,'face intensity           ',5x,e13.6)                           
 9090 format(4x,'at location (x,y,z)      ',5x,3f10.6)                          
      end                                                                       
c                                                                               
c ****************************************************************              
c *                                                              *              
c *     center coordinates of face on an element                 *              
c *                                                              *              
c ****************************************************************              
c                                                                               
      subroutine ccfe( xi, eta, zeta, face, tet_elem )                          
c                                                                               
      implicit none                                                             
c                                                                               
c     parameter declarations                                                    
c     ----------------------                                                    
c                                                                               
      integer face                                                              
      logical tet_elem                                                          
      double precision                                                          
     &     xi, eta, zeta                                                        
c                                                                               
c     declare local variables                                                   
c     -----------------------                                                   
c                                                                               
      double precision                                                          
     &     hex_face_cent(3,6), tet_face_cent(3,4)                               
c                                                                               
c     Set xi, eta, and zeta for hex element face centers                        
c                                                                               
      data hex_face_cent                                                        
     &     / -1.0,  0.0,  0.0,                                                  
     &        1.0,  0.0,  0.0,                                                  
     &        0.0, -1.0,  0.0,                                                  
     &        0.0,  1.0,  0.0,                                                  
     &        0.0,  0.0, -1.0,                                                  
     &        0.0,  0.0,  1.0 /                                                 
c                                                                               
c     Set s2, s3, and s4 for tet element face centers                           
c                                                                               
      data tet_face_cent                                                        
     &     / 0.3333333333, 0.3333333333, 0.0,                                   
     &       0.3333333333, 0.0,          0.3333333333,                          
     &       0.3333333333, 0.3333333333, 0.3333333333,                          
     &       0.0,          0.3333333333, 0.3333333333 /                         
c                                                                               
      if ( tet_elem ) then                                                      
         xi   = tet_face_cent(1,face)                                           
         eta  = tet_face_cent(2,face)                                           
         zeta = tet_face_cent(3,face)                                           
      else                                                                      
         xi   = hex_face_cent(1,face)                                           
         eta  = hex_face_cent(2,face)                                           
         zeta = hex_face_cent(3,face)                                           
      end if                                                                    
c                                                                               
      end                                                                       
c                                                                               
c ****************************************************************              
c *                                                              *              
c *     normal vector for a hex element                          *              
c *                                                              *              
c ****************************************************************              
c                                                                               
      subroutine hex_face_nvec( xi, eta, zeta, ecoord, lvec,                    
     &     element, etype, nnode, face )                                        
c                                                                               
      implicit none                                                             
c                                                                               
c     parameter declarations                                                    
c     ----------------------                                                    
c                                                                               
      integer element, etype, nnode, face                                       
      double precision                                                          
     &     xi, eta, zeta, ecoord(3,*), lvec(3)                                  
c                                                                               
c     declare local variables                                                   
c     -----------------------                                                   
c                                                                               
      integer ierr                                                              
      logical bad, debug                                                        
      double precision                                                          
     &     dsf(32,3), jacob(3,3), jacobi(3,3), det                              
      data debug / .false. /                                                    
      if(debug) write(*,*) '>> In hex_face_nvec'                                
c                                                                               
c              1) evaluate derivatives of shape functions                       
c                 at the node;                                                  
c              2) compute the jacobian at the node;                             
c              3) select vectors veca,vecb that are in the tangent              
c                 plane to the face at the node.  rows of the                   
c                 jacobian are components of the vectors.                       
c              4) construct an inward normal vector to the face                 
c                 (veca x vecb) and normalize to unit length                    
c                 thus producing direction cosines;                             
c              5) reverse the sign of the normal vector to point out            
c                                                                               
      call derivs( etype, xi, eta, zeta, dsf(1,1), dsf(1,2), dsf(1,3) )         
      call eqldjb( dsf, ecoord, nnode, jacob, jacobi, det, ierr )               
      call eqnrmvh( face, jacob, lvec, bad, debug )                             
c                                                                               
      lvec(1) = -lvec(1)                                                        
      lvec(2) = -lvec(2)                                                        
      lvec(3) = -lvec(3)                                                        
      if(debug) write(*,*) '>> Leaving hex_face_nvec'                           
c                                                                               
      end                                                                       
c                                                                               
c ****************************************************************              
c *                                                              *              
c *     normal vector for a tet element                          *              
c *                                                              *              
c ****************************************************************              
c                                                                               
      subroutine tet_face_nvec( xi, eta, zeta, ecoord, nrmvec,                  
     &        element, etype, nnode, face )                                     
c                                                                               
      implicit none                                                             
c                                                                               
c     parameter declarations                                                    
c     ----------------------                                                    
c                                                                               
      integer element, etype, nnode, face                                       
      logical tet_elem                                                          
      double precision                                                          
     &     xi, eta, zeta, ecoord(3,*), nrmvec(3)                                
c                                                                               
c     declare local variables                                                   
c     -----------------------                                                   
c                                                                               
      double precision                                                          
     &     dsf(32,3), temp_dsf(32,3), sf(32), temp_sf(32),                      
     &     darea, mag1, mag2, mag_max, vec1(3), vec2(3), rlen,                  
     &     fcoor_tri6(3,6), fcoor_tri3(3,3), zero                               
      logical debug                                                             
      integer fnodes(10), index, nfnode                                         
      data zero, debug / 0.0, .false. /                                         
c                                                                               
c                                                                               
c                    use shape function derivatives for triangular              
c                    interface elements for tractions on the tet                
c                    faces. this is done to avoid taking derivatives            
c                    with respect to the dependent natural                      
c                    coordinate.                                                
c                                                                               
c                    because we're working in 2D, the isoparametric             
c                    coordinates will be the same regardless of the             
c                    face that is loaded.                                       
c                                                                               
c                                                                               
      data fcoor_tri6                                                           
     & / 0.0, 0.0, 1.0,                                                         
     &   0.5, 0.0, 0.5,                                                         
     &   1.0, 0.0, 0.0,                                                         
     &   0.5, 0.5, 0.0,                                                         
     &   0.0, 1.0, 0.0,                                                         
     &   0.0, 0.5, 0.5 /                                                        
c                                                                               
      data fcoor_tri3                                                           
     & / 0.0, 0.0, 1.0,                                                         
     &   1.0, 0.0, 0.0,                                                         
     &   0.0, 1.0, 0.0 /                                                        
c                                                                               
      if(debug) write(*,*) '>> In tet_face_nvec'                                
c                                                                               
c                     determine if linear or quadratic tet                      
c                                                                               
      if ( etype .eq. 6 )  index = 15                                           
      if ( etype .eq. 13 ) index = 14                                           
c                                                                               
c                    get the tet nodes on the face for mapping                  
c                    derivatives to the correct nodal coordinates.              
c                                                                               
      call tet_get_nodes( fnodes, etype, face, nfnode )                         
c                                                                               
c              1) evaluate derivatives of shape functions                       
c                 at the node;                                                  
c              2) build vec1, vec2                                              
c              3) construct an inward normal vector to the face                 
c                 (vec1 x vec2) and normalize to unit length                    
c                 thus producing direction cosines;                             
c                                                                               
      call derivs( index, xi, eta, zeta, temp_dsf(1,1),                         
     &     temp_dsf(1,2), temp_dsf(1,3) )                                       
      call map_tet_tri( etype, nfnode, fnodes, sf, temp_sf,                     
     &     dsf, temp_dsf )                                                      
      call tet_darea( dsf, nfnode, fnodes, ecoord, darea,                       
     &     vec1, vec2 )                                                         
c                                                                               
c             generate unit normal to surface at node using                     
c             cross product of tangent plane vectors.                           
c             then compute x,y,z components of load.                            
c                                                                               
        nrmvec(1) = vec1(2)*vec2(3) - vec2(2)*vec1(3)                           
        nrmvec(2) = vec2(1)*vec1(3) - vec1(1)*vec2(3)                           
        nrmvec(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)                           
        rlen = sqrt( nrmvec(1)*nrmvec(1) + nrmvec(2)*nrmvec(2) +                
     &               nrmvec(3)*nrmvec(3) )                                      
        nrmvec(1) = nrmvec(1) / rlen                                            
        nrmvec(2) = nrmvec(2) / rlen                                            
        nrmvec(3) = nrmvec(3) / rlen                                            
c                                                                               
      if(debug) write(*,*) '>> Leaving tet_face_nvec'                           
c                                                                               
      end                                                                       
                                                                                
