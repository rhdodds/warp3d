c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine cohes_rot_mat                     *          
c     *                                                              *          
c     *                    written by : aroy                         *          
c     *                 last modified : 02/10/13 rhd                 *          
c     *                                                              *          
c     *     compute rotation matrix bigR[span,3,3] that rotates      *          
c     *     global vector quantities into the local orthogonal       *          
c     *     system at the parametric center of interface elements.   *          
c     *     the local system x-y axes lie in the tangent plane.      *          
c     *     local z is normal to tangent plane                       *          
c     *                                                              *          
c     *     geometrically nonlinear elements: before entering this   *          
c     *     subroutine, the ce vector has been modified such that    *          
c     *     the first nnode/2 points describe the reference surface. *          
c     *     see cohes_ref_surface.                                   *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine cohes_rot_mat(span, felem, nnode, etype, ce, bigR )            
      implicit none                                                             
c                                                                               
      include 'param_def' 
c                                                            
c                  parameters                                                   
c                                                                               
      integer :: span, felem, nnode, etype                                         
      double precision :: ce(mxvl,mxecor), bigR(mxvl,3,3)                                          
c                                                                               
c                   local data                                                  
c                                                                               
      integer :: i, j, imxvl, ispan, innode                                        
      double precision ::                                                          
     & x21,y21,z21, x31,y31,z31, d1, d2, d3, a2, b2, c2,                        
     & a3, b3, c3, x1, y1, z1,                                                  
     & e123_to_ref(3,3),                                           
     & ce_rotated(mxvl,mxecor), len, nzeta(mxndel),                             
     & neta(mxndel), nxi(mxndel), jac(mxvl,3,3), dj(mxvl),                      
     & g_to_e123(mxvl,3,3)      
      double precision, parameter :: zero1 = 0.0d0, zero2 = 0.0d0,
     &                               zero3 = 0.0d0,
     &                               third1 = 1.0d0/3.0d0,
     &                               third2 = 1.0d0/3.0d0                                            
c                                                                               
      logical :: quad_interface, tri_interface    
      logical, parameter :: local_debug = .false.                    
c                                                                               
c                                                                               
c               for geonl these are deformed coordinates of the                 
c               reference surface (most often the deformed                      
c               middle surface. see notes above and in routine                  
c               cohes_ref_surface.                                              
c                                                                               
c               Note: that first nnode/2 terms in                               
c               ce vector describe the reference surface.)                      
c                                                                               
c               shift global element coords so element node 1                   
c               is at origin. really not necessary but makes                    
c               potential debugging easier.                                     
c                                                                               
      x1 = ce(1,1); y1 = ce(1,nnode+1); z1 = ce(1,2*nnode+1)                    
      do j = 1, nnode                                                           
!DIR$ IVDEP                                                                     
         do i = 1, span                                                         
           ce(i,j) = ce(i,j) - x1                                               
           ce(i,nnode+j) = ce(i,nnode+j) - y1                                   
           ce(i,2*nnode+j) = ce(i,2*nnode+j) - z1                               
         end do                                                                 
      end do                                                                    
c                                                                               
c               an element local system based on element nodes                  
c               1, 2,3. local z is normal to plane containing                   
c               element nodes 1, 2, 3. Then rotate the shifted                  
c               global coordinates into this 1-2-3 system.                      
c                                                                               
!DIR$ IVDEP                                                                     
      do i = 1, span                                                            
        x21 = ce(i,2) - ce(i,1)                                                 
        y21 = ce(i,nnode+2) - ce(i,nnode+1)                                     
        z21 = ce(i,2*nnode+2) - ce(i,2*nnode+1)                                 
        x31 = ce(i,3) - ce(i,1)                                                 
        y31 = ce(i,nnode+3) - ce(i,nnode+1)                                     
        z31 = ce(i,2*nnode+3) - ce(i,2*nnode+1)                                 
        d1  = sqrt(x21**2+y21**2+z21**2)                                        
c                                                                               
        g_to_e123(i,1,1) = x21/d1                                               
        g_to_e123(i,1,2) = y21/d1                                               
        g_to_e123(i,1,3) = z21/d1                                               
c                                                                               
        a3 = y21*z31 - y31*z21                                                  
        b3 = x31*z21 - x21*z31                                                  
        c3 = x21*y31 - x31*y21                                                  
        d3 = sqrt(a3**2+b3**2+c3**2)                                            
c                                                                               
        g_to_e123(i,3,1) = a3/d3                                                
        g_to_e123(i,3,2) = b3/d3                                                
        g_to_e123(i,3,3) = c3/d3                                                
c                                                                               
        a2 = g_to_e123(i,3,2)*g_to_e123(i,1,3) -                                
     &       g_to_e123(i,1,2)*g_to_e123(i,3,3)                                  
        b2 = g_to_e123(i,1,1)*g_to_e123(i,3,3) -                                
     &       g_to_e123(i,3,1)*g_to_e123(i,1,3)                                  
        c2 = g_to_e123(i,3,1)*g_to_e123(i,1,2) -                                
     &       g_to_e123(i,1,1)*g_to_e123(i,3,2)                                  
        d2 = sqrt(a2**2+b2**2+c2**2)                                            
c                                                                               
        g_to_e123(i,2,1) = a2/d2                                                
        g_to_e123(i,2,2) = b2/d2                                                
        g_to_e123(i,2,3) = c2/d2                                                
      end do                                                                    
c                                                                               
      imxvl = mxvl; ispan = span; innode = nnode         
      call rotate_cohes_var( imxvl, ispan, innode, g_to_e123,                   
     &                       ce, ce_rotated )                                   
c                                                                               
c              we need the orthogonal system (s1, s2, n) at the                 
c              parametric center of the reference surface, where                
c              s1 and s2 lie in the tangent plane to the reference              
c              surface.                                                         
c                                                                               
c              compute the orthogonal rotation R[3,3] that takes                
c              vector values from the orthogonal 1-2-3 system                   
c              constructed above to s1, s2, n.                                  
c                                                                               
c              for quadrilateral interface elements, use (xi=eta=0)             
c              for triangular interface elements, use s1=s2=1/3                 
c              (s3 is then also 1/3).                                           
c                                                                               
c              use rows 1,2 in computed [J] at parametric center                
c              to define vectors (generally nonorthogonal) in tangent           
c              plane. take cross product to get n-direction.                    
c              use additional cross-products to get orthogonal                  
c              vectors in tangent plane.                                        
c                                                                               
c              final coordinate transformation from global XYZ                  
c              system to interface system on and normal to tangent              
c              plane [bigR] = [e123_to_ref] * [g_to_e123]                       
c                                                                               
c              get derivates of shape functions at parametric center            
c              and coordinate jacobian                                          
c                                                                               
       quad_interface = etype .eq. 12                                           
       tri_interface  = etype .eq. 14  .or. etype .eq. 15                       
c                             
       nxi = zero1; neta = zero1; nzeta = zero1  ! prevents unint failure                                                 
       if( quad_interface )                                                     
     &   call derivs( etype, zero1, zero2, zero3, nxi, neta, nzeta )               
c                                                                               
       if( tri_interface )                                                      
     &   call derivs( etype, third1, third2, zero1, nxi, neta, nzeta )             
c                
       call jacob_cohes( etype, imxvl, ispan, innode, ce_rotated,               
     &                   nxi, neta, nzeta, dj, jac, 2 )                         
c                                                                               
!DIR$ IVDEP                                                                     
       do i = 1, span                                                           
c                                                                               
        d1 = sqrt(jac(i,1,1)*jac(i,1,1) + jac(i,1,2)*jac(i,1,2)                 
     &             +  jac(i,1,3)*jac(i,1,3) )                                   
c                                                                               
        e123_to_ref(1,1) = jac(i,1,1)/d1                                        
        e123_to_ref(1,2) = jac(i,1,2)/d1                                        
        e123_to_ref(1,3) = jac(i,1,3)/d1                                        
c                                                                               
        a3 = jac(i,1,2)*jac(i,2,3) - jac(i,2,2)*jac(i,1,3)                      
        b3 = jac(i,2,1)*jac(i,1,3) - jac(i,1,1)*jac(i,2,3)                      
        c3 = jac(i,1,1)*jac(i,2,2) - jac(i,2,1)*jac(i,1,2)                      
        d3 = sqrt(a3*a3+b3*b3+c3*c3)                                            
c                                                                               
        e123_to_ref(3,1) = a3/d3                                                
        e123_to_ref(3,2) = b3/d3                                                
        e123_to_ref(3,3) = c3/d3                                                
c                                                                               
        a2 = e123_to_ref(3,2)*e123_to_ref(1,3) -                                
     &       e123_to_ref(1,2)*e123_to_ref(3,3)                                  
        b2 = e123_to_ref(1,1)*e123_to_ref(3,3) -                                
     &       e123_to_ref(3,1)*e123_to_ref(1,3)                                  
        c2 = e123_to_ref(3,1)*e123_to_ref(1,2) -                                
     &       e123_to_ref(1,1)*e123_to_ref(3,2)                                  
        d2 = sqrt(a2*a2+b2*b2+c2*c2)                                            
c                                                                               
        e123_to_ref(2,1) = a2/d2                                                
        e123_to_ref(2,2) = b2/d2                                                
        e123_to_ref(2,3) = c2/d2                                                
c                                                                               
        bigR(i,1,1) = e123_to_ref(1,1)*g_to_e123(i,1,1) +                       
     &                e123_to_ref(1,2)*g_to_e123(i,2,1) +                       
     &                e123_to_ref(1,3)*g_to_e123(i,3,1)                         
        bigR(i,1,2) = e123_to_ref(1,1)*g_to_e123(i,1,2) +                       
     &                e123_to_ref(1,2)*g_to_e123(i,2,2) +                       
     &                e123_to_ref(1,3)*g_to_e123(i,3,2)                         
        bigR(i,1,3) = e123_to_ref(1,1)*g_to_e123(i,1,3) +                       
     &                e123_to_ref(1,2)*g_to_e123(i,2,3) +                       
     &                e123_to_ref(1,3)*g_to_e123(i,3,3)                         
        bigR(i,2,1) = e123_to_ref(2,1)*g_to_e123(i,1,1) +                       
     &                e123_to_ref(2,2)*g_to_e123(i,2,1) +                       
     &                e123_to_ref(2,3)*g_to_e123(i,3,1)                         
        bigR(i,2,2) = e123_to_ref(2,1)*g_to_e123(i,1,2) +                       
     &                e123_to_ref(2,2)*g_to_e123(i,2,2) +                       
     &                e123_to_ref(2,3)*g_to_e123(i,3,2)                         
        bigR(i,2,3) = e123_to_ref(2,1)*g_to_e123(i,1,3) +                       
     &                e123_to_ref(2,2)*g_to_e123(i,2,3) +                       
     &                e123_to_ref(2,3)*g_to_e123(i,3,3)                         
        bigR(i,3,1) = e123_to_ref(3,1)*g_to_e123(i,1,1) +                       
     &                e123_to_ref(3,2)*g_to_e123(i,2,1) +                       
     &                e123_to_ref(3,3)*g_to_e123(i,3,1)                         
        bigR(i,3,2) = e123_to_ref(3,1)*g_to_e123(i,1,2) +                       
     &                e123_to_ref(3,2)*g_to_e123(i,2,2) +                       
     &                e123_to_ref(3,3)*g_to_e123(i,3,2)                         
        bigR(i,3,3) = e123_to_ref(3,1)*g_to_e123(i,1,3) +                       
     &                e123_to_ref(3,2)*g_to_e123(i,2,3) +                       
     &                e123_to_ref(3,3)*g_to_e123(i,3,3)                         
      end do                                                                    
c                                                                               
      if( local_debug ) then                                                    
        do i = 1, span                                                          
          write(*,100) bigR(i,1,1), bigR(i,1,2), bigR(i,1,3)                    
          write(*,100) bigR(i,2,1), bigR(i,2,2), bigR(i,2,3)                    
          write(*,100) bigR(i,3,1), bigR(i,3,2), bigR(i,3,3)                    
        end do                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
100   format('Lambda ->',3(2X,f10.6))                                           
101   format('R ->',3(2X,f10.6))                                                
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine rotate_cohes_var                      *          
c     *                                                              *          
c     *                    written by : aroy                         *          
c     *                    modified: 2/10/13 rhd                     *          
c     *                                                              *          
c     *           rotates multiple vectors stored in packed form     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine rotate_cohes_var( mxvl, span, nnode, rot, var, rvar )          
      implicit none                                                             
c                                                                               
c                  parameter declarations                                       
c                                                                               
      integer mxvl, span, nnode                                                 
      double precision                                                          
     & var(mxvl,*), rvar(mxvl,*), rot(mxvl,3,3)                                 
c                                                                               
c                  local declarations                                           
c                                                                               
      integer i, j, jj, k                                                       
c                                                                               
c                                                                               
c              var has nnode * 3 terms. 1-nnode are x components.               
c              then the nnode y compoenents, then the nnode z                   
c              components.                                                      
c                                                                               
      do i = 1, nnode                                                           
         j = i + nnode                                                          
        jj = i + 2*nnode                                                        
!DIR$ IVDEP                                                                     
          do k = 1, span                                                        
             rvar(k,i) = var(k,i)*rot(k,1,1) + var(k,j)*rot(k,1,2)              
     &                   + var(k,jj)*rot(k,1,3)                                 
             rvar(k,j) = var(k,i)*rot(k,2,1) + var(k,j)*rot(k,2,2)              
     &                   + var(k,jj)*rot(k,2,3)                                 
             rvar(k,jj)= var(k,i)*rot(k,3,1) + var(k,j)*rot(k,3,2)              
     &                   + var(k,jj)*rot(k,3,3)                                 
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine jacob_cohes                  *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                       last modified: 2/10/13 rhd             *          
c     *                                                              *          
c     *     compute 3x3 coordinate jacobian matrix & det at a point  *          
c     *     for block of interface elements                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c       uses the formula: [J] = [SFD] * [coords]                                
c         where   [J] = jacobian matrix                                         
c               [SFD] = matrix of shape function derivatives                    
c                       (passed in as nxi, neta, nzeta)                         
c                [coords] = element nodal coordinates                           
c                                                                               
      subroutine jacob_cohes( etype, mxvl, span, nnode, coords,                 
     &                        nxi, neta, nzeta, dj, jac, action )               
      implicit none                                                             
c                                                                               
c                  parameters                                                   
c                                                                               
      integer etype, mxvl, span, nnode, action                                  
c                                                                               
      double precision                                                          
     & coords(mxvl,*), nzeta(*), neta(*), nxi(*), jac(mxvl,3,3), dj(*)          
c                                                                               
c                  locals                                                       
c                                                                               
      integer i, j                                                              
      double precision                                                          
     & zero, j1, j2, j3                                                         
c                                                                               
      data zero / 0.0d0 /                                                       
c                                                                               
!DIR$ IVDEP                                                                     
      do i = 1, span                                                            
         jac(i,1:3,1:3) = zero                                                  
         dj(i)          = zero                                                  
      end do                                                                    
c                                                                               
      do j = 1, nnode                                                           
!DIR$ IVDEP                                                                     
         do i = 1, span                                                         
             jac(i,1,1)= jac(i,1,1)+nxi(j)*coords(i,j)                          
             jac(i,1,2)= jac(i,1,2)+nxi(j)*coords(i,nnode+j)                    
             jac(i,1,3)= jac(i,1,3)+nxi(j)*coords(i,2*nnode+j)                  
             jac(i,2,1)= jac(i,2,1)+neta(j)*coords(i,j)                         
             jac(i,2,2)= jac(i,2,2)+neta(j)*coords(i,nnode+j)                   
             jac(i,2,3)= jac(i,2,3)+neta(j)*coords(i,2*nnode+j)                 
             jac(i,3,1)= jac(i,3,1)+nzeta(j)*coords(i,j)                        
             jac(i,3,2)= jac(i,3,2)+nzeta(j)*coords(i,nnode+j)                  
             jac(i,3,3)= jac(i,3,3)+nzeta(j)*coords(i,2*nnode+j)                
         end do                                                                 
      end do                                                                    
      if ( action .eq. 2 ) return                                               
c                                                                               
c           determinate of the jacobian matrix                                  
c                                                                               
!DIR$ IVDEP                                                                     
      do i = 1, span                                                            
          j1= jac(i,2,2)*jac(i,3,3)-jac(i,2,3)*jac(i,3,2)                       
          j2= jac(i,2,1)*jac(i,3,3)-jac(i,2,3)*jac(i,3,1)                       
          j3= jac(i,2,1)*jac(i,3,2)-jac(i,2,2)*jac(i,3,1)                       
          dj(i)= jac(i,1,1)*j1-jac(i,1,2)*j2+jac(i,1,3)*j3                      
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine cohes_ref_surface                  *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                                                              *          
c     *                   last modified : 12/11/10 rhd comments      *          
c     *                                                              *          
c     *     for the geometrically nonlinear cohesive elements,       *          
c     *     this subroutine computes the coordinates of the          *          
c     *     reference surface for symmetric (top and bottom) and     *          
c     *     non-specialized (middle) configurations.                 *          
c     *     ce_upd is modified here to store the reference           *          
c     *     surface coordinates.                                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c          for cohesive elements, the strain-displacement matrix [B]            
c          is given by [B] = [R][L][N]. for geometrically nonlinear             
c          elements, the [R] matrix is computed using a                         
c          reference surface specified by the user. this                        
c          routine computes the reference surface coordinates                   
c          and stores them in ce_upd.                                           
c                                                                               
c          data structure in cu_upd:                                            
c            on entry, cu_upd stores the current                                
c           (ie deformed) coords of element nodes in following form:            
c                                                                               
c   ce_upd = |        |                                                         
c            |   X    | \                                                       
c            | coor   |   \|    | node   1   of bottom surface      |           
c            |        |   -'    | node   2   of bottom surface      |           
c            |--------|         |    .............                  |           
c            |        |         |    .............                  |           
c            |   Y    |--->     | node 'nnode/2' of bottom surface  |           
c            | coor   |         | node    1    of top    surface    |           
c            |        |         | node    2    of top    surface    |           
c            |--------|   -.    |    .............                  |           
c            |        |   /|    |    .............                  |           
c            |   Z    |  /      | node 'nnode/2' of top surface     |           
c            | coor   | /                                                       
c            |        |                                                         
c                                                                               
c           on exit from this routine, the first nnode/2 locations in           
c           ce_upd store the coordinates of reference surface nodes.            
c           the remaining nnode/2 locations store the reflected surface.        
c           (the reflected surface is so called because at a later              
c           point its coordinates will be reflected upon the                    
c           reference surface to take care of Poisson effect of bulk            
c           material.)                                                          
c                                                                               
c           the definition of reflected surface and the exact structure         
c           of ce_upd on exit depends on reference surface (top, bottom         
c           or middle). these are explained below in each case.                 
c                                                                               
c           for symmetric models, the user must declare the "reference"         
c           surface (top or bottom) as the surface that is connected t          
c           o the solid elements. nodes on the reflected                        
c           surface will be given deformed coordinates                          
c           of the reference surface.                                           
c                                                                               
      subroutine cohes_ref_surface( span, mxvl, mxecor, surf, nnode,            
     &                              totdof, ce, ce_upd, dj )                    
c                                                                               
c              parameters                                                       
c                                                                               
      implicit none                                                             
      integer span, mxvl, mxecor, surf, nnode, totdof                           
      double precision                                                          
     &    ce_upd(mxvl,mxecor), ce(mxvl,mxecor), dj(mxvl)                        
c                                                                               
c              locals. one automatic array.                                     
c                                                                               
      integer i, j, k, shift, xb, xt, yb, yt, zb, zt                            
      double precision                                                          
     &  ce_refsurf(span,mxecor), zero, half                                     
      data zero, half / 0.0d0, 0.5d0 /                                          
      logical local_debug                                                       
      data local_debug /.false./                                                
c                                                                               
      shift = nnode/2                                                           
      xb = 0                                                                    
      xt = xb + shift                                                           
      yb = nnode                                                                
      yt = yb + shift                                                           
      zb = 2*nnode                                                              
      zt = zb + shift                                                           
c                                                                               
      if( local_debug ) then                                                    
        do i = 1, span                                                          
          write(*,*) 'surf = ',surf                                             
        end do                                                                  
      end if                                                                    
c                                                                               
      select case( surf )                                                       
c     ===================                                                       
c                                                                               
      case( 1 )                                                                 
c                                                                               
c                 *** reference surface - top  ***                              
c                                                                               
c      reference surface - top.                                                 
c      reflected surface - bottom.                                              
c      change ce_upd such that on exit it has following structure.              
c                                                                               
c   ce_upd = |        |                                                         
c            |   X    | \                                                       
c            | coor   |   \|   | node   1   of reference (top) surface|         
c            |        |   -'   | node   2   of reference (top) surface|         
c            |--------|        |    .............                     |         
c            |        |        |    .............                     |         
c            |   Y    |---->   |node 'nnode/2' of reference (top) surf|         
c            | coor   |        | node    1    of reflect (bot) surface|         
c            |        |        | node    2    of reflect (bot) surface|         
c            |--------|   -.   |    .............                     |         
c            |        |   /|   |    .............                     |         
c            |   Z    |  /     |node 'nnode/2' of reflect (bot) surf  |         
c            | coor   | /                                                       
c            |        |                                                         
c                                                                               
      do j = 1,shift                                                            
!DIR$ IVDEP                                                                     
        do k = 1, span                                                          
            ce_refsurf(k,xb+j) = ce_upd(k,xt+j)                                 
            ce_refsurf(k,xt+j) = ce_upd(k,xb+j)                                 
            ce_refsurf(k,yb+j) = ce_upd(k,yt+j)                                 
            ce_refsurf(k,yt+j) = ce_upd(k,yb+j)                                 
            ce_refsurf(k,zb+j) = ce_upd(k,zt+j)                                 
            ce_refsurf(k,zt+j) = ce_upd(k,zb+j)                                 
         end do                                                                 
      end do                                                                    
c                                                                               
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
        ce_upd(1:span,j) = ce_refsurf(1:span,j)                                 
      end do                                                                    
c                                                                               
      case( 2 )                                                                 
c                                                                               
c               *** reference surface - middle  ***                             
c                                                                               
c      if penetration has occurred, define reference surface as                 
c      the original (undeformed) (top and bottom) surface of                    
c      cohesive element. the reflect surface is then defined by                 
c      0.5(deformed top+ deformed bottom).                                      
c                                                                               
c      if there is no penetration, define both reference and                    
c      reflect surface by 0.5(deformed top+ deformed bottom).                   
c                                                                               
c                                                                               
c      change ce_upd such that on exit it has following structure.              
c                                                                               
c   ce_upd = |        |                                                         
c            |   X    | \                                                       
c            | coor   |   \|   | node   1   of reference surface|               
c            |        |   -'   | node   2   of reference surface|               
c            |--------|        |    .............               |               
c            |        |        |    .............               |               
c            |   Y    |---->   |node 'nnode/2' of reference surf|               
c            | coor   |        | node    1    of reflect surface|               
c            |        |        | node    2    of reflect surface|               
c            |--------|   -.   |    .............               |               
c            |        |   /|   |    .............               |               
c            |   Z    |  /     |node 'nnode/2' of reflect surf  |               
c            | coor   | /                                                       
c            |        |                                                         
c                                                                               
c                                                                               
c                                                                               
      do k = 1, span                                                            
      if( dj(k) .lt. zero ) then  ! interpenetration                            
             do j = 1, shift                                                    
               ce_refsurf(k,xb+j) =  ce(k,xb+j)                                 
               ce_refsurf(k,xt+j) =                                             
     &         half*( ce_upd(k,xb+j) + ce_upd(k,xt+j) )                         
               ce_refsurf(k,yb+j) =  ce(k,yb+j)                                 
               ce_refsurf(k,yt+j) =                                             
     &         half*( ce_upd(k,yb+j) + ce_upd(k,yt+j) )                         
               ce_refsurf(k,zb+j) =  ce(k,zb+j)                                 
               ce_refsurf(k,zt+j) =                                             
     &         half*( ce_upd(k,zb+j) + ce_upd(k,zt+j) )                         
            end do                                                              
      else                                                                      
             do j = 1, shift                                                    
               ce_refsurf(k,xb+j) =                                             
     &         half*( ce_upd(k,xb+j) + ce_upd(k,xt+j) )                         
               ce_refsurf(k,xt+j) =  ce_refsurf(k,xb+j)                         
               ce_refsurf(k,yb+j) =                                             
     &         half*( ce_upd(k,yb+j) + ce_upd(k,yt+j) )                         
               ce_refsurf(k,yt+j) =  ce_refsurf(k,yb+j)                         
               ce_refsurf(k,zb+j) =                                             
     &         half*( ce_upd(k,zb+j) + ce_upd(k,zt+j) )                         
               ce_refsurf(k,zt+j) =  ce_refsurf(k,zb+j)                         
            end do                                                              
        end if                                                                  
      end do                                                                    
c                                                                               
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
        ce_upd(1:span,j) = ce_refsurf(1:span,j)                                 
      end do                                                                    
c                                                                               
      case( 3 )                                                                 
c                                                                               
c                *** reference surface - bottom  ***                            
c                                                                               
c      reference surface - bottom.                                              
c      reflect surface - top.                                                   
c      ce_upd already has following structure. no change required.              
c                                                                               
c   ce_upd = |        |                                                         
c            |   X    | \                                                       
c            | coor   |   \|   | node   1   of reference (bot) surface|         
c            |        |   -'   | node   2   of reference (bot) surface|         
c            |--------|        |    .............                     |         
c            |        |        |    .............                     |         
c            |   Y    |---->   |node 'nnode/2' of reference (bot) surf|         
c            | coor   |        | node    1    of reflect (top) surface|         
c            |        |        | node    2    of reflect (top) surface|         
c            |--------|   -.   |    .............                     |         
c            |        |   /|   |    .............                     |         
c            |   Z    |  /     |node 'nnode/2' of reflect (top) surf  |         
c            | coor   | /                                                       
c            |        |                                                         
c                                                                               
c                                                                               
      continue                                                                  
c                                                                               
      end select                                                                
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine cohes_mirror_refsurf               *          
c     *                                                              *          
c     *                  written by : aroy                           *          
c     *                  last modified: 2/10/13  rhd                 *          
c     *                                                              *          
c     *     equates the updated nodal coordinates of the reference   *          
c     *     surface to those of the reflected surface. this          *          
c     *     is done only for a symmetric model and the reference     *          
c     *     surface is either bottom or top. (In such cases, all     *          
c     *     dof of the nodes attached to the reference plane must    *          
c     *     be constrained)                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c           the reflected surface deforms with the bulk material.               
c           thus its updated coordinates shows the Poisson's ratio              
c           effect. however, the nodes attached to reference surface            
c           have all dofs constrained. to obtain the true shape of the          
c           deformed cohesive element the inplane coordinates of                
c           reference surface should be equated to those of reflect             
c           surface. here we equate all coordinates of these two                
c           surfaces. since the interface cohesive element is                   
c           treated as a 2-D surface element (2x2 integration),                 
c           the equated normal coodinates will not introduce                    
c           any error. on the other hand, it makes the code                     
c           very simple :)                                                      
c                                                                               
c            NOTE: before entering this subroutine, the ce vector has           
c                been modified such that the first nnode/2 points               
c                describe the reference surface and the remaining               
c                nnode/2 points describe the reflect surface.                   
c                see SR cohes_ref_surface.                                      
c                                                                               
c            when the reference surface is middle this SR is a sham             
c            one since the reference and reflect surface and their              
c            updated coordinates are then the same by definition                
c            (see  SR cohes_ref_surface).                                       
c                                                                               
      subroutine cohes_mirror_refsurf( span, mxvl, totdof,                      
     &    nnode, ce )                                                           
      implicit none                                                             
c                                                                               
c                  parameters                                                   
c                                                                               
      integer span, mxvl, totdof, nnode                                         
      double precision                                                          
     & ce(mxvl,*)                                                               
c                                                                               
c                  locals                                                       
c                                                                               
      integer j, k, shift, xb, xt, yb, yt, zb, zt                               
c                                                                               
      shift = nnode/2                                                           
c                                                                               
      xb = 0                                                                    
      xt = xb + shift                                                           
      yb = nnode                                                                
      yt = yb + shift                                                           
      zb = 2*nnode                                                              
      zt = zb + shift                                                           
c                                                                               
      do k = 1, shift                                                           
       do j = 1, span                                                           
         ce(j,xb+k) = ce(j,xt+k)                                                
         ce(j,yb+k) = ce(j,yt+k)                                                
         ce(j,zb+k) = ce(j,zt+k)                                                
       end do                                                                   
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine chk_cohes_penetrate                   *          
c     *                                                              *          
c     *                    written by : aroy                         *          
c     *                last modified  : 02/10/13 rhd                 *          
c     *                                                              *          
c     *     support for subsequent check interpenetration of         *          
c     *     interface-cohesive elements. checked using jacobian of   *          
c     *     element with updated nodal coordinates. applicable only  *          
c     *     to cohesive zone models which are isotropic under shear. *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c           if penetration occurs, we need to set the                           
c           reference surface as the initial undeformed coordinates.            
c           see SR cohes_ref_surf. coords passed in will be deformed            
c           in current state.                                                   
c                                                                               
      subroutine chk_cohes_penetrate( span, mxvl, felem, mxndel,                
     &                                nnode, etype, coords, dj )                
      implicit none                                                             
c                                                                               
      integer span, mxvl, felem, mxndel, nnode, etype                           
      double precision                                                          
     &    coords(mxvl,*), dj(mxvl)                                              
c                                                                               
c              locally defined                                                  
c                                                                               
      integer  solid_type                                                       
      logical  inter_quad_8, inter_tri_6, inter_tri_12                          
      double precision                                                          
     &    zero, half, third                                                     
c                                                                               
c             automatic arrays                                                  
c                                                                               
      double precision                                                          
     & nzeta(mxndel), neta(mxndel), nxi(mxndel), jac(mxvl,3,3)                  
       data zero, half, third / 0.0d0, 0.5d0, 0.333333333333333333d0 /          
c                                                                               
c         element interpenetration is checked later by the                      
c           sign (+ve or -ve) of the determinant of jacobian.                   
c                                                                               
c         treat the interface elements as solids.                               
c           quadrilateral interface element -> brick                            
c           triangular interface element -> wedge                               
c                                                                               
c        first compute shape function derivatives at element                    
c        parametric center.                                                     
c         for bricks, xi=eta=zeta=0.                                            
c         for wedges, s1=s2=1/3, zeta=0.                                        
c        (see shapef.f and derivs.f for triangle natural                        
c         coordinates s1, s2, s3.)                                              
c                                                                               
c        for quadrilateral cohesive element inter_8,                            
c        call derivs with etype=2 (8-node brick).                               
c        treat triangular elements as wedges.                                   
c                                                                               
      inter_quad_8 = etype .eq. 12                                              
      inter_tri_6  = etype .eq. 14                                              
      inter_tri_12 = etype .eq. 15                                              
c                                                                               
      if( inter_quad_8 ) then                                                   
         solid_type = 2                                                         
         call derivs( solid_type, zero, zero, zero, nxi, neta, nzeta )          
      end if                                                                    
c                                                                               
      if(  inter_tri_6 )                                                        
     &    call deriv_trint6_wedge( etype, third, third, zero,                   
     &                             nxi, neta, nzeta)                            
                                                                                
      if( inter_tri_12 )                                                        
     &     call deriv_trint12_wedge( etype, third, third, zero,                 
     &                               nxi, neta, nzeta)                          
c                                                                               
c           compute coordinate jacobian                                         
c                                                                               
       call jacob_cohes( etype, mxvl, span, nnode, coords,                      
     &                   nxi, neta, nzeta, dj, jac, 1 )                         
       return                                                                   
       end                                                                      
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine deriv_trint6_wedge                 *          
c     *                                                              *          
c     *                    written by : sushovan                     *          
c     *                                 01/05/01                     *          
c     *                                                              *          
c     *     given a point in parametric coordinates, return the      *          
c     *     parametric derivatives for each node shape function,     *          
c     *     for the cohesive element trint6, treating it as wedge.   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c       use triangle natural coordinates (s1,s2,s3) for the wedge base.         
c       (see SR shape14 and deriv14 for triangle natural coordinates).          
c       zeta is along the height of the wedge. the two traingular faces         
c       are at zeta = -1 (nodes 1,2,3) and zeta = +1 (nodes 4,5,6).             
c                                                                               
c       the shape functions for the wedge:                                      
c         q1=0.5*s1*(1-zeta), q2=0.5*s2*(1-zeta), q3=0.5*s3*(1-zeta),           
c         q4=0.5*s1*(1+zeta), q5=0.5*s2*(1+zeta), q6=0.5*s3*(1+zeta).           
c                                                                               
c       s3 is dependent coordinate. replace s3=1-s1-s2 while                    
c       taking the derivatives.                                                 
c                                                                               
      subroutine deriv_trint6_wedge( etype, s1, s2, zeta, qs1, qs2,             
     &                               qzeta )                                    
      implicit none                                                             
c                                                                               
      integer etype                                                             
      double precision                                                          
     &     s1, s2, zeta, qs1(*), qs2(*), qzeta(*)                               
c                                                                               
c             local variables                                                   
c                                                                               
      double precision                                                          
     &     s3, zero, half, one                                                  
      data zero, half, one /0.0d0, 0.5d0, 1.0d0/                                
c                                                                               
      s3 = one - s1 -s2                                                         
c                                                                               
c                                                                               
c          dq_i/ds1, derivatives with respect to s1, the first                  
c          local coordinate : qs1(1) = dq1/ds1,                                 
c                             qs1(2) = dq2/ds1,                                 
c                             qs1(3) = dq3/ds1, etc.                            
c                                                                               
      qs1(1) = half*(1-zeta)                                                    
      qs1(2) = zero                                                             
      qs1(3) = -half*(1-zeta)                                                   
      qs1(4) = half*(1+zeta)                                                    
      qs1(5) = zero                                                             
      qs1(6) = -half*(1+zeta)                                                   
c                                                                               
c          dq_i/ds2, derivatives with respect to s2, the second                 
c          local coordinate: qs2(1) = dq1/ds2,                                  
c                            qs2(2) = dq2/ds2,                                  
c                            qs2(3) = dq3/ds2, etc.                             
c                                                                               
      qs2(1) = zero                                                             
      qs2(2) = half*(1-zeta)                                                    
      qs2(3) = -half*(1-zeta)                                                   
      qs2(4) = zero                                                             
      qs2(5) = half*(1+zeta)                                                    
      qs2(6) = -half*(1+zeta)                                                   
c                                                                               
c          dq_i/d(zeta), derivatives with respect to zeta, the third            
c          local coordinate: qzeta(1) = dq1/d(zeta),                            
c                            qzeta(2) = dq2/d(zeta),                            
c                            qzeta(3) = dq3/d(zeta), etc.                       
c                                                                               
      qzeta(1) = -half*s1                                                       
      qzeta(2) = -half*s2                                                       
      qzeta(3) = -half*s3                                                       
      qzeta(4) = half*s1                                                        
      qzeta(5) = half*s2                                                        
      qzeta(6) = half*s3                                                        
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine deriv_trint12_wedge                *          
c     *                                                              *          
c     *                    written by : sushovan                     *          
c     *                                 01/05/01                     *          
c     *                                                              *          
c     *     given a point in parametric coordinates, return the      *          
c     *     parametric derivatives for each node shape function,     *          
c     *     for the cohesive element trint12, treating it as wedge.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c       use triangle natural coordinates (s1,s2,s3) for the wedge base.         
c       (see SR shape15 and deriv15 for triangle natural coordinates).          
c       zeta is along the height of the wedge. the two traingular faces         
c       are at zeta = -1 (nodes 1-6) and zeta = +1 (nodes 7-12).                
c                                                                               
c       shape functions for the wedge:                                          
c                                                                               
c        q1  = half*(2.0*s1*s1 - s1)(1-zeta)                                    
c        q2  = half*(2.0*s2*s2 - s2)(1-zeta)                                    
c        q3  = half*(2.0*s3*s3 - s3)(1-zeta)                                    
c        q4  = half*four*s1*s2*(1-zeta)                                         
c        q5  = half*four*s2*s3*(1-zeta)                                         
c        q6  = half*four*s3*s1*(1-zeta)                                         
c        q7  = half*(2.0*s1*s1 - s1)(1+zeta)                                    
c        q8  = half*(2.0*s2*s2 - s2)(1+zeta)                                    
c        q9  = half*(2.0*s3*s3 - s3)(1+zeta)                                    
c        q10 = half*four*s1*s2*(1+zeta)                                         
c        q11 = half*four*s2*s3*(1+zeta)                                         
c        q12 = half*four*s3*s1*(1+zeta)                                         
c                                                                               
c       where half=0.5, four=4.0                                                
c                                                                               
c       s3 is dependent coordinate. replace s3=1-s1-s2 while                    
c       taking the derivatives.                                                 
c                                                                               
      subroutine deriv_trint12_wedge( etype, s1, s2, zeta,                      
     &                                qs1, qs2, qzeta )                         
      implicit none                                                             
c                                                                               
      integer etype                                                             
      double precision                                                          
     &     s1, s2, zeta, qs1(*), qs2(*), qzeta(*)                               
c                                                                               
c             local variables                                                   
c                                                                               
      double precision                                                          
     &     s3, zero, half, one, four                                            
      data zero, half, one, four /0.0d0, 0.5d0, 1.0d0, 4.0d0/                   
c                                                                               
      s3 = one - s1 -s2                                                         
c                                                                               
c          dq_i/ds1, derivatives with respect to s1, the first                  
c          local coordinate : qs1(1) = dq1/ds1,                                 
c                             qs1(2) = dq2/ds1,                                 
c                             qs1(3) = dq3/ds1, etc.                            
c                                                                               
      qs1(1)  = half*(four*s1 - one)*(1-zeta)                                   
      qs1(2)  = zero                                                            
      qs1(3)  = half*(one - four*s3)*(1-zeta)                                   
      qs1(4)  = half*four*s2*(1-zeta)                                           
      qs1(5)  = -half*four*s2*(1-zeta)                                          
      qs1(6)  = half*four*(s3 - s1)*(1-zeta)                                    
      qs1(7)  = half*(four*s1 - one)*(1+zeta)                                   
      qs1(8)  = zero                                                            
      qs1(9)  = half*(one - four*s3)*(1+zeta)                                   
      qs1(10) = half*four*s2*(1+zeta)                                           
      qs1(11) = -half*four*s2*(1+zeta)                                          
      qs1(12) = half*four*(s3 - s1)*(1+zeta)                                    
c                                                                               
c          dq_i/ds2, derivatives with respect to s2, the second                 
c          local coordinate: qs2(1) = dq1/ds2,                                  
c                            qs2(2) = dq2/ds2,                                  
c                            qs2(3) = dq3/ds2, etc.                             
c                                                                               
      qs2(1)  = zero                                                            
      qs2(2)  = half*(four*s2 - one)*(1-zeta)                                   
      qs2(3)  = half*(one - four*s3)*(1-zeta)                                   
      qs2(4)  = half*four*s1*(1-zeta)                                           
      qs2(5)  = half*four*(s3 - s2)*(1-zeta)                                    
      qs2(6)  = -half*four*s1*(1-zeta)                                          
      qs2(7)  = zero                                                            
      qs2(8)  = half*(four*s2 - one)*(1+zeta)                                   
      qs2(9)  = half*(one - four*s3)*(1+zeta)                                   
      qs2(10) = half*four*s1*(1+zeta)                                           
      qs2(11) = half*four*(s3 - s2)*(1+zeta)                                    
      qs2(12) = -half*four*s1*(1+zeta)                                          
c                                                                               
c          dq_i/d(zeta), derivatives with respect to zeta, the third            
c          local coordinate: qzeta(1) = dq1/d(zeta),                            
c                            qzeta(2) = dq2/d(zeta),                            
c                            qzeta(3) = dq3/d(zeta), etc.                       
c                                                                               
      qzeta(1)  = -half*(2.0*s1*s1 - s1)                                        
      qzeta(2)  = -half*(2.0*s2*s2 - s2)                                        
      qzeta(3)  = -half*(2.0*s3*s3 - s3)                                        
      qzeta(4)  = -half*four*s1*s2                                              
      qzeta(5)  = -half*four*s2*s3                                              
      qzeta(6)  = -half*four*s3*s1                                              
      qzeta(7)  = half*(2.0*s1*s1 - s1)                                         
      qzeta(8)  = half*(2.0*s2*s2 - s2)                                         
      qzeta(9)  = half*(2.0*s3*s3 - s3)                                         
      qzeta(10) = half*four*s1*s2                                               
      qzeta(11) = half*four*s2*s3                                               
      qzeta(12) = half*four*s3*s1                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine check_cohes_tri                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                   last modified: 2/18/13 rhd                 *          
c     *                                                              *          
c     *     make some sanity checks on the geometry of triangular    *          
c     *     interface elements. uses undeformed coordinates.         *          
c     *     prints warning messages if potential issues are found    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine check_cohes_tri( etype, span, mxvl, felem, nnode,              
     &                            iout,  ce )                                   
      implicit none                                                             
c                                                                               
      integer etype, span, mxvl, felem, nnode, iout                             
      double precision                                                          
     & ce(mxvl,*)                                                               
c                                                                               
c                 local data. some are automatic arrays.                        
c                                                                               
      double precision                                                          
     & xb(nnode), yb(nnode), zb(nnode), xt(nnode), yt(nnode),                   
     & zt(nnode), dxbt, dybt, dzbt, lbt,                                        
     & x(nnode), y(nnode), z(nnode), xm, ym, zm, dx, dy, dz, dmid,              
     & dx12, dy12, dz12, dx23, dy23, dz23, dx13, dy13, dz13,                    
     & lb12, lb23, lb13, lt12, lt23, lt13, l13, tolmin, tol_edge,               
     & long_edge(span), short_edge, mid_edge_tol, top_bott_tol,                 
     & half                                                                     
c                                                                               
      logical flg_bottom, flg_top, collapsed                                    
                                                                                
                                                                                
      integer nsurf, i, j, trint12_list(3,6), n1, n2, n3                        
      data tolmin, tol_edge, mid_edge_tol, top_bott_tol, half                   
     &  / 1.0d-08, 1.0d-03, 0.2d0, 0.05d0, 0.5d0 /                              
      data trint12_list / 1,4,2,  2,5,3,  3,6,1,  7,10,8,  8,11,9,              
     &                    9,12,7 /                                              
c                                                                               
c             1. all 3 sides must have lengths > 0.0 (tolmin)                   
c             2. no side can have less length tol_edge * longest side           
c                                                                               
c             These checks treat both triangle interface elements as            
c             six noded.                                                        
c                                                                               
      nsurf = nnode / 2                                                         
      collapsed = .false.                                                       
      do i = 1, span                                                            
!DIR$ IVDEP                                                                     
         do j = 1, nsurf                                                        
           xb(j) = ce(i,j)                                                      
           yb(j) = ce(i,nnode+j)                                                
           zb(j) = ce(i,2*nnode+j)                                              
           xt(j) = ce(i,nsurf+j)                                                
           yt(j) = ce(i,nnode+nsurf+j)                                          
           zt(j) = ce(i,2*nnode+nsurf+j)                                        
         end do                                                                 
         dx12 = xb(2) - xb(1)                                                   
         dy12 = yb(2) - yb(1)                                                   
         dz12 = zb(2) - zb(1)                                                   
         dx23 = xb(3) - xb(2)                                                   
         dy23 = yb(3) - yb(2)                                                   
         dz23 = zb(3) - zb(2)                                                   
         dx13 = xb(3) - xb(1)                                                   
         dy13 = yb(3) - yb(1)                                                   
         dz13 = zb(3) - zb(1)                                                   
         lb12 = sqrt( dx12*dx12 + dy12*dy12 + dz12*dz12 )                       
         lb23 = sqrt( dx23*dx23 + dy23*dy23 + dz23*dz23 )                       
         lb13 = sqrt( dx13*dx13 + dy13*dy13 + dz13*dz13 )                       
         flg_bottom = lb12 .le. tolmin .or. lb23 .le. tolmin .or.               
     &       lb13 .le. tolmin                                                   
         dx12 = xt(2) - xt(1)                                                   
         dy12 = yt(2) - yt(1)                                                   
         dz12 = zt(2) - zt(1)                                                   
         dx23 = xt(3) - xt(2)                                                   
         dy23 = yt(3) - yt(2)                                                   
         dz23 = zt(3) - zt(2)                                                   
         dx13 = xt(3) - xt(1)                                                   
         dy13 = yt(3) - yt(1)                                                   
         dz13 = zt(3) - zt(1)                                                   
         lt12 = sqrt( dx12*dx12 + dy12*dy12 + dz12*dz12 )                       
         lt23 = sqrt( dx23*dx23 + dy23*dy23 + dz23*dz23 )                       
         lt13 = sqrt( dx13*dx13 + dy13*dy13 + dz13*dz13 )                       
         flg_top = lt12 .le. tolmin .or. lt23 .le. tolmin .or.                  
     &       lt13 .le. tolmin                                                   
         long_edge(i) = max( lb12, lb23, lb13, lt12, lt23, lt13 )               
         short_edge = min( lb12, lb23, lb13, lt12, lt23, lt13 )                 
         if( short_edge .le. tol_edge * long_edge(i) )                          
     &         collapsed = .true.                                               
         if( flg_bottom .or. flg_top .or. collapsed ) then                      
            write(iout,9000) felem+i-1                                          
         end if                                                                 
      end do                                                                    
c                                                                               
c             trint12 elements. check that location of mid-side                 
c             nodes are not too far from mid-point on straight                  
c             line connecting two corners.                                      
                                                                                
      if( nnode .ne. 12 ) go to 300                                             
      do i = 1, span                                                            
!DIR$ IVDEP                                                                     
         do j = 1, 12                                                           
           x(j) = ce(i,j)                                                       
           y(j) = ce(i,12+j)                                                    
           z(j) = ce(i,24+j)                                                    
         end do                                                                 
!DIR$ IVDEP                                                                     
         do j = 1, 6                                                            
           n1 = trint12_list(1,j)                                               
           n2 = trint12_list(2,j)  ! middle node on edge                        
           n3 = trint12_list(3,j)                                               
           dx13 = x(n3) - x(n1)                                                 
           dy13 = y(n3) - y(n1)                                                 
           dz13 = z(n3) - z(n1)                                                 
           l13 = sqrt( dx13*dx13 + dy13*dy13 + dz13*dz13 )                      
           xm = x(n1) + half * dx13                                             
           ym = y(n1) + half * dy13                                             
           zm = z(n1) + half * dz13                                             
           dx = x(n2) - xm                                                      
           dy = y(n2) - ym                                                      
           dz = z(n2) - zm                                                      
           dmid = sqrt( dx*dx + dy*dy + dz*dz )                                 
           if( dmid .gt. mid_edge_tol * l13 ) then                              
             write(iout,9100) felem+i-1, n2                                     
           end if                                                               
         end do                                                                 
      end do                                                                    
c                                                                               
c             check distance between each pair of top and bottom                
c             nodes                                                             
                                                                                
 300  continue                                                                  
      nsurf = nnode / 2                                                         
      do i = 1, span                                                            
!DIR$ IVDEP                                                                     
         do j = 1, nsurf                                                        
           xb(j) = ce(i,j)                                                      
           yb(j) = ce(i,nnode+j)                                                
           zb(j) = ce(i,2*nnode+j)                                              
           xt(j) = ce(i,nsurf+j)                                                
           yt(j) = ce(i,nnode+nsurf+j)                                          
           zt(j) = ce(i,2*nnode+nsurf+j)                                        
           dxbt = xt(j) - xb(j)                                                 
           dybt = yt(j) - yb(j)                                                 
           dzbt = zt(j) - zb(j)                                                 
           lbt  = sqrt( dxbt*dxbt + dybt*dybt + dzbt*dzbt )                     
           if( lbt .gt. top_bott_tol * long_edge(i) ) then                      
            write(iout,9200) felem+i-1, j, j+nsurf, lbt, long_edge(i)           
           end if                                                               
         end do                                                                 
      end do                                                                    
                                                                                
       return                                                                   
c                                                                               
 9000 format(/1x,'>>>>> Warning: interface-cohesive element: ',i8,              
     & /,1x,'               appears to have one or more collapsed',             
     &      ' edges.',//)                                                       
 9100 format(/1x,'>>>>> Warning: interface-cohesive element: ',i8,              
     & /,1x,    '      Element mid-side node: ',i3,' appears to be',            
     &          ' located a significant distance',                              
     & /,1x,    '      from the mid-point on a line connecting the',            
     &          ' two corner nodes.',                                           
     & /,1x,    '      Check coordinates & incidences for errors.',//)          
 9200 format(/1x,'>>>>> Warning: interface-cohesive element: ',i8,              
     & /,1x,    '      Element bottom-top pair node: ',2i3,' are',              
     & /,1x,    '      are a significant distance apart: ',e14.6,               
     & /,1x,    '      compared to longest edge length of element: ',           
     & e14.6,//)                                                                
       end                                                                      
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine check_cohes_quad                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                   last modified: 2/18/13 rhd                 *          
c     *                                                              *          
c     *     make some sanity checks on the geometry of the 8-node    *          
c     *     interface element. uses undeformed coordinates.          *          
c     *     prints warning messages if potential issues are found    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine check_cohes_quad( etype, span, mxvl, felem, nnode,             
     &                             iout, ce )                                   
      implicit none                                                             
c                                                                               
      integer etype, span, mxvl, felem, nnode, iout                             
      double precision                                                          
     & ce(mxvl,*)                                                               
c                                                                               
c                 local data. some are automatic arrays.                        
c                                                                               
      double precision                                                          
     & xb(nnode), yb(nnode), zb(nnode), xt(nnode), yt(nnode),                   
     & zt(nnode), dxbt, dybt, dzbt, lbt,                                        
     & x(nnode), y(nnode), z(nnode), dx41, dy41, dz41,                          
     & dx12, dy12, dz12, dx23, dy23, dz23,                                      
     & dx34, dy34, dz34, lt34, lb34,                                            
     & lb12, lb23, lb41, lt12, lt23, lt41, tolmin, tol_edge,                    
     & long_edge(span), short_edge, top_bott_tol                                
c                                                                               
      logical flg_bottom, flg_top, collapsed                                    
                                                                                
                                                                                
      integer nsurf, i, j                                                       
      data tolmin, tol_edge, top_bott_tol                                       
     &  / 1.0d-08, 1.0d-03, 0.05d0 /                                            
c                                                                               
c             1. all 4 sides must have lengths > 0.0 (tolmin)                   
c             2. no side can have less length tol_edge * longest side           
c                                                                               
      nsurf = nnode / 2                                                         
      collapsed = .false.                                                       
      do i = 1, span                                                            
!DIR$ IVDEP                                                                     
         do j = 1, nsurf                                                        
           xb(j) = ce(i,j)                                                      
           yb(j) = ce(i,nnode+j)                                                
           zb(j) = ce(i,2*nnode+j)                                              
           xt(j) = ce(i,nsurf+j)                                                
           yt(j) = ce(i,nnode+nsurf+j)                                          
           zt(j) = ce(i,2*nnode+nsurf+j)                                        
         end do                                                                 
c                                                                               
         dx12 = xb(2) - xb(1)                                                   
         dy12 = yb(2) - yb(1)                                                   
         dz12 = zb(2) - zb(1)                                                   
c                                                                               
         dx23 = xb(3) - xb(2)                                                   
         dy23 = yb(3) - yb(2)                                                   
         dz23 = zb(3) - zb(2)                                                   
c                                                                               
         dx34 = xb(4) - xb(3)                                                   
         dy34 = yb(4) - yb(3)                                                   
         dz34 = zb(4) - zb(3)                                                   
c                                                                               
         dx41 = xb(1) - xb(4)                                                   
         dy41 = yb(1) - yb(4)                                                   
         dz41 = zb(1) - zb(4)                                                   
c                                                                               
         lb12 = sqrt( dx12*dx12 + dy12*dy12 + dz12*dz12 )                       
         lb23 = sqrt( dx23*dx23 + dy23*dy23 + dz23*dz23 )                       
         lb34 = sqrt( dx34*dx34 + dy34*dy34 + dz34*dz34 )                       
         lb41 = sqrt( dx41*dx41 + dy41*dy41 + dz41*dz41 )                       
         flg_bottom = lb12 .le. tolmin .or. lb23 .le. tolmin .or.               
     &       lb34 .le. tolmin  .or. lb41 .le. tolmin                            
c                                                                               
         dx12 = xt(2) - xt(1)                                                   
         dy12 = yt(2) - yt(1)                                                   
         dz12 = zt(2) - zt(1)                                                   
c                                                                               
         dx23 = xt(3) - xt(2)                                                   
         dy23 = yt(3) - yt(2)                                                   
         dz23 = zt(3) - zt(2)                                                   
c                                                                               
         dx34 = xt(4) - xt(3)                                                   
         dy34 = yt(4) - yt(3)                                                   
         dz34 = zt(4) - zt(3)                                                   
c                                                                               
         dx41 = xt(1) - xt(4)                                                   
         dy41 = yt(1) - yt(4)                                                   
         dz41 = zt(1) - zt(4)                                                   
c                                                                               
         lt12 = sqrt( dx12*dx12 + dy12*dy12 + dz12*dz12 )                       
         lt23 = sqrt( dx23*dx23 + dy23*dy23 + dz23*dz23 )                       
         lt34 = sqrt( dx34*dx34 + dy34*dy34 + dz34*dz34 )                       
         lt41 = sqrt( dx41*dx41 + dy41*dy41 + dz41*dz41 )                       
         flg_top = lt12 .le. tolmin .or. lt23 .le. tolmin .or.                  
     &          lt34 .le. tolmin  .or. lt41 .le. tolmin                         
c                                                                               
         long_edge(i) = max( lb12, lb23, lb34, lb41, lt12, lt23,                
     &                       lt34, lt41)                                        
         short_edge = min( lb12, lb23, lb34, lb41, lt12, lt23,                  
     &                       lt34, lt41 )                                       
         if( short_edge .le. tol_edge * long_edge(i) )                          
     &         collapsed = .true.                                               
         if( flg_bottom .or. flg_top .or. collapsed ) then                      
            write(iout,9000) felem+i-1                                          
         end if                                                                 
      end do                                                                    
c                                                                               
c             check distance between each pair of top and bottom                
c             nodes                                                             
c                                                                               
      nsurf = nnode / 2                                                         
      do i = 1, span                                                            
!DIR$ IVDEP                                                                     
         do j = 1, nsurf                                                        
           xb(j) = ce(i,j)                                                      
           yb(j) = ce(i,nnode+j)                                                
           zb(j) = ce(i,2*nnode+j)                                              
           xt(j) = ce(i,nsurf+j)                                                
           yt(j) = ce(i,nnode+nsurf+j)                                          
           zt(j) = ce(i,2*nnode+nsurf+j)                                        
           dxbt = xt(j) - xb(j)                                                 
           dybt = yt(j) - yb(j)                                                 
           dzbt = zt(j) - zb(j)                                                 
           lbt  = sqrt( dxbt*dxbt + dybt*dybt + dzbt*dzbt )                     
           if( lbt .gt. top_bott_tol * long_edge(i) ) then                      
            write(iout,9200) felem+i-1, j, j+nsurf, lbt, long_edge(i)           
           end if                                                               
         end do                                                                 
      end do                                                                    
                                                                                
       return                                                                   
c                                                                               
 9000 format(/1x,'>>>>> Warning: interface-cohesive element: ',i8,              
     & /,1x,'               appears to have one or more collapsed',             
     &      ' edges.',//)                                                       
 9200 format(/1x,'>>>>> Warning: interface-cohesive element: ',i8,              
     & /,1x,    '      Element bottom-top pair node: ',2i3,' are',              
     & /,1x,    '      are a significant distance apart: ',e14.6,               
     & /,1x,    '      compared to longest edge length of element: ',           
     & e14.6,//)                                                                
c                                                                               
       end                                                                      
                                                                                
