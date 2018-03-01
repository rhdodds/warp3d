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
c *                                                              *              
c ****************************************************************              
c                                                                               
c ****************************************************************              
c                                                                *              
c subroutine to calculate radius "r", and angle "theta" of       *              
c a single point, measured in polar coordinates in the           *              
c local crack-front system.                                      *              
c                                                                *              
c                                 written by: mcw                *              
c                              last modified: 9/18/03            *              
c                                                                *              
c ****************************************************************              
c                                                                               
      subroutine di_calc_r_theta( tag, front_nodes, num_front_nodes,            
     &                         front_coords, domain_type, domain_origin,        
     &                         front_order, point_x, point_y, point_z,          
     &                         elemno, ptno, r1, theta, crack_curvature,        
     &                         debug, out )                                     
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer tag, front_nodes(*), num_front_nodes, domain_type,                
     &        front_order, ptno, p0, p1, out, domain_origin, elemno             
      double precision                                                          
     & front_coords(3,*), point_x, point_y, point_z, r1, theta,                 
     & crack_curvature(*)                                                       
      logical debug                                                             
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, j                                                              
      logical curved_crack                                                      
      double precision                                                          
     & zero, len, pi, x1, y1, z1, x2, y2, z2, dist1, dist2, base,               
     & height, x_local, y_local, half, one, two, x3, y3, z3, a,                 
     & toler, zguess                                                            
c                                                                               
      data zero, half, one, two, pi                                             
     & / 0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.14159265359d0 /                          
c                                                                               
c             determine distance r and angle theta to the integration           
c             point from the straight line connecting the two closest           
c             crack-front nodes.                                                
c                                                                               
c             define: p    = integration point                                  
c                     p0   = node at beginning of element edge                  
c                     p1   = node at end of element edge                        
c                     p2   = projection of integration point 'p' onto           
c                            the local x-z plane                                
c                     r1   = distance from element edge to integration point 'p'
c                     base = distance from p2 to crack front.                   
c                                                                               
      p0 = 0                                                                    
      p1 = 0                                                                    
c                                                                               
c             one linear element: domain types 1, 3, 4.                         
c                                                                               
      if( front_order.eq.1 .and. num_front_nodes .eq. 2 ) then                  
         p0 = 1                                                                 
         p1 = 2                                                                 
      end if                                                                    
c                                                                               
c             one quadratic element: domain types 1, 3, 4.                      
c             measure r and theta from straight line between                    
c             the two front nodes closest to the point.                         
c             use closest end front node to determine                           
c             closest line segment.                                             
c                                                                               
      if( front_order.eq.2 .and. num_front_nodes.eq.3 ) then                    
         x1    = point_x - front_coords(1,1)                                    
         y1    = point_y - front_coords(2,1)                                    
         z1    = point_z - front_coords(3,1)                                    
         dist1 = (x1*x1 + y1*y1 + z1*z1)**half                                  
         x2    = point_x - front_coords(1,3)                                    
         y2    = point_y - front_coords(2,3)                                    
         z2    = point_z - front_coords(3,3)                                    
         dist2 = (x2*x2 + y2*y2 + z2*z2)**half                                  
c                                                                               
         if( dist1.le.dist2 ) then                                              
            p0 = 1                                                              
            p1 = 2                                                              
         else                                                                   
            p0 = 2                                                              
            p1 = 3                                                              
         endif                                                                  
      end if                                                                    
c                                                                               
c             two linear elements: domain_types 2, 4.                           
c             measure r and theta from straight line between                    
c             the two front nodes closest to the point.                         
c             use local z-coordinate to determine closest                       
c             line segment.                                                     
c                                                                               
      if( front_order.eq.1 .and. num_front_nodes.eq.3 ) then                    
c                                                                               
         dist1 = point_z - front_coords(3,2)                                    
c                                                                               
         if( dist1.le.zero ) then                                               
            p0 = 1                                                              
            p1 = 2                                                              
         else                                                                   
            p0 = 2                                                              
            p1 = 3                                                              
         endif                                                                  
      end if                                                                    
c                                                                               
c             two quadratic elements: domain_types 2, 4.                        
c             measure r and theta from straight line between                    
c             the two front nodes closest to the point.                         
c             use local z-coordinate, and then closest end                      
c             front node to determine closest line segment.                     
c                                                                               
      if( front_order.eq.2 .and. num_front_nodes.eq.5 ) then                    
c                                                                               
         dist1 = point_z - front_coords(3,3)                                    
c                                                                               
c             if point is closest to first crack-front segment.                 
c                                                                               
         if( dist1.le.zero ) then                                               
            x1    = point_x - front_coords(1,1)                                 
            y1    = point_y - front_coords(2,1)                                 
            z1    = point_z - front_coords(3,1)                                 
            dist1 = (x1*x1 + y1*y1 + z1*z1)**half                               
            x2    = point_x - front_coords(1,3)                                 
            y2    = point_y - front_coords(2,3)                                 
            z2    = point_z - front_coords(3,3)                                 
            dist2 = (x2*x2 + y2*y2 + z2*z2)**half                               
            if( dist1.le.dist2 ) then                                           
               p0 = 1                                                           
               p1 = 2                                                           
            else                                                                
               p0 = 2                                                           
               p1 = 3                                                           
            endif                                                               
         else                                                                   
c                                                                               
c             if point is closest to second crack-front segment.                
c                                                                               
            x1    = point_x - front_coords(1,3)                                 
            y1    = point_y - front_coords(2,3)                                 
            z1    = point_z - front_coords(3,3)                                 
            dist1 = (x1*x1 + y1*y1 + z1*z1)**half                               
            x2    = point_x - front_coords(1,5)                                 
            y2    = point_y - front_coords(2,5)                                 
            z2    = point_z - front_coords(3,5)                                 
            dist2 = (x2*x2 + y2*y2 + z2*z2)**half                               
            if( dist1.le.dist2 ) then                                           
               p0 = 3                                                           
               p1 = 4                                                           
            else                                                                
               p0 = 4                                                           
               p1 = 5                                                           
            endif                                                               
         end if                                                                 
c                                                                               
      end if                                                                    
c                                                                               
c             if domain_type = 4, type 'd', assume the crack                    
c             front is straight (or nearly straight). this                      
c             is compatible with the local crack-front                          
c             coordinate system, which is defined using the                     
c             first two front nodes in dimrot.f.                                
c                                                                               
      if( domain_type.eq.4 .and. num_front_nodes.gt.5 ) then                    
         p0 = 1                                                                 
         p1 = 2                                                                 
      end if                                                                    
c                                                                               
      if( debug ) write(out,895) p0, p1                                         
c                                                                               
c             compute distance from integration point p to its                  
c             projection onto the crack plane, point p2.                        
c                                                                               
      y_local =  point_y - front_coords(2,p0)                                   
      height  = abs( y_local )                                                  
c                                                                               
c             compute distance from p2 to chord. first compute                  
c             magnitude of cross product of vector p0-p1 and vector             
c             p0-p2, then divide by magnitude of vector p0-p1.                  
c                                                                               
      x1 = front_coords(1,p1) - front_coords(1,p0)                              
      y1 = front_coords(2,p1) - front_coords(2,p0)                              
      z1 = front_coords(3,p1) - front_coords(3,p0)                              
      x2 = point_x - front_coords(1,p0)                                         
      y2 = zero                                                                 
      z2 = point_z - front_coords(3,p0)                                         
c                                                                               
      base = (   (y1*z2 - z1*y2)**two                                           
     &         + (z1*x2 - x1*z2)**two                                           
     &         + (x1*y2 - y1*x2)**two )**half                                   
c                                                                               
      base = base / (x1*x1 + y1*y1 + z1*z1)**half                               
c                                                                               
c             determine if projected point lies ahead of or behind the          
c             crack front. take cross product of vector p0-p1 and vector        
c             p0-p2. if result is positive, projected point lies ahead          
c             of crack front, and if negative, it's behind.                     
c                                                                               
      if( z1*x2 - x1*z2 .ge. zero ) then                                        
         x_local =  one                                                         
      else                                                                      
         x_local = -one                                                         
      end if                                                                    
c                                                                               
c             for curved cracks, call routine to recompute "base."              
c             determine if p2 is ahead of or behind curved front.               
c                                                                               
      if( int( crack_curvature(1) ) .eq. 1 ) then                               
         if( debug ) write(out,894)                                             
         if( num_front_nodes .eq. 3 ) zguess = front_coords(3,2)                
         if( num_front_nodes .ge. 5 ) then                                      
            if( p0 .eq. 1 .or. p0 .eq. 2 ) zguess = front_coords(3,2)           
            if( p0 .eq. 3 .or. p0 .eq. 4 ) zguess = front_coords(3,4)           
         end if                                                                 
         call di_calc_distance( crack_curvature, front_order,                   
     &              num_front_nodes, domain_type, zguess,                       
     &              point_x, point_z, base, x_local, out, debug )               
      end if                                                                    
c                                                                               
      r1 = ( base*base + height*height )**half                                  
c                                                                               
c             calculate angle from the crack plane to the integration           
c             point, measured perpendicular to the crack front segment.         
c             the angle is determined by the quadrant in which the              
c             integration point is located.                                     
c                                                                               
      theta = zero                                                              
c                                                                               
c             quadrant I                                                        
c                                                                               
      if( x_local.gt.zero .and. y_local.gt.zero ) then                          
         theta = acos( abs(base)/r1 )                                           
      end if                                                                    
c                                                                               
c             quadrant II                                                       
c                                                                               
      if( x_local.lt.zero .and. y_local.ge.zero ) then                          
         theta = pi - acos( abs(base)/r1 )                                      
      end if                                                                    
c                                                                               
c             quadrant III                                                      
c                                                                               
      if( x_local.lt.zero .and. y_local.lt.zero ) then                          
         theta = -( pi - acos( abs(base)/r1 ) )                                 
      end if                                                                    
c                                                                               
c             quadrant IV                                                       
c                                                                               
      if( x_local.ge.zero .and. y_local.lt.zero ) then                          
         theta = - acos( abs(base)/r1 )                                         
      end if                                                                    
c                                                                               
      if( debug ) then                                                          
c                                                                               
         if( tag .eq. 1 ) write(out,896) "r, theta for element ",               
     &        elemno, ", node ", ptno                                           
         if( tag .eq. 2 ) write(out,896) "r, theta for element ",               
     &        elemno, ", integration point ", ptno                              
         write(out,897) "domain_origin",                                        
     &        front_nodes(domain_origin),                                       
     &        "domain origin coordinates:",  "x",                               
     &        front_coords(1,domain_origin), "y",                               
     &        front_coords(2,domain_origin), "z",                               
     &        front_coords(3,domain_origin),                                    
     &        " point/node coordinates:", "x", point_x,                         
     &        "y", point_y, "z", point_z                                        
c                                                                               
         write(out,900) front_nodes(p0), front_coords(1,p0),                    
     &                  front_coords(2,p0), front_coords(3,p0)                  
         write(out,905) front_nodes(p1), front_coords(1,p1),                    
     &                  front_coords(2,p1), front_coords(3,p1)                  
         write(out,910) ptno, point_x, point_y, point_z                         
         write(out,915) base*x_local, y_local, r1, theta                        
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 894  format(/,'computing distance for curved crack',/)                         
 895  format(//,5x,'p0 = ',i1,2x,'p1 = ',i1)                                    
 896  format(//,a,i7,a,i6,/)                                                    
 897  format(5x,a27,3x,i7,/,5x,a27,3(/,5x,a27,2x,e13.6),/,5x,a27,               
     &       3(/,5x,a27,2x,e13.6),/)                                            
 900  format(//,5x,'local domain origin',/,                                     
     & 10x,'node: ',i7,2x,'x',2x,e13.6,2x,'y',2x,e13.6,2x,'z',2x,e13.6)         
 905  format(//,5x,'end node of crack segment',/,                               
     & 10x,'node: ',i7,2x,'x',2x,e13.6,2x,'y',2x,e13.6,2x,'z',2x,e13.6)         
 910  format(//,5x,'point/node coordinates',/,                                  
     & 10x,'ptno: ',i7,2x,'x',2x,e13.6,2x,'y',2x,e13.6,2x,'z',2x,e13.6)         
 915  format(//,10x,'x',2x,e13.6,2x,'y',6x,e13.6,                               
     &        /,10x,'r',2x,e13.6,2x,'theta',2x,e13.6,/)                         
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c ****************************************************************              
c                                                                *              
c subroutine to calculate distance from the projection of        *              
c current integration point on the crack plane, to a             *              
c quadratic curve defined by three crack-front nodes.            *              
c use Newton iteration.                                          *              
c                                                                *              
c                                 written by: mcw                *              
c                              last modified: 4/24/04            *              
c                                                                *              
c ****************************************************************              
c                                                                               
      subroutine di_calc_distance( crack_curvature, caseno, zguess,             
     &                point_x, point_z, base, x_local, out, debug )             
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer caseno, out                                                       
      double precision                                                          
     & zguess, point_x, point_z, base, x_local, crack_curvature(*)              
      logical debug                                                             
c                                                                               
c             local variables                                                   
c                                                                               
      integer i                                                                 
      double precision                                                          
     & zero, half, one, two, three, x, fz, fprimez, z, znext,                   
     & d1, d2, dist, x_center, z_center, circle_radius, a, b, c, d, e,          
     & four                                                                     
c                                                                               
      data zero, half, one, two, three, four                                    
     & / 0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0 /                             
c                                                                               
c             for a circular crack, find center and radius of a                 
c             circle defined by three crack-front nodes.                        
c                                                                               
      if( int(crack_curvature(2)) .eq. 1 ) then                                 
         x_center      = crack_curvature(3)                                     
         z_center      = crack_curvature(4)                                     
         circle_radius = crack_curvature(5)                                     
c                                                                               
         d1            = point_x - x_center                                     
         d2            = point_z - z_center                                     
         dist          = (d1*d1 + d2*d2)**half                                  
         if( dist .ge. circle_radius ) then                                     
            base    = dist - circle_radius                                      
            x_local = one                                                       
         else                                                                   
            base    = circle_radius - dist                                      
            x_local = -one                                                      
         end if                                                                 
         return                                                                 
      end if                                                                    
c                                                                               
c             for crack fronts of unknown curvature:                            
c                                                                               
      if( int( crack_curvature(2) ) .eq. 0 ) then                               
c                                                                               
c             for three nodes on front, fit nodes to a                          
c             quadratic polynomial: x = az^2 + bz + c                           
c                                                                               
         if( caseno .eq. 3 ) then                                               
            a = crack_curvature(3)                                              
            b = crack_curvature(4)                                              
            c = crack_curvature(5)                                              
         end if                                                                 
c                                                                               
c             for five nodes on front, fit nodes to a fourth-order              
c             polynomial:  x = az^4 + bz^3 + cz^2 + dz + e                      
c                                                                               
         if( caseno .eq. 5 ) then                                               
            a = crack_curvature(3)                                              
            b = crack_curvature(4)                                              
            c = crack_curvature(5)                                              
            d = crack_curvature(6)                                              
            e = crack_curvature(7)                                              
         end if                                                                 
c                                                                               
c             quadratic polynomial 'g(z)' describes the crack front.            
c             the derivative of the distance from point p2 to g(z) equals       
c             f(z) =  z^3 + 3b/(2a)z^2 + 1/(2a^2)(b^2 + 2ac - 2apoint_x + 1)z   
c                  + 1/(2a^2)(bc - bpoint_x - point_z)                          
c             the second derivative equals                                      
c             f'(z) = 3z^2 + (3b/a)z + 1/(2a^2) (b^2 + 2ac - 2apoint_x + 1)     
c                                                                               
c             a similar procedure generates f(z) and f'(z) for a                
c             fourth-order polynomial                                           
c                                                                               
c             Newton's iteration determines the zero of f(z) as follows:        
c                                                                               
c             z(i+1) = z(i) - f(z(i))/f'(z(i))                                  
c                                                                               
         z = zguess                                                             
         do i=1,5                                                               
            if( caseno .eq. 3 ) then                                            
               fz = z*z*z + three*b*z*z/(two*a)                                 
     &            + z/(two*a*a) * ( b*b + two*a*c - two*a*point_x + 1 )         
     &            + one/(two*a*a) * ( b*c - b*point_x - point_z )               
c                                                                               
               fprimez = three*z*z + three*b*z/a                                
     &           + one/(two*a*a) * ( b*b + two*a*c - two*a*point_x + 1 )        
            end if                                                              
c                                                                               
            if( caseno .eq. 5 ) then                                            
               fz = z - point_z                                                 
     &            + ( a*z*z*z*z + b*z*z*z + c*z*z + d*z + e - point_x )         
     &            * ( four*a*z*z*z + three*b*z*z + two*c*z + d )                
               fprimez = one +(four*a*z*z*z +three*b*z*z +2*c*z +d)**two        
     &                 + ( a*z*z*z*z +b*z*z*z +c*z*z +d*z +e - point_x )        
     &                 * ( three*four*a*z*z + two*three*b*z + two*c )           
            end if                                                              
            znext = z - fz / fprimez                                            
            z     = znext                                                       
         end do                                                                 
c                                                                               
c             use z in the appropriate function to find x.                      
c             use x and z to find distance "base."                              
c                                                                               
         if( caseno .eq. 3 ) x = a*z*z + b*z + c                                
         if( caseno .eq. 5 ) x = a*z*z*z*z + b*z*z*z + c*z*z + d*z + e          
         base = ((point_x - x)**two + (point_z - z)**two)**half                 
c                                                                               
c             determine if the integration point lies ahead of                  
c             or behind the crack front.                                        
c                                                                               
         if( point_x .ge. x ) then                                              
            x_local = one                                                       
         else                                                                   
            x_local = -one                                                      
         end if                                                                 
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c***************************************************************                
c                                                              *                
c    subroutine to calculate constitutive tensors for          *                
c    interaction integral                                      *                
c                                                              *                
c                                 written by: mcw              *                
c                              last modified: 9/18/03          *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_calc_constitutive( dcijkl_x1, sijkl, dsijkl_x1,             
     &                              e, nu, de_x1, dnu_x1, elemno, debug,        
     &                              out )                                       
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer elemno, out                                                       
      double precision                                                          
     &     dcijkl_x1(3), sijkl(3), dsijkl_x1(3), e, nu, de_x1,                  
     &     dnu_x1                                                               
c                                                                               
c             local variables                                                   
c                                                                               
      double precision                                                          
     &     zero, one, two                                                       
      logical debug                                                             
c                                                                               
      data zero, one, two                                                       
     & / 0.d0, 1.d0, 2.d0 /                                                     
c                                                                               
c                                                                               
c             calculate three nonzero components of constitutive tensor         
c             derivative, compliance tensor, and compliance tensor derivative   
c             using e, nu, de_x1 and dnu_x1 at integration points.              
c                                                                               
c             constitutive and compliance relations expressed by (e.g. y.c. fung
c             foundations of solid mechanics, 1965, p. 129):                    
c                                                                               
c                sig_ij = lambda * (eps_11 + eps_22 + eps_33)*kronecker_ij      
c                       + 2 * mu * eps_ij                                       
c                                                                               
c                eps_ij = 1 / (2 * mu) * sig_ij                                 
c                       - nu / e * (sig_11 + sig_22 + sig_33) * kronecker_ij    
c                                                                               
c             derivative of constitutive (c) or compliance (s) tensor in x1 dire
c                d(c,s)ijkl_x1 = d(c,s)_e * de_x1 + d(c,s)_nu * dnu_x1          
c                                                                               
c             the three non-zero components of cijkl are:                       
c                (1) lambda + 2*mu = e(1-nu)/((1+nu)(1-2nu))                    
c                    d(1)/de  = (1-nu)/((1+nu)(1-2nu))                          
c                    d(1)/dnu = 2e*nu(2-nu)/((1+nu)^2*(1-2nu)^2)                
c                (2) lambda = e*nu/((1+nu)(1-2nu))                              
c                    d(2)/de  = nu/((1+nu)(1-2nu))                              
c                    d(2)/dnu = e(1+2nu^2)/((1+nu)^2*(1-2nu)^2)                 
c                (3) 2*mu = e/(1+nu)                                            
c                    d(3)/de  = 1/(1+nu)                                        
c                    d(3)/dnu = -e/(1+nu)^2                                     
c                                                                               
      dcijkl_x1(1:3) = zero                                                     
      sijkl(1:3)     = zero                                                     
      dsijkl_x1(1:3) = zero                                                     
c                                                                               
      dcijkl_x1(1) = ( one-nu ) / ((one+nu)*(one-two*nu)) * de_x1               
     &             + two * e * nu * (two-nu)                                    
     &             / ((one+nu)**two * (one-two*nu)**two)  * dnu_x1              
c                                                                               
      dcijkl_x1(2) = nu/((one+nu)*(one-two*nu))           * de_x1               
     &             + e * (one+two*nu**2)                                        
     &             / ((one+nu)**two * (one-two*nu)**two)  * dnu_x1              
c                                                                               
      dcijkl_x1(3) = one/(one+nu)                         * de_x1               
     &             - e/(one+nu)**two                      * dnu_x1              
c                                                                               
c             the three non-zero components of sijkl are:                       
c                (1) d(1) = 1/e                                                 
c                    d(1)/de = -1/e^2       d(1)/dnu = 0                        
c                (2) d(2) = -nu/e                                               
c                    d(2)/de = nu/e^2       d(2)/dnu = -1/e                     
c                (3) d(3) = (1+nu)/e                                            
c                    d(3)/de = -(1+nu)/e^2  d(3)/dnu = 1/e                      
c                                                                               
      sijkl(1)     =   one/e                                                    
      sijkl(2)     = - nu/e                                                     
      sijkl(3)     =   (one+nu)/e                                               
c                                                                               
      dsijkl_x1(1) = - one/e**two     * de_x1                                   
     &             + zero             * dnu_x1                                  
c                                                                               
      dsijkl_x1(2) = nu/e**two        * de_x1                                   
     &             - one/e            * dnu_x1                                  
c                                                                               
      dsijkl_x1(3) = -(one+nu)/e**two * de_x1                                   
     &             +   one/e          * dnu_x1                                  
c                                                                               
      if( debug ) then                                                          
         write(out,100) "constitutive tensor components"                        
         write(out,110) "e        : ", e,                                       
     &                  "nu       : ", nu,                                      
     &                  "de_x1    : ", de_x1,                                   
     &                  "dnu_x1   : ", dnu_x1,                                  
     &                  "dcijkl(1): ", dcijkl_x1(1),                            
     &                  "dcijkl(2): ", dcijkl_x1(2),                            
     &                  "dcijkl(3): ", dcijkl_x1(3),                            
     &                  "sijkl(1) : ", sijkl(1),                                
     &                  "sijkl(2) : ", sijkl(2),                                
     &                  "sijkl(3) : ", sijkl(3),                                
     &                  "dsijkl(1): ", dsijkl_x1(1),                            
     &                  "dsijkl(2): ", dsijkl_x1(2),                            
     &                  "dsijkl(3): ", dsijkl_x1(3)                             
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 100  format(/,5x,a,/)                                                          
 110  format(13(10x,a11,2x,e13.6,/))                                            
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c***************************************************************                
c                                                              *                
c subroutine to calculate auxiliary stresses and displacements *                
c at element integration points using the expressions obtained *                
c by Williams: On the stress distribution at the base of a     *                
c stationary crack. Journal of Applied Mechanics 24, 109-114,  *                
c 1957.                                                        *                
c                                                              *                
c                                 written by: mcw              *                
c                              last modified: 9/18/03          *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_calc_aux_fields_k( elemno, ptno, r, t, e_front,             
     &                              nu_front, dcijkl_x1, sijkl,                 
     &                              dsijkl_x1, aux_stress, aux_strain,          
     &                              daux_strain_x1, du11, du21, du31,           
     &                              du111, du112, du113,                        
     &                              du211, du212, du213,                        
     &                              du311, du312, du313, out )                  
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer elemno, ptno, out                                                 
      double precision                                                          
     &     r, t, e_front, nu_front,                                             
     &     dcijkl_x1(3), sijkl(3), dsijkl_x1(3),                                
     &     aux_stress(9,8), aux_strain(9,8), daux_strain_x1(9,8),               
     &     du11(8), du21(8), du31(8),                                           
     &     du111(8), du112(8), du113(8),                                        
     &     du211(8), du212(8), du213(8),                                        
     &     du311(8), du312(8), du313(8)                                         
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, j, p, node, snode, enode                                       
      double precision                                                          
     &     mu_front, ucoeff, t2, tt2, kappa, k(3),                              
     &     s11, s22, s33, s12, s13, s23,                                        
     &     ds111, ds121, ds131, ds221, ds231, ds331,                            
     &     du12(8), du22(8), du32(8),                                           
     &     u1(8), u2(8), u3(8),                                                 
     &     e11, e22, e33, e12, e23, e13,                                        
     &     de111, de221, de331, de121, de231, de131,                            
     &     sum1a, sum1b, sum2a, sum2b, sum3a, sum3b,                            
     &     zero, one, two, three, pi, four, five, six, eight,                   
     &     nine, ten, eleven, twelve, fourteen                                  
c                                                                               
      logical debug                                                             
c                                                                               
      data zero, one, two, three, pi, four, five, six, eight, nine, ten,        
     &     eleven, twelve, fourteen                                             
     & / 0.d0, 1.d0, 2.d0, 3.d0, 3.14159265359d0, 4.d0, 5.d0, 6.d0,             
     &   8.d0, 9.d0, 10.d0, 11.d0, 12.d0, 14.d0 /                               
c                                                                               
      debug = .false.                                                           
c                                                                               
c             properties e and nu correspond to the domain origin.              
c                                                                               
      mu_front = e_front / (two * ( one + nu_front ))                           
      ucoeff   = sqrt(r/(two*pi)) / (two*mu_front)                              
      t2       = t / two                                                        
      tt2      = t * three / two                                                
c                                                                               
      if( debug ) then                                                          
         write(out,101) "r", r, "theta", t, "theta/2", t2, "3*theta/2",         
     &                  tt2, "sqrt(r/(2*pi))/(2*mu_front)", ucoeff,             
     &                  "e_front", e_front, "nu_front", nu_front,               
     &                  "mu_front", mu_front                                    
      end if                                                                    
c                                                                               
c             calculate auxiliary fields at integration point including:        
c             stresses, strains, strain derivatives, displacements,             
c             and displacement derivatives. the five cases include:             
c                                                                               
c                j     auxiliary field                                          
c                                                                               
c                1     KI   plane stress                                        
c                2     KI   plane strain                                        
c                3     KII  plane stress                                        
c                4     KII  plane strain                                        
c                5     KIII anti-plane shear                                    
c                                                                               
      aux_stress(1:9,1:8)     = zero                                            
      aux_strain(1:9,1:8)     = zero                                            
      daux_strain_x1(1:9,1:8) = zero                                            
      du11(1:8)               = zero                                            
      du21(1:8)               = zero                                            
      du31(1:8)               = zero                                            
      du111(1:8)              = zero                                            
      du112(1:8)              = zero                                            
      du113(1:8)              = zero                                            
      du211(1:8)              = zero                                            
      du212(1:8)              = zero                                            
      du213(1:8)              = zero                                            
      du311(1:8)              = zero                                            
      du312(1:8)              = zero                                            
      du313(1:8)              = zero                                            
c                                                                               
      do j=1,5                                                                  
         k(1:3) = zero                                                          
c             KI plane stress                                                   
         if( j .eq. 1 ) then                                                    
            kappa = (three - nu_front) / (one + nu_front)                       
            k(1)  = one                                                         
         end if                                                                 
c             KI plane strain                                                   
         if( j .eq. 2 ) then                                                    
            kappa = three - four * nu_front                                     
            k(1)  = one                                                         
         end if                                                                 
c             KII plane stress                                                  
         if( j .eq. 3 ) then                                                    
            kappa = (three - nu_front) / (one + nu_front)                       
            k(2)  = one                                                         
         end if                                                                 
c             KII plane strain                                                  
         if( j .eq. 4 ) then                                                    
            kappa = three - four * nu_front                                     
            k(2)  = one                                                         
         end if                                                                 
c             KIII anti-plane shear                                             
         if( j .eq. 5 ) k(3) = one                                              
c                                                                               
c             calculate auxiliary stresses at current                           
c             integration point, and store in tensor form.                      
c                                                                               
         s11 = one / sqrt(two*pi*r)                                             
     &       * (   k(1) * cos(t2) * (one - sin(t2)*sin(tt2))                    
     &           - k(2) * sin(t2) * (two + cos(t2)*cos(tt2)))                   
c                                                                               
         s22 = one / sqrt(two*pi*r)                                             
     &       * (   k(1) * cos(t2) * (one + sin(t2)*sin(tt2))                    
     &           + k(2) * sin(t2)*cos(t2)*cos(tt2))                             
c                                                                               
         s12 = one / sqrt(two*pi*r)                                             
     &       * (   k(1) * sin(t2)*cos(t2)*cos(tt2)                              
     &           + k(2) * cos(t2) * (one - sin(t2)*sin(tt2)))                   
c                                                                               
         s13 = - one / sqrt(two*pi*r) * k(3) * sin(t2)                          
         s23 =   one / sqrt(two*pi*r) * k(3) * cos(t2)                          
         s33 =   zero                                                           
         if( j.eq.2 .or. j.eq.4 ) s33 = nu_front * (s11 + s22)                  
c                                                                               
         aux_stress(1,j) = s11                                                  
         aux_stress(2,j) = s12                                                  
         aux_stress(3,j) = s13                                                  
         aux_stress(4,j) = s12                                                  
         aux_stress(5,j) = s22                                                  
         aux_stress(6,j) = s23                                                  
         aux_stress(7,j) = s13                                                  
         aux_stress(8,j) = s23                                                  
         aux_stress(9,j) = s33                                                  
c                                                                               
c             calculate x1-derivatives of stress for use                        
c             in calculation of auxiliary strain.                               
c                                                                               
         ds111 = k(1)/(2*r*sqrt(2*pi*r))                                        
     &         * ( - cos(t)*cos(t2) + cos(t)*cos(t2)*sin(t2)*sin(tt2)           
     &             + sin(t)*sin(t2) - sin(t)*sin(tt2)                           
     &             + two*sin(t)*cos(t2)*cos(t2)*sin(tt2)                        
     &             + three*sin(t)*cos(t2)*sin(t2)*cos(tt2) )                    
     &         + k(2)/(2*r*sqrt(2*pi*r))                                        
     &         * (   two*sin(t2)*cos(t)                                         
     &             + cos(t)*sin(t2)*cos(t2)*cos(tt2)                            
     &             + two*sin(t)*cos(t2) - sin(t)*cos(tt2)                       
     &             + two*sin(t)*cos(t2)*cos(t2)*cos(tt2)                        
     &             - three*sin(t)*sin(t2)*cos(t2)*sin(tt2) )                    
c                                                                               
         ds121 = k(1)/(2*r*sqrt(2*pi*r))                                        
     &         * ( - cos(t)*cos(t2)*sin(t2)*cos(tt2)                            
     &             + sin(t)*cos(tt2)                                            
     &             - two*sin(t)*cos(t2)*cos(t2)*cos(tt2)                        
     &             + three*sin(t)*sin(t2)*cos(t2)*sin(tt2) )                    
     &         + k(2)/(2*r*sqrt(2*pi*r))                                        
     &         * ( - cos(t)*cos(t2)                                             
     &             + cos(t)*cos(t2)*sin(t2)*sin(tt2)                            
     &             + sin(t)*sin(t2) - sin(t)*sin(tt2)                           
     &             + two*sin(t)*cos(t2)*cos(t2)*sin(tt2)                        
     &             + three*sin(t)*cos(t2)*sin(t2)*cos(tt2) )                    
c                                                                               
         ds131 = k(3)/(2*r*sqrt(2*pi*r))                                        
     &         * (  sin(t2)*cos(t) + cos(t2)*sin(t) )                           
c                                                                               
         ds221 = k(1)/(2*r*sqrt(2*pi*r))                                        
     &         * ( - cos(t)*cos(t2)                                             
     &             - cos(t)*cos(t2)*sin(t2)*sin(tt2)                            
     &             + sin(t)*sin(t2) + sin(t)*sin(tt2)                           
     &             - two*sin(t)*cos(t2)*cos(t2)*sin(tt2)                        
     &             - three*sin(t)*cos(t2)*sin(t2)*cos(tt2) )                    
     &         + k(2)/(2*r*sqrt(2*pi*r))                                        
     &         * ( - cos(t)*cos(t2)*sin(t2)*cos(tt2)                            
     &             + sin(t)*cos(tt2)                                            
     &             - two*sin(t)*cos(t2)*cos(t2)*cos(tt2)                        
     &             + three*sin(t)*sin(t2)*cos(t2)*sin(tt2) )                    
c                                                                               
         ds231 = k(3)/(2*r*sqrt(2*pi*r))                                        
     &         * ( - cos(t2)*cos(t) + sin(t2)*sin(t) )                          
c                                                                               
         ds331 = zero                                                           
         if( j.eq.2 .or. j.eq.4 ) ds331 = nu_front * (ds111 + ds221)            
c                                                                               
c             calculate auxiliary displacements                                 
c             using material properties at the crack front                      
c                                                                               
         u1(j) = k(1)*ucoeff*cos(t2)*(kappa-one+two*(sin(t2))**two)             
     &         + k(2)*ucoeff*sin(t2)*(kappa+one+two*(cos(t2))**two)             
c                                                                               
         u2(j) = k(1)*ucoeff*sin(t2)*(kappa+one-two*(cos(t2))**two)             
     &         - k(2)*ucoeff*cos(t2)*(kappa-one-two*(sin(t2))**two)             
c                                                                               
         u3(j) = k(3)/mu_front * sqrt(two*r/pi) * sin(t2)                       
c                                                                               
c             calculate auxiliary displacement derivatives. we                  
c             calculate u1,2 and u2,2, and u3,2 only to verify auxiliary strains
c                                                                               
         du11(j) = k(1)/( four * mu_front * sqrt( two * pi * r ))               
     &           * (   cos(t)*cos(t2)*kappa + cos(t)*cos(t2)                    
     &               - two*cos(t)*(cos(t2))**3                                  
     &               + sin(t)*sin(t2)*kappa + sin(t)*sin(t2)                    
     &               - six*sin(t)*sin(t2)*(cos(t2))**two )                      
     &           + k(2)/( four * mu_front * sqrt( two * pi * r ))               
     &           * (   cos(t)*sin(t2)*kappa + cos(t)*sin(t2)                    
     &               + two*cos(t)*sin(t2)*(cos(t2))**two                        
     &               - sin(t)*cos(t2)*kappa + three*sin(t)*cos(t2)              
     &               - six*sin(t)*(cos(t2))**3 )                                
c                                                                               
         du21(j) = k(1)/( four * mu_front * sqrt( two * pi * r ))               
     &           * (   cos(t)*sin(t2)*kappa + cos(t)*sin(t2)                    
     &               - two*cos(t)*sin(t2)*(cos(t2))**two                        
     &               - sin(t)*cos(t2)*kappa - five*sin(t)*cos(t2)               
     &               + six*sin(t)*(cos(t2))**3 )                                
     &           + k(2)/( four * mu_front * sqrt( two * pi * r ))               
     &           * ( - cos(t)*cos(t2)*kappa + three*cos(t)*cos(t2)              
     &               - two*cos(t)*(cos(t2))**3                                  
     &               - sin(t)*sin(t2)*kappa + three*sin(t)*sin(t2)              
     &               - six*sin(t)*sin(t2)*(cos(t2))**two )                      
c                                                                               
         du31(j) = k(3)/( mu_front * sqrt( two * pi * r ))                      
     &           * ( sin(t2)*cos(t) - cos(t2)*sin(t) )                          
c                                                                               
         du12(j) = k(1)/( four * mu_front * sqrt( two * pi * r ))               
     &           * (   cos(t2)*sin(t)*kappa + cos(t2)*sin(t)                    
     &               - two*cos(t2)**3*sin(t) - sin(t2)*cos(t)*kappa             
     &               - sin(t2)*cos(t)                                           
     &               + six*sin(t2)*cos(t)*cos(t2)**two )                        
     &           + k(2)/( four * mu_front * sqrt( two * pi * r ))               
     &           * (   sin(t2)*sin(t)*kappa + sin(t2)*sin(t)                    
     &               + two*sin(t2)*sin(t)*cos(t2)**two                          
     &               + cos(t2)*cos(t)*kappa - three*cos(t2)*cos(t)              
     &               + six*cos(t2)**3*cos(t) )                                  
c                                                                               
         du22(j) = k(1)/( four * mu_front * sqrt( two * pi * r ))               
     &           * (   sin(t2)*sin(t)*kappa + sin(t2)*sin(t)                    
     &               - two*sin(t2)*sin(t)*cos(t2)**two                          
     &               + cos(t2)*cos(t)*kappa + five*cos(t2)*cos(t)               
     &               - six*cos(t2)**3*cos(t) )                                  
     &           + k(2)/( four * mu_front * sqrt( two * pi * r ))               
     &           * ( - cos(t2)*sin(t)*kappa + three*cos(t2)*sin(t)              
     &               - two*cos(t2)**3*sin(t) + sin(t2)*cos(t)*kappa             
     &               - three*sin(t2)*cos(t)                                     
     &               + six*sin(t2)*cos(t)*cos(t2)**two )                        
c                                                                               
         du32(j) = k(3)/( mu_front * sqrt( two * pi * r ))                      
     &           * ( sin(t2)*sin(t) + cos(t2)*cos(t) )                          
c                                                                               
c             calculate and store second derivatives of displacement (uj,1i)    
c                                                                               
         du111(j) = k(1)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * ( - two*cos(t)**two*cos(t2)*kappa                           
     &                + cos(t2)*kappa                                           
     &                + ten*cos(t)**two*cos(t2) - eleven*cos(t2)                
     &                - twelve*cos(t)**two*cos(t2)**3                           
     &                + fourteen*cos(t2)**3                                     
     &                - two*cos(t)*sin(t2)*sin(t)*kappa                         
     &                - two*cos(t)*sin(t2)*sin(t)                               
     &                + twelve*cos(t)*sin(t2)*sin(t)*cos(t2)**two )             
     &            + k(2)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * ( - two*cos(t)**two*sin(t2)*kappa                           
     &                - six*cos(t)**two*sin(t2)                                 
     &                + twelve*cos(t)**two*sin(t2)*cos(t2)**two                 
     &                + two*cos(t)*cos(t2)*sin(t)*kappa                         
     &                - six*cos(t)*cos(t2)*sin(t)                               
     &                + twelve*cos(t)*cos(t2)**3*sin(t)                         
     &                + sin(t2)*kappa + 5*sin(t2)                               
     &                - fourteen*sin(t2)*cos(t2)**two )                         
c                                                                               
         du112(j) = k(1)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * ( - two*cos(t)*cos(t2)*sin(t)*kappa                         
     &                + ten*cos(t)*cos(t2)*sin(t)                               
     &                - twelve*sin(t)*cos(t)*cos(t2)**3                         
     &                - sin(t2)*kappa                                           
     &                + two*cos(t)**two*sin(t2)*kappa                           
     &                - sin(t2) + two*cos(t)**two*sin(t2)                       
     &                + six*sin(t2)*cos(t2)**two                                
     &                - twelve*cos(t)**two*sin(t2)*cos(t2)**two )               
     &            + k(2)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * ( - two*cos(t)*sin(t2)*sin(t)*kappa                         
     &                - six*cos(t)*sin(t2)*sin(t)                               
     &                + twelve*cos(t)*sin(t2)*sin(t)*cos(t2)**two               
     &                + cos(t2)*kappa + six*cos(t)**two*cos(t2)                 
     &                - two*cos(t)**two*cos(t2)*kappa                           
     &                - three*cos(t2) + six*cos(t2)**3                          
     &                - twelve*cos(t)**two*cos(t2)**3 )                         
c                                                                               
         du113(j) = zero                                                        
c                                                                               
         du211(j) = k(1)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * ( - two*cos(t)**two*sin(t2)*kappa                           
     &                + two*cos(t)**two*sin(t2)                                 
     &                - twelve*cos(t)**two*sin(t2)*cos(t2)**two                 
     &                + two*cos(t)*cos(t2)*sin(t)*kappa                         
     &                + ten*cos(t)*cos(t2)*sin(t)                               
     &                - twelve*cos(t)*cos(t2)**3*sin(t)                         
     &                + sin(t2)*kappa - three*sin(t2)                           
     &                + fourteen*sin(t2)*cos(t2)**two )                         
     &            + k(2)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * (   two*cos(t)**two*cos(t2)*kappa                           
     &                + six*cos(t)**two*cos(t2)                                 
     &                - twelve*cos(t)**two*cos(t2)**3                           
     &                + two*cos(t)*sin(t2)*sin(t)*kappa                         
     &                - six*cos(t)*sin(t2)*sin(t)                               
     &                + twelve*cos(t)*sin(t2)*sin(t)*cos(t2)**two               
     &                - cos(t2)*kappa - nine*cos(t2)                            
     &                + fourteen*cos(t2)**3 )                                   
c                                                                               
         du212(j) = k(1)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            *( - two*cos(t)*sin(t2)*sin(t)*kappa                          
     &               + two*cos(t)*sin(t2)*sin(t) - six*cos(t2)**3               
     &               - twelve*cos(t)*sin(t2)*sin(t)*cos(t2)**two                
     &               + cos(t2)*kappa                                            
     &               - two*cos(t)**two*cos(t2)*kappa                            
     &               + five*cos(t2) - ten*cos(t)**two*cos(t2)                   
     &               + twelve*cos(t)**two*cos(t2)**3 )                          
     &            + k(2)/( eight * mu_front * r * sqrt(two*pi*r) )              
     &            * (   two*cos(t)*cos(t2)*sin(t)*kappa                         
     &                + six*cos(t)*cos(t2)*sin(t)                               
     &                + sin(t2)*kappa                                           
     &                - two*cos(t)**two*sin(t2)*kappa                           
     &                - three*sin(t2) + six*cos(t)**two*sin(t2)                 
     &                + six*sin(t2)*cos(t2)**two                                
     &                - twelve*cos(t)**two*sin(t2)*cos(t2)**two                 
     &                - twelve*cos(t)*cos(t2)**3*sin(t) )                       
c                                                                               
         du213(j) = zero                                                        
c                                                                               
         du311(j) = k(3)/( two * mu_front * r * sqrt(two*pi*r) )                
     &            * ( - two*cos(t)**two*sin(t2) + sin(t2)                       
     &                + two*cos(t)*cos(t2)*sin(t) )                             
c                                                                               
         du312(j) = k(3)/( two * mu_front * r * sqrt(two*pi*r) )                
     &            * ( - two*cos(t)*sin(t2)*sin(t) + cos(t2)                     
     &                - two*cos(t)**two*cos(t2) )                               
c                                                                               
         du313(j) = zero                                                        
c                                                                               
c             calculate auxiliary strains at current integration point          
c             (these are tensor strains due to mechanical loads--               
c             the auxiliary field includes no thermal loads).                   
c             store in symmetric tensor form.                                   
c                                                                               
c             eps_ij = sijkl * sig_ij_aux                                       
c                                                                               
         e11 = sijkl(1) * s11                                                   
     &       + sijkl(2) * s22                                                   
     &       + sijkl(2) * s33                                                   
c                                                                               
         e22 = sijkl(2) * s11                                                   
     &       + sijkl(1) * s22                                                   
     &       + sijkl(2) * s33                                                   
c                                                                               
         e33 = sijkl(2) * s11                                                   
     &       + sijkl(2) * s22                                                   
     &       + sijkl(1) * s33                                                   
c                                                                               
         if( j.eq.2 .or. j.eq.4 ) e33 = zero                                    
c                                                                               
         e12 = sijkl(3) * s12                                                   
         e13 = sijkl(3) * s13                                                   
         e23 = sijkl(3) * s23                                                   
c                                                                               
         aux_strain(1,j) = e11                                                  
         aux_strain(2,j) = e12                                                  
         aux_strain(3,j) = e13                                                  
         aux_strain(4,j) = e12                                                  
         aux_strain(5,j) = e22                                                  
         aux_strain(6,j) = e23                                                  
         aux_strain(7,j) = e13                                                  
         aux_strain(8,j) = e23                                                  
         aux_strain(9,j) = e33                                                  
c                                                                               
c             calculate derivatives of auxiliary strain at                      
c             integration points. store in symmetric tensor form.               
c                                                                               
c             eps_ij,1 = sijkl,1 * sig_ij_aux + sijkl * sig_ij,1_aux            
c                                                                               
         de111 = dsijkl_x1(1) * s11                                             
     &         + dsijkl_x1(2) * s22                                             
     &         + dsijkl_x1(2) * s33                                             
     &         + sijkl(1)     * ds111                                           
     &         + sijkl(2)     * ds221                                           
     &         + sijkl(2)     * ds331                                           
c                                                                               
         de221 = dsijkl_x1(2) * s11                                             
     &         + dsijkl_x1(1) * s22                                             
     &         + dsijkl_x1(2) * s33                                             
     &         + sijkl(2)     * ds111                                           
     &         + sijkl(1)     * ds221                                           
     &         + sijkl(2)     * ds331                                           
c                                                                               
c             for plane strain, de331 = 0.                                      
c             for plane stress, de331 is calculated as                          
c                                                                               
c               de331 = dsijkl_x1(2) * s11                                      
c     &               + dsijkl_x1(2) * s22                                      
c     &               + dsijkl_x1(1) * s33                                      
c     &               + sijkl(2)     * ds111                                    
c     &               + sijkl(2)     * ds221                                    
c     &               + sijkl(1)     * ds331                                    
c                                                                               
c             but we set de331 to zero for compatibility                        
c             with uj,1i in the incompatibility formulation                     
c             for fgms because uj,13 = 0.                                       
c                                                                               
         de331 = zero                                                           
c                                                                               
         de121 = dsijkl_x1(3) * s12 + sijkl(3) * ds121                          
         de231 = dsijkl_x1(3) * s23 + sijkl(3) * ds231                          
         de131 = dsijkl_x1(3) * s13 + sijkl(3) * ds131                          
c                                                                               
         daux_strain_x1(1,j) = de111                                            
         daux_strain_x1(2,j) = de121                                            
         daux_strain_x1(3,j) = de131                                            
         daux_strain_x1(4,j) = de121                                            
         daux_strain_x1(5,j) = de221                                            
         daux_strain_x1(6,j) = de231                                            
         daux_strain_x1(7,j) = de131                                            
         daux_strain_x1(8,j) = de231                                            
         daux_strain_x1(9,j) = de331                                            
c                                                                               
         if( debug.and.j.le.5 ) then                                            
            write(out,100) elemno, ptno                                         
            if( j.eq.1 ) write(out,104) "KI = 1, KII = 0, KIII = 0"             
            if( j.eq.2 ) write(out,104) "KI = 1, KII = 0, KIII = 0"             
            if( j.eq.3 ) write(out,104) "KI = 0, KII = 1, KIII = 0"             
            if( j.eq.4 ) write(out,104) "KI = 0, KII = 1, KIII = 0"             
            if( j.eq.5 ) write(out,104) "KI = 0, KII = 0, KIII = 1"             
            write(out,105) (aux_stress(p,j),p=1,9)                              
            write(out,106) ds111, ds121, ds131, ds121, ds221, ds231,            
     &           ds131, ds231, ds331                                            
            write(out,107) (aux_strain(p,j),p=1,9)                              
            write(out,108) (daux_strain_x1(p,j),p=1,9)                          
            write(out,109) u1(j), u2(j), u3(j)                                  
            write(out,110) du11(j), du22(j), du12(j),                           
     &           du21(j), du31(j), du32(j),                                     
     &           0.5*(du21(j) + du12(j)),                                       
     &           0.5*(du31(j)), 0.5*(du32(j))                                   
            write(out,111) du111(j), du112(j), du113(j),                        
     &           du211(j), du212(j), du213(j),                                  
     &           du311(j), du312(j), du313(j)                                   
            sum1a = daux_strain_x1(2,j) + daux_strain_x1(4,j)                   
            sum1b = du112(j)            + du211(j)                              
            sum2a = daux_strain_x1(3,j) + daux_strain_x1(7,j)                   
            sum2b = du113(j)            + du311(j)                              
            sum3a = daux_strain_x1(6,j) + daux_strain_x1(8,j)                   
            sum3b = du213(j)            + du312(j)                              
            write(out,112)                                                      
     &           daux_strain_x1(1,j), du111(j),                                 
     &           daux_strain_x1(5,j), du212(j),                                 
     &           daux_strain_x1(9,j), du313(j),                                 
     &           sum1a, sum1b, sum2a, sum2b, sum3a, sum3b                       
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
 100  format(////,"auxiliary fields for SIFs, element",2x,i7,                   
     &       2x,"point",2x,i2,//)                                               
 101  format(13(5x,a27,2x,e13.6,/))                                             
 104  format(/,10x,a)                                                           
 105  format(/,15x,"tensor of auxiliary stresses",3(/,13x,3(2x,e13.6)))         
 106  format(/,15x,"tensor of auxiliary stress derivatives",                    
     &       3(/,13x,3(2x,e13.6)))                                              
 107  format(/,15x,"tensor of auxiliary strains",3(/,13x,3(2x,e13.6)))          
 108  format(/,15x,"tensor of auxiliary strain derivatives",                    
     &       3(/,13x,3(2x,e13.6)))                                              
 109  format(/,15x,"auxiliary displacements",/,15x,"u1",2x,e13.6,/,15x,         
     &       "u2",2x,e13.6,/,15x,"u3",2x,e13.6)                                 
 110  format(/,15x,"auxiliary displacement uj,i derivatives",/,15x,             
     &       "u1,1",2x,e13.6,/,15x,"u2,2",2x,e13.6,/,15x,"u1,2",2x,             
     &       e13.6,/,15x,"u2,1",2x,e13.6,/,15x,"u3,1",2x,e13.6,/,15x,           
     &       "u3,2",2x,e13.6,/,15x,"0.5(u1,2 + u2,1)",2x,e13.6,/,15x,           
     &       "0.5(u3,1)",9x,e13.6,/,15x,"0.5(u3,2)",9x,e13.6)                   
 111  format(/,15x,"auxiliary displacement uj,1i derivatives",/,15x,            
     &       "u1,11",2x,e13.6,2x,"u1,12",2x,e13.6,2x,"u1,13",2x,                
     &       e13.6,/,15x,"u2,11",2x,e13.6,2x,"u2,12",2x,e13.6,2x,               
     &       "u2,13",2x,e13.6,/,15x,"u3,11",2x,e13.6,2x,"u3,12",2x,             
     &       e13.6,2x,"u3,13",2x,e13.6)                                         
 112  format(/,15x,'****quantities should be equal for homog matl****',         
     &       /,15x,'                    e11,1  u1,11 ', 2(2x,e11.4),            
     &       /,15x,'                    e22,1  u2,12 ', 2(2x,e11.4),            
     &       /,15x,'                    e33,1  u3,13 ', 2(2x,e11.4),            
     &       //,15x,' (e12,1 + e21,1)  (u1,12 + u2,11)', 2(2x,e11.4),           
     &       /,15x, ' (e13,1 + e31,1)  (u1,13 + u3,11)', 2(2x,e11.4),           
     &       /,15x, ' (e23,1 + e32,1)  (u2,13 + u3,12)', 2(2x,e11.4))           
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c***************************************************************                
c                                                              *                
c subroutine to calculate auxiliary stresses and displacements *                
c at integration points for t-stresses using the expressions   *                
c obtained by Michell for a point load on a straight crack:    *                
c (see Timoshenko and Goodier, Theory of Elasticity p. 110.,   *                
c  eq. 72 for stress sigma_r.)                                 *                
c                                                              *                
c                                 written by: mcw              *                
c                              last modified: 9/18/03          *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_calc_aux_fields_t( elemno, ptno, r, t, e_front,             
     &                              nu_front, dcijkl_x1, sijkl,                 
     &                              dsijkl_x1, aux_stress, aux_strain,          
     &                              daux_strain_x1, du11, du21, du31,           
     &                              du111, du112, du113,                        
     &                              du211, du212, du213,                        
     &                              du311, du312, du313, out )                  
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer elemno, ptno, out                                                 
      double precision                                                          
     &     r, t, e_front, nu_front,                                             
     &     dcijkl_x1(3), sijkl(3), dsijkl_x1(3),                                
     &     aux_stress(9,8), aux_strain(9,8), daux_strain_x1(9,8),               
     &     du11(8),  du21(8),  du31(8),                                         
     &     du111(8), du112(8), du113(8),                                        
     &     du211(8), du212(8), du213(8),                                        
     &     du311(8), du312(8), du313(8)                                         
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, j, p, node, snode, enode                                       
      double precision                                                          
     &     two_t, four_t, three_t, kappa, f,                                    
     &     s11, s22, s33, s12, s13, s23,                                        
     &     ds111, ds121, ds131, ds221, ds231, ds331,                            
     &     du12(8), du22(8), du32(8),                                           
     &     e11, e22, e33, e12, e23, e13,                                        
     &     de111, de221, de331, de121, de231, de131,                            
     &     sum1a, sum1b, sum2a, sum2b, sum3a, sum3b,                            
     &     zero, one, two, three, pi, four, five, mu_front                      
c                                                                               
      logical debug                                                             
c                                                                               
      data zero, one, two, three, pi, four, five                                
     & / 0.d0, 1.d0, 2.d0, 3.d0, 3.14159265359d0, 4.d0, 5.0d0 /                 
c                                                                               
      debug = .false.                                                           
c                                                                               
      f        = one                                                            
      two_t    = t * two                                                        
      four_t   = t * four                                                       
      three_t  = t * three                                                      
      mu_front = e_front / (two * (one + nu_front))                             
c                                                                               
c             use properties at domain origin for stresses                      
c             and displacements.                                                
c                                                                               
c             for each integration point, calculate auxiliary fields            
c             for three cases:                                                  
c                                                                               
c                j     auxiliary field                                          
c                                                                               
c                6     T11,T33  plane stress                                    
c                7     T11,T33  plane strain                                    
c                8     T13      anti-plane shear                                
c                                                                               
c             tensor values are stored as follows:                              
c                _11  _12  _13     1  2  3                                      
c                _21  _22  _23  =  4  5  6                                      
c                _31  _32  _33     7  8  9                                      
c                                                                               
c             compute auxiliary stresses at current                             
c             integration point, and store in tensor form.                      
c                                                                               
c                  plane stress, plane strain                                   
c                                                                               
      do j=6,7                                                                  
         aux_stress(1,j) = - f * cos(t)**three / (pi*r)                         
         aux_stress(2,j) = - f * cos(t)**two * sin(t) / (pi*r)                  
         aux_stress(3,j) = zero                                                 
         aux_stress(4,j) = - f * cos(t)**two * sin(t) / (pi*r)                  
         aux_stress(5,j) = - f * cos(t) * sin(t)**two / (pi*r)                  
         aux_stress(6,j) = zero                                                 
         aux_stress(7,j) = zero                                                 
         aux_stress(8,j) = zero                                                 
         aux_stress(9,j) = zero                                                 
      end do                                                                    
      aux_stress(9,7) = -f * nu_front / ( pi * r )                              
     &                * ( cos(t)**three + cos(t)*sin(t)**two )                  
c                                                                               
c                  anti-plane shear                                             
c                                                                               
      aux_stress(1,8) = zero                                                    
      aux_stress(2,8) = zero                                                    
      aux_stress(3,8) = -f * cos(t) / (two * pi * r)                            
      aux_stress(4,8) = zero                                                    
      aux_stress(5,8) = zero                                                    
      aux_stress(6,8) = -f * sin(t) / (two * pi * r)                            
      aux_stress(7,8) = -f * cos(t) / (two * pi * r)                            
      aux_stress(8,8) = -f * sin(t) / (two * pi * r)                            
      aux_stress(9,8) = zero                                                    
c                                                                               
c             compute auxiliary displacement derivatives. we                    
c             compute u1,2 and u2,2 only to verify auxiliary strains.           
c                                                                               
c                  plane stress                                                 
c                                                                               
      du11(6) = -f / (four * pi * e_front * r)                                  
     &        * (   cos(t) * (three - nu_front)                                 
     &            + cos(three_t)*(one + nu_front) )                             
      du21(6) = -f / (four * pi * e_front * r)                                  
     &        * (   sin(t)*(nu_front - three)                                   
     &            + sin(three_t)*(one + nu_front) )                             
      du12(6) = -f / ( four * pi * e_front * r )                                
     &        * (   sin(t)*(five + nu_front)                                    
     &            + sin(three_t) * ( one + nu_front ) )                         
      du22(6) = f / ( four * pi * e_front * r )                                 
     &        * (   cos(t) * ( three*nu_front - one )                           
     &            + cos(three_t)*( one + nu_front ) )                           
      du31(6) = zero                                                            
      du32(6) = zero                                                            
c                                                                               
c                  plane strain                                                 
c                                                                               
      du11(7) = f / (four * pi * e_front * r)                                   
     &        * (   cos(t)*(four*nu_front**two - three + nu_front)              
     &            - cos(three_t)*(one + nu_front) )                             
      du21(7) = - f / (four * pi * e_front * r)                                 
     &        * (   sin(t)*(nu_front - three + four*nu_front**two)              
     &            + sin(three_t)*(one + nu_front) )                             
      du12(7) = f / ( four * pi * e_front * r )                                 
     &        * ( - sin(t)*(five + nu_front -four*nu_front**two)                
     &            - sin(three_t) * ( one + nu_front ) )                         
      du22(7) = f / ( four * pi * e_front * r )                                 
     &        * (   cos(t) * ( three*nu_front                                   
     &            + four*nu_front**two - one )                                  
     &            + cos(three_t)*( one + nu_front ) )                           
      du31(7) = zero                                                            
      du32(7) = zero                                                            
c                                                                               
c                  anti-plane shear                                             
c                                                                               
      du11(8) = zero                                                            
      du21(8) = zero                                                            
      du12(8) = zero                                                            
      du22(8) = zero                                                            
      du31(8) = -f * cos(t) / ( two * r * pi * mu_front )                       
      du32(8) = -f * sin(t) / ( two * r * pi * mu_front )                       
c                                                                               
c             compute second derivatives of displacement (uj,1i)                
c                                                                               
c                  plane stress                                                 
c                                                                               
      du111(6) =  f / (two * pi * e_front * r**two)                             
     &         * (   cos(four_t) * ( one + nu_front )                           
     &             + cos(two_t)  * ( one - nu_front ) )                         
c                                                                               
      du112(6) =  f / (two * pi * e_front * r**two)                             
     &         * (   sin(four_t)*( one + nu_front )                             
     &             + sin(two_t) * two )                                         
c                                                                               
      du211(6) =  f / (two * pi * e_front * r**two)                             
     &         * ( - sin(two_t) * two                                           
     &             + sin(four_t) * ( one + nu_front ) )                         
c                                                                               
      du212(6) = -f / (two * pi * e_front * r**two)                             
     &         * (   cos(two_t) * ( nu_front - one )                            
     &             + cos(four_t) * ( one + nu_front ) )                         
c                                                                               
      du113(6) = zero                                                           
      du213(6) = zero                                                           
      du311(6) = zero                                                           
      du312(6) = zero                                                           
      du313(6) = zero                                                           
c                                                                               
c                  plane strain                                                 
c                                                                               
      du111(7) =  f / (two * pi * e_front * r**two)                             
     &         * (   cos(four_t) * ( one + nu_front )                           
     &             + cos(two_t)  * ( one - nu_front                             
     &             - two * nu_front**two ) )                                    
      du112(7) =  f / (two * pi * e_front * r**two)                             
     &         * (   sin(four_t)*( one + nu_front )                             
     &             + sin(two_t) *( two - two*nu_front**two) )                   
      du211(7) =  f / (two * pi * e_front * r**two)                             
     &         * (   sin(two_t)  * ( two * nu_front**two - two )                
     &             + sin(four_t) * ( one + nu_front ) )                         
      du212(7) = -f / (two * pi * e_front * r**two)                             
     &         * (   cos(two_t) * (   nu_front - one                            
     &                              + two*nu_front**two )                       
     &             + cos(four_t) * ( one + nu_front ) )                         
      du113(7) = zero                                                           
      du213(7) = zero                                                           
      du311(7) = zero                                                           
      du312(7) = zero                                                           
      du313(7) = zero                                                           
c                                                                               
c                  anti-plane shear                                             
c                                                                               
      du111(8) = zero                                                           
      du112(8) = zero                                                           
      du211(8) = zero                                                           
      du212(8) = zero                                                           
      du113(8) = zero                                                           
      du213(8) = zero                                                           
      du311(8) = f * cos(two_t) / (two * r * r * pi * mu_front )                
      du312(8) = f * sin(two_t) / (two * r * r * pi * mu_front )                
      du313(8) = zero                                                           
c                                                                               
c             calculate auxiliary strains at current integration point          
c             (these are tensor strains due to mechanical loads--               
c             the auxiliary field includes no thermal loads).                   
c             store in symmetric tensor form.                                   
c                                                                               
c             eps_ij = sijkl * sig_ij_aux                                       
c                                                                               
c                  plane stress, plane strain, anti-plane shear                 
c                                                                               
      do j=6,8                                                                  
         aux_strain(1,j) = sijkl(1) * aux_stress(1,j)                           
     &                   + sijkl(2) * aux_stress(5,j)                           
     &                   + sijkl(2) * aux_stress(9,j)                           
         aux_strain(2,j) = sijkl(3) * aux_stress(2,j)                           
         aux_strain(3,j) = sijkl(3) * aux_stress(3,j)                           
         aux_strain(4,j) = sijkl(3) * aux_stress(4,j)                           
         aux_strain(5,j) = sijkl(2) * aux_stress(1,j)                           
     &                   + sijkl(1) * aux_stress(5,j)                           
     &                   + sijkl(2) * aux_stress(9,j)                           
         aux_strain(6,j) = sijkl(3) * aux_stress(6,j)                           
         aux_strain(7,j) = sijkl(3) * aux_stress(7,j)                           
         aux_strain(8,j) = sijkl(3) * aux_stress(8,j)                           
         aux_strain(9,j) = sijkl(2) * aux_stress(1,j)                           
     &                   + sijkl(2) * aux_stress(5,j)                           
     &                   + sijkl(1) * aux_stress(9,j)                           
      end do                                                                    
      aux_strain(9,7) = zero                                                    
c                                                                               
c             compute x1-derivatives of stress. then                            
c             compute derivatives of auxiliary strain at                        
c             integration points. store in symmetric tensor form.               
c                                                                               
c             eps_ij,1 = sijkl,1 * sig_ij_aux + sijkl * sig_ij,1_aux            
c                                                                               
      do j=6,8                                                                  
         if( j.eq.6 .or. j.eq.7 ) then                                          
c                                                                               
c                  plane stress, plane strain                                   
c                                                                               
            ds111 = f/(two*pi*r**two) * ( cos(four_t) + cos(two_t) )            
            ds221 = f/(two*pi*r**two) * ( cos(two_t) - cos(four_t) )            
            ds121 = f/(two*pi*r**two) *   sin(four_t)                           
            ds131 = zero                                                        
            ds231 = zero                                                        
            ds331 = zero                                                        
            if( j.eq.7 ) ds331 = nu_front * ( ds111 + ds221 )                   
         end if                                                                 
c                                                                               
c                  anti-plane shear                                             
c                                                                               
         if( j.eq.8 ) then                                                      
            ds111 = zero                                                        
            ds221 = zero                                                        
            ds121 = zero                                                        
            ds131 = f * cos(two_t) / ( two * pi * r * r )                       
            ds231 = f * sin(two_t) / ( two * pi * r * r )                       
            ds331 = zero                                                        
         end if                                                                 
c                                                                               
         de111 = dsijkl_x1(1) * aux_stress(1,j)                                 
     &         + dsijkl_x1(2) * aux_stress(5,j)                                 
     &         + dsijkl_x1(2) * aux_stress(9,j)                                 
     &         + sijkl(1)     * ds111                                           
     &         + sijkl(2)     * ds221                                           
     &         + sijkl(2)     * ds331                                           
c                                                                               
         de221 = dsijkl_x1(2) * aux_stress(1,j)                                 
     &         + dsijkl_x1(1) * aux_stress(5,j)                                 
     &         + dsijkl_x1(2) * aux_stress(9,j)                                 
     &         + sijkl(2)     * ds111                                           
     &         + sijkl(1)     * ds221                                           
     &         + sijkl(2)     * ds331                                           
c                                                                               
c      de331 = dsijkl_x1(2) * s11                                               
c     &      + dsijkl_x1(2) * s22                                               
c     &      + dsijkl_x1(1) * s33                                               
c     &      + sijkl(2)     * ds111                                             
c     &      + sijkl(2)     * ds221                                             
c     &      + sijkl(1)     * ds331                                             
c                                                                               
c                  for plane stress, we set de331 to zero because               
c                  u3,13 equals zero in the incompatibility formulation.        
c                  for plane strain, de331 is zero.                             
c                                                                               
         de331 = zero                                                           
c                                                                               
         de121 = dsijkl_x1(3) * aux_stress(2,j) + sijkl(3) * ds121              
         de231 = dsijkl_x1(3) * aux_stress(6,j) + sijkl(3) * ds231              
         de131 = dsijkl_x1(3) * aux_stress(3,j) + sijkl(3) * ds131              
c                                                                               
         daux_strain_x1(1,j) = de111                                            
         daux_strain_x1(2,j) = de121                                            
         daux_strain_x1(3,j) = de131                                            
         daux_strain_x1(4,j) = de121                                            
         daux_strain_x1(5,j) = de221                                            
         daux_strain_x1(6,j) = de231                                            
         daux_strain_x1(7,j) = de131                                            
         daux_strain_x1(8,j) = de231                                            
         daux_strain_x1(9,j) = de331                                            
c                                                                               
         if( debug.and.j.ge.6 ) then                                            
            if(j.eq.6) write(out,100) elemno, ptno                              
            if(j.eq.7) write(out,101) elemno, ptno                              
            if(j.eq.8) write(out,102) elemno, ptno                              
            write(out,105) (aux_stress(p,j),p=1,9)                              
            write(out,106) ds111, ds121, ds131, ds121, ds221, ds231,            
     &                     ds131, ds231, ds331                                  
            write(out,107) (aux_strain(p,j),p=1,9)                              
            write(out,108) (daux_strain_x1(p,j),p=1,9)                          
            write(out,110) du11(j), du22(j), du12(j),                           
     &                     du21(j), du31(j), du32(j),                           
     &                     0.5*(du21(j) + du12(j)),                             
     &                     0.5*(du31(j)), 0.5*(du32(j))                         
            write(out,111) du111(j), du112(j), du113(j),                        
     &                     du211(j), du212(j), du213(j),                        
     &                     du311(j), du312(j), du313(j)                         
            sum1a = daux_strain_x1(2,j) + daux_strain_x1(4,j)                   
            sum1b = du112(j)            + du211(j)                              
            sum2a = daux_strain_x1(3,j) + daux_strain_x1(7,j)                   
            sum2b = du113(j)            + du311(j)                              
            sum3a = daux_strain_x1(6,j) + daux_strain_x1(8,j)                   
            sum3b = du213(j)            + du312(j)                              
            write(out,112)                                                      
     &           daux_strain_x1(1,j), du111(j),                                 
     &           daux_strain_x1(5,j), du212(j),                                 
     &           daux_strain_x1(9,j), du313(j),                                 
     &           sum1a, sum1b, sum2a, sum2b, sum3a, sum3b                       
         end if                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
 100  format(////,"plane stress aux fields for T11, element",2x,i7,             
     &       2x,"point",2x,i2,//)                                               
 101  format(////,"'plane strain' aux fields for T11, T33, element",            
     &       2x,i7,2x,"point",2x,i2,//)                                         
 102  format(////,"anti-plane shear aux fields for T13, element",2x,i7,         
     &       2x,"point",2x,i2,//)                                               
 105  format(/,15x,"tensor of auxiliary stresses",3(/,13x,3(2x,e13.6)))         
 106  format(/,15x,"tensor of auxiliary stress derivatives",                    
     &       3(/,13x,3(2x,e13.6)))                                              
 107  format(/,15x,"tensor of auxiliary strains",3(/,13x,3(2x,e13.6)))          
 108  format(/,15x,"tensor of auxiliary strain derivatives",                    
     &       3(/,13x,3(2x,e13.6)))                                              
 110  format(/,15x,"auxiliary displacement uj,i derivatives",/,15x,             
     &       "u1,1",2x,e13.6,/,15x,"u2,2",2x,e13.6,/,15x,"u1,2",2x,             
     &       e13.6,/,15x,"u2,1",2x,e13.6,/,15x,"u3,1",2x,e13.6,/,15x,           
     &       "u3,2",2x,e13.6,/,15x,"0.5(u1,2 + u2,1)",2x,e13.6,/,15x,           
     &       "0.5(u3,1)",9x,e13.6,/,15x,"0.5(u3,2)",9x,e13.6)                   
 111  format(/,15x,"auxiliary displacement uj,1i derivatives",/,15x,            
     &       "u1,11",2x,e13.6,2x,"u1,12",2x,e13.6,2x,"u1,13",2x,                
     &       e13.6,/,15x,"u2,11",2x,e13.6,2x,"u2,12",2x,e13.6,2x,               
     &       "u2,13",2x,e13.6,/,15x,"u3,11",2x,e13.6,2x,"u3,12",2x,             
     &       e13.6,2x,"u3,13",2x,e13.6)                                         
 112  format(/,15x,'****quantities should be equal for homog matl****',         
     &       /,15x,'                    e11,1  u1,11 ', 2(2x,e11.4),            
     &       /,15x,'                    e22,1  u2,12 ', 2(2x,e11.4),            
     &       /,15x,'                    e33,1  u3,13 ', 2(2x,e11.4),            
     &       //,15x,' (e12,1 + e21,1)  (u1,12 + u2,11)', 2(2x,e11.4),           
     &       /,15x, ' (e13,1 + e31,1)  (u1,13 + u3,11)', 2(2x,e11.4),           
     &       /,15x, ' (e23,1 + e32,1)  (u2,13 + u3,12)', 2(2x,e11.4))           
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c***************************************************************                
c                                                              *                
c    subroutine to calculate terms of i-integral               *                
c                                                              *                
c                                 written by: mcw              *                
c                              last modified: 9/18/03          *                
c                                                              *                
c***************************************************************                
c                                                                               
      subroutine di_calc_i_terms( ptno, dqx, dqy, dqz, dux, dvx, dwx,           
     &                         dtx, csig, aux_stress, ceps_gp,                  
     &                         aux_strain, dstrain_x1,                          
     &                         daux_strain_x1, dcijkl_x1,                       
     &                         du11_aux,  du21_aux,  du31_aux,                  
     &                         du111_aux, du211_aux, du311_aux,                 
     &                         du112_aux, du212_aux, du312_aux,                 
     &                         du113_aux, du213_aux, du313_aux,                 
     &                         process_temps, elem_alpha, dalpha_x1,            
     &                         point_temp, point_q, weight, elemno,             
     &                         fgm_e, fgm_nu, iterm, out, debug)                
c                                                                               
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer ptno, elemno, out                                                 
      double precision                                                          
     &     dqx, dqy, dqz, dux, dvx, dwx, dtx, csig(10,27),                      
     &     aux_stress(9,8), ceps_gp(9,27), aux_strain(9,8),                     
     &     dstrain_x1(9), daux_strain_x1(9,8), dcijkl_x1(3),                    
     &     du11_aux(8),  du21_aux(8),  du31_aux(8),                             
     &     du111_aux(8), du112_aux(8), du113_aux(8),                            
     &     du211_aux(8), du212_aux(8), du213_aux(8),                            
     &     du311_aux(8), du312_aux(8), du313_aux(8),                            
     &     elem_alpha(6), dalpha_x1(6), point_temp, point_q, weight,            
     &     iterm(8,8)                                                           
      logical process_temps, fgm_e, fgm_nu, debug                               
c                                                                               
c             local variables                                                   
c                                                                               
      integer i, j                                                              
      double precision                                                          
     &     temp1, temp2, temp3, temp4, zero, half, two                          
c                                                                               
      data zero, half, two                                                      
     & / 0.d0, 0.5d0, 2.d0 /                                                    
c                                                                               
c                                                                               
        if( debug ) write(out,1000) elemno, ptno                                
c                                                                               
c              calculate i-integral terms.                                      
c                                                                               
c              term1 = stress * deriv of aux displacement * deriv of q          
c                                                                               
        temp1 = zero                                                            
        temp2 = zero                                                            
        temp3 = zero                                                            
        temp1 = csig(1,ptno)*dqx + csig(4,ptno)*dqy + csig(7,ptno)*dqz          
        temp2 = csig(2,ptno)*dqx + csig(5,ptno)*dqy + csig(8,ptno)*dqz          
        temp3 = csig(3,ptno)*dqx + csig(6,ptno)*dqy + csig(9,ptno)*dqz          
        do j=1,8                                                                
           iterm(1,j) = iterm(1,j)                                              
     &                + weight * (   du11_aux(j) * temp1                        
     &                             + du21_aux(j) * temp2                        
     &                             + du31_aux(j) * temp3 )                      
           if( debug ) then                                                     
              if( j.eq.2.or.j.eq.7 )                                            
     &        write(out,1010) j, du11_aux(j), du21_aux(j),                      
     &                        du31_aux(j), ptno,                                
     &                        (csig(i,ptno),i=1,9),                             
     &                        dqx, dqy, dqz, temp1,                             
     &                        temp2, temp3, iterm(1,j)                          
           end if                                                               
        end do                                                                  
c                                                                               
c              term2 = aux stress * deriv of displacement * deriv of q          
c                                                                               
        if( debug ) write(out,1020)                                             
        do j=1,8                                                                
           temp1 = zero                                                         
           temp2 = zero                                                         
           temp3 = zero                                                         
           temp1 = aux_stress(1,j)*dqx + aux_stress(4,j)*dqy                    
     &           + aux_stress(7,j)*dqz                                          
           temp2 = aux_stress(2,j)*dqx + aux_stress(5,j)*dqy                    
     &           + aux_stress(8,j)*dqz                                          
           temp3 = aux_stress(3,j)*dqx + aux_stress(6,j)*dqy                    
     &           + aux_stress(9,j)*dqz                                          
           iterm(2,j) = iterm(2,j) + weight*( dux * temp1                       
     &                + dvx * temp2 + dwx * temp3 )                             
           if( debug ) then                                                     
              if( j.eq.2.or.j.eq.7 )                                            
     &        write(out,1030) j, ptno,                                          
     &                        (aux_stress(i,j),i=1,9),                          
     &                        dux, dvx, dwx, dqx, dqy, dqz,                     
     &                        weight, temp1, temp2, temp3,                      
     &                        iterm(2,j)                                        
           end if                                                               
        end do                                                                  
c                                                                               
c              term3 = aux stress * strain * deriv of q                         
c                                                                               
        if( debug ) write(out,1040)                                             
        do j=1,8                                                                
           iterm(3,j) = iterm(3,j) - weight * dqx                               
     &                * ( aux_stress(1,j) * ceps_gp(1,ptno)                     
     &                +   aux_stress(2,j) * ceps_gp(2,ptno)                     
     &                +   aux_stress(3,j) * ceps_gp(3,ptno)                     
     &                +   aux_stress(4,j) * ceps_gp(4,ptno)                     
     &                +   aux_stress(5,j) * ceps_gp(5,ptno)                     
     &                +   aux_stress(6,j) * ceps_gp(6,ptno)                     
     &                +   aux_stress(7,j) * ceps_gp(7,ptno)                     
     &                +   aux_stress(8,j) * ceps_gp(8,ptno)                     
     &                +   aux_stress(9,j) * ceps_gp(9,ptno) )                   
           if( debug ) then                                                     
              if( j.eq.2.or.j.eq.7 )                                            
     &        write(out,1050) j, ptno,                                          
     &                        (aux_stress(i,j),i=1,9),                          
     &                        (ceps_gp(i,ptno),i=1,9),                          
     &                         dqx, dqy, dqz, weight, iterm(3,j)                
           end if                                                               
        end do                                                                  
c                                                                               
c             the following is an equivalent calculation of term3 given         
c             by term3 = stress * aux strain * deriv of q.                      
c             this permits a check of the auxiliary strains and the             
c             other version of term3.                                           
c                                                                               
c        if( debug ) write(out,1040)                                            
c        do j=1,8                                                               
c           iterm(3,j) = iterm(3,j) - weight * dqx                              
c     &                * ( csig(1,ptno) * aux_strain(1,j)                       
c     &                +   csig(2,ptno) * aux_strain(2,j)                       
c     &                +   csig(3,ptno) * aux_strain(3,j)                       
c     &                +   csig(4,ptno) * aux_strain(4,j)                       
c     &                +   csig(5,ptno) * aux_strain(5,j)                       
c     &                +   csig(6,ptno) * aux_strain(6,j)                       
c     &                +   csig(7,ptno) * aux_strain(7,j)                       
c     &                +   csig(8,ptno) * aux_strain(8,j)                       
c     &                +   csig(9,ptno) * aux_strain(9,j) )                     
c            if( debug ) then                                                   
c               write(out,1055) j, ptno,                                        
c     &                         (csig(i,ptno),i=1,9),                           
c     &                         (aux_strain(i,j),i=1,9),                        
c     &                         dqx, dqy, dqz, weight, iterm(3,j)               
c            end if                                                             
c        end do                                                                 
c                                                                               
c             terms 4-6 are zero for homogeneous materials                      
c                                                                               
        if( .not. fgm_e .and. .not. fgm_nu ) go to 1111                         
c                                                                               
c             term4 = stress * 2nd deriv of aux displ * q                       
c                                                                               
        if ( debug ) write(out,1060)                                            
        do j=1,8                                                                
           iterm(4,j) = iterm(4,j) + weight * point_q                           
     &                * (   csig(1,ptno) * du111_aux(j)                         
     &                    + csig(2,ptno) * du211_aux(j)                         
     &                    + csig(3,ptno) * du311_aux(j)                         
     &                    + csig(4,ptno) * du112_aux(j)                         
     &                    + csig(5,ptno) * du212_aux(j)                         
     &                    + csig(6,ptno) * du312_aux(j)                         
     &                    + csig(7,ptno) * du113_aux(j)                         
     &                    + csig(8,ptno) * du213_aux(j)                         
     &                    + csig(9,ptno) * du313_aux(j) )                       
           if( debug ) then                                                     
              if( j.eq.2.or.j.eq.7 )                                            
     &        write(out,1070) j, ptno,                                          
     &             (csig(i,ptno),i=1,9),                                        
     &             du111_aux(j), du211_aux(j), du311_aux(j),                    
     &             du112_aux(j), du212_aux(j), du312_aux(j),                    
     &             du113_aux(j), du213_aux(j), du313_aux(j),                    
     &             weight, point_q, iterm(4,j)                                  
           end if                                                               
        end do                                                                  
c                                                                               
c              term5 = stress * x1 deriv of aux strain * q                      
c                                                                               
        if( debug ) write(out,1080)                                             
        do j=1,8                                                                
           iterm(5,j) = iterm(5,j) - weight * point_q                           
     &                * (   csig(1,ptno) * daux_strain_x1(1,j)                  
     &                    + csig(2,ptno) * daux_strain_x1(2,j)                  
     &                    + csig(3,ptno) * daux_strain_x1(3,j)                  
     &                    + csig(4,ptno) * daux_strain_x1(4,j)                  
     &                    + csig(5,ptno) * daux_strain_x1(5,j)                  
     &                    + csig(6,ptno) * daux_strain_x1(6,j)                  
     &                    + csig(7,ptno) * daux_strain_x1(7,j)                  
     &                    + csig(8,ptno) * daux_strain_x1(8,j)                  
     &                    + csig(9,ptno) * daux_strain_x1(9,j) )                
           if( debug ) then                                                     
              if( j.eq.2.or.j.eq.7 )                                            
     &        write(out,1090) j, ptno,                                          
     &                        (csig(i,ptno),i=1,9),                             
     &                        (daux_strain_x1(i,j),i=1,9),                      
     &                        weight, iterm(5,j)                                
           end if                                                               
        end do                                                                  
c                                                                               
c             term6 = dcijkl_x1 * mechanical strain * aux strain * q            
c             an additional term is necessary for thermal loading.              
c                                                                               
        if( debug ) write(out,1100)                                             
        do j=1,8                                                                
           temp1 = zero                                                         
           temp2 = zero                                                         
           temp3 = zero                                                         
           temp4 = zero                                                         
           temp1 = aux_strain(1,j)                                              
     &           * (   dcijkl_x1(1)                                             
     &               * (ceps_gp(1,ptno) - elem_alpha(1)*point_temp)             
     &             +   dcijkl_x1(2)                                             
     &               * (ceps_gp(5,ptno) - elem_alpha(2)*point_temp)             
     &             +   dcijkl_x1(2)                                             
     &               * (ceps_gp(9,ptno) - elem_alpha(3)*point_temp) )           
c                                                                               
           temp2 = aux_strain(5,j)                                              
     &           * (   dcijkl_x1(2)                                             
     &               * (ceps_gp(1,ptno) - elem_alpha(1)*point_temp)             
     &             +   dcijkl_x1(1)                                             
     &               * (ceps_gp(5,ptno) - elem_alpha(2)*point_temp)             
     &             +   dcijkl_x1(2)                                             
     &               * (ceps_gp(9,ptno) - elem_alpha(3)*point_temp) )           
c                                                                               
           temp3 = aux_strain(9,j)                                              
     &           * (   dcijkl_x1(2)                                             
     &               * (ceps_gp(1,ptno) - elem_alpha(1)*point_temp)             
     &             +   dcijkl_x1(2)                                             
     &               * (ceps_gp(5,ptno) - elem_alpha(2)*point_temp)             
     &             +   dcijkl_x1(1)                                             
     &               * (ceps_gp(9,ptno) - elem_alpha(3)*point_temp) )           
c                                                                               
           temp4 = two * dcijkl_x1(3)                                           
     &           * ( aux_strain(2,j) * ceps_gp(2,ptno)                          
     &           +   aux_strain(3,j) * ceps_gp(3,ptno)                          
     &           +   aux_strain(6,j) * ceps_gp(6,ptno) )                        
c                                                                               
           iterm(6,j) = iterm(6,j)                                              
     &                -   weight * point_q                                      
     &                  * ( temp1 + temp2 + temp3 + temp4 )                     
           if( debug ) then                                                     
              if( j.eq.2.or.j.eq.7 )                                            
     &        write(out,1110) j, (dcijkl_x1(i),i=1,3),                          
     &                        (ceps_gp(i,ptno),i=1,9),                          
     &                        (aux_strain(i,j),i=1,9),                          
     &                        temp1, temp2, temp3, temp4,                       
     &                        weight, point_q,                                  
     &                        iterm(6,j)                                        
           end if                                                               
        end do                                                                  
c                                                                               
 1111   continue                                                                
c                                                                               
c             term7 = aux stress * alpha * dtemp_x1 * q                         
c                   + aux stress * dalpha_x1 * temp * q                         
c             (the following code for iterm7 has not been verified.)            
c                                                                               
        if ( .not. process_temps ) goto 2222                                    
c                                                                               
c        if( debug ) write(out,1120)                                            
c        do j=1,8                                                               
c           temp1 = aux_stress(1,j) * elem_alpha(1)                             
c     &           + aux_stress(5,j) * elem_alpha(2)                             
c     &           + aux_stress(9,j) * elem_alpha(3)                             
c     &           + elem_alpha(4) * half * ( aux_stress(2,j)                    
c     &                                    + aux_stress(4,j) )                  
c     &           + elem_alpha(5) * half * ( aux_stress(6,j)                    
c     &                                    + aux_stress(8,j) )                  
c     &           + elem_alpha(6) * half * ( aux_stress(3,j)                    
c     &                                    + aux_stress(7,j) )                  
c           temp1 = temp1 * dtx * weight * point_q                              
c                                                                               
c           temp2 = aux_stress(1,j) * dalpha_x1(1)                              
c     &           + aux_stress(5,j) * dalpha_x1(2)                              
c     &           + aux_stress(9,j) * dalpha_x1(3)                              
c     &           + dalpha_x1(4) * half * ( aux_stress(2,j)                     
c     &                                   + aux_stress(4,j) )                   
c     &           + dalpha_x1(5) * half * ( aux_stress(6,j)                     
c     &                                   + aux_stress(8,j) )                   
c     &           + dalpha_x1(6) * half * ( aux_stress(3,j)                     
c     &                                   + aux_stress(7,j) )                   
c           temp2 = temp2 * point_temp * weight * point_q                       
c                                                                               
c           iterm(7,j) = iterm(7,j) + temp1 + temp2                             
c                                                                               
c           if( debug ) then                                                    
c              write(out,1130) j, ptno,                                         
c     &             (aux_stress(i,j),i=1,9),                                    
c     &             (elem_alpha(i),i=1,6),                                      
c     &             (dalpha_x1(i),i=1,6),                                       
c     &             dtx, point_temp, point_q,                                   
c     &             weight, temp1, temp2,                                       
c     &             iterm(7,j)                                                  
c           end if                                                              
c        end do                                                                 
c                                                                               
 2222   continue                                                                
c                                                                               
        return                                                                  
c                                                                               
 1000 format(////,"interaction integral terms: element",2x,i7,                  
     &       2x,"point",2x,i2,//,'i term 1:')                                   
 1010 format(/,5x,'j = ',i1,/,10x,                                              
     &       'du11_aux du21_aux du31_aux',/,10x,3(e11.4,2x),                    
     &       //,10x,'stress tensor at gp ', i2,/,                               
     &       3(10x,3(e11.4,2x),/),                                              
     &       /,10x,'dqx dqy dqz',/,10x,3(e11.4,2x),                             
     &       //,10x,'temp1 temp2 temp3 iterm1',/,10x,4(e11.4,2x))               
 1020 format(//,'i term 2:')                                                    
 1030 format(/,5x,'j = ',i1,/,10x,'aux stress tensor at gp ',                   
     &       i2,/,3(10x,3(e11.4,2x),/),/,10x,'dux  dvx  dwx',/,10x,             
     &       3(e11.4,2x),//,10x,'dqx  dqy  dqz',/,10x,3(e11.4,2x),//,           
     &       10x,'weight',/,10x,e11.4,//,10x,'temp1 temp2 temp3 iterm2',        
     &       /,10x,4(e13.6,2x))                                                 
 1040 format(//,'i term 3:')                                                    
 1050 format(/,5x,'j = ',i1,/,10x,'aux stress tensor at gp ',                   
     &       i2,/,3(10x,3(e11.4,2x),/),/,10x,'strain tensor',/,                 
     &       3(10x,3(e11.4,2x),/),/,10x,'dqx dqy dqz weight',/,10x,             
     &       4(e11.4,2x),//,10x,'iterm3',/,10x,e11.4)                           
 1055 format(/,5x,'j = ',i1,/,10x,'stress tensor at gp ',                       
     &       i2,/,3(10x,3(e11.4,2x),/),/,10x,'aux strain tensor',/,             
     &       3(10x,3(e11.4,2x),/),/,10x,'dqx dqy dqz weight',/,10x,             
     &       4(e11.4,2x),//,10x,'iterm3',/,10x,e11.4)                           
 1060 format(//,'i term 4:')                                                    
 1070 format(/,5x,'j = ',i1,/,10x,'stress tensor at gp ',                       
     &       i2,/,3(10x,3(e11.4,2x),/),/,10x,                                   
     &       'du111_aux du211_aux du311_aux',/,10x,3(e11.4,2x),//,10x,          
     &       'du112_aux du212_aux du312_aux',/,10x,3(e11.4,2x),//,10x,          
     &       'du113_aux du213_aux du313_aux',/,10x,3(e11.4,2x),//,10x,          
     &       'weight point_q iterm4',/,10x,3(e11.4,2x))                         
 1080 format(//,'i term 5:')                                                    
 1090 format(/,5x,'j = ',i1,/,10x,'stress tensor at gp ',                       
     &       i2,/,3(10x,3(e11.4,2x),/),/,10x,                                   
     &       'daux_strain_x1 tensor',/,3(10x,3(e11.4,2x),/),/,                  
     &       10x,'weight iterm5',/,10x,2(e11.4,2x))                             
 1100 format(//,'i term 6:')                                                    
 1110 format(/,5x,'j = ',i1,/,10x,'dcijkl_x1(1:3)',/,10x,                       
     &       3(e11.4,2x),//,10x,'strain tensor',/,                              
     &       3(10x,3(e11.4,2x),/),/,10x,'aux strain tensor',                    
     &       /,3(10x,3(e11.4,2x),/),/,10x,'temp1 temp2 temp3 temp4',            
     &       /,10x,4(e11.4,2x),//,10x,'weight point_q, iterm6',/,10x,           
     &       3(e11.4,2x))                                                       
 1120 format(//,'i term 7:')                                                    
 1130 format(/,5x,'j = ',i1,/,10x,'aux stress tensor at gp ',                   
     &       i2,/,3(10x,3(e13.6,2x),/),//,10x,'elem_alpha 1-6',                 
     &       /,2(10x,3(e13.6,2x),/),//,10x,'dalpha_x1 1-6',                     
     &       /,2(10x,3(e13.6,2x),/),//,10x,'dtx point_temp point_q',            
     &       /,10x,3(e13.6,2x),//,10x,'weight temp1 temp2',/,10x,               
     &       3(e13.6,2x),//,10x,'iterm7',/,10x,e13.6)                           
 1140 format(//,'i term 9:')                                                    
 1150 format(/,5x,'j = ',i1,/,10x,'aux stress tensor at gp ',                   
     &       i2,/,3(10x,3(e11.4,2x),/),/,10x,'dstrain_x1',                      
     &       /,3(10x,3(e11.4,2x),/),/,10x,'weight point_q, iterm9',/,           
     &       10x,3(e11.4,2x))                                                   
c                                                                               
      end                                                                       
c                                                                               
