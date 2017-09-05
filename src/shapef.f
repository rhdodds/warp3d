c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shapef                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/11/2017 rhd             *          
c     *                                                              *          
c     *     given the element type and a point in parametric         *          
c     *     space, return the element shape functions evaluated      *          
c     *     at that point.                                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine shapef( etype, xi, eta, zeta, n )                              
      implicit integer (a-z)                                                    
      double precision ::  xi, eta, zeta, n(*)                                                  
c                                                                               
c                                                                               
      if ( etype .le. 0 .or. etype .gt. 19 ) then                               
         write(*,9000) etype                                                    
         call die_abort                                                         
      end if                                                                    
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600, 700, 800, 900,                      
     &       1000, 1100, 1200, 1300, 1400, 1500, 1600,
     &       1700, 1800, 1900 ), etype                  
c                                                                               
c                    20-node                                                    
c                                                                               
 100  continue                                                                  
      call shape1( xi, eta, zeta, n )                                           
      return                                                                    
c                                                                               
c                     8-node                                                    
c                                                                               
 200  continue                                                                  
      call shape2( xi, eta, zeta, n )                                           
      return                                                                    
c                                                                               
c                    12-node                                                    
c                                                                               
 300  continue                                                                  
      call shape3( xi, eta, zeta, n )                                           
      return                                                                    
c                                                                               
c                    15-node                                                    
c                                                                               
 400  continue                                                                  
      call shape4( xi, eta, zeta, n )                                           
      return                                                                    
c                                                                               
c                    9-node                                                     
c                                                                               
 500  continue                                                                  
      call shape5( xi, eta, zeta, n )                                           
      return                                                                    
c                                                                               
c                    10-node tetrahedron, "tet10"                               
c                                                                               
 600  continue                                                                  
      call shape6( xi, eta, zeta, n )                                           
      return                                                                    
c                                                                               
c                    15-node wedge, "wedge15"                                   
c                                                                               
 700  continue                                                                  
         write(*,9100) etype                                                    
         call die_abort                                                         
cADD SUBR      call shape7( xi, eta, zeta, n )                                  
      return                                                                    
c                                                                               
c                    6-node planar triangle, "tri6"                             
c                                                                               
 800  continue                                                                  
         write(*,9100) etype                                                    
         call die_abort                                                         
cADD SUBR      call shape8( xi, eta, zeta, n )                                  
      return                                                                    
c                                                                               
c                    8-node planar quadrillateral, "quad8"                      
c                                                                               
 900  continue                                                                  
      call shape9( xi, eta, n )                                                 
      return                                                                    
c                                                                               
c                    8-node axisymmetric quadrillateral, "axiquad8"             
c                                                                               
1000  continue                                                                  
      call shape10( xi, eta, n )                                                
      return                                                                    
c                                                                               
c                    6-node planar quadrillateral, "axitri6"                    
c                                                                               
1100  continue                                                                  
      call shape11( xi, eta, zeta, n )                                          
      return                                                                    
c                                                                               
c                    8-noded interface element, "inter_8"                       
1200  continue                                                                  
      call shape12( xi, eta, n )                                                
      return                                                                    
c                                                                               
c                    4-noded tetrahedron, "tet4"                                
1300  continue                                                                  
      call shape13( xi, eta, zeta, n )                                          
      return                                                                    
c                                                                               
c                    6-noded triangular interface element, "trint6"             
1400  continue                                                                  
      call shape14( xi, eta, n )                                                
      return                                                                    
c                                                                               
c                    12-noded triangular interface element, "trint12"           
1500  continue                                                                  
      call shape15( xi, eta, n )                                                
      return                                                                    
c
c                    4-node quadrillateral element
c
1600  continue
      call shape16( xi, eta, n )
      return
c
c                    3-node line. should not be called in code
c
1700  continue
      write(*,9100) etype                                                    
      call die_abort                                                         

c                                                                               
c                    2-node bar
c                                                                               
1800  continue 
      n(1) = 0.5d0 
      n(2) = 0.5d0                                                                
      return                                                                    
c                                                                               
c                    2-node link
c                                                                               
1900  continue 
      n(1) = 0.5d0 
      n(2) = 0.5d0                                                                
      return                                                                    
c                                                                               
c                                                                               
 9000 format('> FATAL ERROR: shaepf, etype: ',i10,//,                           
     &       '               job aborted' )                                     
 9100 format('> FATAL ERROR: shaepf, etype: ',i10,//,                           
     &       '               element not yet implemented',/,                    
     &       '               shape functions not defined' )                     
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape1                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 19/07/97                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given gauss point    *          
c     *     for 20-node element                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine shape1(r,s,t,n)                                                
      implicit integer (a-z)                                                    
      double precision                                                          
     &     n(*),r,s,t,rp1,rm1,sp1,sm1,tp1,tm1,                                  
     &     one,half,fourth,eighth                                               
c                                                                               
      data one,half,fourth,eighth /1.0,0.5,0.25,0.125/                          
c                                                                               
      rp1 = one+r                                                               
      rm1 = one-r                                                               
      sp1 = one+s                                                               
      sm1 = one-s                                                               
      tp1 = one+t                                                               
      tm1 = one-t                                                               
c                                                                               
      n(1)  = -rm1*sm1*(rp1 + sp1 - t)*tp1*eighth                               
      n(2)  = -rm1*sm1*(rp1 + sp1 + t)*tm1*eighth                               
      n(3)  = -rm1*sp1*(rp1 + sm1 + t)*tm1*eighth                               
      n(4)  = -rm1*sp1*(rp1 + sm1 - t)*tp1*eighth                               
      n(5)  = -rp1*sm1*(rm1 + sp1 - t)*tp1*eighth                               
      n(6)  = -rp1*sm1*(rm1 + sp1 + t)*tm1*eighth                               
      n(7)  = -rp1*sp1*(rm1 + sm1 + t)*tm1*eighth                               
      n(8)  = -rp1*sp1*(rm1 + sm1 - t)*tp1*eighth                               
      n(9)  = rm1*sm1*tm1*tp1*fourth                                            
      n(10) = rm1*sm1*sp1*tm1*fourth                                            
      n(11) = rm1*sp1*tm1*tp1*fourth                                            
      n(12) = rm1*sm1*sp1*tp1*fourth                                            
      n(13) = rp1*sm1*tm1*tp1*fourth                                            
      n(14) = rp1*sm1*sp1*tm1*fourth                                            
      n(15) = rp1*sp1*tm1*tp1*fourth                                            
      n(16) = rp1*sm1*sp1*tp1*fourth                                            
      n(17) = rm1*rp1*sm1*tp1*fourth                                            
      n(18) = rm1*rp1*sm1*tm1*fourth                                            
      n(19) = rm1*rp1*sp1*tm1*fourth                                            
      n(20) = rm1*rp1*sp1*tp1*fourth                                            
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape2                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/18/89                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given gauss point    *          
c     *     for element l3disop                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine shape2( xi, eta, zeta, n )                                     
      implicit integer (a-z)                                                    
      double precision                                                          
     &     xi,eta,zeta,n(*),xp,ep,zp,xm,em,zm,one,eighth                        
      data one, eighth /1.0, 0.125/                                             
c                                                                               
c                       set basic parameters                                    
c            z                                                                  
      xp= one+xi                                                                
      ep= one+eta                                                               
      zp= one+zeta                                                              
      xm= one-xi                                                                
      em= one-eta                                                               
      zm= one-zeta                                                              
c                                                                               
c                       corner nodes                                            
c                                                                               
      n(1)= xm*em*zp*eighth                                                     
      n(2)= xm*em*zm*eighth                                                     
      n(3)= xm*ep*zm*eighth                                                     
      n(4)= xm*ep*zp*eighth                                                     
      n(5)= xp*em*zp*eighth                                                     
      n(6)= xp*em*zm*eighth                                                     
      n(7)= xp*ep*zm*eighth                                                     
      n(8)= xp*ep*zp*eighth                                                     
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape3                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 20/07/97                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given gauss point    *          
c     *     for 12-node transition element                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine shape3(r,s,t,n)                                                
      implicit integer (a-z)                                                    
      double precision                                                          
     &     n(*),r,s,t,rp1,rm1,sp1,sm1,tp1,tm1,spt,smt,                          
     &     one,fourth,eighth                                                    
c                                                                               
      data one,fourth,eighth /1.0,0.25,0.125/                                   
c                                                                               
      rp1 = one+r                                                               
      rm1 = one-r                                                               
      sp1 = one+s                                                               
      sm1 = one-s                                                               
      tp1 = one+t                                                               
      tm1 = one-t                                                               
      spt = s + t                                                               
      smt = s - t                                                               
c                                                                               
      n(1)  = -rm1*sm1*(sp1 + smt*t)*eighth                                     
      n(2)  = -rm1*sm1*(sp1 - spt*t)*eighth                                     
      n(3)  = -rm1*sp1*(sm1 + smt*t)*eighth                                     
      n(4)  = -rm1*sp1*(sm1 - spt*t)*eighth                                     
      n(5)  = rp1*sm1*tp1*eighth                                                
      n(6)  = rp1*sm1*tm1*eighth                                                
      n(7)  = rp1*sp1*tm1*eighth                                                
      n(8)  = rp1*sp1*tp1*eighth                                                
      n(9)  = rm1*sm1*tm1*tp1*fourth                                            
      n(10) = rm1*sm1*sp1*tm1*fourth                                            
      n(11) = rm1*sp1*tm1*tp1*fourth                                            
      n(12) = rm1*sm1*sp1*tp1*fourth                                            
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape4                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 19/07/97                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given gauss point    *          
c     *     for 15-node transition element                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine shape4(r,s,t,n)                                                
      implicit integer (a-z)                                                    
      double precision                                                          
     &     n(*),r,s,t,rp1,rm1,sp1,sm1,tp1,tm1,spt,smt,rpt,rmt,                  
     &     one,fourth,eighth                                                    
c                                                                               
      data one,fourth,eighth /1.0,0.25,0.125/                                   
c                                                                               
      rp1 = one+r                                                               
      rm1 = one-r                                                               
      sp1 = one+s                                                               
      sm1 = one-s                                                               
      tp1 = one+t                                                               
      tm1 = one-t                                                               
      spt = s + t                                                               
      smt = s - t                                                               
      rpt = r + t                                                               
      rmt = r - t                                                               
c                                                                               
      n(1)  = -rm1*sm1*(rp1 + sp1 - t)*tp1*eighth                               
      n(2)  = -rm1*sm1*(rp1 + sp1 + t)*tm1*eighth                               
      n(3)  = -rm1*sp1*(sm1 + smt*t)*eighth                                     
      n(4)  = -rm1*sp1*(sm1 - spt*t)*eighth                                     
      n(5)  = -rp1*sm1*(rm1 - rpt*t)*eighth                                     
      n(6)  = -rp1*sm1*(rm1 + rmt*t)*eighth                                     
      n(7)  = rp1*sp1*tm1*eighth                                                
      n(8)  = rp1*sp1*tp1*eighth                                                
      n(9)  = rm1*sm1*tm1*tp1*fourth                                            
      n(10) = rm1*sm1*sp1*tm1*fourth                                            
      n(11) = rm1*sp1*tm1*tp1*fourth                                            
      n(12) = rm1*sm1*sp1*tp1*fourth                                            
      n(13) = rp1*sm1*tm1*tp1*fourth                                            
      n(14) = rm1*rp1*sm1*tp1*fourth                                            
      n(15) = rm1*rp1*sm1*tm1*fourth                                            
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape5                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 8/22/977                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given gauss point    *          
c     *     for 9-node transition element                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine shape5(r,s,t,n)                                                
      implicit integer (a-z)                                                    
      double precision                                                          
     &     n(*),r,s,t,rp1,rm1,sp1,sm1,tp1,tm1,                                  
     &     one,half,fourth,eighth                                               
c                                                                               
      data one,half,fourth,eighth /1.0,0.5,0.25,0.125/                          
c                                                                               
      rp1 = one+r                                                               
      rm1 = one-r                                                               
      sp1 = one+s                                                               
      sm1 = one-s                                                               
      tp1 = one+t                                                               
      tm1 = one-t                                                               
c                                                                               
      n(1) =  rm1*sm1*tp1*t*eighth                                              
      n(2) = -rm1*sm1*tm1*t*eighth                                              
      n(3) =  rm1*sp1*tm1*eighth                                                
      n(4) =  rm1*sp1*tp1*eighth                                                
      n(5) =  rp1*sm1*tp1*eighth                                                
      n(6) =  rp1*sm1*tm1*eighth                                                
      n(7) =  rp1*sp1*tm1*eighth                                                
      n(8) =  rp1*sp1*tp1*eighth                                                
      n(9) =  rm1*sm1*tm1*tp1*fourth                                            
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape6                       *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/24/98                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the shape         *          
c     *     functions at the isoparametric coordinates of the given  *          
c     *     Gauss integration point for the 10-node tetrahedron.     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         evaluate the 10 shape functions at a given Gauss integration          
c         point for the 10-node tetrahedron volume element.                     
c         See subroutine "derivs" for the derivatives of the shape functions.   
c         The 10-node tetrahedron is element type 6 in Warp3D, use the          
c         subroutine name "shape6" for this element type.                       
c                                                                               
c         Note that the shape functions are with respect to the                 
c         four element local tetrahedron coordinates s1,s2,s3,s4.               
c                                                                               
c         Note that the nodes are numbered at the corners first (nodes 1,2,3    
c         on the face 1 and then node 4 at the peak of the tetrahedron)         
c         and then at the edge nodes on face 1 in a right hand rule sense so    
c         that the surface normal in into the element (nodes 5,6,7), then the   
c         edges from face 1 up to node 4 (edge nodes 8,9,10).                   
c         This should be the usual numbering order and is the same as           
c         used in abaqus.                                                       
c                                                                               
c         Variables:                                                            
c                                                                               
c         s1,s2,s3,s4 = tetrahedron natural coordinates, vary from 0 to 1,      
c                       These coordinates will be at a given Gauss              
c                       integration point.                                      
c         Note that s1 is not passed in but computed locally.                   
c        Properties of coordinates:  s1 + s2 + s3 + s4 = 1                      
c         choose s1 as dependent, then s1 = 1 - s2 - s3 - s4                    
c                               0 <= s1,s2,s3,s4 <= 1                           
c         Note that the Gauss point coordinates s2,s3,s4 must be consistent     
c         with the coordinates given from the "gauss6" subroutine.              
c                                                                               
c         q(numnode) = shape functions evaluated at the given Gauss point,      
c                      evaluated in terms of the 4 tetrahedral coordinates.     
c         (where numnode=10 for the 10-node tetrahedron element)                
c                                                                               
c         Note that since the 4 tetrahedron coordinates are not independent,    
c         s1 = 1 - s2 - s3 - s4 has been substituted into the shape functions;  
c         s1 appears in the shape functions for simplicity of evaluating the    
c         values given a Gauss integration point coordinate set.                
c                                                                               
      subroutine shape6( s2, s3, s4, q )                                        
      implicit none                                                             
      double precision                                                          
     &         s2,s3,s4,q(*)                                                    
c                                                                               
c        s1 = fourth tetrahedron natural coordinate, used here to keep the      
c             simplicity of the functions below                                 
c                                                                               
      double precision                                                          
     &     one, two, four, s1                                                   
c                                                                               
      one  = 1.0D0                                                              
      two  = 2.0D0                                                              
      four = 4.0D0                                                              
      s1 = one - s2 - s3 - s4                                                   
c                                                                               
c        evaluate the shape functions at the current Gauss point,               
c        note that s1 has been used as the dependent coordinate.                
c                                                                               
      q(1) = two*s1*s1 - s1                                                     
      q(2) = two*s2*s2 - s2                                                     
      q(3) = two*s3*s3 - s3                                                     
      q(4) = two*s4*s4 - s4                                                     
      q(5) = four*s1*s2                                                         
      q(6) = four*s2*s3                                                         
      q(7) = four*s1*s3                                                         
      q(8) = four*s1*s4                                                         
      q(9) = four*s2*s4                                                         
      q(10) = four*s3*s4                                                        
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape9                       *          
c     *                                                              *          
c     *                       written by : mcw                       *          
c     *                          (after gvt)                         *          
c     *                                                              *          
c     *                   last modified : 08/21/03                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given  Gauss         *          
c     *     integration point for the 8-node quadrilateral.          *          
c     *                                                              *          
c     *     this subroutine is currently identical to shape10        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine shape9( s1, s2, q )                                            
      implicit none                                                             
      double precision                                                          
     &         s1, s2, q(*)                                                     
c                                                                               
      double precision                                                          
     &         one, two, four                                                   
c                                                                               
      one  = 1.0D0                                                              
      two  = 2.0D0                                                              
      four = 4.0D0                                                              
c                                                                               
c        evaluate the shape functions at the current Gauss point                
c                                                                               
      q(1) = (one - s1)*(one - s2)*(-s1 - s2 - one)/four                        
      q(2) = (one + s1)*(one - s2)*(s1 - s2 - one)/four                         
      q(3) = (one + s1)*(one + s2)*(s1 + s2 - one)/four                         
      q(4) = (one - s1)*(one + s2)*(-s1 + s2 - one)/four                        
      q(5) = (one - s1*s1)*(one - s2)/two                                       
      q(6) = (one + s1)*(one - s2*s2)/two                                       
      q(7) = (one - s1*s1)*(one + s2)/two                                       
      q(8) = (one - s1)*(one - s2*s2)/two                                       
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape10                      *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/24/98                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given  Gauss         *          
c     *     integration point for the 8-node axisymmetric            *          
c     *     quadrillateral.                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine shape10( s1, s2, q )                                           
      implicit none                                                             
      double precision                                                          
     &         s1, s2, q(*)                                                     
c                                                                               
      double precision                                                          
     &         one, two, four                                                   
c                                                                               
      one  = 1.0D0                                                              
      two  = 2.0D0                                                              
      four = 4.0D0                                                              
c                                                                               
c        evaluate the shape functions at the current Gauss point                
c                                                                               
      q(1) = (one - s1)*(one - s2)*(-s1 - s2 - one)/four                        
      q(2) = (one + s1)*(one - s2)*(s1 - s2 - one)/four                         
      q(3) = (one + s1)*(one + s2)*(s1 + s2 - one)/four                         
      q(4) = (one - s1)*(one + s2)*(-s1 + s2 - one)/four                        
      q(5) = (one - s1*s1)*(one - s2)/two                                       
      q(6) = (one + s1)*(one - s2*s2)/two                                       
      q(7) = (one - s1*s1)*(one + s2)/two                                       
      q(8) = (one - s1)*(one - s2*s2)/two                                       
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape11                      *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/24/98                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given Gauss          *          
c     *     integration point for the 6-node axisymmetric triangle.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine shape11( s1, s2, s3, q )                                       
      implicit none                                                             
      double precision                                                          
     &          s1, s2, s3, q(*)                                                
c                                                                               
      double precision                                                          
     &          two, four                                                       
c                                                                               
      two  = 2.0D0                                                              
      four = 4.0D0                                                              
c                                                                               
c          evaluate the shape functions at the current Gauss point              
c                                                                               
      q(1) = two*s1*s1 - s1                                                     
      q(2) = two*s2*s2 - s2                                                     
      q(3) = two*s3*s3 - s3                                                     
      q(4) = four*s1*s2                                                         
      q(5) = four*s2*s3                                                         
      q(6) = four*s1*s3                                                         
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape12                      *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                                                              *          
c     *                   last modified : 5/26/99                    *          
c     *                                                              *          
c     *     this subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given gauss point    *          
c     *     for element inter_8                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine shape12( xi, eta, n )                                          
      implicit integer (a-z)                                                    
      double precision                                                          
     &     xi,eta,n(*),xp,ep,xm,em,one,fourth                                   
      data one, fourth /1.0, 0.25/                                              
c                                                                               
c                 set basic parameters                                          
c     A plane interface element, nodes 1-4 are bottom                           
c     plane nodes, and, 5-8 are top plane nodes                                 
c                                                                               
        xp = one+xi                                                             
        ep = one+eta                                                            
        xm = one-xi                                                             
        em = one-eta                                                            
c                                                                               
c                       corner nodes                                            
c                                                                               
      n(1) = xm*em*fourth                                                       
      n(2) = xp*em*fourth                                                       
      n(3) = xp*ep*fourth                                                       
      n(4) = xm*ep*fourth                                                       
      n(5) = n(1)                                                               
      n(6) = n(2)                                                               
      n(7) = n(3)                                                               
      n(8) = n(4)                                                               
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape13                      *          
c     *                                                              *          
c     *                       written by : rau                       *          
c     *                                                              *          
c     *                   last modified : 12/21/00                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the shape         *          
c     *     functions at the isoparametric coordinates of the given  *          
c     *     Gauss integration point for the 4-node tetrahedron.      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         Note that the shape functions are with respect to the                 
c         four element local tetrahedron coordinates s1,s2,s3,s4.               
c                                                                               
c         Variables:                                                            
c                                                                               
c         s1,s2,s3,s4 = tetrahedron natural coordinates, vary from 0 to 1,      
c                       These coordinates will be at a given Gauss              
c                       integration point.                                      
c         Note that s1 is not passed in but computed locally.                   
c        Properties of coordinates:  s1 + s2 + s3 + s4 = 1                      
c         choose s1 as dependent, then s1 = 1 - s2 - s3 - s4                    
c                               0 <= s1,s2,s3,s4 <= 1                           
c         Note that the Gauss point coordinates s2,s3,s4 must be consistent     
c         with the coordinates given from the "gauss6" subroutine.              
c                                                                               
c         q(numnode) = shape functions evaluated at the given Gauss point,      
c                      evaluated in terms of the 4 tetrahedral coordinates.     
c         (where numnode=4 for the 4-node tetrahedron element)                  
c                                                                               
c         Note that since the 4 tetrahedron coordinates are not independent,    
c         s1 = 1 - s2 - s3 - s4 has been substituted into the shape functions;  
c         s1 appears in the shape functions for simplicity of evaluating the    
c         values given a Gauss integration point coordinate set.                
c                                                                               
      subroutine shape13( s2, s3, s4, q )                                       
      implicit none                                                             
      double precision                                                          
     &         s2,s3,s4,q(*)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     one, s1                                                              
c                                                                               
      one  = 1.0D0                                                              
      s1 = one - s2 - s3 - s4                                                   
c                                                                               
c        evaluate the shape functions at the current Gauss point,               
c        note that s1 has been used as the dependent coordinate.                
c                                                                               
      q(1) = s1                                                                 
      q(2) = s2                                                                 
      q(3) = s3                                                                 
      q(4) = s4                                                                 
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape14                      *          
c     *                                                              *          
c     *                       written by : sushovan                  *          
c     *                                                              *          
c     *                   last modified : 2/8/13 rhd                 *          
c     *                                                              *          
c     *     computes shape functions at a given parametric point for *          
c     *     6-noded triangular interface element trint6              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c   Shape functions are in terms of triangle natural coordinates.               
c                                                                               
c   3-noded bottom surface of the 6-noded triangular element:                   
c                                                                               
c                                                                               
c        3                                                                      
c         o                                                                     
c         |  \                                                                  
c         |    \                                                                
c         |      \                                                              
c         |        \                                                            
c         |          \                                                          
c         |            \                                                        
c         |              \                                                      
c         o - - - - - - -  o                                                    
c        1                  2                                                   
c                                                                               
c                                                                               
c      Variables:                                                               
c                                                                               
c       s1,s2,s3 = triangle natural coordinates, vary from 0 to 1,              
c                  example: s1=0 on the triangle side 2-3.                      
c                             =1 at corner 1                                    
c             Properties of coordinates:  s1 + s2 + s3 = 1                      
c                           and 0 <= s1,s2,s3 <= 1                              
c             Note: Only s1 and s2 are input to this subroutine.                
c                   use s3 = 1 - s1 - s2.                                       
c                                                                               
      subroutine shape14( s1, s2, q )                                           
      implicit none                                                             
      double precision                                                          
     &          s1, s2, s3, q(*)                                                
c                                                                               
      double precision                                                          
     &          one, two, four                                                  
c                                                                               
      one  = 1.0d0                                                              
c                                                                               
      s3 = one - s1 - s2                                                        
c                                                                               
c          evaluate the shape functions at the current point                    
c                                                                               
      q(1) = s1                                                                 
      q(2) = s2                                                                 
      q(3) = s3                                                                 
      q(4) = q(1)                                                               
      q(5) = q(2)                                                               
      q(6) = q(3)                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape15                      *          
c     *                                                              *          
c     *                       written by : sushovan                  *          
c     *                                                              *          
c     *                   last modified : 2/8/13  rhd                *          
c     *                                                              *          
c     *     computes values of the shape at a given parametric point *          
c     *     for the 12-noded triangular interface element trint12    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c   Shape functions are in terms of triangle natural coordinates.               
c                                                                               
c   6-noded bottom surface of the 12-noded triangular element:                  
c                                                                               
c       3                                                                       
c         o                                                                     
c         |  \                                                                  
c         |    \                                                                
c         |      \  5                                                           
c       6 o        o                                                            
c         |          \                                                          
c         |            \                                                        
c         |              \                                                      
c         o - - - o - - - o                                                     
c        1        4        2                                                    
c                                                                               
c                                                                               
c      Variables:                                                               
c                                                                               
c       s1,s2,s3 = triangle natural coordinates, vary from 0 to 1,              
c                  example: s1=0 on the triangle side 2-3                       
c                             =0.5 at nodes 4 and 6                             
c                             =1 at corner 1                                    
c             Properties of coordinates:  s1 + s2 + s3 = 1                      
c                           and 0 <= s1,s2,s3 <= 1                              
c             Note: Only s1 and s2 are input to this subroutine.                
c                   use s3 = 1 - s1 - s2.                                       
c                                                                               
                                                                                
      subroutine shape15( s1, s2, q )                                           
      implicit none                                                             
      double precision                                                          
     &          s1, s2, s3, q(*)                                                
c                                                                               
      double precision                                                          
     &          one, two, four                                                  
c                                                                               
      one  = 1.0d0                                                              
      two  = 2.0d0                                                              
      four = 4.0d0                                                              
c                                                                               
      s3 = one - s1 - s2                                                        
c                                                                               
c          evaluate the shape functions at the current point                    
c                                                                               
      q(1) = two*s1*s1 - s1                                                     
      q(2) = two*s2*s2 - s2                                                     
      q(3) = two*s3*s3 - s3                                                     
      q(4) = four*s1*s2                                                         
      q(5) = four*s2*s3                                                         
      q(6) = four*s3*s1                                                         
      q(7) = q(1)                                                               
      q(8) = q(2)                                                               
      q(9) = q(3)                                                               
      q(10) = q(4)                                                              
      q(11) = q(5)                                                              
      q(12) = q(6)                                                              
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine shape16                      *          
c     *                                                              *          
c     *                              by : mcw                        *          
c     *                            (after gvt)                       *          
c     *                                                              *          
c     *                   last modified : 08/21/03                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the shape         *          
c     *     functions at the coordinates of the given  Gauss         *          
c     *     integration point for the 4-node quadrillateral.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine shape16( s1, s2, q )                                           
      implicit none                                                             
      double precision                                                          
     &         s1, s2, q(*)                                                     
c                                                                               
      double precision                                                          
     &         one, four                                                        
c                                                                               
      one  = 1.0D0                                                              
      four = 4.0D0                                                              
c                                                                               
c        evaluate the shape functions at the current Gauss point                
c                                                                               
      q(1) = (one - s1)*(one - s2)/four                                         
      q(2) = (one + s1)*(one - s2)/four                                         
      q(3) = (one + s1)*(one + s2)/four                                         
      q(4) = (one - s1)*(one + s2)/four                                         
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
