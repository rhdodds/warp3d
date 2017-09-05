c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine derivs                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/11/2017 rhd             *          
c     *                                                              *          
c     *     given the element type, integration point coordinates    *          
c     *     return the parametric derivatives for each node          *          
c     *     shape function                                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine derivs( etype, xi, eta, zeta, nxi, neta, nzeta )   
c                  
      implicit none                                                             
      double precision :: xi, eta, zeta, nxi(*), neta(*), nzeta(*)                             
      integer :: etype  
      double precision, parameter :: zero = 0.0d0                                                           
c                                                                               
      if ( etype .le. 0 .or. etype .gt. 19 ) then                               
         write(*,9000) etype                                                    
         call die_abort                                                         
      end if                                                                    
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600, 700, 800, 900,                      
     &       1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
     &       1800, 1900 ), etype                  
c                                                                               
c                    20-node brick                                              
c                                                                               
 100  continue                                                                  
      call deriv1( xi, eta, zeta, nxi, neta, nzeta )                            
      return                                                                    
c                                                                               
c                     8-node brick                                              
c                                                                               
 200  continue                                                                  
      call deriv2( xi, eta, zeta, nxi, neta, nzeta )                            
      return                                                                    
c                                                                               
c                    12-node brick                                              
c                                                                               
 300  continue                                                                  
      call deriv3( xi, eta, zeta, nxi, neta, nzeta )                            
      return                                                                    
c                                                                               
c                    15-node brick                                              
c                                                                               
 400  continue                                                                  
      call deriv4( xi, eta, zeta, nxi, neta, nzeta )                            
      return                                                                    
c                                                                               
c                    9-node brick                                               
c                                                                               
 500  continue                                                                  
      call deriv5( xi, eta, zeta, nxi, neta, nzeta )                            
      return                                                                    
c                                                                               
c                    10-node tetrahedron, "tet10"                               
c                                                                               
 600  continue                                                                  
      call deriv6( xi, eta, zeta, nxi, neta, nzeta )                            
      return                                                                    
c                                                                               
c                    15-node wedge, "wedge15" -- not implemented                
c                                                                               
 700  continue                                                                  
      write(*,9100) etype                                                       
      call die_abort                                                            
c                                                                               
c                    6-node planar triangle, "tri6" -- not implemented          
c                                                                               
 800  continue                                                                  
      write(*,9100) etype                                                       
      call die_abort                                                            
c                                                                               
c                    8-node planar quadrilateral, "quad8"                       
c                    (currently identical to deriv10)                           
c                                                                               
 900  continue                                                                  
      call deriv9( xi, eta, nxi, neta )                                         
      return                                                                    
c                                                                               
c                    8-node axisymmetric quadrilateral, "axiquad8"              
c                                                                               
1000  continue                                                                  
      call deriv10( xi, eta, nxi, neta )                                        
      return                                                                    
c                                                                               
c                    6-node planar quadrilateral, "axitri6"                     
c                                                                               
1100  continue                                                                  
      call deriv11( xi, eta, zeta, nxi, neta )                                  
      return                                                                    
c                                                                               
c                     8-node interface element, "inter_8"                       
c                                                                               
1200  continue                                                                  
      call deriv12( xi, eta, nxi, neta )                                        
      return                                                                    
c                                                                               
c                     4-node tetrahedron, "tet4"                                
c                                                                               
1300  continue                                                                  
      call deriv13( xi, eta, zeta, nxi, neta, nzeta )                           
      return                                                                    
c                                                                               
c                     6-node triangular interface element, "trint6"             
c                                                                               
1400  continue                                                                  
      call deriv14( xi, eta, nxi, neta )                                        
      return                                                                    
c                                                                               
c                     12-node triangular interface element, "trint12"           
c                                                                               
1500  continue                                                                  
      call deriv15( xi, eta, nxi, neta )                                        
      return                                                                    
c                                                                               
c                    4-node planar quadrilateral, "quad4"  -not 
c                    really and element. used in various places                      
c                                                                               
1600  continue                                                                  
      call deriv16( xi, eta, nxi, neta )                                        
      return                                                                    
c                                                                               
c                    3-node line.  not implemented for derivs.          
c                                                                               
1700  continue                                                                  
      write(*,9100) etype                                                       
      call die_abort                                                            
c                                                                               
c                     2-node bar element. just dummy values since never
c                                         used in code
c                                                                               
1800  continue                                                                  
      nxi(1) = zero
      nxi(2) = zero
      neta(1) = zero
      neta(2) = zero
      nzeta(1) = zero
      nzeta(2) = zero
      return                                                                    
c                     2-node link element. just dummy values since never
c                                         used in code
c                                                                               
1900  continue                                                                  
      nxi(1) = zero
      nxi(2) = zero
      neta(1) = zero
      neta(2) = zero
      nzeta(1) = zero
      nzeta(2) = zero
      return                                                                    
c                                                                               
 9000 format('> FATAL ERROR: derivs, etype: ',i10,//,                           
     &       '               job aborted' )                                     
 9100 format('> FATAL ERROR: derivs, etype: ',i10,//,                           
     &       '               element not yet implemented',/,                    
     &       '               shape function derivatives not defined' )          
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv1                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 02/09/13 rhd               *          
c     *                                                              *          
c     *     derivatives of shape functions at the coordinates        *          
c     *     of integration point for 20-node element                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine deriv1(r,s,t,nr,ns,nt)                                         
      implicit integer (a-z)                                                    
      double precision                                                          
     &     nr(*),ns(*),nt(*),rp1,rm1,sp1,sm1,tp1,tm1,rps,rms,spt,smt,           
     &     r,s,t,one,half,fourth,eighth                                         
c                                                                               
      data one,half,fourth,eighth /1.0d0,0.5d0,0.25d0,0.125d0/                  
c                                                                               
      rp1 = one+r                                                               
      rm1 = one-r                                                               
      sp1 = one+s                                                               
      sm1 = one-s                                                               
      tp1 = one+t                                                               
      tm1 = one-t                                                               
      rps = r + s                                                               
      rms = r - s                                                               
      spt = s + t                                                               
      smt = s - t                                                               
c                                                                               
      nr(1)  =  sm1*(rp1 + rps - t)*tp1*eighth                                  
      nr(2)  =  sm1*(rp1 + rps + t)*tm1*eighth                                  
      nr(3)  =  sp1*(rp1 + rms + t)*tm1*eighth                                  
      nr(4)  =  sp1*(rp1 + rms - t)*tp1*eighth                                  
      nr(5)  = -sm1*(rm1 - rms - t)*tp1*eighth                                  
      nr(6)  = -sm1*(rm1 - rms + t)*tm1*eighth                                  
      nr(7)  = -sp1*(rm1 - rps + t)*tm1*eighth                                  
      nr(8)  = -sp1*(rm1 - rps - t)*tp1*eighth                                  
      nr(9)  = -sm1*tm1*tp1*fourth                                              
      nr(10) = -sm1*sp1*tm1*fourth                                              
      nr(11) = -sp1*tm1*tp1*fourth                                              
      nr(12) = -sm1*sp1*tp1*fourth                                              
      nr(13) =  sm1*tm1*tp1*fourth                                              
      nr(14) =  sm1*sp1*tm1*fourth                                              
      nr(15) =  sp1*tm1*tp1*fourth                                              
      nr(16) =  sm1*sp1*tp1*fourth                                              
      nr(17) =  -r*sm1*tp1*half                                                 
      nr(18) =  -r*sm1*tm1*half                                                 
      nr(19) =  -r*sp1*tm1*half                                                 
      nr(20) =  -r*sp1*tp1*half                                                 
c                                                                               
      ns(1)  =  rm1*(rp1 + s + smt)*tp1*eighth                                  
      ns(2)  =  rm1*(rp1 + s + spt)*tm1*eighth                                  
      ns(3)  = -rm1*(rp1 - s - smt)*tm1*eighth                                  
      ns(4)  = -rm1*(rp1 - s - spt)*tp1*eighth                                  
      ns(5)  =  rp1*(rm1 + s + smt)*tp1*eighth                                  
      ns(6)  =  rp1*(rm1 + s + spt)*tm1*eighth                                  
      ns(7)  = -rp1*(rm1 - s - smt)*tm1*eighth                                  
      ns(8)  = -rp1*(rm1 - s - spt)*tp1*eighth                                  
      ns(9)  = -rm1*tm1*tp1*fourth                                              
      ns(10) = -rm1*s*tm1*half                                                  
      ns(11) =  rm1*tm1*tp1*fourth                                              
      ns(12) = -rm1*s*tp1*half                                                  
      ns(13) = -rp1*tm1*tp1*fourth                                              
      ns(14) = -rp1*s*tm1*half                                                  
      ns(15) =  rp1*tm1*tp1*fourth                                              
      ns(16) = -rp1*s*tp1*half                                                  
      ns(17) = -rm1*rp1*tp1*fourth                                              
      ns(18) = -rm1*rp1*tm1*fourth                                              
      ns(19) =  rm1*rp1*tm1*fourth                                              
      ns(20) =  rm1*rp1*tp1*fourth                                              
c                                                                               
      nt(1)  = -rm1*sm1*(rp1 + smt - t)*eighth                                  
      nt(2)  =  rm1*sm1*(rp1 + spt + t)*eighth                                  
      nt(3)  =  rm1*sp1*(rp1 - smt + t)*eighth                                  
      nt(4)  = -rm1*sp1*(rp1 - spt - t)*eighth                                  
      nt(5)  = -rp1*sm1*(rm1 + smt - t)*eighth                                  
      nt(6)  =  rp1*sm1*(rm1 + spt + t)*eighth                                  
      nt(7)  =  rp1*sp1*(rm1 - smt + t)*eighth                                  
      nt(8)  = -rp1*sp1*(rm1 - spt - t)*eighth                                  
      nt(9)  = -rm1*sm1*t*half                                                  
      nt(10) = -rm1*sm1*sp1*fourth                                              
      nt(11) = -rm1*sp1*t*half                                                  
      nt(12) = rm1*sm1*sp1*fourth                                               
      nt(13) = -rp1*sm1*t*half                                                  
      nt(14) = -rp1*sm1*sp1*fourth                                              
      nt(15) = -rp1*sp1*t*half                                                  
      nt(16) = rp1*sm1*sp1*fourth                                               
      nt(17) = rm1*rp1*sm1*fourth                                               
      nt(18) = -rm1*rp1*sm1*fourth                                              
      nt(19) = -rm1*rp1*sp1*fourth                                              
      nt(20) = rm1*rp1*sp1*fourth                                               
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv2                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/18/89                   *          
c     *                                 : 02/08/94                   *          
c     *                                                              *          
c     *     this subroutine computes the derivatives of the          *          
c     *     shape functions at the coordinates of the given          *          
c     *     gauss point for element l3disop                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine deriv2( xi, eta, zeta, nxi, neta, nzeta )                      
      implicit integer (a-z)                                                    
      double precision                                                          
     & xi,eta,zeta,nxi(*),neta(*),nzeta(*),xp,ep,zp,xm,em,zm                    
      double precision                                                          
     &  one, one25                                                              
      data one, one25 / 1.0d0, 0.125d0 /                                        
c                                                                               
c                       set basic parameters                                    
c                                                                               
      xp= one+xi                                                                
      ep= one+eta                                                               
      zp= one+zeta                                                              
      xm= one-xi                                                                
      em= one-eta                                                               
      zm= one-zeta                                                              
c                                                                               
c                       corner nodes                                            
c                                                                               
      nxi(1)=   -em*zp*one25                                                    
      neta(1)=  -xm*zp*one25                                                    
      nzeta(1)=  xm*em*one25                                                    
      nxi(2)=   -em*zm*one25                                                    
      neta(2)=  -xm*zm*one25                                                    
      nzeta(2)= -xm*em*one25                                                    
      nxi(3)=   -ep*zm*one25                                                    
      neta(3)=   xm*zm*one25                                                    
      nzeta(3)= -xm*ep*one25                                                    
      nxi(4)=   -ep*zp*one25                                                    
      neta(4)=   xm*zp*one25                                                    
      nzeta(4)=  xm*ep*one25                                                    
      nxi(5)=    em*zp*one25                                                    
      neta(5)=  -xp*zp*one25                                                    
      nzeta(5)=  xp*em*one25                                                    
      nxi(6)=    em*zm*one25                                                    
      neta(6)=  -xp*zm*one25                                                    
      nzeta(6)= -xp*em*one25                                                    
      nxi(7)=    ep*zm*one25                                                    
      neta(7)=   xp*zm*one25                                                    
      nzeta(7)= -xp*ep*one25                                                    
      nxi(8)=    ep*zp*one25                                                    
      neta(8)=   xp*zp*one25                                                    
      nzeta(8)=  xp*ep*one25                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv3                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 20/07/97                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     gauss point for 12-node transition element               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine deriv3(r,s,t,nr,ns,nt)                                         
      implicit integer (a-z)                                                    
      double precision                                                          
     &     nr(*),ns(*),nt(*),rp1,rm1,sp1,sm1,tp1,tm1,spt,smt,                   
     &     r,s,t,one,half,fourth,eighth                                         
c                                                                               
      data one,half,fourth,eighth /1.0d0,0.5d0,0.25d0,0.125d0/                  
                                                                                
      rp1 = one + r                                                             
      rm1 = one - r                                                             
      sp1 = one + s                                                             
      sm1 = one - s                                                             
      tp1 = one + t                                                             
      tm1 = one - t                                                             
      spt = s + t                                                               
      smt = s - t                                                               
c                                                                               
      nr(1)  =  sm1*(sp1 + smt*t)*eighth                                        
      nr(2)  =  sm1*(sp1 - spt*t)*eighth                                        
      nr(3)  =  sp1*(sm1 + smt*t)*eighth                                        
      nr(4)  =  sp1*(sm1 - spt*t)*eighth                                        
      nr(5)  =  sm1*tp1*eighth                                                  
      nr(6)  =  sm1*tm1*eighth                                                  
      nr(7)  =  sp1*tm1*eighth                                                  
      nr(8)  =  sp1*tp1*eighth                                                  
      nr(9)  = -sm1*tm1*tp1*fourth                                              
      nr(10) = -sm1*sp1*tm1*fourth                                              
      nr(11) = -sp1*tm1*tp1*fourth                                              
      nr(12) = -sm1*sp1*tp1*fourth                                              
c                                                                               
      ns(1)  =  rm1*(s + smt)*tp1*eighth                                        
      ns(2)  =  rm1*(s + spt)*tm1*eighth                                        
      ns(3)  =  rm1*(s + smt)*tm1*eighth                                        
      ns(4)  =  rm1*(s + spt)*tp1*eighth                                        
      ns(5)  = -rp1*tp1*eighth                                                  
      ns(6)  = -rp1*tm1*eighth                                                  
      ns(7)  =  rp1*tm1*eighth                                                  
      ns(8)  =  rp1*tp1*eighth                                                  
      ns(9)  = -rm1*tm1*tp1*fourth                                              
      ns(10) = -rm1*s*tm1*half                                                  
      ns(11) =  rm1*tm1*tp1*fourth                                              
      ns(12) = -rm1*s*tp1*half                                                  
c                                                                               
      nt(1)  = -rm1*sm1*(smt - t)*eighth                                        
      nt(2)  =  rm1*sm1*(spt + t)*eighth                                        
      nt(3)  = -rm1*sp1*(smt - t)*eighth                                        
      nt(4)  =  rm1*sp1*(spt + t)*eighth                                        
      nt(5)  =  rp1*sm1*eighth                                                  
      nt(6)  = -rp1*sm1*eighth                                                  
      nt(7)  = -rp1*sp1*eighth                                                  
      nt(8)  =  rp1*sp1*eighth                                                  
      nt(9)  = -rm1*sm1*t*half                                                  
      nt(10) = -rm1*sm1*sp1*fourth                                              
      nt(11) = -rm1*sp1*t*half                                                  
      nt(12) =  rm1*sm1*sp1*fourth                                              
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv4                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 20/07/97                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     gauss point for 15-node transition element               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine deriv4(r,s,t,nr,ns,nt)                                         
      implicit integer (a-z)                                                    
      double precision                                                          
     &     nr(*),ns(*),nt(*),rp1,rm1,sp1,sm1,tp1,tm1,rpt,rmt,spt,smt,           
     &     r,s,t,one,half,fourth,eighth                                         
c                                                                               
      data one,half,fourth,eighth /1.0d0,0.5d0,0.25d0,0.125d0/                  
c                                                                               
      rp1 = one + r                                                             
      rm1 = one - r                                                             
      sp1 = one + s                                                             
      sm1 = one - s                                                             
      tp1 = one + t                                                             
      tm1 = one - t                                                             
      rpt = r + t                                                               
      rmt = r - t                                                               
      spt = s + t                                                               
      smt = s - t                                                               
c                                                                               
      nr(1)  =  sm1*(rp1 + r + smt)*tp1*eighth                                  
      nr(2)  =  sm1*(rp1 + r + spt)*tm1*eighth                                  
      nr(3)  =  sp1*(sm1 + smt*t)*eighth                                        
      nr(4)  =  sp1*(sm1 - spt*t)*eighth                                        
      nr(5)  =  sm1*(r + rpt)*tp1*eighth                                        
      nr(6)  =  sm1*(r + rmt)*tm1*eighth                                        
      nr(7)  =  sp1*tm1*eighth                                                  
      nr(8)  =  sp1*tp1*eighth                                                  
      nr(9)  = -sm1*tm1*tp1*fourth                                              
      nr(10) = -sm1*sp1*tm1*fourth                                              
      nr(11) = -sp1*tm1*tp1*fourth                                              
      nr(12) = -sm1*sp1*tp1*fourth                                              
      nr(13) =  sm1*tm1*tp1*fourth                                              
      nr(14) = -r*sm1*tp1*half                                                  
      nr(15) = -r*sm1*tm1*half                                                  
c                                                                               
      ns(1)  =  rm1*(rp1 + s + smt)*tp1*eighth                                  
      ns(2)  =  rm1*(rp1 + s + spt)*tm1*eighth                                  
      ns(3)  =  rm1*(s + smt)*tm1*eighth                                        
      ns(4)  =  rm1*(s + spt)*tp1*eighth                                        
      ns(5)  =  rp1*(rm1 - rpt*t)*eighth                                        
      ns(6)  =  rp1*(rm1 + rmt*t)*eighth                                        
      ns(7)  =  rp1*tm1*eighth                                                  
      ns(8)  =  rp1*tp1*eighth                                                  
      ns(9)  = -rm1*tm1*tp1*fourth                                              
      ns(10) = -rm1*s*tm1*half                                                  
      ns(11) =  rm1*tm1*tp1*fourth                                              
      ns(12) = -rm1*s*tp1*half                                                  
      ns(13) = -rp1*tm1*tp1*fourth                                              
      ns(14) = -rm1*rp1*tp1*fourth                                              
      ns(15) = -rm1*rp1*tm1*fourth                                              
c                                                                               
      nt(1)  = -rm1*sm1*(rp1 + smt - t)*eighth                                  
      nt(2)  =  rm1*sm1*(rp1 + spt + t)*eighth                                  
      nt(3)  = -rm1*sp1*(smt - t)*eighth                                        
      nt(4)  =  rm1*sp1*(spt + t)*eighth                                        
      nt(5)  =  rp1*sm1*(rpt + t)*eighth                                        
      nt(6)  = -rp1*sm1*(rmt - t)*eighth                                        
      nt(7)  = -rp1*sp1*eighth                                                  
      nt(8)  =  rp1*sp1*eighth                                                  
      nt(9)  = -rm1*sm1*t*half                                                  
      nt(10) = -rm1*sm1*sp1*fourth                                              
      nt(11) = -rm1*sp1*t*half                                                  
      nt(12) =  rm1*sm1*sp1*fourth                                              
      nt(13) = -rp1*sm1*t*half                                                  
      nt(14) =  rm1*rp1*sm1*fourth                                              
      nt(15) = -rm1*rp1*sm1*fourth                                              
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv5                       *          
c     *                                                              *          
c     *                       written by : cm                        *          
c     *                                                              *          
c     *                   last modified : 08/22/97                   *          
c     *                                                              *          
c     *     this subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     gauss point for 9-node transition  element               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine deriv5(r,s,t,nr,ns,nt)                                         
      implicit integer (a-z)                                                    
      double precision                                                          
     &     nr(*),ns(*),nt(*),rp1,rm1,sp1,sm1,tp1,tm1,                           
     &     r,s,t,one,half,fourth,eighth                                         
c                                                                               
      data one,half,fourth,eighth /1.0d0,0.5d0,0.25d0,0.125d0/                  
c                                                                               
      rp1 = one+r                                                               
      rm1 = one-r                                                               
      sp1 = one+s                                                               
      sm1 = one-s                                                               
      tp1 = one+t                                                               
      tm1 = one-t                                                               
c                                                                               
      nr(1) = -sm1*tp1*t*eighth                                                 
      nr(2) =  sm1*tm1*t*eighth                                                 
      nr(3) = -sp1*tm1*eighth                                                   
      nr(4) = -sp1*tp1*eighth                                                   
      nr(5) =  sm1*tp1*eighth                                                   
      nr(6) =  sm1*tm1*eighth                                                   
      nr(7) =  sp1*tm1*eighth                                                   
      nr(8) =  sp1*tp1*eighth                                                   
      nr(9) = -sm1*tm1*tp1*fourth                                               
c                                                                               
      ns(1) = -rm1*tp1*t*eighth                                                 
      ns(2) =  rm1*tm1*t*eighth                                                 
      ns(3) =  rm1*tm1*eighth                                                   
      ns(4) =  rm1*tp1*eighth                                                   
      ns(5) = -rp1*tp1*eighth                                                   
      ns(6) = -rp1*tm1*eighth                                                   
      ns(7) =  rp1*tm1*eighth                                                   
      ns(8) =  rp1*tp1*eighth                                                   
      ns(9) = -rm1*tm1*tp1*fourth                                               
c                                                                               
      nt(1) =  rm1*sm1*( tp1 + t)*eighth                                        
      nt(2) =  rm1*sm1*(-tm1 + t)*eighth                                        
      nt(3) = -rm1*sp1*eighth                                                   
      nt(4) =  rm1*sp1*eighth                                                   
      nt(5) =  rp1*sm1*eighth                                                   
      nt(6) = -rp1*sm1*eighth                                                   
      nt(7) = -rp1*sp1*eighth                                                   
      nt(8) =  rp1*sp1*eighth                                                   
      nt(9) = -rm1*sm1*t*half                                                   
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv6                       *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/19/98                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     integration point for the 10-node tetrahedron.           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c       Evaluate the 10 shape function derivatives at a given                   
c       integration point for the 10-node tetrahedron volume element.           
c       The 10-node tetrahedron is element type 6 in Warp3D, use the            
c       subroutine name "deriv6" for this element type.                         
c       See subroutine "shape6" in "shapef" for the shape functions.            
c                                                                               
c       Note that the shape function derivatives are with respect to the        
c       four element local tetrahedron coordinates s1,s2,s3,s4.                 
c                                                                               
c       Note that the nodes are numbered at the corners first (nodes 1,2,3      
c       on the face 1 and then node 4 at the peak of the tetrahedron)           
c       and then at the edge nodes on face 1 in a right hand rule sense so      
c       that the surface normal in into the element (nodes 5,6,7), then the     
c       edges from face 1 up to node 4 (edge nodes 8,9,10).                     
c       This should be the usual numbering order and is the same as used in     
c       abaqus.                                                                 
c                                                                               
c       Variables:                                                              
c                                                                               
c       s1,s2,s3,s4 = tetrahedron natural coordinates, vary from 0 to 1,        
c                     These coordinates will be at a given Gauss                
c                     integration point.                                        
c        Note that s1 is not passed in but computed locally.                    
c       Properties of coordinates:  s1 + s2 + s3 + s4 = 1                       
c       choose s1 as dependent, then s1 = 1 - s2 - s3 - s4                      
c                               0 <= s1,s2,s3,s4 <= 1                           
c       Note that the Gauss point coordinates s2,s3,s4 must be consistent       
c       with the coordinates given from the "gauss6" subroutine.                
c                                                                               
c       qs2(numnode) = shape function derivatives with respect to triangle      
c                      coordinate s2, evaluated at the given Gauss point        
c       qs3(numnode) = shape function derivatives with respect to triangle      
c                      coordinate s3, evaluated at the given Gauss point        
c       qs4(numnode) = shape function derivatives with respect to triangle      
c                      coordinate s4, evaluated at the given Gauss point        
c       (where numnode=10 for the 10-node tetrahedron element)                  
c                                                                               
c       Note that since the 4 tetrahedron coordinates are not independent,      
c       s1 = 1 - s2 - s3 - s4 has been substituted into the shape functions     
c       when taking derivatives.  s1 appears in the shape functions and         
c       derivatives for simplicity of evaluating the values given a Gauss       
c       integration point coordinate set.                                       
c                                                                               
      subroutine deriv6( s2, s3, s4, qs2, qs3, qs4 )                            
      implicit none                                                             
      double precision                                                          
     &         s2, s3, s4, qs2(*), qs3(*), qs4(*)                               
c                                                                               
c           s1 = fourth tetrahedron natural coordinate, used here               
c                to keep the simplicity of the functions below                  
c                                                                               
      double precision                                                          
     &     zero, one, four, s1                                                  
c                                                                               
      zero = 0.0d0; one  = 1.0d0; four = 4.0d0                                  
      s1 = one - s2 - s3 - s4                                                   
c                                                                               
c           evaluate the shape function derivatives with respect to the         
c           tetrahedron coordinates at the current Gauss point.                 
c           (s1 = 1 - s2 - s3 - s4 was substituted when taking derivatives,     
c           but appears below for simplicity of evaluating values)              
c                                                                               
c           dq_i/ds2, derivatives with respect to s2, the first independent     
c           local coordinate (substituted for s1)                               
c           ex: qs2(1) = dq1/s2, qs2(2) = dq2/s2, qs2(3) = dq3/s2, etc.         
c                                                                               
      qs2(1)  = one - four*s1                                                   
      qs2(2)  = four*s2 - one                                                   
      qs2(3)  = zero                                                            
      qs2(4)  = zero                                                            
      qs2(5)  = four*(s1 - s2)                                                  
      qs2(6)  = four*s3                                                         
      qs2(7)  = -four*s3                                                        
      qs2(8)  = -four*s4                                                        
      qs2(9)  = four*s4                                                         
      qs2(10) = zero                                                            
c                                                                               
c           dq_i/ds3, derivatives with respect to s3, the second independent    
c           local coordinate (substituted for s1)                               
c           ex: qs3(1) = dq1/s3, qs3(2) = dq2/s3, qs3(3) = dq3/s3, etc.         
c                                                                               
      qs3(1)  = one - four*s1                                                   
      qs3(2)  = zero                                                            
      qs3(3)  = four*s3 - one                                                   
      qs3(4)  = zero                                                            
      qs3(5)  = -four*s2                                                        
      qs3(6)  = four*s2                                                         
      qs3(7)  = four*(s1 - s3)                                                  
      qs3(8)  = -four*s4                                                        
      qs3(9)  = zero                                                            
      qs3(10) = four*s4                                                         
c                                                                               
c           dq_i/ds4, derivatives with respect to s4, the third independent     
c           local coordinate (substituted for s1)                               
c           ex: qs3(1) = dq1/s4, qs3(2) = dq2/s4, qs3(3) = dq3/s4, etc.         
c                                                                               
      qs4(1)  = one - four*s1                                                   
      qs4(2)  = zero                                                            
      qs4(3)  = zero                                                            
      qs4(4)  = four*s4 - one                                                   
      qs4(5)  = -four*s2                                                        
      qs4(6)  = zero                                                            
      qs4(7)  = -four*s3                                                        
      qs4(8)  = four*(s1 - s4)                                                  
      qs4(9)  = four*s2                                                         
      qs4(10) = four*s3                                                         
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv9                       *          
c     *                      currently identical to deriv10          *          
c     *                       written by : mcw                       *          
c     *                       (from gvt)                             *          
c     *                   last modified : 08/30/03                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     Gauss integration point for the 8-node quadrilateral.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         evaluate the 8 derivatives of shape functions at a given Gauss        
c         integration point for the 8-node quadrilateral element.               
c         In Warp3D the 8-node quadrilateral element is                         
c         element type 9, use "deriv9" for this element type.                   
c         See "shape9" in "shapef" for the shape functions.                     
c                                                                               
c         Note that the shape function derivatives are with respect to the      
c         two element local quadrilateral natural coordinates s1,s2 (Xi,Eta).   
c                                                                               
c         Note from the picture below that the nodes are numbered at the        
c         corners and then at the edges in a counter-clockwise order.           
c                                                                               
c         8-node quadrilateral element (shown as a rectangle,                   
c         but the geometry can of course be more general):                      
c                                                                               
c          4            7            3                                          
c           o - - - - - o - - - - - o                                           
c           |                       |                                           
c           |           ^s2         |                                           
c           |           |           |                                           
c         8 o           +-->s1      o 6                                         
c           |                       |                                           
c           |                       |                                           
c           |                       |                                           
c           o - - - - - o - - - - - o                                           
c          1            5             2                                         
c                                                                               
c                                                                               
c       Variables:                                                              
c                                                                               
c        s1,s2, = quadrilateral natural coordinates, vary from -1 to 1,         
c                 These coordinates will be at a given Gauss integration        
c                 point. (or think of s1=Xi, s2=Eta)                            
c        qs1(numnode) = shape function derivatives with respect to the          
c                       quadrilateral coordinate s1, evaluated at               
c                       the given Gauss point                                   
c        qs2(numnode) = shape function derivatives with respect to the          
c                       quadrilateral coordinate s2, evaluated at the           
c                       given Gauss point                                       
c        (where numnode = 8 for the 8-node quadrilateral element)               
c                                                                               
      subroutine deriv9( s1, s2, qs1, qs2 )                                     
      implicit none                                                             
      double precision                                                          
     &         s1, s2, qs1(*), qs2(*)                                           
c                                                                               
      double precision                                                          
     &         one, two, four                                                   
c                                                                               
      one  = 1.0d0; two  = 2.0d0; four = 4.0d0                                  
c                                                                               
c         dq_i/ds1, derivatives with respect to s1, the first local             
c         coordinate (or dq_i/dXi with Xi as the first local coordinate)        
c         ex: qs1(1) = dq1/s1, qs1(2) = dq2/s1, qs1(3) = dq3/s1, etc.           
c                                                                               
      qs1(1) = (two*s1 + s2)*(one - s2)/four                                    
      qs1(2) = (two*s1 - s2)*(one - s2)/four                                    
      qs1(3) = (two*s1 + s2)*(one + s2)/four                                    
      qs1(4) = (two*s1 - s2)*(one + s2)/four                                    
      qs1(5) = -s1*(one - s2)                                                   
      qs1(6) = (one - s2*s2)/two                                                
      qs1(7) = -s1*(one + s2)                                                   
      qs1(8) = -(one - s2*s2)/two                                               
c                                                                               
c         dq_i/ds2, derivatives with respect to s2, the second local            
c         coordinate (or dq_i/dEta with Eta as the second local coordinate)     
c         ex: qs2(1) = dq1/s2, qs2(2) = dq2/s2, qs2(3) = dq3/s2, etc.           
c                                                                               
      qs2(1) = (one - s1)*(two*s2 + s1)/four                                    
      qs2(2) = (one + s1)*(two*s2 - s1)/four                                    
      qs2(3) = (one + s1)*(two*s2 + s1)/four                                    
      qs2(4) = (one - s1)*(two*s2 - s1)/four                                    
      qs2(5) = -(one - s1*s1)/two                                               
      qs2(6) = -s2*(one + s1)                                                   
      qs2(7) = (one - s1*s1)/two                                                
      qs2(8) = -s2*(one - s1)                                                   
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv10                      *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/20/98                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     Gauss integration point for the 8-node axisymmetric      *          
c     *     quadrilateral.                                          *           
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         evaluate the 8 derivatives of shape functions at a given Gauss        
c         integration point for the 8-node axisymmetric quadrilateral           
c         element (no center node).                                             
c         In Warp3D the 8-node axisymmetric quadrilateral element is            
c         element type 10, use "deriv10" for this element type.                 
c         See "shape10" in "shapef" for the shape functions.                    
c                                                                               
c         Note that the shape function derivatives are with respect to the      
c         two element local quadrilateral natural coordinates s1,s2 (Xi,Eta).   
c                                                                               
c         Note from the picture below that the nodes are numbered at the        
c         corners and then at the edges in a counter-clockwise order.           
c                                                                               
c         8-node axisymmetric quadrilateral element (shown as a rectangle,      
c         but the geometry can of course be more general):                      
c                                                                               
c          4            7            3                                          
c           o - - - - - o - - - - - o                                           
c           |                       |                                           
c           |           ^s2         |                                           
c           |           |           |                                           
c         8 o           +-->s1      o 6                                         
c           |                       |                                           
c           |                       |                                           
c           |                       |                                           
c           o - - - - - o - - - - - o                                           
c          1            5             2                                         
c                                                                               
c                                                                               
c       Variables:                                                              
c                                                                               
c        s1,s2, = quadrilateral natural coordinates, vary from -1 to 1,         
c                 These coordinates will be at a given Gauss integration        
c                 point. (or think of s1=Xi, s2=Eta)                            
c        qs1(numnode) = shape function derivatives with respect to the          
c                       quadrilateral coordinate s1, evaluated at               
c                       the given Gauss point                                   
c        qs2(numnode) = shape function derivatives with respect to the          
c                       quadrilateral coordinate s2, evaluated at the           
c                       given Gauss point                                       
c        (where numnode = 8 for the 8-node axisymmetric quadrilateral           
c                         element)                                              
c                                                                               
      subroutine deriv10( s1, s2, qs1, qs2 )                                    
      implicit none                                                             
      double precision                                                          
     &         s1, s2, qs1(*), qs2(*)                                           
c                                                                               
      double precision                                                          
     &         one, two, four                                                   
c                                                                               
      one  = 1.0d0; two = 2.0d0; four = 4.0d0                                   
c                                                                               
c         dq_i/ds1, derivatives with respect to s1, the first local             
c         coordinate (or dq_i/dXi with Xi as the first local coordinate)        
c         ex: qs1(1) = dq1/s1, qs1(2) = dq2/s1, qs1(3) = dq3/s1, etc.           
c                                                                               
      qs1(1) = (two*s1 + s2)*(one - s2)/four                                    
      qs1(2) = (two*s1 - s2)*(one - s2)/four                                    
      qs1(3) = (two*s1 + s2)*(one + s2)/four                                    
      qs1(4) = (two*s1 - s2)*(one + s2)/four                                    
      qs1(5) = -s1*(one - s2)                                                   
      qs1(6) = (one - s2*s2)/two                                                
      qs1(7) = -s1*(one + s2)                                                   
      qs1(8) = -(one - s2*s2)/two                                               
c                                                                               
c         dq_i/ds2, derivatives with respect to s2, the second local            
c         coordinate (or dq_i/dEta with Eta as the second local coordinate)     
c         ex: qs2(1) = dq1/s2, qs2(2) = dq2/s2, qs2(3) = dq3/s2, etc.           
c                                                                               
      qs2(1) = (one - s1)*(two*s2 + s1)/four                                    
      qs2(2) = (one + s1)*(two*s2 - s1)/four                                    
      qs2(3) = (one + s1)*(two*s2 + s1)/four                                    
      qs2(4) = (one - s1)*(two*s2 - s1)/four                                    
      qs2(5) = -(one - s1*s1)/two                                               
      qs2(6) = -s2*(one + s1)                                                   
      qs2(7) = (one - s1*s1)/two                                                
      qs2(8) = -s2*(one - s1)                                                   
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv11                      *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/20/98                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     Gauss integration point for the 6-node axisymmetric      *          
c     *     triangle.                                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c           evaluate the 6 derivatives of shape functions at a given Gauss      
c           integration point for the 6-node axisymmetric triangle element.     
c           In Warp3D the axisymmetric 6-node triangle is element               
c           type 11, use "deriv11" for the Gauss quadrature integration         
c           for this element type.                                              
c           See "shape11" in "shapef" for the shape functions.                  
c                                                                               
c           Note that the shape function derivatives are with respect           
c           to the first two of the three element local triangle                
c           coordinates s1,s2,s3; derivatives with respect to s1 and s2         
c           (s3 is dependent on s1 and s2).                                     
c                                                                               
c           Note from the picture below that the nodes are numbered at the      
c           corners and then at the edges in a counter-clockwise order.         
c                                                                               
c 6-node triangle element:                                                      
c                                                                               
c                  3                                                            
c                 o                                                             
c                /  \                                                           
c               /     \                                                         
c              /        \                                                       
c           6 o           o 5                                                   
c            /              \                                                   
c           /                 \                                                 
c          /                    \                                               
c         o - - - - - o - - - - - o                                             
c        1            4             2                                           
c                                                                               
c                                                                               
c      Variables:                                                               
c                                                                               
c       s1,s2,s3 = triangle natural coordinates, vary from 0 to 1,              
c                 example: s1=0 on the triangle side across from                
c                  corner 1 (through nodes 2, 5, and 3), s1=0.5 along           
c                  line that bisects the other two sides parallel to the        
c                  far side and through the middle the triangle (through        
c                  nodes 4 and 6), and s1=1 at corner 1.                        
c                  These coordinates will be at a given Gauss                   
c                  integration point.                                           
c      Properties of coordinates:  s1 + s2 + s3 = 1                             
c                           use s3 = 1 - s1 - s2 here                           
c                           and 0 <= s1,s2,s3 <= 1                              
c                                                                               
c       Note: pass in all three coordinates even through s3 is dependent        
c             on s1 and s2                                                      
c                                                                               
c       qs1(numnode) = shape function derivatives with respect to triangle      
c                     coordinate s1, evaluated at the given Gauss point         
c       qs2(numnode) = shape function derivatives with respect to triangle      
c                     coordinate s2, evaluated at the given Gauss point         
c       (where numnode = 6 for the 6-node triangle element)                     
c                                                                               
      subroutine deriv11( s1, s2, s3, qs1, qs2 )                                
      implicit none                                                             
      double precision                                                          
     &          s1, s2, s3, qs1(*), qs2(*)                                      
c                                                                               
      double precision                                                          
     &          zero, one, four                                                 
c                                                                               
      zero = 0.0d0; one = 1.0d0;  four = 4.0d0                                  
c                                                                               
c          evaluate the shape function derivatives with respect to the          
c          triangle coordinates at the current Gauss point.                     
c          Note that before differentiating the shape functions,                
c          the dependent coordinate s3 was substituted into the shape           
c          functions.  After the differentiation is complete,                   
c          s3 can be substituted back into the functions for simplicity.        
c                                                                               
c          dq_i/ds1, derivatives with respect to s1, the first                  
c          local coordinate ex: qs1(1) = dq1/s1, qs1(2) =                       
c          dq2/s1, qs1(3) = dq3/s1, etc.                                        
c                                                                               
      qs1(1) = four*s1 - one                                                    
      qs1(2) = zero                                                             
      qs1(3) = one - four*s3                                                    
      qs1(4) = four*s2                                                          
      qs1(5) = -four*s2                                                         
      qs1(6) = four*(s3 - s1)                                                   
c                                                                               
c          dq_i/ds2, derivatives with respect to s2, the second                 
c          local coordinatec ex: qs2(1) = dq1/s2, qs2(2) = dq2/s2,              
c          qs2(3) = dq3/s2, etc.                                                
c                                                                               
      qs2(1) = zero                                                             
      qs2(2) = four*s2 - one                                                    
      qs2(3) = one - four*s3                                                    
      qs2(4) = four*s1                                                          
      qs2(5) = four*(s3 - s2)                                                   
      qs2(6) = -four*s1                                                         
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv12                      *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                                                              *          
c     *                   last modified : 02/9/13 rhd                *          
c     *                                                              *          
c     *     derivatives of shape functions at integration point      *          
c     *     for element inter_8                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine deriv12( xi, eta, nxi, neta )                                  
      implicit integer (a-z)                                                    
      double precision                                                          
     & xi,eta,nxi(*),neta(*),xp,ep,xm,em, zero                                  
      double precision                                                          
     &  one, one4                                                               
      data one, one4, zero / 1.0d0, 0.25d0, 0.0d0 /                             
c                                                                               
c            nodes 1-4 are bottom; 5-8 top surface                              
c                                                                               
      xp = one+xi                                                               
      ep = one+eta                                                              
      xm = one-xi                                                               
      em = one-eta                                                              
c                                                                               
c                       corner nodes                                            
c                                                                               
      nxi(1)  = -em*one4                                                        
      neta(1) = -xm*one4                                                        
      nxi(2)  =  em*one4                                                        
      neta(2) = -xp*one4                                                        
      nxi(3)  =  ep*one4                                                        
      neta(3) =  xp*one4                                                        
      nxi(4)  = -ep*one4                                                        
      neta(4) =  xm*one4                                                        
      nxi(5:8)  =  zero                                                         
      neta(5:8) =  zero                                                         
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv13                      *          
c     *                                                              *          
c     *                       written by : N. Rau                    *          
c     *                                                              *          
c     *                   last modified : 02/09/13 rhd               *          
c     *                                                              *          
c     *     derivatives of shape functions at point for the 4-node   *          
c     *     tetrahedron.                                             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c       Shape function derivatives are with respect to the                      
c       four element local tetrahedron coordinates s1,s2,s3,s4.                 
c                                                                               
c       Integration point coordinates s2,s3,s4 must be consistent               
c       with the coordinates given from the "gauss6" subroutine.                
c                                                                               
c       This element has linear shape functions, so the derivatives             
c       will be constant in all cases.                                          
c                                                                               
c       qs2(numnode) = shape function derivatives with respect to triangle      
c                      coordinate s2, evaluated at the given Gauss point        
c       qs3(numnode) = shape function derivatives with respect to triangle      
c                      coordinate s3, evaluated at the given Gauss point        
c       qs4(numnode) = shape function derivatives with respect to triangle      
c                      coordinate s4, evaluated at the given Gauss point        
c       (where numnode=4 for the 4-node tetrahedron element)                    
c                                                                               
      subroutine deriv13( s2, s3, s4, qs2, qs3, qs4 )                           
      implicit none                                                             
      double precision                                                          
     &         s2, s3, s4, qs2(*), qs3(*), qs4(*)                               
      double precision                                                          
     &     zero, one, minusone, s1                                              
c                                                                               
      zero = 0.0d0; one = 1.0d0; minusone = -1.0d0                              
c                                                                               
c           s1 = fourth tetrahedron natural coordinate, used here               
c                to keep the simplicity of the functions below                  
c                                                                               
      s1 = one - s2 - s3 - s4                                                   
c                                                                               
c           dq_i/ds2, derivatives with respect to s2, the first                 
c           independent local coordinate (substituted for s1)                   
c           ex: qs2(1) = dq1/s2, qs2(2) = dq2/s2, qs2(3) = dq3/s2, etc.         
c                                                                               
      qs2(1)  = minusone                                                        
      qs2(2)  = one                                                             
      qs2(3)  = zero                                                            
      qs2(4)  = zero                                                            
c                                                                               
c           dq_i/ds3, derivatives with respect to s3, the second                
c           independent local coordinate (substituted for s1)                   
c           ex: qs3(1) = dq1/s3, qs3(2) = dq2/s3, qs3(3) = dq3/s3, etc.         
c                                                                               
      qs3(1)  = minusone                                                        
      qs3(2)  = zero                                                            
      qs3(3)  = one                                                             
      qs3(4)  = zero                                                            
c                                                                               
c           dq_i/ds4, derivatives with respect to s4, the third                 
c           independent local coordinate (substituted for s1)                   
c           ex: qs3(1) = dq1/s4, qs3(2) = dq2/s4, qs3(3) = dq3/s4, etc.         
c                                                                               
      qs4(1)  = minusone                                                        
      qs4(2)  = zero                                                            
      qs4(3)  = zero                                                            
      qs4(4)  = one                                                             
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv14                      *          
c     *                                                              *          
c     *                       written by : sushovan                  *          
c     *                                                              *          
c     *                   last modified : 02/09/13 rhd               *          
c     *                                                              *          
c     *                                                              *          
c     *     shape function derivatives at a given parametric point   *          
c     *     for the 6-node triangular interface element trint6.      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         evaluate the derivatives of 6 shape functions at a given              
c         parametric point for the 6-node triangular interface                  
c         element trin6 (elem type 14)                                          
c         See "shape14" in "shapef" for the shape functions.                    
c                                                                               
c         Shape function derivatives are with respect                           
c         to the first two (s1,s2) of the three element local triangle          
c         coordinates s1,s2,s3                                                  
c                                                                               
c         For the 6-node interface element, nodes 1-3 define                    
c         the bottom surface and 4-6 define the top surface.                    
c         Points inside the element are defined only on the top and             
c         bottom surface.                                                       
c                                                                               
c         3-node bottom surface of the 6-node interface element:                
c                                                                               
c        3                                                                      
c         o                                                                     
c         |  \                                                                  
c         |    \                                                                
c         |      \                                                              
c         |        \   corresponding nodes 4-6 on top surface                   
c         |          \                                                          
c         |            \                                                        
c         |              \                                                      
c         o - - - - - - -  o                                                    
c        1                  2                                                   
c                                                                               
c                                                                               
c                                                                               
c      Variables:                                                               
c                                                                               
c       s1,s2,s3 = triangle natural coordinates, vary from 0 to 1,              
c                     example: s1 = 0 on the triangle side 2-3                  
c                              = 0.5 at nodes 4 and 6                           
c                              = 1 at corner 1                                  
c          s1 + s2 + s3 = 1. Use s3 = 1 - s1 - s2 as convenient                 
c                                                                               
c       Only s1 and s2 are input to this subroutine.                            
c                                                                               
c       qs1(num_enodes) = derivatives wrt s1                                    
c       qs2(num_enodes) = derivatives wrt s2                                    
c       (where num_enodes = 6)                                                  
c                                                                               
      subroutine deriv14( s1, s2, qs1, qs2 )                                    
      implicit none                                                             
      double precision                                                          
     &          s1, s2, s3, qs1(*), qs2(*)                                      
c                                                                               
      double precision                                                          
     &          zero, one                                                       
c                                                                               
      zero = 0.0d0; one  = 1.0d0                                                
c                                                                               
c          before differentiating the shape functions,                          
c          the dependent coordinate s3 was substituted into the shape           
c          functions.  After the differentiation is complete,                   
c          s3 can be substituted back into the derivatives for simplicity.      
c                                                                               
      s3 = one - s1 -s2                                                         
      qs1(1) = one                                                              
      qs1(2) = zero                                                             
      qs1(3) =  -one                                                            
      qs1(4) = zero                                                             
      qs1(5) = zero                                                             
      qs1(6) = zero                                                             
c                                                                               
      qs2(1) = zero                                                             
      qs2(2) = one                                                              
      qs2(3) = -one                                                             
      qs2(4) = zero                                                             
      qs2(5) = zero                                                             
      qs2(6) = zero                                                             
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv15                      *          
c     *                                                              *          
c     *                       written by : sushovan                  *          
c     *                                                              *          
c     *                   last modified : 02/9/13 rhd                *          
c     *                                                              *          
c     *     shape function derivatives at a given parametric point   *          
c     *     for the 12-node triangular interface element trint12.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         evaluate the derivatives of 12 shape functions at a given             
c         parametric point for the 12-node triangular interface                 
c         element trint12 (elem type 15)                                        
c         See "shape15" in "shapef" for the shape functions.                    
c                                                                               
c         Shape function derivatives are with respect                           
c         to the first two (s1,s2) of the three element local triangle          
c         coordinates s1,s2,s3                                                  
c                                                                               
c         For the 12-node interface element, nodes 1-6 define                   
c         the bottom surface, 7-12 define the top surface.                      
c         Points inside the element are defined only on the top and             
c         bottom surface.                                                       
c                                                                               
c         6-node bottom surface of the 12-node interface element:               
c                                                                               
c       3                                                                       
c         o                                                                     
c         |  \                                                                  
c         |    \                                                                
c         |      \  5                                                           
c       6 o        o          corresponding nodes 7-12 on top surface           
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
c                     example: s1 = 0 on the triangle side 2-3                  
c                              = 0.5 at nodes 4 and 6                           
c                              = 1 at corner 1                                  
c          s1 + s2 + s3 = 1. Use s3 = 1 - s1 - s2 as convenient                 
c                                                                               
c       Only s1 and s2 are input to this subroutine.                            
c                                                                               
c       qs1(num_enodes) = derivatives wrt s1                                    
c       qs2(num_enodes) = derivatives wrt s2                                    
c       (where num_enodes = 12)                                                 
c                                                                               
      subroutine deriv15( s1, s2, qs1, qs2 )                                    
      implicit none                                                             
      double precision                                                          
     &          s1, s2, qs1(*), qs2(*)                                          
c                                                                               
      double precision                                                          
     &          zero, one, four, s3                                             
      data zero, one, four                                                      
     &  / 0.0d0, 1.0d0, 4.0d0 /                                                 
c                                                                               
c          before differentiating the shape functions,                          
c          the dependent coordinate s3 was substituted into the shape           
c          functions.  After the differentiation is complete,                   
c          s3 can be substituted back into the derivatives for simplicity.      
c                                                                               
c          qs1(7-12) and qs2(7-12)  = 0 for subsequent computations             
c                                                                               
      s3 = one - s1 - s2                                                        
      qs1(1) = four*s1 - one                                                    
      qs1(2) = zero                                                             
      qs1(3) = one - four*s3                                                    
      qs1(4) = four*s2                                                          
      qs1(5) = -four*s2                                                         
      qs1(6) = four*(s3 - s1)                                                   
      qs1(7:12) = zero                                                          
c                                                                               
      qs2(1) = zero                                                             
      qs2(2) = four*s2 - one                                                    
      qs2(3) = one - four*s3                                                    
      qs2(4) = four*s1                                                          
      qs2(5) = four*(s3 - s2)                                                   
      qs2(6) = -four*s1                                                         
      qs2(7:12) = zero                                                          
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine deriv16                      *          
c     *                                                              *          
c     *                       written by : mcw                       *          
c     *                       (from gvt)                             *          
c     *                   last modified : 08/30/03                   *          
c     *                                                              *          
c     *     This subroutine computes the values of the derivatives   *          
c     *     of the shape functions at the coordinates of the given   *          
c     *     Gauss integration point for the 4-node quadrilateral.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         evaluate the 4 derivatives of shape functions at a given Gauss        
c         integration point for the 4-node quadrilateral element.               
c         In Warp3D the 4-node quadrilateral element is                         
c         element type 16?, use "deriv16" for this element type.                
c         See "shape16" in "shapef" for the shape functions.                    
c                                                                               
c         Note that the shape function derivatives are with respect to the      
c         two element local quadrilateral natural coordinates s1,s2 (Xi,Eta).   
c                                                                               
c         Note from the picture below that the nodes are numbered at the        
c         corners and then at the edges in a counter-clockwise order.           
c                                                                               
c         4-node quadrilateral element (shown as a rectangle,                   
c         but the geometry can of course be more general):                      
c                                                                               
c          4                         3                                          
c           o - - - - - - - - - - - o                                           
c           |                       |                                           
c           |           ^s2         |                                           
c           |           |           |                                           
c           |           +-->s1      |                                           
c           |                       |                                           
c           |                       |                                           
c           |                       |                                           
c           o - - - - - - - - - - - o                                           
c          1                         2                                          
c                                                                               
c                                                                               
c       Variables:                                                              
c                                                                               
c        s1,s2, = quadrilateral natural coordinates, vary from -1 to 1,         
c                 These coordinates will be at a given Gauss integration        
c                 point. (or think of s1=Xi, s2=Eta)                            
c        qs1(numnode) = shape function derivatives with respect to the          
c                       quadrilateral coordinate s1, evaluated at               
c                       the given Gauss point                                   
c        qs2(numnode) = shape function derivatives with respect to the          
c                       quadrilateral coordinate s2, evaluated at the           
c                       given Gauss point                                       
c        (where numnode = 4 for the 4-node quadrilateral element)               
c                                                                               
      subroutine deriv16( s1, s2, qs1, qs2 )                                    
      implicit none                                                             
      double precision                                                          
     &         s1, s2, qs1(*), qs2(*)                                           
c                                                                               
      double precision                                                          
     &         one, two, four                                                   
c                                                                               
      one  = 1.0D0                                                              
      two  = 2.0D0                                                              
      four = 4.0D0                                                              
c                                                                               
c         dq_i/ds1, derivatives with respect to s1, the first local             
c         coordinate (or dq_i/dXi with Xi as the first local coordinate)        
c         ex: qs1(1) = dq1/s1, qs1(2) = dq2/s1, qs1(3) = dq3/s1, etc.           
c                                                                               
      qs1(1) = -(one - s2)/four                                                 
      qs1(2) =  (one - s2)/four                                                 
      qs1(3) =  (one + s2)/four                                                 
      qs1(4) = -(one + s2)/four                                                 
c                                                                               
c         dq_i/ds2, derivatives with respect to s2, the second local            
c         coordinate (or dq_i/dEta with Eta as the second local coordinate)     
c         ex: qs2(1) = dq1/s2, qs2(2) = dq2/s2, qs2(3) = dq3/s2, etc.           
c                                                                               
      qs2(1) = -(one - s1)/four                                                 
      qs2(2) = -(one + s1)/four                                                 
      qs2(3) =  (one + s1)/four                                                 
      qs2(4) =  (one - s1)/four                                                 
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
