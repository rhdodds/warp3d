c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulgf                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/22/97                   *          
c     *                                                              *          
c     *     given the element type, order of integration and the     *          
c     *     integration point number, return the isoparametric       *          
c     *     coordinates for the point and the weight value for       *          
c     *     integration                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oulgf( etype, xi, eta, zeta, sf, iorder )                      
      implicit integer (a-z)                                                    
      double precision                                                          
     &     xi, eta, zeta, sf(*)                                                 
c                                                                               
      if ( etype .le. 0 .or. etype .gt. 5 ) then                                
         write(*,9000) etype                                                    
         call die_abort                                                         
         stop                                                                   
      end if                                                                    
c                                                                               
      go to ( 100, 200, 300, 400, 500 ), etype                                  
c                                                                               
c                    20-node                                                    
c                                                                               
 100  continue                                                                  
      call  oulgr1( xi, eta, zeta, sf, iorder )                                 
      return                                                                    
c                                                                               
c                     8-node                                                    
c                                                                               
 200  continue                                                                  
      call  oulgr2( xi, eta, zeta, sf, iorder )                                 
      return                                                                    
c                                                                               
c                    12-node                                                    
c                                                                               
 300  continue                                                                  
      call  oulgr1( xi, eta, zeta, sf, iorder )                                 
      return                                                                    
c                                                                               
c                    15-node                                                    
c                                                                               
 400  continue                                                                  
      call  oulgr1( xi, eta, zeta, sf, iorder )                                 
      return                                                                    
c                                                                               
c                    9-node                                                     
c                                                                               
 500  continue                                                                  
      call  oulgr1( xi, eta, zeta, sf, iorder )                                 
      return                                                                    
c                                                                               
 9000 format('> FATAL ERROR: oulgf, etype: ',i10,//,                            
     &       '               job aborted' )                                     
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulgr1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 05/25/90                   *          
c     *                   last modified : 07/22/97                   *          
c     *                                                              *          
c     *     this subroutine computes the lagrange polynomials        *          
c     *     at each gauss point for a single point in a 3dqisop      *          
c     *     element given the parametric coordinates of that point.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine oulgr1( r, s, t, lg, ord )                                     
      implicit integer (a-z)                                                    
      double precision                                                          
     &     r,s,t,lg(*),p,q,p3,p6,pq4,p2q2,rmp,rpp,smp,spp,tmp,tpp,              
     &     rmq,rpq,smq,spq,tmq,tpq,                                             
     &     p774, eight, four, two, zero, p577, p758                             
      data p774, eight, four, two, zero, p577, p758                             
     &     /0.774596669241483, 8.0, 4.0, 2.0, 0.0,                              
     &     0.577350269189626, 0.758786911/                                      
c                                                                               
c                                                                               
c                       branch on the order of quadrature. possible             
c                       integration schemes are :                               
c                                                                               
c                          order = 1 --  27 point rule ( 3x3x3 )                
c                          order = 2 --   9 point rule                          
c                          order = 3 --  not used                               
c                          order = 4 --  not used                               
c                          order = 5 --  not used                               
c                          order = 6 --  not used                               
c                          order = 7 --  not used                               
c                          order = 8 --   8 point rule ( 2x2x2 )                
c                          order = 9 --  14 point rule                          
c                                                                               
c             extrapolations to nodes based on 9 and 14 pt                      
c             rules is not supported.                                           
c                                                                               
c                                                                               
      if( ord .eq. 1 ) then                                                     
c                                                                               
c                       27 point rule                                           
c                                                                               
         p= p774                                                                
c                                                                               
         rmp= r-p                                                               
         rpp= r+p                                                               
         smp= s-p                                                               
         spp= s+p                                                               
         tmp= t-p                                                               
         tpp= t+p                                                               
         p6= p*p*p*p*p*p                                                        
c                                                                               
         lg(1) = (r*s*t*rmp*smp*tmp)/(eight*p6)                                 
         lg(2) = -(r*t*rmp*smp*spp*tmp)/(four*p6)                               
         lg(3) = (r*s*t*rmp*spp*tmp)/(eight*p6)                                 
         lg(4) = -(s*t*rmp*rpp*smp*tmp)/(four*p6)                               
         lg(5) = (t*rmp*rpp*smp*spp*tmp)/(two*p6)                               
         lg(6) = -(s*t*rmp*rpp*spp*tmp)/(four*p6)                               
         lg(7) = (r*s*t*rpp*smp*tmp)/(eight*p6)                                 
         lg(8) = -(r*t*rpp*smp*spp*tmp)/(four*p6)                               
         lg(9) = (r*s*t*rpp*spp*tmp)/(eight*p6)                                 
         lg(10)= -(r*s*rmp*smp*tmp*tpp)/(four*p6)                               
         lg(11)= (r*rmp*smp*spp*tmp*tpp)/(two*p6)                               
         lg(12)= -(r*s*rmp*spp*tmp*tpp)/(four*p6)                               
         lg(13)= (s*rmp*rpp*smp*tmp*tpp)/(two*p6)                               
         lg(14)= -(rmp*rpp*smp*spp*tmp*tpp)/p6                                  
         lg(15)= (s*rmp*rpp*spp*tmp*tpp)/(two*p6)                               
         lg(16)= -(r*s*rpp*smp*tmp*tpp)/(four*p6)                               
         lg(17)= (r*rpp*smp*spp*tmp*tpp)/(two*p6)                               
         lg(18)= -(r*s*rpp*spp*tmp*tpp)/(four*p6)                               
         lg(19)= (r*s*t*rmp*smp*tpp)/(eight*p6)                                 
         lg(20)= -(r*t*rmp*smp*spp*tpp)/(four*p6)                               
         lg(21)= (r*s*t*rmp*spp*tpp)/(eight*p6)                                 
         lg(22)= -(s*t*rmp*rpp*smp*tpp)/(four*p6)                               
         lg(23)= (t*rmp*rpp*smp*spp*tpp)/(two*p6)                               
         lg(24)= -(s*t*rmp*rpp*spp*tpp)/(four*p6)                               
         lg(25)= (r*s*t*rpp*smp*tpp)/(eight*p6)                                 
         lg(26)= -(r*t*rpp*smp*spp*tpp)/(four*p6)                               
         lg(27)= (r*s*t*rpp*spp*tpp)/(eight*p6)                                 
c                                                                               
      else if( ord .eq. 8 ) then                                                
c                                                                               
c                       8 point rule                                            
c                                                                               
         p= p577                                                                
c                                                                               
         rmp= r-p                                                               
         rpp= r+p                                                               
         smp= s-p                                                               
         spp= s+p                                                               
         tmp= t-p                                                               
         tpp= t+p                                                               
         p3= p*p*p                                                              
c                                                                               
         lg(1) = -(rmp*smp*tmp)/(eight*p3)                                      
         lg(2) = (rmp*spp*tmp)/(eight*p3)                                       
         lg(3) = (rpp*smp*tmp)/(eight*p3)                                       
         lg(4) = -(rpp*spp*tmp)/(eight*p3)                                      
         lg(5) = (rmp*smp*tpp)/(eight*p3)                                       
         lg(6) = -(rmp*spp*tpp)/(eight*p3)                                      
         lg(7) = -(rpp*smp*tpp)/(eight*p3)                                      
         lg(8) = (rpp*spp*tpp)/(eight*p3)                                       
c                                                                               
      else if( ord .eq. 9 ) then                                                
c                                                                               
c                       14 point rule.   not verified.                          
c                                                                               
         p= p758                                                                
c                                                                               
         rmp= r-p                                                               
         rpp= r+p                                                               
         smp= s-p                                                               
         spp= s+p                                                               
         tmp= t-p                                                               
         tpp= t+p                                                               
         p3= p*p*p                                                              
c                                                                               
         lg(1)  =  zero                                                         
         lg(2)  =  zero                                                         
         lg(3)  =  zero                                                         
         lg(4)  =  zero                                                         
         lg(5)  =  zero                                                         
         lg(6)  =  zero                                                         
         lg(7)  =  (rmp*smp*tpp)/(eight*p3)                                     
         lg(8)  = -(rmp*smp*tmp)/(eight*p3)                                     
         lg(9)  =  (rmp*spp*tmp)/(eight*p3)                                     
         lg(10) = -(rmp*spp*tpp)/(eight*p3)                                     
         lg(11) = -(rpp*smp*tpp)/(eight*p3)                                     
         lg(12) =  (rpp*smp*tmp)/(eight*p3)                                     
         lg(13) = -(rpp*spp*tmp)/(eight*p3)                                     
         lg(14) =  (rpp*spp*tpp)/(eight*p3)                                     
c                                                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulgr2                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/18/89                   *          
c     *                                                              *          
c     *     this subroutine computes the lagrange polynomials        *          
c     *     at each gauss point for a single point in a l3disop      *          
c     *     element given the parametric coordinates of that point.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine oulgr2( r, s, t, lg, ord )                                     
      implicit integer (a-z)                                                    
      double precision                                                          
     &     r,s,t,lg(*),p,p3,rmp,rpp,smp,spp,tmp,tpp,                            
     &     p557, eight, p795, two                                               
      data p557, eight, p795, two                                               
     &     /0.577350269189626, 8.0, 0.795822426, 2.0/                           
c                                                                               
c                                                                               
c                       branch on the order of quadrature. possible             
c                       integration schemes are :                               
c                                                                               
c                          ord = 1 --  8 point rule ( 2x2x2 )                   
c                          ord = 2 --  6 point rule                             
c                                                                               
c                                                                               
c                                                                               
      if(ord.eq.1) then                                                         
c                                                                               
c                       8 point rule                                            
c                                                                               
         p= p557                                                                
c                                                                               
         rmp= r-p                                                               
         rpp= r+p                                                               
         smp= s-p                                                               
         spp= s+p                                                               
         tmp= t-p                                                               
         tpp= t+p                                                               
         p3= p*p*p                                                              
c                                                                               
         lg(1) = -(rmp*smp*tmp)/(eight*p3)                                      
         lg(2) = (rmp*spp*tmp)/(eight*p3)                                       
         lg(3) = (rpp*smp*tmp)/(eight*p3)                                       
         lg(4) = -(rpp*spp*tmp)/(eight*p3)                                      
         lg(5) = (rmp*smp*tpp)/(eight*p3)                                       
         lg(6) = -(rmp*spp*tpp)/(eight*p3)                                      
         lg(7) = -(rpp*smp*tpp)/(eight*p3)                                      
         lg(8) = (rpp*spp*tpp)/(eight*p3)                                       
c                                                                               
      else if(ord.eq.2) then                                                    
c                                                                               
c                       6 point rule.                                           
c                                                                               
         p= p795                                                                
c                                                                               
         rmp= r-p                                                               
         rpp= r+p                                                               
         smp= s-p                                                               
         spp= s+p                                                               
         tmp= t-p                                                               
         tpp= t+p                                                               
         p2= p*p                                                                
c                                                                               
         lg(1) = -rmp/(two*p2)                                                  
         lg(2) =  rpp/(two*p2)                                                  
         lg(3) = -smp/(two*p2)                                                  
         lg(4) =  spp/(two*p2)                                                  
         lg(5) = -tmp/(two*p2)                                                  
         lg(6) =  tpp/(two*p2)                                                  
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulgr3                       *          
c     *                                                              *          
c     *                       written by : mcw                       *          
c     *                                                              *          
c     *                     last modified: 8/21/03 by mcw            *          
c     *                                                              *          
c     *     this subroutine computes the lagrange polynomials        *          
c     *     at each gauss point for a single point in a 2D quad      *          
c     *     element given the parametric coordinates of that point.  *          
c     *     this is for extrapolating gauss-point values to nodes.   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine oulgr3( r, s, lg, ord )                                        
      implicit integer (a-z)                                                    
      double precision                                                          
     &     r, s, lg(*), p, p2, p4, ppr, pmr, pps, pms, p557, p774,              
     &     one, four, eight                                                     
      data p557, p774, one, two, four, eight                                    
     &     / 0.577350269189626, 0.774596669241483, 1.0, 2.0, 4.0, 8.0 /         
c                                                                               
c                       branch on the order of quadrature. possible             
c                       integration schemes are :                               
c                                                                               
c                       ord = 2 --  4 point rule ( 2x2 )                        
c                       ord = 3 --  9 point rule ( 3x3 )                        
c                                                                               
      if(ord.eq.2) then                                                         
c                                                                               
c                       4 point rule (2x2)                                      
c                       these are the same as for co-planar gauss points        
c                       in an 8-point rule for brick elements.                  
c                                                                               
         p   = p557                                                             
c                                                                               
         ppr = p+r                                                              
         pmr = p-r                                                              
         pps = p+s                                                              
         pms = p-s                                                              
         p2  = p*p                                                              
c                                                                               
         lg(1) = (pmr*pms)/(four*p2)                                            
         lg(2) = (pmr*pps)/(four*p2)                                            
         lg(3) = (ppr*pms)/(four*p2)                                            
         lg(4) = (ppr*pps)/(four*p2)                                            
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
      if(ord.eq.3) then                                                         
c                                                                               
c                       9 point rule (3x3)                                      
c                                                                               
         p   = p774                                                             
c                                                                               
         pmr = p-r                                                              
         ppr = p+r                                                              
         pms = p-s                                                              
         pps = p+s                                                              
         p4  = p*p*p*p                                                          
c                                                                               
         lg(1) =  (r*s*pmr*pms)     / (four*p4)                                 
         lg(2) = -(r*pmr*pms*pps)   / (two*p4)                                  
         lg(3) = -(r*s*pmr*pps)     / (four*p4)                                 
         lg(4) = -(s*ppr*pmr*pms)   / (two*p4)                                  
         lg(5) =  (ppr*pmr*pps*pms) / p4                                        
         lg(6) =  (s*ppr*pmr*pps)   / (two*p4)                                  
         lg(7) = -(r*s*ppr*pms)     / (four*p4)                                 
         lg(8) =  (r*ppr*pps*pms)   / (two*p4)                                  
         lg(9) =  (r*s*ppr*pps)     / (four*p4)                                 
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
