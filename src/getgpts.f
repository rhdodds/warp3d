c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine getgpts                      *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/11/2017 rhd             *          
c     *                                                              *          
c     *     given the element type, order of integration and the     *          
c     *     integration point number, return the parametric          *          
c     *     coordinates for the point and the weight value for       *          
c     *     integration                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getgpts( etype, order, gpn, xi, eta, zeta, weight )            
      implicit none                                                             
      double precision :: xi, eta, zeta, weight                                                
      integer etype, order, gpn                                                 
c                                                                               
      if ( etype .le. 0 .or. etype .gt. 19 ) then                               
         write(*,9000) etype                                                    
         call die_abort                                                         
      end if 
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,                
     &       1100, 1200, 1300, 1400, 1500, 1600, 1700, 
     &       1800, 1900 ), etype                  
c                                                                               
c                20-node brick                                                  
c                                                                               
 100  continue                                                                  
      call quad1( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                8-node brick                                                   
c                                                                               
 200  continue                                                                  
      call quad2( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                12-node brick                                                  
c                                                                               
 300  continue                                                                  
      call quad1( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                15-node brick                                                  
c                                                                               
 400  continue                                                                  
      call quad1( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                9-node brick                                                   
c                                                                               
 500  continue                                                                  
      call quad1( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                10-node tetrahedron, "tet10"                                   
c                                                                               
 600  continue                                                                  
      call quad6( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                15-node wedge, "wedge15" -- not implemented                    
c                                                                               
 700  continue                                                                  
      write(*,9100) etype                                                       
      call die_abort                                                            
c                                                                               
c                6-node planar triangle, "tri6"                                 
c                                                                               
 800  continue                                                                  
      call quad11( order, gpn, xi, eta, zeta, weight )                          
      return                                                                    
c                                                                               
c                8-node planar quadrilateral, "quad8"                           
c                                                                               
 900  continue                                                                  
      call quad10( order, gpn, xi, eta, zeta, weight )                          
      return                                                                    
c                                                                               
c                8-node axisymmetric quadrilateral, "axiquad8"                  
c                                                                               
1000  continue                                                                  
      call quad10( order, gpn, xi, eta, zeta, weight )                          
      return                                                                    
c                                                                               
c                6-node axisymmetric triangle, "axitri6"                        
c                                                                               
1100  continue                                                                  
      call quad11( order, gpn, xi, eta, zeta, weight )                          
      return                                                                    
c                                                                               
c                8-node interface element                                       
c                                                                               
1200  continue                                                                  
      call inter_gp( order, gpn, xi, eta, zeta, weight )                        
      return                                                                    
c                                                                               
c                4-node tetrahedron, "tet4"                                     
c                                                                               
1300  continue                                                                  
      call quad6( order, gpn, xi, eta, zeta, weight )                           
      return                                                                    
c                                                                               
c                6-node triangular interface element, "trint6"                  
c                                                                               
1400  continue                                                                  
      call trint6_gp( order, gpn, xi, eta, zeta, weight )                       
      return                                                                    
c                                                                               
c                12-node triangular interface element, "trint12"                
c                                                                               
1500  continue                                                                  
      call trint12_gp( order, gpn, xi, eta, zeta, weight )                      
      return                                                                    
c
c                4-node quadrilateral - used for various
c                purposes - not a real element
c
1600  continue
      call quad10( order, gpn, xi, eta, zeta, weight )
      return
c
c                3-node line -- - used for various
c                purposes - not a real element
c
1700  continue
      call line3( order, gpn, xi, eta, zeta, weight )
      return                                                                               
c                                                                                                                                                            
c                2-node bar element                                   
c                                                                               
1800  continue 
      xi = 0.5d0
      eta = 0.0d0
      zeta = 0.0d0
      weight = 1.0d0
      return
c                                                                                                                                                            
c                2-node link element
c                                                                               
1900  continue 
      xi = 0.5d0
      eta = 0.0d0
      zeta = 0.0d0
      weight = 1.0d0
      return
c
                                                                     
 9000 format('> FATAL ERROR: getgpts, etype: ',i10,//,                          
     &       '               job aborted' )                                     
 9100 format('> FATAL ERROR: getgpts, etype: ',i10,//,                          
     &       '               element not yet implemented',/,                    
     &       '               Gauss integration not defined' )                   
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine quad1                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 02/09/13 rhd               *          
c     *                                                              *          
c     *     given the gauss point number and the order of            *          
c     *     integration, this subroutine assigns the coordinates     *          
c     *     of the gauss point and its itegration weight for         *          
c     *     element 20,15,12 node elements                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine quad1( order, gpn, xi, eta, zeta, w )                          
      implicit double precision (a-z)                                           
      integer order, gpn                                                        
      logical valid(9)                                                          
      data valid / .true., .true., .false., .false., .false.,                   
     &             .false., .false., .true., .true. /                           
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
c                                                                               
      if ( order .lt. 1 .or. order .gt. 9 ) then                                
        write(*,9000) order                                                     
        call die_abort                                                          
      end if                                                                    
      if ( .not. valid(order) ) then                                            
        write(*,9000) order     
        call die_abort                                                 
      end if                                                                    
c                                                                               
      select case( order )                                                      
c     ====================                                                      
      case( 1 )                                                                 
c                                                                               
c                       27 point rule  (3x3x3)                                  
c                                                                               
      l1 = 0.774596669241483d0                                                  
      l2 = 0.0d0                                                                
      w1 = 0.1714677640603571d0                                                 
      w2 = 0.2743484224965712d0                                                 
      w3 = 0.4389574759945135d0                                                 
      w4 = 0.7023319615912210d0                                                 
c                                                                               
c                       branch on gauss point number to get the proper          
c                       information.                                            
c                                                                               
      go to (1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,            
     &          1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,         
     &          1230,1240,1250,1260,1270), gpn                                  
c                                                                               
 1010    xi= -l1                                                                
         eta= -l1                                                               
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1020    xi= -l1                                                                
         eta= l2                                                                
         zeta= -l1                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1030    xi= -l1                                                                
         eta= l1                                                                
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1040    xi= l2                                                                 
         eta= -l1                                                               
         zeta= -l1                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1050    xi= l2                                                                 
         eta= l2                                                                
         zeta= -l1                                                              
         w= w3                                                                  
         go to 9999                                                             
c                                                                               
 1060    xi= l2                                                                 
         eta= l1                                                                
         zeta= -l1                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1070    xi= l1                                                                 
         eta= -l1                                                               
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1080    xi= l1                                                                 
         eta= l2                                                                
         zeta= -l1                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1090    xi= l1                                                                 
         eta= l1                                                                
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1100    xi= -l1                                                                
         eta= -l1                                                               
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1110    xi= -l1                                                                
         eta= l2                                                                
         zeta= l2                                                               
         w= w3                                                                  
         go to 9999                                                             
c                                                                               
 1120    xi= -l1                                                                
         eta= l1                                                                
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1130    xi= l2                                                                 
         eta= -l1                                                               
         zeta= l2                                                               
         w= w3                                                                  
         go to 9999                                                             
c                                                                               
 1140    xi= l2                                                                 
         eta= l2                                                                
         zeta= l2                                                               
         w= w4                                                                  
         go to 9999                                                             
c                                                                               
 1150    xi= l2                                                                 
         eta= l1                                                                
         zeta= l2                                                               
         w= w3                                                                  
         go to 9999                                                             
c                                                                               
 1160    xi= l1                                                                 
         eta= -l1                                                               
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1170    xi= l1                                                                 
         eta= l2                                                                
         zeta= l2                                                               
         w= w3                                                                  
         go to 9999                                                             
c                                                                               
 1180    xi= l1                                                                 
         eta= l1                                                                
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1190    xi= -l1                                                                
         eta= -l1                                                               
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1200    xi= -l1                                                                
         eta= l2                                                                
         zeta= l1                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1210    xi= -l1                                                                
         eta= l1                                                                
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1220    xi= l2                                                                 
         eta= -l1                                                               
         zeta= l1                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1230    xi= l2                                                                 
         eta= l2                                                                
         zeta= l1                                                               
         w= w3                                                                  
         go to 9999                                                             
c                                                                               
 1240    xi= l2                                                                 
         eta= l1                                                                
         zeta= l1                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1250    xi= l1                                                                 
         eta= -l1                                                               
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 1260    xi= l1                                                                 
         eta= l2                                                                
         zeta= l1                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 1270    xi= l1                                                                 
         eta= l1                                                                
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
c                                                                               
      case( 2 )                                                                 
c                                                                               
c                       9 point rule (2x2x2 + 1)                                
c                                                                               
         g1   = 0.774596669241483d0                                             
         w1   = 0.555555555555556d0                                             
         w2   = 3.555555555555556d0                                             
         zero = 0.0d0                                                           
         go to (2010,2020,2030,2040,2050,2060,2070,2080,2090), gpn              
c                                                                               
 2010       xi = -g1                                                            
            eta = -g1                                                           
            zeta = -g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2020       xi = -g1                                                            
            eta =  g1                                                           
            zeta = -g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2030       xi =  g1                                                            
            eta = -g1                                                           
            zeta = -g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2040       xi =  g1                                                            
            eta =  g1                                                           
            zeta = -g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2050       xi = -g1                                                            
            eta = -g1                                                           
            zeta =  g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2060       xi = -g1                                                            
            eta =  g1                                                           
            zeta =  g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2070       xi =  g1                                                            
            eta = -g1                                                           
            zeta =  g1                                                          
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2080       xi = g1                                                             
            eta = g1                                                            
            zeta = g1                                                           
            w = w1                                                              
            go to 9999                                                          
c                                                                               
 2090       xi = zero                                                           
            eta = zero                                                          
            zeta = zero                                                         
            w = w2                                                              
            go to 9999                                                          
c                                                                               
      case( 8 )                                                                 
c                                                                               
c                       8 point rule (2x2x2)                                    
c                                                                               
         l1= 0.577350269189626d0                                                
         w1= 1.0d0                                                              
         go to (8010,8020,8030,8040,8050,8060,8070,8080), gpn                   
c                                                                               
 8010    xi= -l1                                                                
         eta= -l1                                                               
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8020    xi= -l1                                                                
         eta= l1                                                                
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8030    xi= l1                                                                 
         eta= -l1                                                               
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8040    xi= l1                                                                 
         eta= l1                                                                
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8050    xi= -l1                                                                
         eta= -l1                                                               
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8060    xi= -l1                                                                
         eta= l1                                                                
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8070    xi= l1                                                                 
         eta= -l1                                                               
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 8080    xi= l1                                                                 
         eta= l1                                                                
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
c                                                                               
      case( 9 )                                                                 
c                                                                               
c                       14 point rule.                                          
c                                                                               
         l1= 0.795822426d0                                                      
         l2= 0.758786911d0                                                      
         l3= 0.0d0                                                              
         w1= 0.886426593d0                                                      
         w2= 0.335180055d0                                                      
c                                                                               
c                       branch on gauss point number to get the proper          
c                       information.                                            
c                                                                               
         go to (9010,9020,9030,9040,9050,9060,9070,9080,9090,9100,9110,         
     &          9120,9130,9140), gpn                                            
c                                                                               
 9010    xi= -l1                                                                
         eta= l3                                                                
         zeta= l3                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 9020    xi= l1                                                                 
         eta= l3                                                                
         zeta= l3                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 9030    xi= l3                                                                 
         eta= -l1                                                               
         zeta= l3                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
c                                                                               
 9040    xi= l3                                                                 
         eta= l1                                                                
         zeta= l3                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 9050    xi= l3                                                                 
         eta= l3                                                                
         zeta= -l1                                                              
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 9060    xi= l3                                                                 
         eta= l3                                                                
         zeta= l1                                                               
         w= w1                                                                  
         go to 9999                                                             
c                                                                               
 9070    xi= -l2                                                                
         eta= -l2                                                               
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9080    xi= -l2                                                                
         eta= -l2                                                               
         zeta= -l2                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9090    xi= -l2                                                                
         eta= l2                                                                
         zeta= -l2                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9100    xi= -l2                                                                
         eta= l2                                                                
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9110    xi= l2                                                                 
         eta= -l2                                                               
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9120    xi= l2                                                                 
         eta= -l2                                                               
         zeta= -l2                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9130    xi= l2                                                                 
         eta= l2                                                                
         zeta= -l2                                                              
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
 9140    xi= l2                                                                 
         eta= l2                                                                
         zeta= l2                                                               
         w= w2                                                                  
         go to 9999                                                             
c                                                                               
      case default                                                              
        write(*,9000) order                                                     
c                                                                               
      end select                                                                
c     ==========                                                                
c                                                                               
 9999 return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                20, 15 or 12-node element. order= ',i9,           
     & /,    '                job terminated.',//)                              
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine quad2                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 11/14/89                   *          
c     *                                                              *          
c     *     given the gauss point number and the order of            *          
c     *     integration, this subroutine assigns the coordinates     *          
c     *     of the gauss point and its itegration weight for         *          
c     *     element l3disop                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine quad2( order, gpn, xi, eta, zeta, w )                          
      implicit integer (a-z)                                                    
      double precision                                                          
     &  xi,eta,zeta,w,l1,w1,l2                                                  
      double precision                                                          
     &  zero, one, three, four, root33                                          
      data zero, one, three, four, root33                                       
     & / 0.0d0, 1.0d0, 3.0d0, 4.0d0, 0.5773502691896257645d0 /                  
c                                                                               
c                                                                               
c                       branch on the order of quadrature. possible             
c                       integration schemes are :                               
c                                                                               
c                          order = 1 --  8 point rule ( 2x2x2 )                 
c                                                                               
      if ( order .ne. 1 ) then                                                  
        write(*,9000) order                                                     
        call die_abort                                                          
      end if                                                                    
c                                                                               
c                    8 point rule                                               
c                                                                               
      l1 = root33                                                               
      w1 = one                                                                  
c                                                                               
c                    branch on gauss point number to get the proper             
c                    information.                                               
c                                                                               
      go to (10,20,30,40,50,60,70,80) gpn                                       
c                                                                               
 10   xi= -l1                                                                   
      eta= -l1                                                                  
      zeta= -l1                                                                 
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 20   xi= -l1                                                                   
      eta= l1                                                                   
      zeta= -l1                                                                 
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 30   xi= l1                                                                    
      eta= -l1                                                                  
      zeta= -l1                                                                 
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 40   xi= l1                                                                    
      eta= l1                                                                   
      zeta= -l1                                                                 
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 50   xi= -l1                                                                   
      eta= -l1                                                                  
      zeta= l1                                                                  
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 60   xi= -l1                                                                   
      eta= l1                                                                   
      zeta= l1                                                                  
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 70   xi= l1                                                                    
      eta= -l1                                                                  
      zeta= l1                                                                  
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 80   xi= l1                                                                    
      eta= l1                                                                   
      zeta= l1                                                                  
      w= w1                                                                     
c                                                                               
 9999 continue                                                                  
      return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                8-node element. order= ',i3,                      
     & /,    '                job terminated.',//)                              
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine quad6                        *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/19/98                   *          
c     *                                                              *          
c     *     Given the current gauss point number and the order of    *          
c     *     integration (the total number of Gauss points), this     *          
c     *     subroutine assigns the coordinates of the gauss point    *          
c     *     and its integration weight for element "tet10".          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c        Gauss point coordinates and integration weight for tetrahedron         
c        elements. The 10-node tetrahedron is element type 6 in Warp3D;         
c        use the subroutine name "quad6" for the Gauss quadrature               
c        integration for this element type.                                     
c                                                                               
c        Input the integration order (the total number of Gauss points)         
c        and the current point number to get the tetrahedron coordinates        
c        and integration weight at that Gauss integration point.                
c                                                                               
c        Variables:                                                             
c                                                                               
c          order = specifies the integration order, and is the number of        
c                  Gauss integration points:  number of Gauss points = 1,4,5.   
c                  Note: if the value for order does not match a known value,   
c                  return all zero values as a default.  Also write an error    
c                  message and stop the analysis.                               
c          point = current integration point, 1 <= point <= numGauss            
c          s2,s3,s4 = natural tetrahedron coordinates returned, note that       
c                     the s1 coordinate is dependant on s2,s3,s4 and            
c                     doesn't need to be returned; if needed it can be          
c                     computed elsewhere.                                       
c                     property of coordinates:  s1 + s2 + s3 + s4 = 1           
c                     and 0 <= s1,s2,s3,s4 <= 1                                 
c                                                                               
c          NOTE that the coordinates returned must be compatible with the       
c          coordinates expected in deriv6 to compute shape function             
c          derivatives for the tetrahedron.                                     
c                                                                               
c          weight = Gauss integration weight at the current point returned      
c                                                                               
c                                                                               
      subroutine quad6( order, point, s2, s3, s4, weight )                      
      implicit none                                                             
      integer order,point                                                       
      double precision                                                          
     &         s2,s3,s4,weight                                                  
c                                                                               
      double precision                                                          
     &        a,b,w,s1, r, r1, t1, b1, r2, t2, b2, u, v, c,sqrt15               
      logical valid                                                             
c                                                                               
c         check the requested integration order for permitted values.           
c                                                                               
      valid = order .eq. 1 .or. order .eq. 4 .or. order .eq. 5 .or.             
     &        order .eq. 15                                                     
c                                                                               
      if( .not. valid ) then                                                    
        write(*,9000) order                                                     
        stop                                                                    
      end if                                                                    
c                                                                               
c         get the Gauss point coordinates and weights.                          
c                                                                               
      if( order.eq.1 ) then                                                     
c                                                                               
c             1-point rule, linear function, point at the center of             
c             the tetrahedron                                                   
c                                                                               
        a = 0.25D0                                                              
        w = 1.0D0                                                               
        s1 = a                                                                  
        s2 = a                                                                  
        s3 = a                                                                  
        s4 = a                                                                  
        weight = w                                                              
c                                                                               
      else if( order.eq.4 ) then                                                
c                                                                               
c             4-point rule, quadratic functions, equal weight for all           
c             points points just inside the tetrahedron corners                 
c                                                                               
        a = 0.58541020D0                                                        
        b = 0.13819660D0                                                        
        w = 0.25D0                                                              
        weight = w                                                              
c                                                                               
        if( point.eq.1 ) then                                                   
            s1 = a                                                              
            s2 = b                                                              
            s3 = b                                                              
            s4 = b                                                              
        else if( point.eq.2 ) then                                              
            s1 = b                                                              
            s2 = a                                                              
            s3 = b                                                              
            s4 = b                                                              
        else if( point.eq.3 ) then                                              
            s1 = b                                                              
            s2 = b                                                              
            s3 = a                                                              
            s4 = b                                                              
        else if( point.eq.4 ) then                                              
            s1 = b                                                              
            s2 = b                                                              
            s3 = b                                                              
            s4 = a                                                              
        end if                                                                  
c                                                                               
      else if( order.eq.5 ) then                                                
c                                                                               
c              5-point rule, cubic functions, point 1 at the center of          
c              the tetrahedron, note the negative weight value at point 1       
c              points 2,3,4,5 inside the tetrahedron corners                    
c                                                                               
        if( point.eq.1 ) then                                                   
            a = 0.25D0                                                          
            w = -4.0D0/5.0D0                                                    
        else                                                                    
            a = 1.0D0/2.0D0                                                     
            b = 1.0D0/6.0D0                                                     
            w = 9.0D0/20.0D0                                                    
        end if                                                                  
c                                                                               
        if( point.eq.1 ) then                                                   
            s1 = a                                                              
            s2 = a                                                              
            s3 = a                                                              
            s4 = a                                                              
            weight = w                                                          
        else if( point.eq.2 ) then                                              
            s1 = a                                                              
            s2 = b                                                              
            s3 = b                                                              
            s4 = b                                                              
            weight = w                                                          
        else if( point.eq.3 ) then                                              
            s1 = b                                                              
            s2 = a                                                              
            s3 = b                                                              
            s4 = b                                                              
            weight = w                                                          
        else if( point.eq.4 ) then                                              
            s1 = b                                                              
            s2 = b                                                              
            s3 = a                                                              
            s4 = b                                                              
            weight = w                                                          
        else if( point.eq.5 ) then                                              
            s1 = b                                                              
            s2 = b                                                              
            s3 = b                                                              
            s4 = a                                                              
            weight = w                                                          
        end if                                                                  
c                                                                               
      else if ( order.eq.15 ) then                                              
c                                                                               
c                                                                               
c             fifteen point rule. this is the scheme referenced                 
c             in the abaqus theory manual for quadratic tet                     
c             mass matrix computation.                                          
c                                                                               
c             ref: "approximate calculation of multiple integrals"              
c                   - a. h. stroud, 1971, prentice-hall inc.                    
c                     pp. 315, section 5-1                                      
c                                                                               
c                                                                               
          sqrt15 = dsqrt(15.0D0)                                                
c                                                                               
        if( point .eq. 1) then                                                  
          r = 0.25D0                                                            
          a = 16.D0/135.D0                                                      
        else if( point .ge. 2 .and. point .le. 5 ) then                         
          r1 = (7.0D0 - sqrt15)/34.0D0                                          
          t1 = (13.0D0 + 3.D0*sqrt15)/34.0D0                                    
          b1 = (2665.0D0 + 14.D0*sqrt15)/37800.0D0                              
        else if( point .ge. 6 .and. point .le. 9) then                          
          r2 = (7.0D0 + sqrt15)/34.0D0                                          
          t2 = (13.0D0 - 3.D0*sqrt15)/34.0D0                                    
          b2 = (2665.0D0 - 14.D0*sqrt15)/37800.0D0                              
        else if( point .ge. 10 .and. point .le. 15) then                        
          u = (10.0D0 - 2.D0*sqrt15)/40.0D0                                     
          v = (10.0D0 + 2.D0*sqrt15)/40.0D0                                     
          c = 20.D0/378.0D0                                                     
        end if                                                                  
c                                                                               
c                                                                               
        if( point .eq. 1 ) then                                                 
              s1 = r                                                            
              s2 = r                                                            
              s3 = r                                                            
              s4 = r                                                            
              weight = a                                                        
        else if( point .eq. 2 ) then                                            
              s1 = r1                                                           
              s2 = r1                                                           
              s3 = r1                                                           
              s4 = t1                                                           
              weight = b1                                                       
        else if( point .eq. 3 ) then                                            
              s1 = r1                                                           
              s2 = r1                                                           
              s3 = t1                                                           
              s4 = r1                                                           
              weight = b1                                                       
        else if( point .eq. 4 ) then                                            
              s1 = r1                                                           
              s2 = t1                                                           
              s3 = r1                                                           
              s4 = r1                                                           
              weight = b1                                                       
        else if( point .eq. 5 ) then                                            
              s1 = t1                                                           
              s2 = r1                                                           
              s3 = r1                                                           
              s4 = r1                                                           
              weight = b1                                                       
        else if( point .eq. 6 ) then                                            
              s1 = r2                                                           
              s2 = r2                                                           
              s3 = r2                                                           
              s4 = t2                                                           
              weight = b2                                                       
        else if( point .eq. 7 ) then                                            
              s1 = r2                                                           
              s2 = r2                                                           
              s3 = t2                                                           
              s4 = r2                                                           
              weight = b2                                                       
        else if( point .eq. 8 ) then                                            
              s1 = r2                                                           
              s2 = t2                                                           
              s3 = r2                                                           
              s4 = r2                                                           
              weight = b2                                                       
        else if( point .eq. 9 ) then                                            
              s1 = t2                                                           
              s2 = r2                                                           
              s3 = r2                                                           
              s4 = r2                                                           
              weight = b2                                                       
        else if( point .eq. 10 ) then                                           
              s1 = u                                                            
              s2 = u                                                            
              s3 = v                                                            
              s4 = v                                                            
              weight = c                                                        
        else if( point .eq. 11 ) then                                           
              s1 = u                                                            
              s2 = v                                                            
              s3 = v                                                            
              s4 = u                                                            
              weight = c                                                        
        else if( point .eq. 12 ) then                                           
              s1 = v                                                            
              s2 = u                                                            
              s3 = v                                                            
              s4 = u                                                            
              weight = c                                                        
        else if( point .eq. 13 ) then                                           
              s1 = v                                                            
              s2 = v                                                            
              s3 = u                                                            
              s4 = u                                                            
              weight = c                                                        
        else if( point .eq. 14 ) then                                           
              s1 = u                                                            
              s2 = v                                                            
              s3 = u                                                            
              s4 = v                                                            
              weight = c                                                        
        else if( point .eq. 15 ) then                                           
              s1 = v                                                            
              s2 = u                                                            
              s3 = u                                                            
              s4 = v                                                            
              weight = c                                                        
        end if                                                                  
c                                                                               
      end if                                                                    
c                                                                               
      return                                                                    
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                10-node tetrahedron element. order= ',i3,         
     & /,    '                (order = 4 is recommended)',                      
     & /,    '                job terminated.',//)                              
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine quad10                       *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 10/24/03                   *          
c     *                              by : mcw                        *          
c     *                                                              *          
c     *     Given the current gauss point number and the order of    *          
c     *     integration (the total number of Gauss points), this     *          
c     *     subroutine assigns the coordinates of the gauss point    *          
c     *     and its itegration weight for elements "quad8" or        *          
c     *     "axiquad8".                                              *          
c     *                                                              *          
c     *     This subroutine also works for 4-noded quads.            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c        gauss point coordinates and integration weight for quadrilateral       
c        elements.  Use the 1-D Gauss integration twice.                        
c        In Warp3D the 8-node axisymmetric quadrilateral element is             
c        element type 10, use "quad10" for this element type.                   
c                                                                               
c        Input the integration order (the total number of Gauss points)         
c        and the current point number to get the triangular coordinates         
c        and integration weight at that Gauss integration point.                
c                                                                               
c        Variables:                                                             
c                                                                               
c         order = specifies the integration order, and is the number of         
c                 Gauss integration points in each direction                    
c         order = 1,2,3,4 for the                                               
c                number of Gauss points = 1,4 for 2x2,9 for 3x3,16 for 4x4      
c         Note: if the value for numGauss does not match a known value,         
c         return all zero values as a default.                                  
c         point = current integration point, 1 <= point <= numGauss             
c         s1,s2 = Xi,Eta = natural quadrilateral coordinates returned,          
c                          -1 <= s1,s2 <= 1                                     
c         s3 = dummy coordinate for "Zeta", set to zero                         
c         weight = Gauss integration weight at the current point returned       
c                                                                               
      subroutine quad10( order, point, s1, s2, s3, weight )                     
      implicit none                                                             
      integer order,point                                                       
      double precision                                                          
     &         s1, s2, s3, weight                                               
c                                                                               
c           local variables                                                     
c                                                                               
      double precision                                                          
     &        w1, w2, zero                                                      
      integer p1, p2                                                            
      logical flag, debug                                                       
c                                                                               
      data zero                                                                 
     & / 0.0d0 /                                                                
c                                                                               
      flag  = .false.                                                           
      debug = .false.                                                           
c                                                                               
c           get the Gauss point coordinates and weights,                        
c           for quadrilaterals, use the 1-D Gauss quadrature twice,             
c           one for the s1 (Xi) direction and another for the s2                
c           (Eta) direction.                                                    
c                                                                               
      if( order.eq.1 ) then                                                     
        p1 = 1                                                                  
        call gauss1d( order, p1, s1, w1 )                                       
        call gauss1d( order, p1, s2, w2 )                                       
        weight = w1 * w2                                                        
        flag   = .true.                                                         
c                                                                               
      else if( order.eq.2 ) then                                                
c                                                                               
c           2x2 gauss quadrature, 4 integration points total                    
c           set the two point indexes to call the 1-D Gauss rule                
c                                                                               
        if( point.eq.1 ) then                                                   
          p1 = 1                                                                
          p2 = 1                                                                
        else if( point.eq.2 ) then                                              
          p1 = 1                                                                
          p2 = 2                                                                
        else if( point.eq.3 ) then                                              
          p1 = 2                                                                
          p2 = 1                                                                
        else if( point.eq.4 ) then                                              
          p1 = 2                                                                
          p2 = 2                                                                
        end if                                                                  
        call gauss1d( order, p1, s1, w1 )                                       
        call gauss1d( order, p2, s2, w2 )                                       
        weight = w1 * w2                                                        
        flag   = .true.                                                         
c                                                                               
      else if( order.eq.3 ) then                                                
c                                                                               
c           3x3 gauss quadrature, 9 integration points total                    
c           set the two point indexes to call the 1-D Gauss rule                
c                                                                               
        if( point.eq.1 ) then                                                   
          p1 = 1                                                                
          p2 = 1                                                                
        else if( point.eq.2 ) then                                              
          p1 = 1                                                                
          p2 = 2                                                                
        else if( point.eq.3 ) then                                              
          p1 = 1                                                                
          p2 = 3                                                                
        else if( point.eq.4 ) then                                              
          p1 = 2                                                                
          p2 = 1                                                                
        else if( point.eq.5 ) then                                              
          p1 = 2                                                                
          p2 = 2                                                                
        else if( point.eq.6 ) then                                              
          p1 = 2                                                                
          p2 = 3                                                                
        else if( point.eq.7 ) then                                              
          p1 = 3                                                                
          p2 = 1                                                                
        else if( point.eq.8 ) then                                              
          p1 = 3                                                                
          p2 = 2                                                                
        else if( point.eq.9 ) then                                              
          p1 = 3                                                                
          p2 = 3                                                                
        end if                                                                  
        call gauss1d( order, p1, s1, w1 )                                       
        call gauss1d( order, p2, s2, w2 )                                       
        weight = w1 * w2                                                        
        flag   = .true.                                                         
c                                                                               
      else if( order.eq.4 ) then                                                
c                                                                               
c          4x4 gauss quadrature, 16 integration points total                    
c          set the two point indexes to call the 1-D Gauss rule                 
c                                                                               
        if( point.eq.1 ) then                                                   
          p1 = 1                                                                
          p2 = 1                                                                
        else if( point.eq.2 ) then                                              
          p1 = 1                                                                
          p2 = 2                                                                
        else if( point.eq.3 ) then                                              
          p1 = 1                                                                
          p2 = 3                                                                
        else if( point.eq.4 ) then                                              
          p1 = 1                                                                
          p2 = 4                                                                
        else if( point.eq.5 ) then                                              
          p1 = 2                                                                
          p2 = 1                                                                
        else if( point.eq.6 ) then                                              
          p1 = 2                                                                
          p2 = 2                                                                
        else if( point.eq.7 ) then                                              
          p1 = 2                                                                
          p2 = 3                                                                
        else if( point.eq.8 ) then                                              
          p1 = 2                                                                
          p2 = 4                                                                
        else if( point.eq.9 ) then                                              
          p1 = 3                                                                
          p2 = 1                                                                
        else if( point.eq.10 ) then                                             
          p1 = 3                                                                
          p2 = 2                                                                
        else if( point.eq.11 ) then                                             
          p1 = 3                                                                
          p2 = 3                                                                
        else if( point.eq.12 ) then                                             
          p1 = 3                                                                
          p2 = 4                                                                
        else if( point.eq.13 ) then                                             
          p1 = 4                                                                
          p2 = 1                                                                
        else if( point.eq.14 ) then                                             
          p1 = 4                                                                
          p2 = 2                                                                
        else if( point.eq.15 ) then                                             
          p1 = 4                                                                
          p2 = 3                                                                
        else if( point.eq.16 ) then                                             
          p1 = 4                                                                
          p2 = 4                                                                
        end if                                                                  
        call gauss1d( order, p1, s1, w1 )                                       
        call gauss1d( order, p2, s2, w2 )                                       
        weight = w1 * w2                                                        
        flag   = .true.                                                         
c                                                                               
      end if                                                                    
c                                                                               
      if( debug ) write(*,*) point, order                                       
      if( debug ) write(*,*) s1, s2, w1, w2, weight                             
c                                                                               
c           set s3=0 (dummy coordinate for planar or axisymmetric elements).    
c                                                                               
      s3 = zero                                                                 
c                                                                               
c           check the requested integration order for permitted values.         
c                                                                               
      if( .not. flag ) then                                                     
         write(*,9000) order                                                    
         call die_abort                                                         
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '             axisymmetric or planar          '                    
     & /,    '             quadrilateral element:  order=  ',i3,                
     & /,    '             (order = 3 is recommended)',                         
     & /,    '             job terminated.',//)                                 
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gauss1d                      *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/19/98                   *          
c     *                                                              *          
c     *     1-D Gauss integration coordinates and weights.           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c            1-d gauss quadrature integration coordinate and weight             
c            values for a selected integration order (1,2,3,4 points)           
c            at a given point.                                                  
c                                                                               
c            Variables:                                                         
c                                                                               
c            num = number of Gauss integration points, the integration order    
c            p   = current integration point number, 1 <= p <= num              
c            s   = 1-D coordinate at the current integration point,             
c                   -1 <= s <= 1                                                
c            w = integration weight at the current point                        
c                                                                               
      subroutine gauss1d( num, p, s, w )                                        
      implicit none                                                             
      integer num,p                                                             
      double precision                                                          
     &         s, w                                                             
c                                                                               
      double precision                                                          
     &         zero, one, two, a1, a2, w1, w2                                   
c                                                                               
      zero = 0.0D0                                                              
      one  = 1.0D0                                                              
      two  = 2.0D0                                                              
c                                                                               
      if( num.eq.1 ) then                                                       
        s = zero                                                                
        w = two                                                                 
      else if( num.eq.2 ) then                                                  
        a1 = 0.5773502692D0                                                     
        w = one                                                                 
        if( p.eq.1 ) then                                                       
          s = -a1                                                               
        else if( p.eq.2 ) then                                                  
          s = a1                                                                
        end if                                                                  
      else if( num.eq.3 ) then                                                  
        a1 = 0.7745966692D0                                                     
        w1 = 0.5555555556D0                                                     
        w2 = 0.8888888889D0                                                     
        if( p.eq.1 ) then                                                       
          s = -a1                                                               
          w = w1                                                                
        else if( p.eq.2 ) then                                                  
          s = zero                                                              
          w = w2                                                                
        else if( p.eq.3 ) then                                                  
          s = a1                                                                
          w = w1                                                                
        end if                                                                  
      else if( num.eq.4 ) then                                                  
        a1 = 0.8611363116D0                                                     
        a2 = 0.3399810436D0                                                     
        w1 = 0.3478548451D0                                                     
        w2 = 0.6521451549D0                                                     
        if( p.eq.1 ) then                                                       
          s = -a1                                                               
          w = w1                                                                
        else if( p.eq.2 ) then                                                  
          s = -a2                                                               
          w = w2                                                                
        else if( p.eq.3 ) then                                                  
          s = a2                                                                
          w = w2                                                                
        else if( p.eq.4 ) then                                                  
          s = a1                                                                
          w = w1                                                                
        end if                                                                  
      else                                                                      
c                                                                               
c          assume a 1-point rule for an unknown number of points                
c                                                                               
        s = zero                                                                
        w = two                                                                 
c                                                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine quad11                       *          
c     *                                                              *          
c     *                       written by : gvt                       *          
c     *                                                              *          
c     *                   last modified : 08/19/98                   *          
c     *                                                              *          
c     *     Given the current gauss point number and the order of    *          
c     *     integration (the total number of Gauss points), this     *          
c     *     subroutine assigns the coordinates of the gauss point    *          
c     *     and its itegration weight for elements "tri6" or         *          
c     *     "axitri6".                                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         gauss point coordinates and integration weight for                    
c         triangle elements. in Warp3D the axisymmetric 6-node                  
c         triangle is element type 11, use "quad11" for the Gauss               
c         quadrature integration for this element type.                         
c                                                                               
c         Input the integration order (the total number of Gauss points)        
c         and the current point number to get the triangular coordinates        
c         and integration weight at that Gauss integration point.               
c                                                                               
c         Variables:                                                            
c                                                                               
c         order = specifies the integration order, and is the number of         
c                 Gauss integration points                                      
c                 number of Gauss points = 1,3,3,4,6,7                          
c                 (two possible 3-point rules)                                  
c         Note: if the value for order does not match a known value,            
c               write an error message and stop the analysis.                   
c         point = current integration point, 1 <= point <= numGauss             
c         s1,s2,s3 = natural triangular coordinates returned,                   
c         property of coordinates:  s1 + s2 + s3 = 1                            
c                                   0 <= s1,s2,s3 <= 1                          
c         weight = Gauss integration weight at the current point returned       
c                                                                               
c                                                                               
      subroutine quad11( order, point, s1, s2, s3, weight )                     
      implicit none                                                             
      integer order,point                                                       
      double precision                                                          
     &         s1, s2, s3, weight                                               
c                                                                               
c        Local Variables:                                                       
c                                                                               
c         subcase = integer flag set in this subroutine to choose on of the     
c                   two 3-point integration rules. For axisymmetric             
c                   elements use the second choice with the points              
c                   just inside the triangle corners.                           
c                                                                               
      integer subcase                                                           
      double precision                                                          
     &        sqrt15,a1,a2,b1,b2,w1,w2,w3                                       
      logical valid(7)                                                          
      data valid / .true.,.false.,.true.,.true.,.false.,.true.,.true. /         
c                                                                               
c          check the requested integration order for permitted values.          
c                                                                               
      if( order.lt.1 .or. order.gt.7 ) then                                     
        write(*,9000) order                                                     
        call die_abort                                                          
      end if                                                                    
      if( .not.valid(order) ) then                                              
        write(*,9000)                                                           
        call die_abort                                                          
      end if                                                                    
c                                                                               
c          get the gauss point coordinates and weights                          
c                                                                               
      if( order.eq.1 ) then                                                     
c                                                                               
c         1-point rule, linear function, point at the center of the triangle    
c                                                                               
        a1 = 1.0D0/3.0D0                                                        
        w1 = 1.0D0                                                              
        s1 = a1                                                                 
        s2 = a1                                                                 
        s3 = a1                                                                 
        weight = w1                                                             
c                                                                               
      else if( order.eq.3 ) then                                                
c                                                                               
c         3-point rules, quadratic functions, equal weight at all               
c         three points choose one of the two 3-point integration rules          
        subcase = 2                                                             
        w1 = 1.0D0/3.0D0                                                        
        weight = w1                                                             
        if( subcase.eq.1 ) then                                                 
c                                                                               
c         3-point midpoint rule, points on the element edge midpoints           
c                                                                               
          a1 = 0.5D0                                                            
          b1 = 0.0D0                                                            
c                                                                               
          if( point.eq.1 ) then                                                 
            s1 = a1                                                             
            s2 = a1                                                             
            s3 = b1                                                             
          else if( point.eq.2 ) then                                            
            s1 = b1                                                             
            s2 = a1                                                             
            s3 = a1                                                             
          else                                                                  
            s1 = a1                                                             
            s2 = b1                                                             
            s3 = a1                                                             
          end if                                                                
        else if( subcase.eq.2 ) then                                            
c                                                                               
c          3-point rule, points just inside the triangle corners,               
c          note: this is the one to use with axisymmetric triangles             
          a1 = 1.0D0/6.0D0                                                      
          b1 = 2.0D0/3.0D0                                                      
c                                                                               
          if( point.eq.1 ) then                                                 
            s1 = b1                                                             
            s2 = a1                                                             
            s3 = a1                                                             
          else if( point.eq.2 ) then                                            
            s1 = a1                                                             
            s2 = b1                                                             
            s3 = a1                                                             
          else                                                                  
            s1 = a1                                                             
            s2 = a1                                                             
            s3 = b1                                                             
          end if                                                                
      end if                                                                    
                                                                                
      else if( order.eq.4 ) then                                                
c                                                                               
c           4-point rule, cubic functions, note the negative weight             
c           for point 1 at the triangle center point 1 at center,               
c           points 2,3,4 just inside triangle corners                           
c                                                                               
        if( point.eq.1 ) then                                                   
          a1 = 1.0D0/3.0D0                                                      
          w1 = -27.0D0/48.0D0                                                   
        else                                                                    
          a1 = 0.6D0                                                            
          b1 = 0.2D0                                                            
          w2 = 25.0D0/48.0D0                                                    
        end if                                                                  
c                                                                               
        if( point.eq.1 ) then                                                   
          s1 = a1                                                               
          s2 = a1                                                               
          s3 = a1                                                               
          weight = w1                                                           
        else if( point.eq.2 ) then                                              
          s1 = a1                                                               
          s2 = b1                                                               
          s3 = b1                                                               
          weight = w2                                                           
        else if( point.eq.3 ) then                                              
          s1 = b1                                                               
          s2 = a1                                                               
          s3 = b1                                                               
          weight = w2                                                           
        else if( point.eq.4 ) then                                              
          s1 = b1                                                               
          s2 = b1                                                               
          s3 = a1                                                               
          weight = w2                                                           
        end if                                                                  
c                                                                               
      else if( order.eq.6 ) then                                                
c                                                                               
c         6-point rule, quartic functions                                       
c         points 1,2,3 just inside the triangle corners                         
c         points 4,5,6 just inside the element edge mid points                  
c                                                                               
        if( point.le.3 ) then                                                   
          a1 = 0.8168475730D0                                                   
          b1 = 0.0915762135D0                                                   
          w1 = 0.1099517437D0                                                   
        else                                                                    
          a2 = 0.1081030182D0                                                   
          b2 = 0.4459484909D0                                                   
          w2 = 0.2233815897D0                                                   
        end if                                                                  
c                                                                               
        if( point.eq.1 ) then                                                   
          s1 = a1                                                               
          s2 = b1                                                               
          s3 = b1                                                               
          weight = w1                                                           
        else if( point.eq.2 ) then                                              
          s1 = b1                                                               
          s2 = a1                                                               
          s3 = b1                                                               
          weight = w1                                                           
        else if( point.eq.3 ) then                                              
          s1 = b1                                                               
          s2 = b1                                                               
          s3 = a1                                                               
          weight = w1                                                           
        else if( point.eq.4 ) then                                              
          s1 = a2                                                               
          s2 = b2                                                               
          s3 = b2                                                               
          weight = w2                                                           
        else if( point.eq.5 ) then                                              
          s1 = b2                                                               
          s2 = a2                                                               
          s3 = b2                                                               
          weight = w2                                                           
        else if( point.eq.6 ) then                                              
          s1 = b2                                                               
          s2 = b2                                                               
          s3 = a2                                                               
          weight = w2                                                           
        end if                                                                  
c                                                                               
      else if( order.eq.7 ) then                                                
c                                                                               
c            7-point rule                                                       
c            point 1 at center of triangle                                      
c            points 3,4,5 just inside the triangle corners                      
c            points 6,7,8 just inside the element edge mid points               
c                                                                               
        if( point.gt.1 .and. point.le.4 ) then                                  
          sqrt15 = dsqrt(15.0D0)                                                
          a1 = (9.0D0 + 2.0D0*sqrt15)/21.0D0                                    
          b1 = (6.0D0 - sqrt15)/21.0D0                                          
          w2 = (155.0D0 - sqrt15)/1200.0D0                                      
        else if( point.ge.5 ) then                                              
          sqrt15 = dsqrt(15.0D0)                                                
          a2 = (9.0D0 - 2.0D0*sqrt15)/21.0D0                                    
          b2 = (6.0D0 + sqrt15)/21.0D0                                          
          w3 = (31.0D0/120.0D0) - (155.0D0 - sqrt15)/1200.0D0                   
        else                                                                    
          a1 = 1.0D0/3.0D0                                                      
          w1 = 9.0D0/40.0D0                                                     
        end if                                                                  
c                                                                               
        if( point.eq.1 ) then                                                   
          s1 = a1                                                               
          s2 = a1                                                               
          s3 = a1                                                               
          weight = w1                                                           
        else if( point.eq.2 ) then                                              
          s1 = a1                                                               
          s2 = b1                                                               
          s3 = b1                                                               
          weight = w2                                                           
        else if( point.eq.3 ) then                                              
          s1 = b1                                                               
          s2 = a1                                                               
          s3 = b1                                                               
          weight = w2                                                           
        else if( point.eq.4 ) then                                              
          s1 = b1                                                               
          s2 = b1                                                               
          s3 = a1                                                               
          weight = w2                                                           
        else if( point.eq.5 ) then                                              
          s1 = a2                                                               
          s2 = b2                                                               
          s3 = b2                                                               
          weight = w3                                                           
        else if( point.eq.6 ) then                                              
          s1 = b2                                                               
          s2 = a2                                                               
          s3 = b2                                                               
          weight = w3                                                           
        else if( point.eq.7 ) then                                              
          s1 = b2                                                               
          s2 = b2                                                               
          s3 = a2                                                               
          weight = w3                                                           
        end if                                                                  
c                                                                               
      end if                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                6-node axisymmetric or planar triangle'           
     & /,    '                element:  order= ',i3,                            
     & /,    '                (order = 3 is recommended)',                      
     & /,    '                job terminated.',//)                              
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine inter_gp                     *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                                                              *          
c     *                   last modified : 05/26/99                   *          
c     *                                                              *          
c     *     given the gauss point number and the order of            *          
c     *     integration, this subroutine assigns the coordinates     *          
c     *     of the gauss point and its itegration weight for         *          
c     *     element inter_8                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine inter_gp( order, gpn, xi, eta, zeta, w )                       
      implicit integer (a-z)                                                    
      double precision                                                          
     &  xi,eta,zeta,w,l1,w1                                                     
      double precision                                                          
     &  zero, one, three, four, root33                                          
      data zero, one, three, four, root33                                       
     & / 0.0d0, 1.0d0, 3.0d0, 4.0d0, 0.5773502691896257645d0 /                  
c                                                                               
c                       branch on the order of quadrature. possible             
c                       integration schemes are :                               
c                                                                               
c         order = 4    2x2 nodal lumping                                        
c         order = 5    2x2 gauss                                                
c         order = 6    1 point gauss                                            
c                                                                               
      if ( order .ne. 4 .and. order .ne. 5 .and. order .ne. 6 ) then            
        write(*,9000) order                                                     
        call die_abort                                                          
      end if                                                                    
c                                                                               
      zeta = zero                                                               
      if ( order .eq. 6 ) then                                                  
        l1 = zero                                                               
        w1 = four                                                               
      else if ( order .eq. 5 ) then                                             
        l1 = root33                                                             
        w1 = one                                                                
      else                                                                      
        l1 = one                                                                
        w1 = one                                                                
      end if                                                                    
c                                                                               
c                    branch on gauss point number to get the proper             
c                    information.                                               
c                                                                               
      go to (10,20,30,40), gpn                                                  
c                                                                               
 10   xi= -l1                                                                   
      eta= -l1                                                                  
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 20   xi=  l1                                                                   
      eta= -l1                                                                  
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 30   xi= l1                                                                    
      eta= l1                                                                   
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 40   xi= -l1                                                                   
      eta= l1                                                                   
      w= w1                                                                     
      go to 9999                                                                
c                                                                               
 9999 continue                                                                  
      return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                8-node  cohesive element. order= ',i3,            
     & /,    '                job terminated.',//)                              
c                                                                               
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *   subroutine trint6_gp    (6 node triangular interface)      *          
c     *                                                              *          
c     *                 written by : sushovan                        *          
c     *                                                              *          
c     *            last modified : 02/09/13 rhd                      *          
c     *                                                              *          
c     *     given the integration point number and the order of      *          
c     *     integration, return parametric coordinates of            *          
c     *     integration point and integration weight                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         Variables:                                                            
c                                                                               
c         order = integration scheme = 1, 2, 3                                  
c                 number of points   = 1, 3, 3                                  
c                 (two possible 3-point rules)                                  
c                                                                               
c         point = integration point                                             
c         s1,s2,s3 = natural triangular coordinates.                            
c         property of coordinates:  s1 + s2 + s3 = 1                            
c         For interface element s1,s2,zeta(=0) returned,                        
c         weight = integration weight at point returned                         
c                                                                               
c         integration schemes:                                                  
c                                                                               
c          order = 1    1 point rule, point at triangle centroid                
c                          (exact for linear functions)                         
c          order = 2    3 point rule, points just inside triangle corners       
c                          (quadratic)                                          
c          order = 3    3 point midpoint rule, points on edge midpoints         
c                          (quadratic)                                          
c                                                                               
c         references:                                                           
c              1.  The Finite Element Method, vol. 1                            
c                    by Zienkiewicz and Taylor                                  
c              2.  Concepts and Application of FEM                              
c                    by Cook, Malkus and Plesha                                 
c                                                                               
      subroutine trint6_gp( order, point, s1, s2, zeta, weight )                
      implicit none                                                             
      integer order, point                                                      
      double precision                                                          
     &         s1, s2, s3, zeta, weight                                         
c                                                                               
c     Local variables                                                           
c                                                                               
      double precision                                                          
     &  zero, one, half, third, sixth, two_third                                
c                                                                               
      data zero, one, half                                                      
     & /0.0d0, 1.0d0, 0.5d0/                                                    
c                                                                               
      third = 1.0d0/3.0d0; sixth = 1.0d0/6.0d0                                  
      two_third = 2.0d0/3.0d0                                                   
c                                                                               
      if( (order.lt.1) .or. (order.gt.3) ) then                                 
        write(*,9000) order                                                     
        call die_abort                                                          
      end if                                                                    
c                                                                               
      zeta = zero                                                               
c                                                                               
      corder: select case( order )                                              
c     ============================                                              
                                                                                
      case( 1 ) ! 1 pt at triangle centroid                                     
        if( point .ne. 1 ) then                                                 
           write(*,9100) order, point                                           
           call die_abort                                                       
        end if                                                                  
        s1 = third                                                              
        s2 = third                                                              
        s3 = third                                                              
        weight = one                                                            
c                                                                               
      case( 2 ) ! 3 pts just inside triangle corners                            
        cpt1: select case( point )                                              
        case( 1)                                                                
          s1 = two_third                                                        
          s2 = sixth                                                            
          s3 = sixth                                                            
          weight = third                                                        
        case( 2 )                                                               
          s1 = sixth                                                            
          s2 = two_third                                                        
          s3 = sixth                                                            
          weight = third                                                        
        case( 3 )                                                               
          s1 = sixth                                                            
          s2 = sixth                                                            
          s3 = two_third                                                        
          weight = third                                                        
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt1                                                         
c                                                                               
      case( 3 ) ! 3 pts on edge midpoints                                       
        cpt2: select case( point )                                              
        case( 1 )                                                               
          s1 = half                                                             
          s2 = half                                                             
          s3 = zero                                                             
          weight = third                                                        
        case( 2 )                                                               
          s1 = zero                                                             
          s2 = half                                                             
          s3 = half                                                             
          weight = third                                                        
        case( 3 )                                                               
          s1 = half                                                             
          s2 = zero                                                             
          s3 = half                                                             
          weight = third                                                        
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt2                                                         
c                                                                               
      end select corder                                                         
c     =================                                                         
c                                                                               
c           integration of f(s1,s2) over area is:                               
c               sum [0.5*f(i)*w(i)], i = 1 ,num_int_pts                         
c           the 0.5 comes area of triangle = 0.5 * |J|.                         
c           here we include it in the weight for int point                      
c                                                                               
      weight = half * weight                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                6-node triangular interface element.',            
     & /,    '                order: ',i8,                                      
     & /,    '                job terminated.',//)                              
 9100 format('>> FATAL ERROR: invalid integration point for',                   
     & /,    '                6-node triangular interface element.',            
     & /,    '                order: ',i3,' point: ',i8,                        
     & /,    '                job terminated.',//)                              
c                                                                               
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                                                              *          
c     *   subroutine trint12_gp    (12 node triangular interface)    *          
c     *                                                              *          
c     *                 written by : sushovan                        *          
c     *                                                              *          
c     *            last modified : 02/09/13 rhd                      *          
c     *                                                              *          
c     *     given the integration point number and the order of      *          
c     *     integration, return parametric coordinates of            *          
c     *     integration point and integration weight                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         Variables:                                                            
c                                                                               
c         order = integration scheme  = 1, 2, 3, 4, 6, 7                        
c                 number of points    = 1, 3, 3, 4, 6, 7                        
c                                                                               
c         point = integration point                                             
c         s1,s2,s3 = natural triangular coordinates.                            
c         property of coordinates:  s1 + s2 + s3 = 1                            
c         For interface element s1,s2,zeta(=0) returned,                        
c         weight = integration weight at point returned                         
c                                                                               
c        integration schemes:                                                   
c                                                                               
c          order = 1  1 point rule, point at triangle centroid                  
c                       (exact for linear functions)                            
c          order = 2  3 point rule, points just inside triangle corners         
c                       (quadratic)                                             
c          order = 3  3 point midpoint rule, points on edge midpoints           
c                       (quadratic)                                             
c          order = 4  4 point rule, 1 point at triangle centroid,               
c                       3 points just inside triangle corners (cubic)           
c                       (involves negative weight)                              
c          order = 6  6 point rule, 3 points just inside triangle corners,      
c                       3 points just inside edge midpoints (quartic,4)         
c          order = 7  7 point rule, 1 point at triangle centroid,               
c                       3 points just inside triangle corners,                  
c                       3 points just inside edge midpoints (quintic,5)         
c                                                                               
c        References:                                                            
c              1.  The Finite Element Method, vol. 1                            
c                    by Zienkiewicz and Taylor                                  
c              2.  Concepts and Application of FEM                              
c                    by Cook, Malkus and Plesha                                 
c                                                                               
      subroutine trint12_gp( order, point, s1, s2, zeta, weight )               
      implicit none                                                             
      integer order,point                                                       
      double precision                                                          
     &         s1, s2, s3, zeta, weight                                         
c                                                                               
c        Local Variables:                                                       
c                                                                               
      double precision                                                          
     & half, zero, one, two, three, four, five, six,                            
     & sqrt15, a1, a2, b1, b2, w1, w2, w3, pt6, pt2, third                      
      logical valid(7)                                                          
c                                                                               
      data half, zero, one, two, three, four, five, six                         
     &  / 0.5d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /              
      data  pt6, pt2, sqrt15, third                                             
     &   /  0.6d0, 0.2d0, 3.872983346207417d0, 0.33333333333333333d0 /          
      data valid / .true.,.true.,.true.,.true.,.false.,.true.,.true. /          
c                                                                               
      if( (order.lt.1) .or. (order.gt.7) .or. (order.eq.5) ) then               
        write(*,9000) order                                                     
        call die_abort                                                          
      end if                                                                    
      if( .not.valid(order) ) then                                              
        write(*,9000)                                                           
        call die_abort                                                          
      end if                                                                    
c                                                                               
      zeta = zero                                                               
c                                                                               
      corder: select case( order )                                              
c     ============================                                              
c                                                                               
      case( 1 )  ! 1 pt at the center                                           
        a1 = third                                                              
        s1 = a1                                                                 
        s2 = a1                                                                 
        s3 = a1                                                                 
        weight = one                                                            
c                                                                               
      case( 2 )  ! 3 pts just inside triangle corners                           
        a1 = one / six                                                          
        b1 = two / three                                                        
        w1 = third                                                              
        cpt2: select case( point )                                              
        case( 1)                                                                
            s1 = b1                                                             
            s2 = a1                                                             
            s3 = a1                                                             
            weight = w1                                                         
        case( 2 )                                                               
            s1 = a1                                                             
            s2 = b1                                                             
            s3 = a1                                                             
            weight = w1                                                         
        case( 3 )                                                               
            s1 = a1                                                             
            s2 = a1                                                             
            s3 = b1                                                             
            weight = w1                                                         
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt2                                                         
c                                                                               
      case( 3 )  ! 3 pts on  edge midpoints                                     
        a1 = half                                                               
        b1 = zero                                                               
        w1 = third                                                              
        cpt3: select case( point )                                              
        case( 1 )                                                               
            s1 = a1                                                             
            s2 = a1                                                             
            s3 = b1                                                             
            weight = w1                                                         
        case( 2 )                                                               
            s1 = b1                                                             
            s2 = a1                                                             
            s3 = a1                                                             
            weight = w1                                                         
          case( 3 )                                                             
            s1 = a1                                                             
            s2 = b1                                                             
            s3 = a1                                                             
            weight = w1                                                         
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt3                                                         
c                                                                               
      case( 4 )  ! 4 pts w/ 1 @ center                                          
        if( point .eq. 1 ) then                                                 
           a1 = third                                                           
           w1 = -27.0d0/48.0d0  ! note < 0 weight                               
        else                                                                    
           a2 = pt6                                                             
           b2 = pt2                                                             
           w2 = 25.0d0/48.0d0                                                   
        end if                                                                  
        cpt4: select case( point )                                              
        case( 1 )                                                               
          s1 = a1                                                               
          s2 = a1                                                               
          s3 = a1                                                               
          weight = w1                                                           
        case( 2 )                                                               
          s1 = a2                                                               
          s2 = b2                                                               
          s3 = b2                                                               
          weight = w2                                                           
        case( 3 )                                                               
          s1 = b2                                                               
          s2 = a2                                                               
          s3 = b2                                                               
          weight = w2                                                           
        case( 4 )                                                               
          s1 = b2                                                               
          s2 = b2                                                               
          s3 = a2                                                               
          weight = w2                                                           
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt4                                                         
c                                                                               
      case( 6 ) ! 6  pts                                                        
c                 1-3 just inside triangle corners                              
c                 4-6 just edge mid points                                      
        if( point.le.3 ) then                                                   
          a1 = 0.816847572980459d0                                              
          b1 = 0.091576213509771d0                                              
          w1 = 0.109951743655322d0                                              
        else                                                                    
          a2 = 0.108103018168070d0                                              
          b2 = 0.445948490915965d0                                              
          w2 = 0.223381589678011d0                                              
        end if                                                                  
        cpt6: select case( point )                                              
        case( 1 )                                                               
          s1 = a1                                                               
          s2 = b1                                                               
          s3 = b1                                                               
          weight = w1                                                           
        case( 2 )                                                               
          s1 = b1                                                               
          s2 = a1                                                               
          s3 = b1                                                               
          weight = w1                                                           
        case( 3 )                                                               
          s1 = b1                                                               
          s2 = b1                                                               
          s3 = a1                                                               
          weight = w1                                                           
        case( 4 )                                                               
          s1 = a2                                                               
          s2 = b2                                                               
          s3 = b2                                                               
          weight = w2                                                           
        case( 5 )                                                               
          s1 = b2                                                               
          s2 = a2                                                               
          s3 = b2                                                               
          weight = w2                                                           
        case( 6 )                                                               
          s1 = b2                                                               
          s2 = b2                                                               
          s3 = a2                                                               
          weight = w2                                                           
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt6                                                         
c                                                                               
      case( 7 )  ! 7 pts                                                        
c                1 @ center of triangle                                         
c                2-4 just inside triangle corners                               
c                5-7 just inside triangle edge mid-points                       
        if( point.ge.2  .and.  point.le.4 ) then                                
          a1 = (9.0d0 + two*sqrt15)/21.0d0                                      
          b1 = (six - sqrt15)/21.0d0                                            
          w2 = (155.0d0 - sqrt15)/1200.0d0                                      
        else if( point.ge.5 ) then                                              
          a2 = (9.0d0 - two*sqrt15)/21.0d0                                      
          b2 = (six + sqrt15)/21.0d0                                            
          w3 = (31.0d0/120.0d0) - (155.0d0 - sqrt15)/1200.0d0                   
        else                                                                    
          a1 = third                                                            
          w1 = 0.225d0   ! 9/40                                                 
        end if                                                                  
c                                                                               
        cpt7: select case( point )                                              
        case( 1 )                                                               
          s1 = a1                                                               
          s2 = a1                                                               
          s3 = a1                                                               
          weight = w1                                                           
        case( 2 )                                                               
          s1 = a1                                                               
          s2 = b1                                                               
          s3 = b1                                                               
          weight = w2                                                           
        case( 3 )                                                               
          s1 = b1                                                               
          s2 = a1                                                               
          s3 = b1                                                               
          weight = w2                                                           
        case( 4 )                                                               
          s1 = b1                                                               
          s2 = b1                                                               
          s3 = a1                                                               
          weight = w2                                                           
        case( 5 )                                                               
          s1 = a2                                                               
          s2 = b2                                                               
          s3 = b2                                                               
          weight = w3                                                           
        case( 6 )                                                               
          s1 = b2                                                               
          s2 = a2                                                               
          s3 = b2                                                               
          weight = w3                                                           
        case( 7 )                                                               
          s1 = b2                                                               
          s2 = b2                                                               
          s3 = a2                                                               
          weight = w3                                                           
        case default                                                            
           write(*,9100) order, point                                           
           call die_abort                                                       
        end select cpt7                                                         
c                                                                               
      end select corder                                                         
c     =================                                                         
c                                                                               
c           integration of f(s1,s2) over area is:                               
c               sum [0.5*f(i)*w(i)], i = 1 ,num_int_pts                         
c           the 0.5 comes area of triangle = 0.5 * |J|.                         
c           here we include it in the weight for int point                      
c                                                                               
      weight = half * weight                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                12-node triangular interface element;',           
     & /,    '                order: ',i4,                                      
     & /,    '                job terminated.',//)                              
 9100 format('>> FATAL ERROR: invalid integration point for',                   
     & /,    '                12-node triangular interface element.',           
     & /,    '                order: ',i4,' point: ',i8,                        
     & /,    '                job terminated.',//)                              
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine line3                        *          
c     *                                                              *          
c     *                       written by : mcw                       *          
c     *                                                              *          
c     *                   last modified : 10/10/03                   *          
c     *                                                              *          
c     *     given the gauss point number and the order of            *          
c     *     integration, this subroutine assigns the coordinates     *          
c     *     of the gauss point and its itegration weight for         *          
c     *     a quadratic line element                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine line3( order, gpn, xi, eta, zeta, w )                          
      implicit none                                                             
c                                                                               
c             dummy variables                                                   
c                                                                               
      integer order, gpn                                                        
      double precision                                                          
     &     xi, eta, zeta, w                                                     
c                                                                               
c             local variables                                                   
c                                                                               
      double precision                                                          
     &     op2, p6, zero, one, two, three, four, five, six, seven,              
     &     eight, nine                                                          
      data op2, p6, zero, one, two, three, four, five, six, seven,              
     &     eight, nine                                                          
     & / 0.2d0, 0.6d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0,                
     &   7.d0, 8.d0, 9.d0 /                                                     
c                                                                               
      if ( order.lt.1 .or. order.gt.4 ) then                                    
         write(*,9000) order                                                    
         stop                                                                   
      end if                                                                    
c                                                                               
c             branch on the order of quadrature. possible                       
c             integration schemes are :                                         
c                                                                               
c                          order = 1 -- 1 point rule                            
c                          order = 2 -- 2 point rule                            
c                          order = 3 -- 3 point rule                            
c                          order = 4 -- 4 point rule                            
c                                                                               
c             1 point rule                                                      
c                                                                               
      if( order .eq. 1 ) then                                                   
         xi = zero                                                              
         w  = two                                                               
         go to 9999                                                             
      end if                                                                    
c                                                                               
c             2 point rule                                                      
c                                                                               
      if( order .eq. 2 ) then                                                   
c                                                                               
c             branch on gauss point number to get the proper                    
c             information.                                                      
c                                                                               
         go to ( 210, 220 ) gpn                                                 
c                                                                               
 210     xi = - one/sqrt(three)                                                 
         w  = one                                                               
         go to 9999                                                             
c                                                                               
 220     xi = one/sqrt(three)                                                   
         w  = one                                                               
         go to 9999                                                             
c                                                                               
      end if                                                                    
c                                                                               
c             3 point rule                                                      
c                                                                               
      if( order .eq. 3 ) then                                                   
c                                                                               
c             branch on gauss point number to get the proper                    
c             information.                                                      
c                                                                               
         go to ( 310, 320, 330 ) gpn                                            
c                                                                               
 310     xi = - sqrt(p6)                                                        
         w  = five / nine                                                       
         go to 9999                                                             
c                                                                               
 320     xi = zero                                                              
         w  = eight / nine                                                      
         go to 9999                                                             
c                                                                               
 330     xi = sqrt(p6)                                                          
         w  = five / nine                                                       
         go to 9999                                                             
c                                                                               
      end if                                                                    
c                                                                               
c             4 point rule                                                      
c                                                                               
      if( order .eq. 4 ) then                                                   
c                                                                               
c             branch on gauss point number to get the proper                    
c             information.                                                      
c                                                                               
         go to ( 410, 420, 430, 440 ) gpn                                       
c                                                                               
 410     xi = - sqrt(three + two * sqrt(op2) / seven)                           
         w  = one / two - one / (six * sqrt(op2))                               
         go to 9999                                                             
c                                                                               
 420     xi = - sqrt(three - two * sqrt(op2) / seven)                           
         w  = one / two + one / (six * sqrt(op2))                               
         go to 9999                                                             
c                                                                               
 430     xi = sqrt(three - two * sqrt(op2) / seven)                             
         w  = one / two + one / (six * sqrt(op2))                               
         go to 9999                                                             
c                                                                               
 440     xi = sqrt(three + two * sqrt(op2) / seven)                             
         w  = one / two - one / (six * sqrt(op2))                               
         go to 9999                                                             
c                                                                               
      end if                                                                    
c                                                                               
 9999 continue                                                                  
      return                                                                    
c                                                                               
 9000 format('>> FATAL ERROR: invalid integration order for',                   
     & /,    '                3-node line element. order= ',i3,                 
     & /,    '                job terminated.',//)                              
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
