c *******************************************************************           
c *                                                                 *           
c *    subroutine di1dsf ---- calculates shape function values and  *           
c *                           derivatives for 1-dimension           *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
      subroutine di1dsf( xsi, dsf, sf, nlnode )                                 
c                                                                               
c              parameter declarations                                           
c                                                                               
      double precision                                                          
     & xsi, dsf(*), sf(*)                                                       
c                                                                               
c              local declarations                                               
c                                                                               
      double precision                                                          
     & xsisqr, half, one, two                                                   
      data half, one, two / 0.5, 1.0, 2.0 /                                     
c                                                                               
      go to ( 200, 200, 300 ), nlnode                                           
c                                                                               
c             linear variation.                                                 
c                                                                               
 200  continue                                                                  
      sf( 1 )  = half * ( one - xsi )                                           
      sf( 2 )  = half * ( one + xsi )                                           
      dsf( 1 ) = -half                                                          
      dsf( 2 ) =  half                                                          
      return                                                                    
c                                                                               
c             quadratic variation                                               
c                                                                               
 300  continue                                                                  
      xsisqr = xsi * xsi                                                        
      sf( 1 )  = half * ( xsisqr - xsi )                                        
      sf( 2 )  = one - xsisqr                                                   
      sf( 3 )  = half * ( xsisqr + xsi )                                        
      dsf( 1 ) = xsi - half                                                     
      dsf( 2 ) = -two*xsi                                                       
      dsf( 3 ) = xsi + half                                                     
      return                                                                    
c                                                                               
      end                                                                       
