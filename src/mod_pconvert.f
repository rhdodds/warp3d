c                                                                               
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *     module: pvars                                                  *        
c *                                                                    *        
c *     contains global variables                                      *        
c *                                                                    *        
c *     Use saved_dp1, saved_dp2, saved_dp3, saved_step, saved_int1,   *        
c *     saved_int2, and saved_int3, when data from one packet          *        
c *     subroutine is needed in following packet subroutines.          *        
c *     If more than 3 double precision or integer variables are       *        
c *     needed, add additional variables to this module.               *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      module pvars                                                              
c                                                                               
      double precision saved_dp1, saved_dp2, saved_dp3                          
      integer saved_step, saved_int1, saved_int2, saved_int3                    
      logical debug                                                             
c                                                                               
      end module pvars                                                          
