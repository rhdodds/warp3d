c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine bmod                         *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 03/15/2017 rhd             *          
c     *                                                              *          
c     *     modifiy the linear-displacement [b] matrix of            *          
c     *     the l3disop element to implement b-bar.                  *          
c     *     include stability term to safeguard against              *          
c     *     spurious modes                                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine bmod ( b, vol, span, mxvl, eps_stab, mxedof )                  
      implicit none                                                             
c                                                                               
      integer :: span, mxvl, mxedof                                             
      double precision :: b(mxvl,mxedof,*), vol(mxvl,8,*), eps_stab             
c                                                                               
      integer :: i                                                              
      double precision :: alpha, beta,                                          
     &   dummy_1, dummy_2, dummy_3, dummy_4,                                    
     &   dummy_5, dummy_6, dummy_7, dummy_8                                     
      double precision, parameter :: one = 1.0d0, two = 2.0d0,                  
     &                               third = 1.0d0/3.0d0                        
c                                                                               
c                                                                               
c                  the vol data structure contains volume averaged              
c                  interpolating function derivatives for elements              
c                  in block. see vol_terms.f for complete description.          
c                  here we implement the b-bar modifications to the             
c                  usual b matrix. our b-bar dilatation terms follow            
c                  the mean dilation as suggested by Nagtegaal.                 
c                  we take the usual [B] matrix, subtract out the               
c                  dilational part and add back the mean-dilatation             
c                  for the element. the result is an element that does          
c                  not lock and has the same mean stress at each of             
c                  the 2 x 2 x 2 Gauss points. our procedure follows            
c                  that described by Highes book, p. 233-235. see also          
c                  the warp user guide for a summary of these eqns.             
c                                                                               
        alpha = (two + eps_stab)*third                                          
        beta  = (one - eps_stab)*third                                          
c                                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
        do i = 1, span                                                          
c                                                                               
c                  save value for use in later calculations                     
c                                                                               
          dummy_1  = b (i,1,1)                                                  
          dummy_2  = b (i,2,1)                                                  
          dummy_3  = b (i,3,1)                                                  
          dummy_4  = b (i,4,1)                                                  
          dummy_5  = b (i,5,1)                                                  
          dummy_6  = b (i,6,1)                                                  
          dummy_7  = b (i,7,1)                                                  
          dummy_8  = b (i,8,1)                                                  
c                                                                               
c                  modify first column of b matrix                              
c                                                                               
c        b (i,j,1) = (two * dummy + vol (i,j, 1)) * third                       
c                                                                               
        b (i,1,1) = alpha * dummy_1 + beta * vol (i,1, 1)                       
        b (i,2,1) = alpha * dummy_2 + beta * vol (i,2, 1)                       
        b (i,3,1) = alpha * dummy_3 + beta * vol (i,3, 1)                       
        b (i,4,1) = alpha * dummy_4 + beta * vol (i,4, 1)                       
        b (i,5,1) = alpha * dummy_5 + beta * vol (i,5, 1)                       
        b (i,6,1) = alpha * dummy_6 + beta * vol (i,6, 1)                       
        b (i,7,1) = alpha * dummy_7 + beta * vol (i,7, 1)                       
        b (i,8,1) = alpha * dummy_8 + beta * vol (i,8, 1)                       
c                                                                               
c                                                                               
c                 modify second column of b matrix                              
c                                                                               
c        b (i,j,2) =  (vol (i,j,1) - dummy) * third                             
c                                                                               
        b (i,1,2) = beta * ( vol (i,1,1) - dummy_1 )                            
        b (i,2,2) = beta * ( vol (i,2,1) - dummy_2 )                            
        b (i,3,2) = beta * ( vol (i,3,1) - dummy_3 )                            
        b (i,4,2) = beta * ( vol (i,4,1) - dummy_4 )                            
        b (i,5,2) = beta * ( vol (i,5,1) - dummy_5 )                            
        b (i,6,2) = beta * ( vol (i,6,1) - dummy_6 )                            
        b (i,7,2) = beta * ( vol (i,7,1) - dummy_7 )                            
        b (i,8,2) = beta * ( vol (i,8,1) - dummy_8 )                            
      end do                                                                    
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span                                                            
c                                                                               
c        b (i,bpos1+j,1) = (vol (i,j, 2) - b (i,bpos1+j,2))                     
c     &                      * third                                            
c                                                                               
c                                                                               
        b (i,9,1)  = beta * ( vol (i,1, 2) - b (i, 9,2) )                       
        b (i,10,1) = beta * ( vol (i,2, 2) - b (i,10,2) )                       
        b (i,11,1) = beta * ( vol (i,3, 2) - b (i,11,2) )                       
        b (i,12,1) = beta * ( vol (i,4, 2) - b (i,12,2) )                       
        b (i,13,1) = beta * ( vol (i,5, 2) - b (i,13,2) )                       
        b (i,14,1) = beta * ( vol (i,6, 2) - b (i,14,2) )                       
        b (i,15,1) = beta * ( vol (i,7, 2) - b (i,15,2) )                       
        b (i,16,1) = beta * ( vol (i,8, 2) - b (i,16,2) )                       
c                                                                               
c                                                                               
c        b (i,bpos2+j,1) = (vol (i,j, 3) - b (i,bpos2+j,3))                     
c     &                      * third                                            
c                                                                               
        b (i,17,1) = beta * ( vol (i,1,3) - b (i,17,3) )                        
        b (i,18,1) = beta * ( vol (i,2,3) - b (i,18,3) )                        
        b (i,19,1) = beta * ( vol (i,3,3) - b (i,19,3) )                        
        b (i,20,1) = beta * ( vol (i,4,3) - b (i,20,3) )                        
        b (i,21,1) = beta * ( vol (i,5,3) - b (i,21,3) )                        
        b (i,22,1) = beta * ( vol (i,6,3) - b (i,22,3) )                        
        b (i,23,1) = beta * ( vol (i,7,3) - b (i,23,3) )                        
        b (i,24,1) = beta * ( vol (i,8,3) - b (i,24,3) )                        
                                                                                
                                                                                
c                                                                               
c        b (i,bpos1+j,2) = (two * b (i,bpos1+j,2) + vol (i,j,2))                
c     &                      * third                                            
c                                                                               
        b (i,9,2)  =   alpha * b (i,9,2)  + beta * vol (i,1,2)                  
        b (i,10,2) =   alpha * b (i,10,2) + beta * vol (i,2,2)                  
        b (i,11,2) =   alpha * b (i,11,2) + beta * vol (i,3,2)                  
        b (i,12,2) =   alpha * b (i,12,2) + beta * vol (i,4,2)                  
        b (i,13,2) =   alpha * b (i,13,2) + beta * vol (i,5,2)                  
        b (i,14,2) =   alpha * b (i,14,2) + beta * vol (i,6,2)                  
        b (i,15,2) =   alpha * b (i,15,2) + beta * vol (i,7,2)                  
        b (i,16,2) =   alpha * b (i,16,2) + beta * vol (i,8,2)                  
                                                                                
c                                                                               
c        b (i,bpos2+j,2) =  b (i,bpos2+j,1)                                     
c                                                                               
       b (i,17,2) =  b (i,17,1)                                                 
       b (i,18,2) =  b (i,18,1)                                                 
       b (i,19,2) =  b (i,19,1)                                                 
       b (i,20,2) =  b (i,20,1)                                                 
       b (i,21,2) =  b (i,21,1)                                                 
       b (i,22,2) =  b (i,22,1)                                                 
       b (i,23,2) =  b (i,23,1)                                                 
       b (i,24,2) =  b (i,24,1)                                                 
c                                                                               
c           modify third column of b matrix                                     
c                                                                               
c        b (i,j,3) =  b (i,j,2)                                                 
c                                                                               
        b (i,1,3) =  b (i,1,2)                                                  
        b (i,2,3) =  b (i,2,2)                                                  
        b (i,3,3) =  b (i,3,2)                                                  
        b (i,4,3) =  b (i,4,2)                                                  
        b (i,5,3) =  b (i,5,2)                                                  
        b (i,6,3) =  b (i,6,2)                                                  
        b (i,7,3) =  b (i,7,2)                                                  
        b (i,8,3) =  b (i,8,2)                                                  
c                                                                               
c        b (i,bpos1+j,3) =  b (i,bpos2+j,1)                                     
c                                                                               
        b (i,9,3)  =  b (i,9,1)                                                 
        b (i,10,3) =  b (i,10,1)                                                
        b (i,11,3) =  b (i,11,1)                                                
        b (i,12,3) =  b (i,12,1)                                                
        b (i,13,3) =  b (i,13,1)                                                
        b (i,14,3) =  b (i,14,1)                                                
        b (i,15,3) =  b (i,15,1)                                                
        b (i,16,3) =  b (i,16,1)                                                
c                                                                               
c        b (i,bpos2+j,3) =  (two * b (i,bpos2+j,3) + vol(i,j,3))                
c     &                     * third                                             
c                                                                               
        b (i,17,3) =  alpha * b (i,17,3) + beta * vol(i,1,3)                    
        b (i,18,3) =  alpha * b (i,18,3) + beta * vol(i,2,3)                    
        b (i,19,3) =  alpha * b (i,19,3) + beta * vol(i,3,3)                    
        b (i,20,3) =  alpha * b (i,20,3) + beta * vol(i,4,3)                    
        b (i,21,3) =  alpha * b (i,21,3) + beta * vol(i,5,3)                    
        b (i,22,3) =  alpha * b (i,22,3) + beta * vol(i,6,3)                    
        b (i,23,3) =  alpha * b (i,23,3) + beta * vol(i,7,3)                    
        b (i,24,3) =  alpha * b (i,24,3) + beta * vol(i,8,3)                    
                                                                                
c                                                                               
c           remaining columns need not be modified because they                 
c           are deviatoric, not volumetric                                      
c                                                                               
        end do                                                                  
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
