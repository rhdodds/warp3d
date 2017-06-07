c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine zero_vector                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/28/94                   *          
c     *                                                              *          
c     *     zero a vector of specified length w/ floating zero       *          
c     *     signle or double based on this port                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine zero_vector( vec, n )                                          
      double precision                                                          
     &  vec(*), zero                                                            
      data zero / 0.0d00 /                                                      
c                                                                               
!DIR$ IVDEP                                                                     
      vec(1:n) = zero                                                           
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine vec_ops                      *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/27/96                   *          
c     *                                                              *          
c     *     perform simple vector-vector operations on single or     *          
c     *     double precision vectors (based on the port)             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine vec_ops( veca, vecb, vecc, n, opcode )                         
      double precision                                                          
     &  veca(*), vecb(*), vecc(*), zero, const                                  
      integer opcode                                                            
      data zero / 0.0d0/                                                        
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600, 700 ), opcode                       
c                                                                               
c            opcode 1:   c = a * b                                              
c                                                                               
 100  continue                                                                  
!DIR$ IVDEP                                                                     
      vecc(1:n) = veca(1:n) * vecb(1:n)                                         
      return                                                                    
c                                                                               
c            opcode 2:   b = b * a                                              
c                                                                               
 200  continue                                                                  
!DIR$ IVDEP                                                                     
      vecb(1:n) = vecb(1:n) * veca(1:n)                                         
      return                                                                    
c                                                                               
c            opcode 3:   c = a / b                                              
c                                                                               
 300  continue                                                                  
!DIR$ IVDEP                                                                     
      vecc(1:n) = veca(1:n) / vecb(1:n)                                         
      return                                                                    
c                                                                               
c            opcode v:   c = zero                                               
c                                                                               
 400  continue                                                                  
!DIR$ IVDEP                                                                     
      vecc(1:n) = zero                                                          
      return                                                                    
c                                                                               
c            opcode v:   a = b                                                  
c                                                                               
 500  continue                                                                  
!DIR$ IVDEP                                                                     
      veca(1:n) = vecb(1:n)                                                     
      return                                                                    
c                                                                               
c            opcode v:   a = a + b                                              
c                                                                               
 600  continue                                                                  
!DIR$ IVDEP                                                                     
      veca(1:n) = veca(1:n) + vecb(1:n)                                         
      return                                                                    
c                                                                               
c            opcode v:   a = const * b ; const = vecc(1)                        
c                                                                               
 700  continue                                                                  
      const = vecc(1)                                                           
!DIR$ IVDEP                                                                     
      veca(1:n) = const * vecb(1:n)                                             
      return                                                                    
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine warp3d_sort_float            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/28/00                   *          
c     *                                                              *          
c     *    sort a vector of floating point numbers and an associated *          
c     *    index vector into ascending order                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine warp3d_sort_float( n, dvec, index_vec )                        
      double precision                                                          
     &  dvec(*), a                                                              
      integer index_vec(*), b                                                   
c                                                                               
      do j = 2, n                                                               
        a = dvec(j)                                                             
        b = index_vec(j)                                                        
        do i = j-1, 1, -1                                                       
          if ( dvec(i) .le. a ) go to 10                                        
          dvec(i+1) = dvec(i)                                                   
          index_vec(i+1) = index_vec(i)                                         
        end do                                                                  
        i = 0                                                                   
 10     continue                                                                
        dvec(i+1) = a                                                           
        index_vec(i+1) = b                                                      
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
