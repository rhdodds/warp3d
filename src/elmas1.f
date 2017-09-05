c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine elmas1                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 8/11/2017 rhd              *          
c     *                                                              *          
c     *     this subroutine scales the appropriate mass              *          
c     *     matrices and parses them to all dof for a block of       *          
c     *     similar, non-conflicting elements in                     *          
c     *     uniform global coordinates.                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine elmas1( span, nnode, emass, volume, mel, totdof )              
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                      parameter declarations                                   
c                                                                               
      integer ::  span, nnode, totdof                                           
      double precision :: volume(*), mel(totdof,*), emass(*)                                       
c                                                                               
c                      locals                                                   
c                                                                               
      integer :: i, j                                                           
      double precision :: scale(mxvl), trace(mxvl)
      double precision, parameter :: zero=0.0d0, mass_tol=1.0d-20                                
      logical, parameter :: local_debug = .false.                                                    
c                                                                               
c                       compute the element diagonal mass matrix.               
c                       compute the trace of the consistent diagonal            
c                       mass matrix.                                            
c                                                                               
      trace = zero                                                              
      scale = zero                                                              
      do j = 1, nnode                                                           
       do i = 1, span                                                           
         trace(i) = trace(i) + mel(j,i)                                         
       end do                                                                   
      end do                                                                    
c                                                                               
c                    compute the scale factor used to make the                  
c                    terms of the lumped mass matrix sum to the                 
c                    total masses for each of the three nodal                   
c                    degrees of freedom (total masses =                         
c                    3 * total mass of elements).                               
c                                                                               
      do i = 1, span                                                            
       if ( abs(trace(i)) .gt. mass_tol )                                       
     &        scale(i) =  emass(i) / trace(i)                                   
      end do                                                                    
                                                                                
      if ( local_debug ) then                                                   
        write(*,*) '>> emass:'                                                  
        write(*,*) (i,emass(i),i=1,span)                                        
      end if                                                                    
c                                                                               
c                    fill in the remaining terms of the                         
c                    lumped mass matrix corresponding to the                    
c                    remaining two global dof. multiply in                      
c                    the scale factor computed above.                           
c                                                                               
      do j = 1, nnode                                                           
       do i = 1, span                                                           
         mel(j,i) = mel(j,i) * scale(i)                                         
       end do                                                                   
      end do                                                                    
c                                                                               
      do j = 1, nnode                                                           
       do i = 1, span                                                           
         mel(j+nnode,i)   = mel(j,i)                                            
         mel(j+2*nnode,i) = mel(j,i)                                            
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
