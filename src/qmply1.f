c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine qmply1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 03/17/2017 rhd             *          
c     *                                                              *          
c     *      performs the rotation or transformation                 *          
c     *      of 3x3 tensors expressed in 6x1 vector form             *          
c     *      for a block of elements                                 *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine qmply1( span, mxvl, nstr, q, m1, m2 )                          
      implicit none                                                             
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer :: span, mxvl, nstr                                               
      double precision :: q(mxvl,nstr,*), m1(mxvl,*), m2(mxvl,*)                
c                                                                               
      integer :: i, j                                                           
c                                                                               
c                      {m2} = [q] * {m1}   (6x1 vectors and 6x6 q]              
c                                                                               
      if( nstr .eq. 6 ) then                                                    
        do j = 1, 6                                                             
!DIR$ VECTOR ALIGNED                                                            
           do i = 1, span                                                       
              m2(i,j)= q(i,j,1)*m1(i,1)+q(i,j,2)*m1(i,2)+                       
     &                 q(i,j,3)*m1(i,3)+q(i,j,4)*m1(i,4)+                       
     &                 q(i,j,5)*m1(i,5)+q(i,j,6)*m1(i,6)                        
           end do                                                               
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
      if( nstr .eq. 3 ) then   ! possible future 2-d elements                   
         do j = 1, 3                                                            
!DIR$ VECTOR ALIGNED                                                            
            do i = 1, span                                                      
               m2(i,j)= q(i,j,1)*m1(i,1)+q(i,j,2)*m1(i,2)+                      
     &                 q(i,j,3)*m1(i,3)+q(i,j,4)*m1(i,4)+                       
     &                 q(i,j,5)*m1(i,5)+q(i,j,6)*m1(i,6)                        
           end do                                                               
         end do                                                                 
         return                                                                 
      endif                                                                     
c                                                                               
      end                                                                       
