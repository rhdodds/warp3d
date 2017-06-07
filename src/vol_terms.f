c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine vol_terms                    *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 03/16/2017 rhd             *          
c     *                                                              *          
c     *          include this gauss point contributions to volume    *          
c     *          terms for bbar at all elements in block             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine vol_terms ( in_jacob, dj, vol, nxi, neta, nzeta,               
     &                       volume, span, mxvl )                               
      implicit none                                                             
c                                                                               
      integer :: span, mxvl                                                     
      double precision :: in_jacob(mxvl,3,*), dj(*), vol(mxvl,8,*),             
     &  nxi(*), neta(*), nzeta(*), volume(*)                                    
c                                                                               
      integer :: i                                                              
      double precision :: a, b, c, d                                            
c                                                                               
c             compute b-bar terms for the 8-node 3-D                            
c             element using the "mean dilatation" scheme.                       
c             vol is a span x 8 nodes x 3 derivatives array.                    
c             this routine executes once for each gauss point                   
c             in a block of elements. vol array is zeroed before                
c             computations start for block.                                     
c                                                                               
c             to examine values for an element look in an 'i'                   
c             plane of the array. at i,j,1 we compute                           
c             and store partial N sub j / partial x,                            
c             at i,j,2 we compute and store partial N sub j / partial y,        
c             and at i,j,3 we compute and store                                 
c             partial N sub j / partial z. when this routine has been           
c             called 8 times (2x2x2) for elements in the block,                 
c             i,j,1 has the summed partial N sub j / partial x                  
c             derivatives for element i. same for i,j,2 and i,j,3               
c             (y and z derivatives). when each of these values                  
c             is divided by the corresponding element volume,                   
c             the result is the "mean" derivatives needed for the               
c             mean-dilatation b-bar. the averaging is done the                  
c             later routine vol_avg.                                            
c                                                                               
c             loop over all elements in the block                               
c                                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span                                                            
c                                                                               
c                       calculate 1st term for all nodes:                       
c                       partial Nj / partial x                                  
c                                                                               
c                                                                               
      a = in_jacob(i,1,1)                                                       
      b = in_jacob(i,1,2)                                                       
      c = in_jacob(i,1,3)                                                       
      d = dj(i)                                                                 
      vol(i,1,1) = vol(i,1,1) + (a*nxi(1) + b*neta(1) + c*nzeta(1)) * d         
      vol(i,2,1) = vol(i,2,1) + (a*nxi(2) + b*neta(2) + c*nzeta(2)) * d         
      vol(i,3,1) = vol(i,3,1) + (a*nxi(3) + b*neta(3) + c*nzeta(3)) * d         
      vol(i,4,1) = vol(i,4,1) + (a*nxi(4) + b*neta(4) + c*nzeta(4)) * d         
      vol(i,5,1) = vol(i,5,1) + (a*nxi(5) + b*neta(5) + c*nzeta(5)) * d         
      vol(i,6,1) = vol(i,6,1) + (a*nxi(6) + b*neta(6) + c*nzeta(6)) * d         
      vol(i,7,1) = vol(i,7,1) + (a*nxi(7) + b*neta(7) + c*nzeta(7)) * d         
      vol(i,8,1) = vol(i,8,1) + (a*nxi(8) + b*neta(8) + c*nzeta(8)) * d         
c                                                                               
c                       calculate 2nd term for all nodes                        
c                       partial Nj / partial y                                  
c                                                                               
      a = in_jacob(i,2,1)                                                       
      b = in_jacob(i,2,2)                                                       
      c = in_jacob(i,2,3)                                                       
      vol(i,1,2) = vol(i,1,2) + (a*nxi(1) + b*neta(1) + c*nzeta(1)) * d         
      vol(i,2,2) = vol(i,2,2) + (a*nxi(2) + b*neta(2) + c*nzeta(2)) * d         
      vol(i,3,2) = vol(i,3,2) + (a*nxi(3) + b*neta(3) + c*nzeta(3)) * d         
      vol(i,4,2) = vol(i,4,2) + (a*nxi(4) + b*neta(4) + c*nzeta(4)) * d         
      vol(i,5,2) = vol(i,5,2) + (a*nxi(5) + b*neta(5) + c*nzeta(5)) * d         
      vol(i,6,2) = vol(i,6,2) + (a*nxi(6) + b*neta(6) + c*nzeta(6)) * d         
      vol(i,7,2) = vol(i,7,2) + (a*nxi(7) + b*neta(7) + c*nzeta(7)) * d         
      vol(i,8,2) = vol(i,8,2) + (a*nxi(8) + b*neta(8) + c*nzeta(8)) * d         
c                                                                               
c                       calculate 3rd term for all nodes                        
c                       partial Nj / partial z                                  
c                                                                               
      a = in_jacob(i,3,1)                                                       
      b = in_jacob(i,3,2)                                                       
      c = in_jacob(i,3,3)                                                       
      vol(i,1,3) = vol(i,1,3) + (a*nxi(1) + b*neta(1) + c*nzeta(1)) * d         
      vol(i,2,3) = vol(i,2,3) + (a*nxi(2) + b*neta(2) + c*nzeta(2)) * d         
      vol(i,3,3) = vol(i,3,3) + (a*nxi(3) + b*neta(3) + c*nzeta(3)) * d         
      vol(i,4,3) = vol(i,4,3) + (a*nxi(4) + b*neta(4) + c*nzeta(4)) * d         
      vol(i,5,3) = vol(i,5,3) + (a*nxi(5) + b*neta(5) + c*nzeta(5)) * d         
      vol(i,6,3) = vol(i,6,3) + (a*nxi(6) + b*neta(6) + c*nzeta(6)) * d         
      vol(i,7,3) = vol(i,7,3) + (a*nxi(7) + b*neta(7) + c*nzeta(7)) * d         
      vol(i,8,3) = vol(i,8,3) + (a*nxi(8) + b*neta(8) + c*nzeta(8)) * d         
c                                                                               
c                       calculate total volume of element one gp at a time      
c                                                                               
      volume(i) = volume(i) + d                                                 
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
