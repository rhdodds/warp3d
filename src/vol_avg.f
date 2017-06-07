c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine vol_avg                      *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 03/16/2017 rhd             *          
c     *                                                              *          
c     *     compute averaged volume terms for bbar                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine vol_avg( vol, volume, span, mxvl )                             
      implicit none                                                             
c                                                                               
      integer :: span, mxvl                                                     
      double precision :: vol(mxvl,8,*), volume(*)                              
c                                                                               
      integer :: i                                                              
      double precision, parameter :: one = 1.0d0                                
      double precision :: voli                                                  
c                                                                               
c               refer to extensive comments in vol_terms.f                      
c                                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span                                                            
c                                                                               
       voli = one / volume(i)                                                   
                                                                                
c                                                                               
c                       term 1: partial Nj / partial x                          
c                                                                               
                vol(i,1,1 ) = vol(i,1,1) * voli                                 
                vol(i,2,1 ) = vol(i,2,1) * voli                                 
                vol(i,3,1 ) = vol(i,3,1) * voli                                 
                vol(i,4,1 ) = vol(i,4,1) * voli                                 
                vol(i,5,1 ) = vol(i,5,1) * voli                                 
                vol(i,6,1 ) = vol(i,6,1) * voli                                 
                vol(i,7,1 ) = vol(i,7,1) * voli                                 
                vol(i,8,1 ) = vol(i,8,1) * voli                                 
c                                                                               
c                     term 2: partial Nj / partial y                            
c                                                                               
                vol(i,1,2 ) = vol(i,1,2) * voli                                 
                vol(i,2,2 ) = vol(i,2,2) * voli                                 
                vol(i,3,2 ) = vol(i,3,2) * voli                                 
                vol(i,4,2 ) = vol(i,4,2) * voli                                 
                vol(i,5,2 ) = vol(i,5,2) * voli                                 
                vol(i,6,2 ) = vol(i,6,2) * voli                                 
                vol(i,7,2 ) = vol(i,7,2) * voli                                 
                vol(i,8,2 ) = vol(i,8,2) * voli                                 
c                                                                               
c                     term 3: partial Nj / partial z                            
c                                                                               
                vol(i,1,3 ) = vol(i,1,3) * voli                                 
                vol(i,2,3 ) = vol(i,2,3) * voli                                 
                vol(i,3,3 ) = vol(i,3,3) * voli                                 
                vol(i,4,3 ) = vol(i,4,3) * voli                                 
                vol(i,5,3 ) = vol(i,5,3) * voli                                 
                vol(i,6,3 ) = vol(i,6,3) * voli                                 
                vol(i,7,3 ) = vol(i,7,3) * voli                                 
                vol(i,8,3 ) = vol(i,8,3) * voli                                 
c                                                                               
      end do                                                                    
      return                                                                    
      end                                                                       
                                                                                
