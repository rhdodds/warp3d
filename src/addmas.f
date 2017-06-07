c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine addmas                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/02/91                   *          
c     *                                                              *          
c     *     assembles diagonal mass matrix for each block element    *          
c     *     into structure diagonal mass matrix                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine addmas( span, bedst, totdof, mdiag, mel )                      
      implicit none                                                             
c                                                                               
      integer :: span, totdof, bedst(totdof,*)                                  
      double precision ::                                                       
     &   mdiag(*), mel(totdof,*)                                                
c                                                                               
      integer :: i, j                                                           
c                                                                               
      do j = 1, totdof                                                          
         do i = 1, span                                                         
            mdiag(bedst(j,i)) = mdiag(bedst(j,i)) + mel(j,i)                    
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
