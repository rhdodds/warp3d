c ************************************************************                  
c *                                                          *                  
c *         isoparametric coordinates of element nodes       *                  
c *         by face number or element node number            *                  
c *                                                          *                  
c *                 modified by mcw on 8-21-03               *                  
c *                                                          *                  
c ************************************************************                  
c                                                                               
c                                                                               
      subroutine ndpts1( fnodes, nfnode, fcoor, etype, enode,                   
     &                   xi, eta, zeta )                                        
      implicit double precision (a-h,o-z)                                       
c                                                                               
c                                                                               
c              fill the fcoor table for element etype,                          
c              face "face" with isoparametric coordinates of                    
c              the face nodes.                                                  
c                                                                               
c                                                                               
      dimension    fnodes(*), fcoor(3,*)                                        
      real         coor20(3,20), coor12(3,12), coor15(3,15),                    
     &             coor9(3,9), coorq4(2,4), coorq8(2,8)                         
      integer      etype, fnodes, enode                                         
c                                                                               
      data  coor20  /                                                           
     &  -1.0,-1.0, 1.0,                                                         
     &  -1.0,-1.0,-1.0,                                                         
     &  -1.0, 1.0,-1.0,                                                         
     &  -1.0, 1.0, 1.0,                                                         
     &   1.0,-1.0, 1.0,                                                         
     &   1.0,-1.0,-1.0,                                                         
     &   1.0, 1.0,-1.0,                                                         
     &   1.0, 1.0, 1.0,                                                         
     &  -1.0,-1.0, 0.0,                                                         
     &  -1.0, 0.0,-1.0,                                                         
     &  -1.0, 1.0, 0.0,                                                         
     &  -1.0, 0.0, 1.0,                                                         
     &   1.0,-1.0, 0.0,                                                         
     &   1.0, 0.0,-1.0,                                                         
     &   1.0, 1.0, 0.0,                                                         
     &   1.0, 0.0, 1.0,                                                         
     &   0.0,-1.0, 1.0,                                                         
     &   0.0,-1.0,-1.0,                                                         
     &   0.0, 1.0,-1.0,                                                         
     &   0.0, 1.0, 1.0 /                                                        
c                                                                               
      data  coor9  /                                                            
     &  -1.0,-1.0, 1.0,                                                         
     &  -1.0,-1.0,-1.0,                                                         
     &  -1.0, 1.0,-1.0,                                                         
     &  -1.0, 1.0, 1.0,                                                         
     &   1.0,-1.0, 1.0,                                                         
     &   1.0,-1.0,-1.0,                                                         
     &   1.0, 1.0,-1.0,                                                         
     &   1.0, 1.0, 1.0,                                                         
     &  -1.0,-1.0, 0.0 /                                                        
c                                                                               
      data  coor12  /                                                           
     &  -1.0,-1.0, 1.0,                                                         
     &  -1.0,-1.0,-1.0,                                                         
     &  -1.0, 1.0,-1.0,                                                         
     &  -1.0, 1.0, 1.0,                                                         
     &   1.0,-1.0, 1.0,                                                         
     &   1.0,-1.0,-1.0,                                                         
     &   1.0, 1.0,-1.0,                                                         
     &   1.0, 1.0, 1.0,                                                         
     &  -1.0,-1.0, 0.0,                                                         
     &  -1.0, 0.0,-1.0,                                                         
     &  -1.0, 1.0, 0.0,                                                         
     &  -1.0, 0.0, 1.0 /                                                        
c                                                                               
      data   coor15  /                                                          
     &  -1.0,-1.0, 1.0,                                                         
     &  -1.0,-1.0,-1.0,                                                         
     &  -1.0, 1.0,-1.0,                                                         
     &  -1.0, 1.0, 1.0,                                                         
     &   1.0,-1.0, 1.0,                                                         
     &   1.0,-1.0,-1.0,                                                         
     &   1.0, 1.0,-1.0,                                                         
     &   1.0, 1.0, 1.0,                                                         
     &  -1.0,-1.0, 0.0,                                                         
     &  -1.0, 0.0,-1.0,                                                         
     &  -1.0, 1.0, 0.0,                                                         
     &  -1.0, 0.0, 1.0,                                                         
     &   1.0,-1.0, 0.0,                                                         
     &   0.0,-1.0, 1.0,                                                         
     &   0.0,-1.0,-1.0 /                                                        
c                                                                               
      data  coorq4 /                                                            
     &  -1.0,-1.0,                                                              
     &   1.0,-1.0,                                                              
     &   1.0, 1.0,                                                              
     &  -1.0, 1.0 /                                                             
c                                                                               
      data  coorq8 /                                                            
     &  -1.0,-1.0,                                                              
     &   1.0,-1.0,                                                              
     &   1.0, 1.0,                                                              
     &  -1.0, 1.0,                                                              
     &   0.0,-1.0,                                                              
     &   1.0, 0.0,                                                              
     &   0.0, 1.0,                                                              
     &  -1.0, 0.0 /                                                             
c                                                                               
c                                                                               
c                                                                               
      if( etype .eq. 9  ) etype = 6                                             
      if( etype .eq. 16 ) etype = 7                                             
c                                                                               
      go to ( 100, 200, 300, 400, 500, 900, 1600 ), etype                       
c                                                                               
c             20 node brick                                                     
c                                                                               
 100  continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coor20(1,enode)                                                  
        eta  = coor20(2,enode)                                                  
        zeta = coor20(3,enode)                                                  
        return                                                                  
      end if                                                                    
      do  inode = 1, 8                                                          
        j = fnodes(inode)                                                       
        fcoor(1,inode) = coor20(1,j)                                            
        fcoor(2,inode) = coor20(2,j)                                            
        fcoor(3,inode) = coor20(3,j)                                            
      end do                                                                    
c                                                                               
c             8 node brick                                                      
c                                                                               
 200  continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coor20(1,enode)                                                  
        eta  = coor20(2,enode)                                                  
        zeta = coor20(3,enode)                                                  
        return                                                                  
      end if                                                                    
      do  inode = 1, 4                                                          
        j = fnodes(inode)                                                       
        fcoor(1,inode) = coor20(1,j)                                            
        fcoor(2,inode) = coor20(2,j)                                            
        fcoor(3,inode) = coor20(3,j)                                            
      end do                                                                    
      return                                                                    
c                                                                               
c             12 node brick                                                     
c                                                                               
 300  continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coor12(1,enode)                                                  
        eta  = coor12(2,enode)                                                  
        zeta = coor12(3,enode)                                                  
        return                                                                  
      end if                                                                    
      do  inode = 1, nfnode                                                     
        j = fnodes(inode)                                                       
        fcoor(1,inode) = coor12(1,j)                                            
        fcoor(2,inode) = coor12(2,j)                                            
        fcoor(3,inode) = coor12(3,j)                                            
      end do                                                                    
      return                                                                    
c                                                                               
c             15 node brick                                                     
c                                                                               
 400  continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coor15(1,enode)                                                  
        eta  = coor15(2,enode)                                                  
        zeta = coor15(3,enode)                                                  
        return                                                                  
      end if                                                                    
      do  inode = 1, nfnode                                                     
        j = fnodes(inode)                                                       
        fcoor(1,inode) = coor15(1,j)                                            
        fcoor(2,inode) = coor15(2,j)                                            
        fcoor(3,inode) = coor15(3,j)                                            
      end do                                                                    
      return                                                                    
c                                                                               
c             9 node brick                                                      
c                                                                               
 500  continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coor9(1,enode)                                                   
        eta  = coor9(2,enode)                                                   
        zeta = coor9(3,enode)                                                   
        return                                                                  
      end if                                                                    
      do  inode = 1, nfnode                                                     
        j = fnodes(inode)                                                       
        fcoor(1,inode) = coor9(1,j)                                             
        fcoor(2,inode) = coor9(2,j)                                             
        fcoor(3,inode) = coor9(3,j)                                             
      end do                                                                    
      return                                                                    
c                                                                               
c                                                                               
c             8 node quad                                                       
c                                                                               
 900  continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coorq8(1,enode)                                                  
        eta  = coorq8(2,enode)                                                  
        zeta = 0.0                                                              
        return                                                                  
      end if                                                                    
      do inode = 1, 8                                                           
        fcoor(1,inode) = coorq8(1,inode)                                        
        fcoor(2,inode) = coorq8(2,inode)                                        
        fcoor(3,inode) = 0.0                                                    
      end do                                                                    
      return                                                                    
c                                                                               
c             4 node quad                                                       
c                                                                               
 1600 continue                                                                  
      if ( enode .gt. 0 ) then                                                  
        xi   = coorq4(1,enode)                                                  
        eta  = coorq4(2,enode)                                                  
        zeta = 0.0                                                              
        return                                                                  
      end if                                                                    
      do inode = 1, 4                                                           
        fcoor(1,inode) = coorq4(1,inode)                                        
        fcoor(2,inode) = coorq4(2,inode)                                        
        fcoor(3,inode) = 0.0                                                    
      end do                                                                    
      return                                                                    
c                                                                               
      end                                                                       
