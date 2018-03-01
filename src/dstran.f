c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dstran                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/02/90                   *          
c     *                                                              *          
c     *     this subroutine displays the transformation matrix from  *          
c     *     uniform global corrdinates to constraint compatable glo- *          
c     *     bal coordinates for the given node.                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dstran( node )                                                 
      use global_data ! old common.main
c                                                                               
      use main_data, only : trnmat, trn, inverse_incidences                     
c                                                                               
      implicit integer (a-z)                                                    
      real dumr                                                                 
      double precision                                                          
     &     dumd                                                                 
      character :: dums                                                         
c                                                                               
c                                                                               
      ndof = iprops(4,inverse_incidences(node)%element_list(1))                 
c                                                                               
c                       print title                                             
c                                                                               
      write(out,9000) node                                                      
c                                                                               
c                       print transformation matrix for given node,             
c                       if it is indexed.                                       
c                                                                               
      if ( .not. trn(node) ) then                                               
         call errmsg(80,node,dums,dumr,dumd)                                    
         go to 9999                                                             
      end if                                                                    
c                                                                               
      do row = 1, ndof                                                          
         write(out,9010) row,(trnmat(node)%mat(row,col),col=1,ndof)             
      end do                                                                    
c                                                                               
      write(out,9005)                                                           
c                                                                               
c                                                                               
 9000 format(////20x,'transformation matrix for node ',i7/20x,37('-')//)        
c                                                                               
 9005 format(//1x,' ')                                                          
c                                                                               
 9010 format(9x,'row',i2,':',3(5x,f12.6)/14x,3(5x,f12.6)//)                     
c                                                                               
c                                                                               
 9999 return                                                                    
      end                                                                       
