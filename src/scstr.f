c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine scstr                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 03/18/04 rhd               *          
c     *                                                              *          
c     *     this subroutine scatters element stresses to the global  *          
c     *     stress data structure from a block of similar, non-      *          
c     *     conflicting elements for all gauss points.               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine scstr( ml, mg, ngp, nprm, span )                               
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision                                                          
     &     ml(mxvl,nprm,*),mg(nprm,ngp,*)                                       
c                                                                               
c                                                                               
      if ( ngp .ne. 8 ) then                                                    
        do k = 1, ngp                                                           
           do j = 1, nprm                                                       
              do i = 1, span                                                    
                 mg(j,k,i) = ml(i,j,k)                                          
              end do                                                            
           end do                                                               
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c                       number of gauss points = 8                              
c                                                                               
      do j = 1, nprm                                                            
        do i = 1, span                                                          
          mg(j,1,i) = ml(i,j,1)                                                 
          mg(j,2,i) = ml(i,j,2)                                                 
          mg(j,3,i) = ml(i,j,3)                                                 
          mg(j,4,i) = ml(i,j,4)                                                 
          mg(j,5,i) = ml(i,j,5)                                                 
          mg(j,6,i) = ml(i,j,6)                                                 
          mg(j,7,i) = ml(i,j,7)                                                 
          mg(j,8,i) = ml(i,j,8)                                                 
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine scstr_history                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/18/04 rhd               *          
c     *                                                              *          
c     *     this subroutine scatters element stresses to the global  *          
c     *     stress data structure from a block of similar, non-      *          
c     *     conflicting elements for all gauss points.               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine scstr_history( local_history, global_history,                  
     &                          ngp, hist_size, span )                          
      implicit integer (a-z)                                                    
      double precision                                                          
     &      local_history(span,hist_size,ngp),                                  
     &      global_history(hist_size,ngp,span)                                  
c                                                                               
c                                                                               
      if ( ngp .ne. 8 ) then                                                    
        do k = 1, ngp                                                           
           do j = 1, hist_size                                                  
              do i = 1, span                                                    
                 global_history(j,k,i) = local_history(i,j,k)                   
              end do                                                            
           end do                                                               
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c                       number of gauss points = 8                              
c                                                                               
      do j = 1, hist_size                                                       
        do i = 1, span                                                          
          global_history(j,1,i) = local_history(i,j,1)                          
          global_history(j,2,i) = local_history(i,j,2)                          
          global_history(j,3,i) = local_history(i,j,3)                          
          global_history(j,4,i) = local_history(i,j,4)                          
          global_history(j,5,i) = local_history(i,j,5)                          
          global_history(j,6,i) = local_history(i,j,6)                          
          global_history(j,7,i) = local_history(i,j,7)                          
          global_history(j,8,i) = local_history(i,j,8)                          
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
