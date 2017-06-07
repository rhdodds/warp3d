c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouext1                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 3/11/04 (rhd)              *          
c     *                                                              *          
c     *     this subroutine computes derived stress/strain values    *          
c     *     at gauss points or node points for later output.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouext1( span, stress, nnode, ngp, nodpts,                      
     &                   num_short_stress, num_short_strain,                    
     &                   op_code, mxvl )                                        
      use elblk_data, only : elestr                                             
      implicit integer (a-z)                                                    
      logical stress, nodpts                                                    
c                                                                               
c                       find the number of strain points.                       
c                                                                               
      nstrpt = ngp                                                              
      if ( nodpts ) nstrpt = nnode                                              
c                                                                               
c                       compute additional data for output at                   
c                       each of the element nodes or gauss points               
c                       of the element block.                                   
c                                                                               
c                       op_code = 1 : compute mises equiv. stress or            
c                                     equiv. strain.                            
c                       op_code = 2 : compute invariants, principal             
c                                     values.                                   
c                                                                               
      if ( op_code .eq. 1 ) then                                                
        do strpt = 1, nstrpt                                                    
         if( stress) then                                                       
            call ouyld1( span, elestr(1,1,strpt), elestr(1,8,strpt),            
     &                   mxvl )                                                 
         else                                                                   
            call oueff1( span, elestr(1,1,strpt), elestr(1,7,strpt),            
     &                   mxvl )                                                 
         end if                                                                 
       end do                                                                   
       return                                                                   
      end if                                                                    
c                                                                               
c                       compute invariants, principal values.                   
c                                                                               
      loc = num_short_strain + 1                                                
      if ( stress ) loc = num_short_stress + 1                                  
c                                                                               
      do strpt = 1, nstrpt                                                      
c                                                                               
c                       compute the principal invariants of desired             
c                       stress or strain.                                       
c                                                                               
         call ouinv1( span, elestr(1,1,strpt), elestr(1,loc,strpt),             
     &                stress )                                                  
c                                                                               
c                       compute the desired principal stresses or               
c                       principal strains and the direction                     
c                       cosines of their corresponding normals.                 
c                                                                               
         call oupri1( span, elestr(1,1,strpt), elestr(1,loc+3,strpt),           
     &                elestr(1,loc+6,strpt), stress )                           
c                                                                               
      end do                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouyld1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 1/20/2017 rhd              *          
c     *                                                              *          
c     *     this subroutine computes the mises equiv. stress at a    *          
c     *     given strain point for a block of elements               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouyld1( span, gpstr, yf, mxvl )                                
      implicit none                                                             
c                                                                               
      integer :: span, mxvl                                                     
      double precision :: gpstr(mxvl,*), yf(*)                                  
c                                                                               
      integer :: i                                                              
      double precision, parameter :: iroot2=0.70711d0, six=6.0                  
c                                                                               
c                       compute the von-mises stress.                           
c                                                                               
      do i = 1, span                                                            
         yf(i) = sqrt( (gpstr(i,1)-gpstr(i,2))**2+                              
     &                 (gpstr(i,2)-gpstr(i,3))**2+                              
     &                 (gpstr(i,1)-gpstr(i,3))**2+                              
     &             six*(gpstr(i,4)**2+gpstr(i,5)**2+                            
     &                  gpstr(i,6)**2) )*iroot2                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oueff1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 04/18/90                   *          
c     *                                                              *          
c     *     this subroutine computes the effective strain measure    *          
c     *     for a block of similar, non-conflicting q3disop or       *          
c     *     l3disop elements.                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oueff1( span, str, efeps, mxvl )                               
      implicit integer (a-z)                                                    
      double precision                                                          
     &     str(mxvl,*), efeps(*), root23, onep5                                 
      data root23, onep5 / 0.471404, 1.5 /                                      
c                                                                               
c                       compute the effective strain measure.                   
c                                                                               
      do i = 1, span                                                            
       efeps(i) = root23 *                                                      
     &      sqrt( (str(i,1)-str(i,2))**2+(str(i,2)-str(i,3))**2+                
     &            (str(i,1)-str(i,3))**2+                                       
     &            onep5*(str(i,4)**2+str(i,5)**2+str(i,6)**2) )                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouinv1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/16/89                   *          
c     *                                                              *          
c     *     this subroutine computes the invariants of the cauchy    *          
c     *     stress or the almansi strain at a strain point, either a *          
c     *     gauss point or a node point, for a block of similar,     *          
c     *     non-conflicting q3disop or l3disop elements.             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouinv1( span, str, inv, stress )                               
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision                                                          
     &     str(mxvl,*), inv(mxvl,*)                                             
      logical stress                                                            
c                                                                               
c                    locally allocated                                          
c                                                                               
      double precision                                                          
     &     tstr(mxvl,nstr), factor, one, half                                   
      data one, half / 1.0, 0.5 /                                               
c                                                                               
      factor = one                                                              
      if ( .not. stress ) factor = half                                         
      do i = 1, span                                                            
         tstr(i,1) = str(i,1)                                                   
         tstr(i,2) = str(i,2)                                                   
         tstr(i,3) = str(i,3)                                                   
         tstr(i,4) = str(i,4) * factor                                          
         tstr(i,5) = str(i,5) * factor                                          
         tstr(i,6) = str(i,6) * factor                                          
      enddo                                                                     
c                                                                               
c                       compute the principal invariants.                       
c                                                                               
      do i = 1, span                                                            
       inv(i,1) = tstr(i,1)+tstr(i,2)+tstr(i,3)                                 
       inv(i,2) = tstr(i,4)*tstr(i,4)+tstr(i,5)*tstr(i,5) +                     
     &            tstr(i,6)*tstr(i,6)-tstr(i,1)*tstr(i,2) -                     
     &            tstr(i,2)*tstr(i,3)-tstr(i,1)*tstr(i,3)                       
       inv(i,3)=                                                                
     &   tstr(i,1)*( tstr(i,2)*tstr(i,3)-tstr(i,5)*tstr(i,5) ) -                
     &   tstr(i,4)*( tstr(i,4)*tstr(i,3)-tstr(i,5)*tstr(i,6) ) +                
     &   tstr(i,6)*( tstr(i,4)*tstr(i,5)-tstr(i,2)*tstr(i,6) )                  
      end do                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupri1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/04/91                   *          
c     *                                                              *          
c     *     this subroutine computes the principal cauchy stresses   *          
c     *     or almansi strains and the direction cosines of their    *          
c     *     corresponding normals at a strain point, either a gauss  *          
c     *     point or a node point, for a block of similar, non-      *          
c     *     conflicting q3disop or l3disop elements.                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupri1( span, str, prstr, angles, stress )                     
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision                                                          
     &     str(mxvl,*), prstr(mxvl,*), angles(mxvl,ndim,*)                      
      logical stress                                                            
c                                                                               
c                    locally allocated                                          
c                                                                               
      double precision                                                          
     &     tstr(nstr,mxvl), wk(ndim,mxvl), ev(nstr,mxvl),                       
     &     evec(ndim,ndim,mxvl), factor, one, half                              
      data one, half / 1.0, 0.5 /                                               
c                                                                               
      factor = one                                                              
      if ( .not. stress ) factor = half                                         
      do i = 1, span                                                            
         tstr(1,i) = str(i,1)                                                   
         tstr(2,i) = str(i,4) * factor                                          
         tstr(3,i) = str(i,2)                                                   
         tstr(4,i) = str(i,6) * factor                                          
         tstr(5,i) = str(i,5) * factor                                          
         tstr(6,i) = str(i,3)                                                   
      end do                                                                    
c                                                                               
      do i = 1, span                                                            
         call ou3dpr( tstr(1,i), ndim, 1, ev(1,i), evec(1,1,i),                 
     &                ndim, wk(1,i), ier )                                      
      end do                                                                    
c                                                                               
c                                                                               
      do i = 1, span                                                            
         prstr(i,1)    = ev(1,i)                                                
         prstr(i,2)    = ev(2,i)                                                
         prstr(i,3)    = ev(3,i)                                                
         angles(i,1,1) = evec(1,1,i)                                            
         angles(i,1,2) = evec(1,2,i)                                            
         angles(i,1,3) = evec(1,3,i)                                            
         angles(i,2,1) = evec(2,1,i)                                            
         angles(i,2,2) = evec(2,2,i)                                            
         angles(i,2,3) = evec(2,3,i)                                            
         angles(i,3,1) = evec(3,1,i)                                            
         angles(i,3,2) = evec(3,2,i)                                            
         angles(i,3,3) = evec(3,3,i)                                            
      end do                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
