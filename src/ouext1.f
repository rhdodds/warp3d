c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouext1                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/21/20 (rhd)             *          
c     *                                                              *          
c     *     computes derived stress/strain values                    *          
c     *     at gauss points or node points for later output.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouext1( span, stresses, nnode, ngp, nodpts,                      
     &                   num_short_stress, num_short_strain,                    
     &                   op_code, mxvl )  
c                                      
      use elblk_data, only : elestr 
c                                    
      implicit none
c
      integer :: span, nnode, ngp, num_short_stress, num_short_strain,                    
     &           op_code, mxvl                                                   
      logical :: stresses, nodpts          
c
      integer :: nstrpt, strpt, invar_col, prin_col, evec_col                                        
c                                                                               
c                       find the number of strain points.                       
c                                                                               
      nstrpt = ngp                                                              
      if( nodpts ) nstrpt = nnode                                              
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
      if( op_code .eq. 1 ) then                                                
       do strpt = 1, nstrpt                                                    
         if( stresses) then                                                       
           call get_mises( span, elestr(1,1,strpt), elestr(1,8,strpt),            
     &                     mxvl )                                                 
         else                                                                   
            call get_equiv_strain( span, elestr(1,1,strpt), 
     &                             elestr(1,7,strpt), mxvl )                                                 
         end if                                                                 
       end do                                                                   
       return                                                                   
      end if                                                                    
c                                                                               
c                       compute invariants, principal values.                   
c                                                                               
      invar_col = num_short_strain + 1                                                
      if( stresses ) invar_col = num_short_stress + 1  
      prin_col = invar_col + 3   ! 1st principal value
      evec_col = prin_col + 3    ! 1st direction cosines
c                                                                               
      do strpt = 1, nstrpt                                                      
c                                                                               
c                       compute the principal invariants of desired             
c                       stress or strain.                                       
c                                                                               
       call get_invariants( span, mxvl, elestr(1,1,strpt), 
     &                      elestr(1,invar_col,strpt), stresses )                                                  
c                                                                               
c                       compute the principal stresses or               
c                       principal strains and the direction                     
c                       cosines of their corresponding normals.     
c                       descending order            
c                                                                               
       call get_princ_vals( span, mxvl, elestr(1,1,strpt), 
     &                      elestr(1,prin_col,strpt),           
     &                      elestr(1,evec_col,strpt), stresses )                           
c                                                                               
      end do                                                                    
c                                                                                                                                                              
      return     
c
      contains   !  just to hide/protect the routine names
c     ========
                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_mises                    *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 1/20/2017 rhd              *          
c     *                                                              *          
c     *     mises equiv. stress                                      *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine get_mises( span, gpstr, mises, mxvl )  
      use constants                              
      implicit none                                                             
c                                                                               
      integer :: span, mxvl                                                     
      double precision :: gpstr(mxvl,*), mises(*)                                  
c                                                                               
      integer :: i                                                              
c                                                                               
c                       compute the von-mises stress.                           
c                                                                               
      do i = 1, span                                                            
        mises(i) = sqrt( (gpstr(i,1)-gpstr(i,2))**2+                              
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
c     *                      subroutine get_equiv_strain             *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/21/20 rhd               *          
c     *                                                              *          
c     *     computes effective strain measure  for block             *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine get_equiv_strain( span, str, efeps, mxvl ) 
c
      use constants
c                              
      implicit none
c
      integer :: span, mxvl
      double precision :: str(mxvl,*), efeps(*)  
c
      integer :: i     
      double precision :: exx, eyy, ezz, gamxy, gamyz, gamxz                       
c                                                                               
c                       compute the effective strain measure.                   
c
      do i = 1, span
       exx = str(i,1)                                               
       eyy = str(i,2)                                               
       ezz = str(i,3)                                               
       gamxy = str(i,4)                                          
       gamyz = str(i,5)                                          
       gamxz = str(i,6)                                           
       efeps(i) = (root2/three) *
     &      sqrt( (exx-eyy)**2 + (exx-ezz)**2 + (eyy-ezz)**2 +                                       
     &            onept5*( gamxy**2 + gamyz**2 + gamxz**2 ) )
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_invariants               *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/21/20 rhd               *          
c     *                                                              *          
c     *     invariants of the stress or strain at a strain point     *
c     *                                                              *          
c     ****************************************************************          
c                                                                                                                                                             
      subroutine get_invariants( span, mxvl, str, inv, stresses ) 
c
      use constants
c                              
      implicit none                                                    
c
      integer :: span, mxvl                                                    
      double precision :: str(mxvl,*), inv(mxvl,*)                                             
      logical :: stresses    ! true if doing stresses     
c
c                    locals                                          
c         
      integer :: i                                                                      
      double precision :: xx, yy, zz, xy, yz, xz, t1, t2, t3, t4, t5,
     &                    factor
c                                                                               
      factor = half
      if( stresses ) factor = one
c
      do i = 1, span                                                            
        xx = str(i,1)                                               
        yy = str(i,2)                                               
        zz = str(i,3)                                               
        xy = str(i,4) * factor                                          
        yz = str(i,5) * factor                                         
        xz = str(i,6) * factor
        inv(i,1) = xx + yy + zz
        inv(i,2) = -xy**2 - yz**2 - xz**2 + xx*yy + yy*zz + zz*xx   
        t1 = xx*yy*zz
        t2 = two * xy*yz*xz
        t3 = xy*xy*zz
        t4 = yz*yz*xx     
        t5 = xz*xz*yy
        inv(i,3) = t1 + t2 - t3 - t4 - t5
      end do   
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_princ_vals               *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     computes the principal values and direction cosines      *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine get_princ_vals( span, mxvl, str, prstr, angles, 
     &                           stresses ) 
c  
      use constants
c                  
      implicit none                                                    
c
      integer :: span, mxvl
      logical :: stresses                                                   
      double precision :: str(mxvl,*), prstr(mxvl,*), angles(mxvl,*)                      
c                                                                               
c                    locals
c                
      integer :: i
      double precision :: factor, values(6), evalues(3), evectors(3,3)                                                               
c
      factor = half                                                             
      if( stresses ) factor = one
c
      do i = 1, span
        values(1) =  str(i,1)                                       
        values(2) =  str(i,2)                                       
        values(3) =  str(i,3)      
        values(4) =  str(i,4) * factor                                       
        values(5) =  str(i,5) * factor                                       
        values(6) =  str(i,6) * factor        
        call principal_values_ext1( values,  evectors, evalues )                               
        prstr(i,1:3)= evalues                                                
        angles(i,1:3) = evectors(1:3,1)
        angles(i,4:6) = evectors(1:3,2)
        angles(i,7:9) = evectors(1:3,3)
      end do
c
      return                                                                    
      end                                                                       
c                                                                                
      end subroutine ouext1  


c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine principal_values_ext          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                last modified : 12/18/2020 rhd                *          
c     *                                                              *          
c     *  eigenvalues of symmetric 3x3 returned in descending order   *
c     *                                                              *          
c     ****************************************************************          
c
      subroutine principal_values_ext1( vector, evectors, evals)
c
      use global_data, only : out
      use constants
      implicit none
c
c              calculates the eigenvalues and normalized eigenvectors 
c              of a symmetric 3x3 matrix using the jacobi algorithm. 
c              not the fastest method but the most stable and accurate
c              especially for ill-conditioned input.
c
c              input vector not altered
c
c              vector: xx, yy, zz, xy, yz, xz
c
c              reference:  Numerical diagonalization of 3x3 matrcies
c                          Copyright (C) 2006  Joachim Kopp
c                          GNU Lesser General Public License
c                          his name for routine: dsyevj3
c
      double precision :: vector(6), evectors(3,3), evals(3)
c
      integer :: i, x, y, r
      double precision :: q(3,3), a(3,3), w(3), sd, so, s, c, t, g, h,
     &                    z, theta, thresh
c
      a(1,1) = vector(1)
      a(1,2) = vector(4)
      a(1,3) = vector(6)
      a(2,2) = vector(2)
      a(2,3) = vector(5)
      a(3,3) = vector(3)
c
      do x = 1, 3
        q(x,x) = one
        do y = 1, x-1
          q(x, y) = zero
          q(y, x) = zero
        end do
      end do
c
      do x = 1, 3
        w(x) = a(x,x)
      end do
c
      sd = zero
      do x = 1, 3
        sd = sd + abs(w(x))
      end do
      sd = sd**2
c 
c              main iteration loop
c
      do i = 1, 50  ! max iterations
        so = zero
        do x = 1, 3
          do y = x+1, 3
            so = so + abs(a(x,y))
          end do
        end do
        if( so == zero ) then ! converged
           call principal_values_reorder_ext1
           return
        end if
c
        if( i < 4) then ! looser initial threshold
          thresh = 0.2d0 * so / 3**2
        else
          thresh = zero
        end if
c
c                do  a sweep
c
        do x = 1, 3
          do y = x+1, 3
            g = hundred * ( abs(a(x,y)) )
            if( i > 4 .and. abs(w(x)) + g == abs(w(x))
     &           .and. abs(w(y)) + g == abs(w(y)) ) then
              a(x, y) = zero
            else if( abs(a(x, y)) > thresh ) then
c                calculate jacobi transformation
              h = w(y) - w(x)
              if( abs(h) + g == abs(h) ) then
                t = a(x, y) / h
              else
                theta = half * h / a(x, y)
                if( theta < zero ) then
                  t = -one / (sqrt(one + theta**2) - theta)
                else
                  t = one / (sqrt(one + theta**2) + theta)
                end if
              end if
c
              c = one / sqrt( one + t**2 )
              s = t * c
              z = t * a(x, y)
c              
c                apply jacobi transformation
c
              a(x, y) = zero
              w(x)    = w(x) - z
              w(y)    = w(y) + z
              do r = 1, x-1
                t       = a(r, x)
                a(r, x) = c * t - s * a(r, y)
                a(r, y) = s * t + c * a(r, y)
              end do
              do r = x+1, y-1
                t       = a(x, r)
                a(x, r) = c * t - s * a(r, y)
                a(r, y) = s * t + c * a(r, y)
              end do
              do r = y+1, 3
                t       = a(x, r)
                a(x, r) = c * t - s * a(y, r)
                a(y, r) = s * t + c * a(y, r)
              end do
c
c             update eigenvectors
c
              do r = 1, 3
                t       = q(r, x)
                q(r, x) = c * t - s * q(r, y)
                q(r, y) = s * t + c * q(r, y)
              end do
            end if
          end do ! on y
        end do ! on x
      end do ! on i

      write(out,9000)
      call die_abort
 9000 format(/1x,'>>>>> FATAL ERROR: routine principal_values.',
     &          ' Jacobi iterations failed to converge',/
     & /,'.....  analysis terminated  .....',// )
c
      contains     
c     ========       
c  
      subroutine principal_values_reorder_ext1
      implicit none
c
c              put eigenvalues in descending order.
c              make eigenvector columns match ordering.
c      
      integer :: k, i
c
      do i = 1, 3 
        k = maxloc( w,dim=1 )
        evals(i) = w(k)
        evectors(1:3,i) = q(1:3,k)
        w(k) = -1.0d20
      end do
c      
      return     
      end subroutine principal_values_reorder_ext1
c      
      end subroutine principal_values_ext1
   
                                                                              
