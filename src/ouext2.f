c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouext2                       *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/21/20 rhd               *          
c     *                                                              *          
c     *     this subroutine computes derived stress/strain values    *          
c     *     after the primary values have been averaged              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2( str, nrowd, numele, do_stresses )                  
      implicit none                                                             
c                                                                               
      integer :: nrowd, numele                                                   
      double precision :: str(nrowd,*)                                      
      logical :: do_stresses                                                    
c                                                                               
      if( do_stresses ) then                                                    
         call get_stress_invariants( str, nrowd, numele )                        
         call get_princ_stresses( str, nrowd, numele )                            
         call get_mises( str, nrowd, numele )                           
      else                                                                      
         call get_strain_invariants( str, nrowd, numele )                        
         call get_princ_strains( str, nrowd, numele )                            
         call get_equiv_strain( str, nrowd, numele )                             
      end if                                                                    
c                                                                               
      return   
c
      contains   ! just to hide these routines
c     ========                                                                 
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine get_strain_invariants               *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     strain invariants                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_strain_invariants( str, nrowd, numele ) 
c
      use constants
c
      implicit none
c
      integer :: nrowd, numele
      double precision :: str(nrowd,*)
c
c                    locals                                          
c         
      integer :: i                                                                      
      double precision :: xx, yy, zz, xy, yz, xz, t1, t2, t3, t4, t5
c  
      do i = 1, numele                                                                             
       xx = str(i,1)                                               
       yy = str(i,2)                                               
       zz = str(i,3)                                               
       xy = str(i,4) * half
       yz = str(i,5) * half
       xz = str(i,6) * half
       str(i,8) = xx + yy + zz
       str(i,9) = -xy**2 - yz**2 - xz**2 + xx*yy + yy*zz + zz*xx                   
       t1 = xx*yy*zz
       t2 = two * xy*yz*xz
       t3 = xy*xy*zz
       t4 = yz*yz*xx     
       t5 = xz*xz*yy
       str(i,10) =  t1 + t2 - t3 - t4 - t5                                                       
       end do                                                                   
c                                                                               
       return                                                                   
       end                                                                      
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine get_stress_invariants            *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     stress invariants                                        *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_stress_invariants( str, nrowd, numele ) 
c
      use constants
c                        
      implicit none
c
      integer :: nrowd, numele                                                    
      double precision :: str(nrowd,*)
c
      integer :: i                                                      
      double precision :: xx, yy, zz, xy, yz, xz, t1, t2, t3, t4, t5
c                                                                                
       do i = 1, numele                
        xx = str(i,1)                                               
        yy = str(i,2)                                               
        zz = str(i,3)                                               
        xy = str(i,4)                                        
        yz = str(i,5)                                       
        xz = str(i,6)
        str(i,12) = xx + yy + zz
        str(i,13) = -xy**2 - yz**2 - xz**2 + xx*yy + yy*zz + zz*xx
        t1 = xx*yy*zz
        t2 = two * xy*yz*xz
        t3 = xy*xy*zz
        t4 = yz*yz*xx     
        t5 = xz*xz*yy
        str(i,14) = t1 + t2 - t3 - t4 - t5
       end do                                                                   
c                                                                               
       return                                                                   
       end                                                                      
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_princ_strains            *         
c     *                        (descending order)                    *
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     principal strains                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_princ_strains( str, nrowd, numele ) 
c
      use constants
c                           
      implicit none
c
      integer :: nrowd, numele 
      double precision :: str(nrowd,*)                                                     
c                                                                               
c                    locals
c   
      integer :: i                                                                            
      double precision :: values(6), evalues(3), evectors(3,3)
c                                                                                                                                                              
      do i = 1, numele
        values(1) = str(i,1)                                       
        values(2) = str(i,2)                                       
        values(3) = str(i,3)      
        values(4) = str(i,4) * half
        values(5) = str(i,5) * half
        values(6) = str(i,6) * half
        call principal_values_ext2( values,  evectors, evalues )                               
        str(i,11) = evalues(1)                                                 
        str(i,12) = evalues(2)                                                 
        str(i,13) = evalues(3)                                                 
        str(i,14) = evectors(1,1)                                             
        str(i,15) = evectors(2,1)                                             
        str(i,16) = evectors(3,1)                                             
        str(i,17) = evectors(1,2)                                             
        str(i,18) = evectors(2,2)                                             
        str(i,19) = evectors(3,2)                                             
        str(i,20) = evectors(1,3)                                             
        str(i,21) = evectors(2,3)                                             
        str(i,22) = evectors(3,3)                                             
      end do     
c                                                               
      return                                                                    
      end                                                                       

c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_princ_stresses           *          
c     *                        (descending order)                    *
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     computes principal stress values at the nodes/elem       *
c     *     after the primary values have been averaged              *                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_princ_stresses( str, nrowd, numele ) 
      implicit none
c
      integer :: nrowd, numele                                                 
      double precision :: str(nrowd,*)                                                     
c                                                                               
c                    locally allocated                                          
c         
      integer :: i                                                                      
      double precision :: temp_stress(6), ev(3), evec(3,3)   
c                                                                               
      do i = 1, numele                                                            
          temp_stress(1:6) = str(i,1:6)                                         
          call principal_values_ext2( temp_stress, evec, ev )         
          str(i,15) = ev(1)                                                 
          str(i,16) = ev(2)                                                 
          str(i,17) = ev(3)    
          str(i,18) = evec(1,1)                                             
          str(i,19) = evec(2,1)                                             
          str(i,20) = evec(3,1)                                             
          str(i,21) = evec(1,2)                                             
          str(i,22) = evec(2,2)                                             
          str(i,23) = evec(3,2)                                             
          str(i,24) = evec(1,3)                                             
          str(i,25) = evec(2,3)                                             
          str(i,26) = evec(3,3)                                             
      end do    
c                                                                
      return                                                                    
      end                                                                       

c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine get_mises                       *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     mises equivalent stress                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_mises( str, nrowd, numele )  
c
      use constants
c                      
      implicit none
c
      integer :: nrowd, numele
      double precision :: str(nrowd,*)
c
      integer :: i
      double precision :: xx, yy, zz, xy, yz, xz                                            
c                                                                               
      do i = 1, numele  
        xx = str(i,1)                                               
        yy = str(i,2)                                               
        zz = str(i,3)                                               
        xy = str(i,4)                                        
        yz = str(i,5)                                       
        xz = str(i,6)
        str(i,8) = sqrt( (xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2  + 
     &                   six*(xy**2 + yz**2 + xz**2) ) * iroot2                                                                              
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine get_equiv_strain                 *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/21/20 rhd               *          
c     *                                                              *          
c     *     equivalent strain measure                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_equiv_strain( str, nrowd, numele )      
c
      use constants  
c                  
      implicit none  
c                                                  
      integer :: nrowd, numele
      double precision :: str(nrowd,*)
c
      integer :: i
      double precision :: exx, eyy, ezz, gamxy, gamyz, gamxz                       
c                                                                               
c                       compute the effective strain measure.                   
c
      do i = 1, numele
       exx = str(i,1)                                               
       eyy = str(i,2)                                               
       ezz = str(i,3)                                               
       gamxy = str(i,4)                                          
       gamyz = str(i,5)                                          
       gamxz = str(i,6)                                           
       str(i,7) = (root2/three) *
     &      sqrt( (exx-eyy)**2 + (exx-ezz)**2 + (eyy-ezz)**2 +                                       
     &            onept5*( gamxy**2 + gamyz**2 + gamxz**2 ) )
      end do                                                                    
c                                                                                
      return                                                                    
      end                                                                       
c
      end subroutine ouext2


c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine principal_values_ext2          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                last modified : 12/18/2020 rhd                *          
c     *                                                              *          
c     *  eigenvalues of symmetric 3x3 returned in descending order   *
c     *                                                              *          
c     ****************************************************************          
c
      subroutine principal_values_ext2( vector, evectors, evals)
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
           call principal_values_reorder_ext2
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
      subroutine principal_values_reorder_ext2
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
      end subroutine principal_values_reorder_ext2
c      
      end subroutine principal_values_ext2
   
