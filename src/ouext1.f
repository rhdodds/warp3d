c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouext1                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/17/2025 rhd              *
c     *                                                              *          
c     *     computes derived stress/strain values                    *          
c     *     at gauss points or node points for later output.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouext1( span, do_stresses, nnode, ngp, nodpts,                      
     &                   op_code, mxvl )  
c  
      use global_data, only : out                                    
      use elblk_data, only : elestr 
      use output_value_indexes
c                                    
      implicit none
c
      integer :: span, nnode, ngp, op_code, mxvl, Lode_col                                                   
      logical :: do_stresses, nodpts          
c
      integer :: nstrpt, strpt, invar_col_strains, prin_col, evec_col,
     &           mises_col, equiv_eps_col, prin_col_strains, 
     &           prin_col_stresses
      logical :: do_strains                                       
c                              
c                       elestr is size maxvl x max possible output 
c                       values allowed in code
c                                                 
      do_strains = .not. do_stresses                                                                            
      nstrpt = ngp                                                              
      if( nodpts ) nstrpt = nnode    
      mises_col     = index_sig_vm
      equiv_eps_col = index_eps_eff
c                                                                               
c                       compute additional data for output at                   
c                       each of the element nodes or gauss points               
c                       of the element block.                                   
c                                                                               
c                       op_code = 1 : compute mises equiv. stress or            
c                                     equiv. strain. both from
c                                     extrapolated/averaged components                          
c                       op_code = 2 : compute strain invariants, Lode 
c                                     angle parms, principal             
c                                     values.                                   
c                                                                               
c                                                    
      if( op_code == 1 .and.  do_stresses ) then
           do strpt = 1, nstrpt                                                    
             call get_mises( span, elestr(1,1,strpt), 
     &                       elestr(1,mises_col,strpt), mxvl )    
           end do
           return 
      end if 
c                                                            
      if( op_code == 1  .and.  do_strains ) then
           do strpt = 1, nstrpt                                                    
             call get_equiv_strain( span, elestr(1,1,strpt), 
     &                           elestr(1,equiv_eps_col,strpt), mxvl )                                                 
           end do  
           return                                                                 
      end if   
c      
      if( op_code .ne. 2 ) then
        write(out,9995); write(out,9996)
        call die_abort
      end if  
c                                                                               
c                       compute strain invariants, stress values for 
c                       Lode angles, and principal values for strains
c                       or stresses.                   
c                                                                               
      invar_col_strains = index_eps_I1
      if( do_stresses ) Lode_col = index_sig_Lode_xi
c      
      prin_col_strains  = index_eps_1   ! 1st principal value
      prin_col_stresses = index_sig_1   ! 1st principal value
      if( do_strains )  prin_col = prin_col_strains
      if( do_stresses ) prin_col = prin_col_stresses
      evec_col = prin_col + 3    ! 1st direction cosines
c                                                                               
c                       compute the principal stresses or               
c                       principal strains and the direction                     
c                       cosines of their corresponding normals.     
c                       descending order. lower level threaded over
c                       span since the eigenvalue problem at each point
c                       can be done in parallel           
c                                                                               
      do strpt = 1, nstrpt     
       call get_princ_vals_ext( span, mxvl, elestr(1,1,strpt), 
     &                      elestr(1,prin_col,strpt),           
     &                      elestr(1,evec_col,strpt), do_stresses )                           
      end do
c                                                                               
c                       for strains, compute the 3 invariants I1, I2, I3.
c                       for stresses, comoute the mean stress, stress 
c                       triaxiality and Lode angle parameter. With 
c                       values from above, the strain invariants and Lode
c                       parameters have simpler forms
c                                                                               
c     
      if( do_strains ) then                                                                          
        do strpt = 1, nstrpt    
          call get_strain_invariants( span, mxvl, elestr(1,1,strpt), 
     &              elestr(1,invar_col_strains,strpt) )  
        end do
        return
      end if 
c
      if( do_stresses ) then                                                                          
        do strpt = 1, nstrpt    
          call  get_Lode_values( span, mxvl, elestr(1,1,strpt),
     &                           elestr(1,index_sig_mean,strpt) ) 
        end do   
        return 
      end if                                                                
c                                                                                                                                                              
      return 
c      
 9995 format(/1x,'>>>>> Fatal Error: routine ouext1. Invalid op_code')
 9996 format(/1x,'>>>>>              Job terminated    .....',//)
          
c
      end
c
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_mises                    *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 1/20/2017 rhd              *          
c     *                                                              *          
c     *                mises equiv. stress all elements in block     *
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
!$omp simd
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
c     *   computes effective strain measure. all elements in block   *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine get_equiv_strain( span, str, efeps, mxvl ) 
c
      use constants
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
!$omp simd
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
c     *                      subroutine get_Lode_values              *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 9/16/2025 rhd              *          
c     *                                                              *          
c     *     mean stress, triaxiality, Lode parameter (-1<=L<=1)      *
c     *                                                              *          
c     ****************************************************************          
c                                                                                                                                                             
      subroutine get_Lode_values( span, mxvl, stresses, values )
c
      use constants 
      use output_value_indexes                             
      implicit none                                                    
c
      integer :: span, mxvl                                                    
      double precision :: stresses(mxvl,*), values(mxvl,*)                                             
c
c                    locals                                          
c         
      integer :: i                                                                      
      double precision :: sig1, sig2, sig3, mean, mises, triaxiality,
     &                    t1, t2, t3, r, xi
      double precision, parameter :: mises_zero_tol = 1.0d-10
c
!$omp simd
      do i = 1, span                                                            
        sig1 = stresses(i,index_sig_1)  ! principal values 1>2>3
        sig2 = stresses(i,index_sig_2)
        sig3 = stresses(i,index_sig_3)
        mean = ( sig1 + sig2 + sig3 ) * third
        mises = max( stresses(i,index_sig_vm), mises_zero_tol )
        triaxiality = mean / mises
        t1 = sig1 - mean
        t2 = sig2 - mean
        t3 = sig3 - mean
        r = thirteenpt5 * t1 * t2 * t3  ! except for cube root
        xi = r / mises**three
        values(i,1) = mean  ! mean
        values(i,2) = triaxiality  ! triaxiality
        values(i,3) = xi ! Lode parameter -1 <= xi <= 1
      end do   
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_strain_invariants        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 8/7/2025 rhd               *          
c     *                                                              *          
c     *       invariants of strain all elements in block             *
c     *                                                              *          
c     ****************************************************************          
c                                                                                                                                                             
      subroutine get_strain_invariants( span, mxvl, strains, inv ) 
c
      use constants
      implicit none                                                    
c
      integer :: span, mxvl                                                    
      double precision :: strains(mxvl,*), inv(mxvl,*)                                             
c
c                    locals                                          
c         
      integer :: i                                                                      
      double precision :: eps1, eps2, eps3
c                                                                               
!$omp simd
      do i = 1, span        
        eps1 = strains(i,11)
        eps2 = strains(i,12)
        eps3 = strains(i,13)
        inv(i,1) = eps1 + eps2 + eps3
        inv(i,2) = eps1*eps2 + eps2*eps3 + eps3*eps1
        inv(i,3) = eps1 * eps2 * eps3
      end do   
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine get_princ_vals_ext               *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 8/7/25 rhd                 *          
c     *                                                              *          
c     *     computes the principal values and direction cosines      *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine get_princ_vals_ext( span, mxvl, str, prstr, angles, 
     &                               do_stresses ) 
c  
      use constants
      use omp_lib
c                  
      implicit none                                                    
c
      integer :: span, mxvl
      logical :: do_stresses                                                   
      double precision :: str(mxvl,*), prstr(mxvl,*), angles(mxvl,*)                      
c                                                                               
c                    locals
c                
      integer :: i
      double precision :: factor, values(6,mxvl), evalues(3,mxvl), 
     &                    evectors(9,mxvl)
      logical :: in_parallel 
c
c              we can run the loop to call eigenvalue routine with threads
c              unless this routine is called froma loop already running 
c              threads      
c
      in_parallel = omp_in_parallel()
      factor = half                                                             
      if( do_stresses ) factor = one
c      
!$omp simd
      do i = 1, span
        values(1,i) =  str(i,1)                                       
        values(2,i) =  str(i,2)                                       
        values(3,i) =  str(i,3)      
        values(4,i) =  str(i,4) * factor                                       
        values(5,i) =  str(i,5) * factor                                       
        values(6,i) =  str(i,6) * factor    
      end do             
c
      if( in_parallel ) then  ! don't use threads
        do i = 1, span
          call principal_values_ext( values(1,i), evectors(1,i),
     &                              evalues(1,i) )   
         end do
      end if
c
      if( .not. in_parallel ) then  ! ok to run thread parallel
c$OMP PARALLEL DO  SCHEDULE( static, 64 ) PRIVATE( i )
        do i = 1, span
          call principal_values_ext( values(1,i), evectors(1,i),
     &                               evalues(1,i) )   
        end do
c$OMP END PARALLEL DO    
      end if  
c                                          
!$omp simd
      do i = 1, span
        prstr(i,1:3)  = evalues(1:3,i)                                                
        angles(i,1:9) = evectors(1:9,i)
      end do
c
      return                                                                    
      end                                                                       
c                                                                                


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
      subroutine principal_values_ext( vector, evectors, evals)
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
           call principal_values_reorder_ext
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
      subroutine principal_values_reorder_ext
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
      end subroutine principal_values_reorder_ext
c      
      end subroutine principal_values_ext
   
                                                                              
