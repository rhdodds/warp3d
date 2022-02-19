c                                                                               
c           user_routines_support.f   Distribution version                      
c                                                                               
c           Updated:  12/21/20 rhd                                            
c                                                                               
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *         integer function warp3d_get_device_number()          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/26/12                    *          
c     *                                                              *          
c     *     find a non-connected Fortran device number               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      integer function warp3d_get_device_number()                               
      implicit none                                                             
      integer :: iunit                                                             
      logical :: connected                                                         
c                                                                               
c        1.  find an available device number, open the neutral                  
c            file  
c        2. alternatively use the newunit feature of modern fortran                                                             
c                                                                               
      warp3d_get_device_number = -1                                             
      do iunit = 11, 5000                                                       
        inquire( unit=iunit, opened=connected )                                 
        if ( connected ) cycle                                                  
        warp3d_get_device_number = iunit                                        
        return                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c *********************************************************************         
c **                                                                 **         
c **              Abaqus Compatible Support Routines                 **         
c **                                                                 **         
c *********************************************************************         
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine xit  (called by UMATs)            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 3/22/12                     *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine xit  
      use global_data, only : out                                                          
      implicit integer (a-z)                                                    
c     
      write(out,*) ' '                                                                          
      write(out,*) '>>> UMAT called to abort execution'        
      write(out,*) ' '                                                                          
c
      call die_abort                                                            
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine rotsig                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 05/31/19 vsp                   *          
c     *                                                              *          
c     *     tensor rotation routine for Abaqus compatible UMAT       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine rotsig( in_vec, drot, out_vec, type, nrow, ncol )  
      use constants            
      implicit none                                                             
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision, intent(in) :: in_vec(6), drot(3,3)
      double precision, intent(out) :: out_vec(6)
      integer, intent(in) :: type, nrow, ncol                                                  
c                                                                               
c                     locally defined arrays-variables                          
c                                                                               
      double precision :: factor,  a(3,3), t(3,3), c(3,3)                          
      logical, parameter :: local_debug = .true.                               
c                                                                               
c                                                                               
c                     out = drot * in * trans(drot)                             
c 
c                     type = 1 input tensor has stress                          
c                     type = 2 input tensor has engr strain                     
c 
c                     Abaqus umat ordering: x,y,z,xy,xz,yz                      
c 
c                     put input tensor into 3x3 matrix form.                    
c                                                                               
      factor = one                                                              
      if( type .eq. 2 ) factor = half                                           
c                                                                               
      a(1,1) = in_vec(1)                                                        
      a(2,1) = in_vec(4) * factor                                               
      a(3,1) = in_vec(5) * factor                                               
      a(1,2) = a(2,1)                                                           
      a(2,2) = in_vec(2)                                                        
      a(3,2) = in_vec(6) * factor                                               
      a(1,3) = a(3,1)                                                           
      a(2,3) = a(3,2)                                                           
      a(3,3) = in_vec(3)                                                        
c                                                                               
c                     t = a * trans(drot)                                       
c                                                                               
      t(1,1) = a(1,1)*drot(1,1) + a(1,2)*drot(1,2)  + a(1,3)*drot(1,3)          
      t(2,1) = a(2,1)*drot(1,1) + a(2,2)*drot(1,2)  + a(2,3)*drot(1,3)          
      t(3,1) = a(3,1)*drot(1,1) + a(3,2)*drot(1,2)  + a(3,3)*drot(1,3)          
c                                                                               
      t(1,2) = a(1,1)*drot(2,1) + a(1,2)*drot(2,2)  + a(1,3)*drot(2,3)          
      t(2,2) = a(2,1)*drot(2,1) + a(2,2)*drot(2,2)  + a(2,3)*drot(2,3)          
      t(3,2) = a(3,1)*drot(2,1) + a(3,2)*drot(2,2)  + a(3,3)*drot(2,3)          
c                                                                               
      t(1,3) = a(1,1)*drot(3,1) + a(1,2)*drot(3,2)  + a(1,3)*drot(3,3)          
      t(2,3) = a(2,1)*drot(3,1) + a(2,2)*drot(3,2)  + a(2,3)*drot(3,3)          
      t(3,3) = a(3,1)*drot(3,1) + a(3,2)*drot(3,2)  + a(3,3)*drot(3,3)          
c                                                                               
c                     c = drot * t                                              
c                                                                               
      c(1,1) = drot(1,1)*t(1,1) + drot(1,2)*t(2,1) + drot(1,3)*t(3,1)           
      c(2,1) = drot(2,1)*t(1,1) + drot(2,2)*t(2,1) + drot(2,3)*t(3,1)           
      c(3,1) = drot(3,1)*t(1,1) + drot(3,2)*t(2,1) + drot(3,3)*t(3,1)           
c                                                                               
      c(1,2) = drot(1,1)*t(1,2) + drot(1,2)*t(2,2) + drot(1,3)*t(3,2)           
      c(2,2) = drot(2,1)*t(1,2) + drot(2,2)*t(2,2) + drot(2,3)*t(3,2)           
      c(3,2) = drot(3,1)*t(1,2) + drot(3,2)*t(2,2) + drot(3,3)*t(3,2)           
c                                                                               
      c(1,3) = drot(1,1)*t(1,3) + drot(1,2)*t(2,3) + drot(1,3)*t(3,3)           
      c(2,3) = drot(2,1)*t(1,3) + drot(2,2)*t(2,3) + drot(2,3)*t(3,3)           
      c(3,3) = drot(3,1)*t(1,3) + drot(3,2)*t(2,3) + drot(3,3)*t(3,3)           
c                                                                               
      factor = one                                                              
      if( type .eq. 2 ) factor = two                                            
      out_vec(1) = c(1,1)                                                       
      out_vec(2) = c(2,2)                                                       
      out_vec(3) = c(3,3)                                                       
      out_vec(4) = c(2,1) * factor                                                 
      out_vec(5) = c(3,1) * factor                                                 
      out_vec(6) = c(3,2) * factor                                                 
c                                                                               
      return                                                                    
      end                                                                       
                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine sinv                         *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 05/14/12                       *          
c     *                                                              *          
c     *              invariants 1 & 2 of 3D stress tensor            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine sinv( sig, sinv1, sinv2, ndi, nshr )
      use constants                           
      implicit none                                                             
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision :: sig(*), sinv1, sinv2                                                    
      integer :: ndi, nshr                                                         
c                                                                               
c                     locally defined arrays-variables                          
c                                                                               
      double precision :: sig_dev(6), t1, t2                           
c                                                                               
      sinv1 = third * ( sig(1) + sig(2) + sig(3) )                          
c                                                                               
      sig_dev(1) = sig(1) - sinv1                                               
      sig_dev(2) = sig(2) - sinv1                                               
      sig_dev(3) = sig(3) - sinv1                                               
      sig_dev(4) = sig(4)                                                       
      sig_dev(5) = sig(5)                                                       
      sig_dev(6) = sig(6)                                                       
c                                                                               
      t1 = sig_dev(1)**2 + sig_dev(2)**2 + sig_dev(3)**2                        
      t2 = sig_dev(4)**2 + sig_dev(5)**2 + sig_dev(6)**2                        
      sinv2 = sqrt( onept5 * (t1 + two * t2 ) )                              
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine sprinc                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/20/20 rhd               *          
c     *                                                              *          
c     *     compute principal strain or stresses compatible with     *          
c     *     UMAT specifications                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine sprinc( s, ps, lstr, ndi, nshr ) 
c
      use constants
c                              
      implicit none  
      integer :: lstr, ndi, nshr, ier                                                 
      double precision :: s(6), ps(3)                                                          
c                                                                               
c                    locally allocated                                          
c                                                                               
      double precision :: temp(6), evectors(3,3), evals(3), factor                            
c                                                                               
c        calculate the principal strains or stresses for UMAT support.          
c        Abaqus stress/strain ordering: xx, yy, zz, xy, xz, yz                  
c        lstr = 1 (stresses), lstr = 2 (strains)                                
c                                                                               
      factor = one                                                              
      if( lstr .eq. 2 ) factor = half  ! strains                                         
c                                                                               
      temp(1) = s(1)                                                            
      temp(2) = s(2)                                                 
      temp(3) = s(3)                                                            
      temp(4) = s(4) * factor                                                   
      temp(5) = s(5) * factor                                                   
      temp(6) = s(6) * factor                                                              
      call principal_values( temp, evectors, evals )       
      ps = evals
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine sprind                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/14/12                    *          
c     *                                                              *          
c     *     compute principal strain or stresses and direction       *
c     *     cosines following UMAT specifications                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine sprind( s, ps, an, lstr, ndi, nshr )   
c
      use constants
c                        
      implicit none
      integer :: lstr, ndi, nshr, ier                                                    
      double precision :: s(6), ps(3), an(3,3)                                                 
c                                                                               
c                    locals
c
      double precision :: temp(6), evectors(3,3), evals(3), factor                            
c                                                                                                                                                             
c        calculate the principal strains or stresses and eigenvectors           
c        for UMAT support.                                                      
c        Abaqus stress/strain ordering: xx, yy, zz, xy, xz, yz                  
c        lstr = 1 (stresses), lstr = 2 (strains)                                
c                                                                               
      factor = one                                                              
      if( lstr .eq. 2 ) factor = half                                           
c 
      temp(1) = s(1)                                                            
      temp(2) = s(2)                                                 
      temp(3) = s(3)                                                            
      temp(4) = s(4) * factor                                                   
      temp(5) = s(5) * factor                                                   
      temp(6) = s(6) * factor                                                              
      call principal_values( temp, evectors, evals )       
      ps = evals
      an = evectors
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine getnumcpus                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/13/2017                  *          
c     *                                                              *          
c     *     return the number of ranks during MPI execution          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getnumcpus( nranks )                                           
      use global_data, only : use_mpi, numprocs
      implicit none
      integer :: nranks                                                   
c                                                                               
      nranks = 1                                                                
      if( use_mpi ) nranks = numprocs                                           
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine getrank                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/13/2017                  *          
c     *                                                              *          
c     *     return the MPI rank number for this process              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getrank( thisrank )                                            
      use global_data, only : use_mpi, myid
      implicit none
      integer :: thisrank                                                   
c                                                                               
      thisrank = 0                                                              
      if( use_mpi ) thisrank = myid                                             
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine getmodelsizes                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/13/2017                  *          
c     *                                                              *          
c     *     return number of nodes and elements                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getmodelsizes( num_model_nodes, num_model_elements )           
      use global_data,only : nonode, noelem
      implicit none
      integer :: num_model_nodes, num_model_elements                                                    
c                                                                               
      num_model_nodes    = nonode                                               
      num_model_elements = noelem                                               
c                                                                               
      return                                                                    
      end                                                                       
          

c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine principal_values               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                last modified : 2/19/22 rhd                   *          
c     *                                                              *          
c     *  eigenvalues of symmetric 3x3 returned in descending order.  *
c     *  eigenvectors stored to match.                               *                 
c     *                                                              *          
c     ****************************************************************          
c
      subroutine principal_values( vector, evectors, evals)
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
           call principal_values_reorder
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
      subroutine principal_values_reorder
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
      end subroutine principal_values_reorder
c      
      end subroutine principal_values
      

                                                                      
