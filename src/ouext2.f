c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouext2                       *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *     compute derived stress/strain values                     *          
c     *     after the primary values have been averaged              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2( values, nrowd, num_elems_nodes, do_stresses )
c
      use output_value_indexes
c      
      implicit none                                                             
c                                                                               
      integer :: nrowd, num_elems_nodes
      double precision :: values(nrowd,*)    ! strains or stresses                                  
      logical :: do_stresses
c      
      logical :: do_strains      
      integer :: prin_col, evec_col                                              
c  
c              numele is number of elements in block. 
c              this routines also used to compute values
c
c              for a single node -- then numele = 1.
c              really need a better design.
c  
      do_strains = .not. do_stresses
      if( do_strains )  prin_col = index_eps_1
      if( do_stresses ) prin_col = index_sig_1
      evec_col = prin_col + 3    ! 1st direction cosines
c
      call get_princ_vals_ext( num_elems_nodes, nrowd, values, 
     &                     values(1,prin_col),
     &                     values(1,evec_col), do_stresses )
c                                                                                         
      if( do_stresses ) then      
         call ouext2_get_mises( values, nrowd, num_elems_nodes )                           
         call ouext2_get_Lode_parms( values, nrowd, num_elems_nodes )
      end if 
c
      if( do_strains ) then       
         call ouext2_get_strain_invariants( values, nrowd,
     &          num_elems_nodes )                        
         call ouext2_get_equiv_strain( values, nrowd, 
     &          num_elems_nodes )                             
      end if                                                                    
c                                                                               
      return   
c
      end
c     ===                                                                 

c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine ouext2_get_mises                *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *                    mises equivalent stress                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2_get_mises( stresses, nrowd, count )  
c
      use constants
      use output_value_indexes, only : index_sig_vm
c                      
      implicit none
c
      integer :: nrowd, count
      double precision :: stresses(nrowd,count)
c
      integer :: i
      double precision :: xx, yy, zz, xy, yz, xz                                            
c                                                                               
!$omp simd
      do i = 1, count 
        xx = stresses(i,1)                                               
        yy = stresses(i,2)                                               
        zz = stresses(i,3)                                               
        xy = stresses(i,4)                                        
        yz = stresses(i,5)                                       
        xz = stresses(i,6)
        stresses(i,index_sig_vm) = sqrt( (xx-yy)**2 + (yy-zz)**2 +
     &                   (zz-xx)**2  + 
     &                   six*(xy**2 + yz**2 + xz**2) ) * iroot2                                                                              
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *           subroutine ouext2_get_equiv_strain                 *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *                 total equivalent strain measure              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2_get_equiv_strain( strains, nrowd, count )      
c
      use constants  
      use output_value_indexes, only : index_eps_eff
c                  
      implicit none  
c                                                  
      integer :: nrowd, count
      double precision :: strains(nrowd,count)
c
      integer :: i
      double precision :: exx, eyy, ezz, gamxy, gamyz, gamxz                       
c                                                                               
c                       compute the effective strain measure.                   
c
!$omp simd
      do i = 1, count
       exx = strains(i,1)                                               
       eyy = strains(i,2)                                               
       ezz = strains(i,3)                                               
       gamxy = strains(i,4)                                          
       gamyz = strains(i,5)                                          
       gamxz = strains(i,6)                                           
       strains(i,index_eps_eff) = (root2/three) *
     &      sqrt( (exx-eyy)**2 + (exx-ezz)**2 + (eyy-ezz)**2 +                                       
     &            onept5*( gamxy**2 + gamyz**2 + gamxz**2 ) )
      end do                                                                    
c                                                                                
      return                                                                    
      end                                                                       
c

c     ****************************************************************          
c     *                                                              *          
c     *               subroutine ouext2_get_strain_invariants        *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *                         strain invariants                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2_get_strain_invariants( strains, nrowd, count ) 
c
      use constants
      use output_value_indexes
c
      implicit none
c
      integer :: nrowd, count
      double precision :: strains(nrowd,count)
c
c                    locals                                          
c         
      integer :: i                                                                      
      double precision :: eps1, eps2, eps3
      
c  
!$omp simd
      do i = 1, count  
        eps1 = strains(i,index_eps_1)
        eps2 = strains(i,index_eps_2)
        eps3 = strains(i,index_eps_3)
        strains(i,index_eps_I1) = eps1 + eps2 + eps3
        strains(i,index_eps_I2) = eps1*eps2 + eps2*eps3 + eps3*eps1
        strains(i,index_eps_I3) = eps1 * eps2 * eps3
       end do                                                                   
c                                                                               
      return                                                                   
      end                                                                      
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine ouext2_get_Lode_parms                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *     mean stress, triaxiality, Lode parameter (-1<=L<=1)      *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2_get_Lode_parms( stresses, nrowd, count ) 
c
      use constants
      use output_value_indexes
c                        
      implicit none
c
      integer :: nrowd, count
      double precision :: stresses(nrowd,count)
c
      integer :: i 
      double precision :: sig1, sig2, sig3, mean, mises, triaxiality,
     &                    t1, t2, t3, r, xi
      double precision, parameter :: mises_zero_tol = 1.0d-10
c
!$omp simd
      do i = 1, count                                                            
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
        stresses(i,index_sig_mean)        = mean  ! mean
        stresses(i,index_sig_triaxiality) = triaxiality  ! triaxiality
        stresses(i,index_sig_Lode_xi)     = xi ! Lode parameter -1 <= xi <= 1
       end do                                                                   
c                                                                               
       return                                                                   
       end                                                                      
 