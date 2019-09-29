c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine energy                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/28/2019 rhd              *          
c     *                                                              *          
c     *     following convergence of a load step. drive the          *          
c     *     computation of the current kinetic energy and work of    *          
c     *     the applied loads. the data values are appended to a     *          
c     *     'energy' file for the problem.                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine energy ( step, adapt_result, mdiag )                           
      use global_data ! old common.main
      implicit none   
c
      integer :: step, adapt_result
      double precision :: mdiag(*)
                                            
      integer :: kout
      double precision :: k_energy, total_int_energy,
     &                    total_pls_energy                   
      logical :: there                                                             
c                                                                               
c           open the energy file in append mode. compute and                    
c           write one record. on step 1, delete the energy file                 
c           if it exists.                                                       
c                                                                               
c                                                                               
      inquire (file = 'energy', exist = there)                                  
      if( there .and. (step .eq.1) )  then                                      
         open( newunit=kout, file='energy', form='formatted',                      
     &         status='old', access='sequential' )                               
         close(unit=kout,status='delete')                                       
      end if                                                                    
c                                                                               
      if( step .eq. 1 .or. .not. there ) then                                  
        open( newunit=kout, file='energy', form='formatted',                       
     &        status='new', access='sequential',                                
     &        position='append' )                                               
      else                                                                      
        open( newunit=kout, file='energy', form='formatted',                       
     &        status='old', access='sequential',                                
     &        position='append' )                                               
      end if                                                                    
c                                                                               
c           calculate kinetic energy and external work                          
c                                                                               
      if( step .eq. 1 .or. .not. there ) then                                  
        write(kout,*) ' '                                                       
        write(kout,*) ' '                                                       
        write(kout,8900)                                                        
        write(kout,*) ' '                                                       
        write(kout,*) ' '                                                       
        write(kout,9000)                                                        
      end if                                                                    
c                                                                               
      call kinetic_energy( mdiag, v, nodof, k_energy )                          
c                                                                               
c           only output results if adaptive solution is at the end of           
c           the given step                                                      
c                                                                               
      if( adapt_result .ne. 2 ) then                                           
          total_int_energy = internal_energy + killed_ele_int_work              
          total_pls_energy = plastic_work +  killed_ele_pls_work                
          write(kout,1000)                                                      
     &               step, k_energy,                                            
     &               total_int_energy,                                          
     &               total_int_energy + k_energy,                               
     &               total_pls_energy,                                          
     &               killed_ele_int_work, killed_ele_pls_work                   
      end if                                                                    
      close( kout )                                                            
c                                                                               
      return                                                                    
c                                                                               
 1000 format( 3x, i4, 6(3x,e16.8) )                                             
 8900 format(2x,                                                                
     &   'internal: total of all internal work composed of'                     
     & /,12x,'elastic + plastic work of all active elements',                   
     & /,12x,'work of killed elements (cells, cohesive)',//,                    
     & 2x,                                                                      
     & 'plastic:  total of all plastic contributions to internal work',         
     & /,12x,'plastic work of all active elements',                             
     & /,12x,'plastic work of killed elements (cells, cohesive)',//,
     & 2x,'Note: internal energy is not correct for *large* elastic ',
     &    'deformations computed with',
     & /,8x,'the hypoelastic material models (bilinear, ',
     &      'mises, Gurson, ..).',
     & /,8x,'WARP3D uses the Green-Nagdhi stress-rate. Abaqus uses ',
     &       'Jaumann rate. Both codes compute',
     &  /,8x,'incorrect internal work for large elastic strains ',
     & 'with hypoelastic material models.'
     & /8x,' ',
     & 'w/ v. small elastic strains.')            
 9000 format(                                                                   
     &'   step       kinetic          internal',                                
     &'         internal+kinetic      plastic work',                            
     &'   killed element work',                                                 
     &'  killed element plastic work')                                          
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine kinetic_energy                      *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 9/28/2019 rhd              *          
c     *                                                              *          
c     *     compute current kinetic energy for structure using       *          
c     *     lumped mass and nodal velocities                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine kinetic_energy( mdiag, v, nodof, k_energy )  
      implicit none
c  
      integer :: nodof, i                
      double precision :: mdiag(*), v(*), k_energy, zero, half                                   
      data zero, half / 0.0d0, 0.5d0 /                                              
c                                                                               
c           start energy calculations at zero for each time step                
c           because kinetic energy is for point in time, not history            
c                                                                               
      k_energy = zero                                                           
      do i = 1, nodof                                                           
         k_energy = k_energy + half*mdiag(i)*v(i)**2                            
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
