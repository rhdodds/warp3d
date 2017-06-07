c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine inclmass                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : rhd 1/20/2015              *          
c     *                                                              *          
c     *     includes the element lumped mass matrices into the       *          
c     *     element stiffness matrices using newmark factors         *          
c     *                                                              *          
c     *     now includes modified for asymmetric assembly            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine inclmass                                                       
      use global_data ! old common.main
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      integer :: blk, now_thread                                                
c                                                                               
c                                                                               
c             MPI:                                                              
c               tell the slave processors to enter this routine and             
c               add the mass contribution to the element stiffnesses            
c               they own.                                                       
c             Non-MPI:                                                          
c               these are dummy routines which return immediately               
c                                                                               
      call wmpi_alert_slaves ( 6 )                                              
c                                                                               
c             update diagonal terms of element stiffness                        
c             matrices to include nodal mass scaled by                          
c             time increment and newmark beta factor.                           
c             do the update by blocks. elblks(2,blk) holds                      
c             which processor owns the block (=0 non-MPI).                      
c                                                                               
c$OMP PARALLEL DO PRIVATE( blk, now_thread )                                    
c$OMP&            SHARED( nelblk, elblks, myid )                                
      do blk = 1, nelblk                                                        
         if( elblks(2,blk) .ne. myid ) cycle                                    
         now_thread = omp_get_thread_num() + 1                                  
         call inclmass_blk( blk )                                               
      end do ! over blocks                                                      
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine inclmass_blk                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : rhd 3/16/2017 rhd         *           
c     *                                                              *          
c     *     process a block for updating element stiffness with      *          
c     *     effective mass for newmark integration                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine inclmass_blk( blk )                                            
      use global_data ! old common.main
c                                                                               
      use elem_block_data, only: mass_blocks, estiff_blocks                     
      use main_data, only: asymmetric_assembly                                  
      implicit none                                                             
c                                                                               
      integer :: blk                                                            
c                                                                               
c             locals. use pointers to simplify indexing within                  
c             innermost loop.                                                   
c                                                                               
      integer :: felem, num_enode, totdof, span,i, j, k, l                      
      double precision :: nfac, mel(mxvl,mxedof)                                
      double precision, parameter :: one = 1.0d0                                
      double precision, dimension(:,:), pointer :: emat, mmat                   
      logical :: symmetric_assembly                                             
c                                                                               
c             newmark multiplication factor.                                    
c                                                                               
      nfac = one / (nbeta*dt*dt)                                                
c                                                                               
c             update diagonal terms of element stiffness                        
c             matrices to include nodal mass scaled by                          
c             time increment and newmark beta factor.                           
c                                                                               
c             expand the compact element masses in block to length              
c             3*num_enode from num_enode to simplify updating loops.            
c             flip array so the processing loops can access element             
c             mass by down the column                                           
c                                                                               
      symmetric_assembly = .not. asymmetric_assembly                            
      felem     = elblks(1,blk)                                                 
      num_enode = iprops(2,felem)                                               
      totdof    = 3 * num_enode                                                 
      span      = elblks(0,blk)                                                 
      emat      => estiff_blocks(blk)%ptr                                       
      mmat      => mass_blocks(blk)%ptr                                         
c                                                                               
      do i = 1, num_enode                                                       
        k = i + num_enode                                                       
        l = i + num_enode + num_enode                                           
!DIR$ VECTOR ALIGNED                                                            
        do j = 1, span                                                          
           mel(j,i) =  mmat(i,j) * nfac                                         
           mel(j,k) =  mmat(i,j) * nfac                                         
           mel(j,l) =  mmat(i,j) * nfac                                         
        end do                                                                  
      end do                                                                    
c                                                                               
      if( symmetric_assembly ) then                                             
           do i = 1, totdof                                                     
            k = dcp(i)                                                          
!DIR$ VECTOR ALIGNED                                                            
            do j = 1, span                                                      
              emat(k,j) = emat(k,j) + mel(j,i)                                  
            end do                                                              
           end do                                                               
      else                                                                      
           do i = 1, totdof                                                     
            k = (i-1)*totdof + i ! location of diagonal                         
!DIR$ VECTOR ALIGNED                                                            
            do j = 1, span                                                      
               emat(k,j) = emat(k,j) + mel(j,i)                                 
            end do                                                              
           end do                                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
