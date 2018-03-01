c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine update                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 2/17/2017 rhd              *          
c     *                                                              *          
c     *     various updates of vectors required after the iterative  *          
c     *     solution procedure for a step has been completed.        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine update                                                         
      use global_data ! old common.main
c                                                                               
      use main_data, only :  temper_nodes, temper_elems,                        
     &                       dtemp_nodes, dtemp_elems, mdiag,                   
     &                       nonlocal_analysis, pbar                            
      use elem_block_data, only : history_blocks, history1_blocks,              
     &                            eps_n_blocks, eps_n1_blocks,                  
     &                            urcs_n_blocks, urcs_n1_blocks,                
     &                            history_blk_list, eps_blk_list,               
     &                            urcs_blk_list,                                
     &                            nonlocal_flags, nonlocal_data_n,              
     &                            nonlocal_data_n1                              
      use stiffness_data, only : total_lagrange_forces,                         
     &                           d_lagrange_forces                              
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
                                                                                
      double precision ::                                                       
     &    dlf, tlf                                                              
c                                                                               
      integer :: blk, felem, span, ngp, hblock_size, eblock_size,               
     &           ublock_size, k, n, i, dof                                      
      logical :: chk, update_lag_forces                                         
      logical, parameter :: local_debug = .false.                               
c                                                                               
      if( local_debug ) write(out,9300)                                         
c                                                                               
c                       tell those slaves to get in here.                       
c                                                                               
      call wmpi_alert_slaves (17)                                               
c                                                                               
c                       update displacements                                    
c                                                                               
      u(1:nodof) = u(1:nodof) + du(1:nodof)                                     
c                                                                               
c                       update the velocities and accelerations                 
c                       to n+1 using Newmark beta method, velocities            
c                       and accelerations at n and the converged                
c                       displacement increment (du) over the step.              
c                       skip if we are on worker ranks                          
c                                                                               
      if( root_processor ) call newmrk( nodof, nbeta, dt, du, v, a )            
c                                                                               
c                       compute the equivalent loads for the con-               
c                       strained dof.                                           
c                                                                               
      dof = csthed                                                              
      do while ( dof .ne. -1 )                                                  
         load(dof) = mdiag(dof)*a(dof) + ifv(dof)                               
         dof       = cstmap(dof)                                                
      end do                                                                    
c                                                                               
c                      update the structural vectors so that states             
c                      at n and n+1 are identical                               
c                        1) element histories                                   
c                        2) element strains                                     
c                        3) element stresses                                    
c                                                                               
c$OMP PARALLEL DO PRIVATE( blk, felem, span, ngp, hblock_size,                  
c$OMP&                     eblock_size, ublock_size )                           
      do blk = 1, nelblk                                                        
         if( myid .ne. elblks(2,blk) ) cycle                                    
         felem        = elblks(1,blk)                                           
         span         = elblks(0,blk)                                           
         ngp          = iprops(6,felem)                                         
         hblock_size  = span * ngp * history_blk_list(blk)                      
         eblock_size  = span * ngp * nstr                                       
         ublock_size  = span * ngp * nstrs                                      
c                                                                               
         if( hblock_size .gt. 0 )                                               
     &      call update_copy( history_blocks(blk)%ptr(1),                       
     &                  history1_blocks(blk)%ptr(1), hblock_size )              
         if( eps_blk_list(blk) .eq. 1 )                                         
     &      call update_copy( eps_n_blocks(blk)%ptr(1),                         
     &                  eps_n1_blocks(blk)%ptr(1), eblock_size)                 
         if( urcs_blk_list(blk) .eq. 1 )                                        
     &      call update_copy( urcs_n_blocks(blk)%ptr(1),                        
     &                  urcs_n1_blocks(blk)%ptr(1), ublock_size )               
      end do ! on blk                                                           
c$OMP END PARALLEL DO                                                           
c                                                                               
c                      update nonlocal shared state values                      
c                      if they exist                                            
c                                                                               
      if( nonlocal_analysis ) then                                              
       n = nonlocal_shared_state_size                                           
       do i = 1, noelem                                                         
         if( .not. nonlocal_flags(i) ) cycle                                    
         chk = allocated( nonlocal_data_n(i)%state_values ) .and.               
     &         allocated( nonlocal_data_n1(i)%state_values )                    
        if( chk ) then                                                          
!DIR$ VECTOR ALIGNED                                                            
           nonlocal_data_n(i)%state_values(1:n) =                               
     &              nonlocal_data_n1(i)%state_values(1:n)                       
        else                                                                    
           write(out,9100) i                                                    
           call die_abort                                                       
        end if                                                                  
       end do                                                                   
      end if                                                                    
c                                                                               
c                       update the nodal and element temperatures               
c                                                                               
      if( temperatures ) then                                                   
!DIR$ VECTOR ALIGNED                                                            
        temper_nodes(1:nonode) = temper_nodes(1:nonode) +                       
     &                           dtemp_nodes(1:nonode)                          
!DIR$ VECTOR ALIGNED                                                            
        temper_elems(1:noelem) = temper_elems(1:noelem) +                       
     &                           dtemp_elems(1:noelem)                          
      end if                                                                    
c                                                                               
c                       update contact geometry terms if needed                 
c                                                                               
      call updt_contact ! does not change contact forces                        
c                                                                               
c                       update lagrange nodal forces for MPCs.                  
c                       not supported on MPI                                    
c                                                                               
      update_lag_forces = allocated(total_lagrange_forces) .and.                
     &                    allocated(d_lagrange_forces)                          
      if( update_lag_forces ) then                                              
!DIR$ VECTOR ALIGNED                                                            
       total_lagrange_forces(1:nodof) = total_lagrange_forces(1:nodof)          
     &                                 + d_lagrange_forces(1:nodof)             
       if( local_debug ) then                                                   
          write(out,9205); write(out,9210)                                      
          do i = 1, nodof                                                       
            dlf = d_lagrange_forces(i)                                          
            tlf = total_lagrange_forces(i)                                      
            write(out,9220) i, tlf, tlf-dlf, dlf                                
          end do                                                                
       end if                                                                   
      end if                                                                    
c                                                                               
      if( local_debug ) write(out,9305)                                         
c                                                                               
      return                                                                    
c                                                                               
 9200 format(10x,i5,3f15.5)                                                     
 9190 format(15x,"pbar",10x,"bob_lag",10x,"dlag")                               
      return                                                                    
c                                                                               
 9000 format(3x,i5, 2e14.6)                                                     
 9100 format(">>>>> FATAL ERROR. update. nonlocal. elem: ",i8,                  
     &      /,"      Job terminated." )                                         
 9205 format(5x,"... updating Lagrange forces to n+1 ...")                      
 9210 format(t5,"sdof", t16, "tlf (new)", t31, "tlf(old)",                      
     &  t50, "dlf" )                                                            
 9220 format(2x,i7,3f16.6)                                                      
 9300 format(1x,"--- entering update ---" )                                     
 9305 format(1x,"--- leaving update ---" )                                      
                                                                                
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine update_copy                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/20/2015 rhd             *          
c     *                                                              *          
c     *     copy vector a = b. using this routine hides the blocked  *          
c     *     structures of input arrays from pointers. should get     *          
c     *     inlined                                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine update_copy( a, b, n )                                         
      implicit none                                                             
      integer :: n                                                              
      double precision :: a(n), b(n)                                            
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      a = b                                                                     
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine newmrk                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 04/20/2015 rhd             *          
c     *                                                              *          
c     *     computation of the final velocities and accelerations    *          
c     *     at n+1 from newmark beta method.                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine newmrk( nodof, nbeta, dt, du, velocity,                        
     &                   acceleration )                                         
      implicit none                                                             
c                                                                               
      integer :: nodof                                                          
      double precision ::                                                       
     & du(nodof), velocity(nodof),  acceleration(nodof), dt, nbeta              
                                                                                
      integer :: i                                                              
      double precision ::                                                       
     &     bdt, bdt2, exp1, exp2, exp3, exp4, zero, lt1, one,                   
     &     two, three, four, veln, acceln                                       
c                                                                               
      data zero, lt1, one, two, three, four / 0.0d00,                           
     &     0.999999999999999d00, 1.0d00, 2.0d00, 3.0d00, 4.0d00 /               
c                                                                               
      bdt  = nbeta * dt                                                         
      bdt2 = bdt * dt                                                           
c                                                                               
      if( nbeta .ge. lt1 ) then                                                 
c                                                                               
c                       code kluge to eliminate accelerations and               
c                       velocity using numerical damping with                   
c                       gamma = 1.5, best damaping with nbeta = 1.0             
c                                                                               
         exp1 = three / (two*bdt)                                               
         exp2 = (three-two*nbeta) / (two*nbeta)                                 
         exp3 = (three-four*nbeta) / (four*nbeta)                               
         exp4 = (one-two*nbeta) / (two*nbeta)                                   
c                                                                               
      else                                                                      
c                                                                               
c                       perform the computations based on newmark's             
c                       beta method.                                            
c                                                                               
         exp1 = one / (two*bdt)                                                 
         exp2 = (one-two*nbeta) / (two*nbeta)                                   
         exp3 = (one-four*nbeta) / (four*nbeta)                                 
         exp4 = (one-two*nbeta) / (two*nbeta)                                   
c                                                                               
c                                                                               
      end if                                                                    
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, nodof                                                           
         veln            = velocity(i)                                          
         acceln          = acceleration(i)                                      
         velocity(i)     = du(i)*exp1 - veln*exp2 - acceln*exp3*dt              
         acceleration(i) = du(i)/bdt2 - veln/bdt  - acceln*exp4                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
