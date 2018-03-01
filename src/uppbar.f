c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine uppbar                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 04/19/2015 rhd             *          
c     *                                                              *          
c     *        set the starting effective load vector for step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine uppbar( pbar, mdiag, accel, veloc, load, dload,                
     &                   nbeta, dt, adapt_load_fact, nodof, out )               
c                                                                               
      use stiffness_data, only : total_lagrange_forces                          
      implicit none                                                             
c                                                                               
c          parameters                                                           
c                                                                               
      integer :: nodof, out                                                     
      double precision ::                                                       
     &  pbar(nodof), mdiag(nodof), accel(nodof), veloc(nodof),                  
     &  load(nodof), dload(nodof),                                              
     &  nbeta, dt, adapt_load_fact                                              
c                                                                               
c          locals                                                               
c                                                                               
      integer :: i                                                              
      double precision ::                                                       
     & nfac1, nfac2, one, two, adfact, force_lag, zero, inertia                 
      logical :: have_mpcs                                                      
      logical, parameter :: local_debug = .false.                               
      data zero, one, two / 0.0d00, 1.0d00, 2.0d00 /                            
c                                                                               
c          compute first part of estimated total applied forces                 
c          on nodes at n+1. Section 1.6.2 of manual (Feb 2015).                 
c                                                                               
c          pbar = P_{n+1}^d in Section 1.6                                      
c                                                                               
c          non-zero imposed displacements, multi-point constraints,             
c          forces from temperatures are treated subsequently.                   
c                                                                               
c          this routine is called once at the start of a load (time)            
c          step (before iteration "0") and again at the start                   
c          of each adaptive substep when adaptive treatment of the              
c          current time step is in effect.                                      
c                                                                               
c          static:                                                              
c                                                                               
c          p-effective is the set of user                                       
c          defined loads applied to nodes at n+1. includes                      
c          increment of directly applied nodal forces, increment                
c          of equivalent nodal forces from distributed loads/                   
c                                                                               
c          the adaptive                                                         
c          solution procedure adjusts the increment of load over                
c          the step in each increment as required. the adaptive                 
c          scheme uses a scale from 0.0->1.0 over the step                      
c          to subdivide total load into smaller increments.                     
c                                                                               
c          When adaptive is NOT in effect, adapt_load_fact = 1.0.               
c          load here already has the full dload included                        
c          (see modify_load called from stpdrv).                                
c                                                                               
c          This seemingly unusual implementation does make                      
c          adaptive processing simpler.                                         
c                                                                               
c          total_lagr... nodal forces at time n that enforce                    
c          all MPCs. these are self-equilibrating.                              
c                                                                               
c          dynamic analysis:                                                    
c                                                                               
c          include the mass times the                                           
c          estimated acceleraion at n+1. the actual acceleration                
c          at n+1 as modified by displacement increment over the                
c          step is accounted for in the upres routine.                          
c                                                                               
c                                                                               
c          newmark beta method factors. accel and veloc are                     
c          at time n - start of step.                                           
c                                                                               
      if( local_debug ) write(out,9010)                                         
                                                                                
      nfac1 = (one-two*nbeta)/(two*nbeta)                                       
      nfac2 = one/(nbeta*dt)                                                    
c                                                                               
      adfact =  adapt_load_fact - one                                           
      have_mpcs = allocated( total_lagrange_forces )                            
c                                                                               
      if( have_mpcs ) then                                                      
        do i = 1, nodof                                                         
          pbar(i) =                                                             
     &       mdiag(i) * ( nfac1*accel(i) + nfac2*veloc(i) ) +                   
     &       load(i) + ! includes full dload for step                           
     &       adfact * dload(i) ! adfact=1 for no adaptive                       
     &       + total_lagrange_forces(i)                                         
        end do                                                                  
      else                                                                      
        do i = 1, nodof                                                         
          pbar(i) =                                                             
     &       mdiag(i) * ( nfac1*accel(i) + nfac2*veloc(i) ) +                   
     &       load(i) + ! includes full dload for step                           
     &       adfact * dload(i) ! adfact=1 for no adaptive                       
        end do                                                                  
      end if                                                                    
c                                                                               
      if( .not. local_debug ) return                                            
c                                                                               
      write(out,9020) nodof, nfac1, nfac2, adapt_load_fact                      
      write(out,9030)                                                           
      do i = 1, nodof                                                           
        force_lag = zero                                                        
        inertia   = mdiag(i) * ( nfac1*accel(i) + nfac2*veloc(i) )              
        if( have_mpcs ) force_lag = total_lagrange_forces(i)                    
        write(out,9040) i, pbar(i), inertia, load(i),                           
     &                  dload(i), force_lag                                     
      end do                                                                    
c                                                                               
      write(out,9015)                                                           
c                                                                               
      return                                                                    
c                                                                               
 9010 format(/,1x,"--- entering uppbar")                                        
 9015 format(1x,"--- leaving uppbar ---" )                                      
 9020 format(5x,"... nodof, nfac1, nfac2, adapt_load_fact: ",i8,                
     &       3e14.6)                                                            
 9030 format(t5,"sdof", t16, "pbar @ n", t32,"inertia @ n", t49,                
     &  "load @ n",                                                             
     &  t66, "dload",                                                           
     &  t78, "lagrange @ n" )                                                   
 9040 format(2x,i7,5f16.6)                                                      
c                                                                               
      end                                                                       
