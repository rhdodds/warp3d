c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine upres_iter_0                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/19/2015                 *          
c     *                                                              *          
c     *     residual nodal force vector for iteration "0" of step    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine upres_iter_0( norm_load, nbeta, dt, nodof, iout,               
     &                         cstmap, pbar, ifv, mdiag, du, res,               
     &                         dstmap )                                         
      use contact, only : contact_force                                         
c                                                                               
      implicit none                                                             
c                                                                               
c                    parameters                                                 
c                                                                               
      integer ::  nodof, iout, cstmap(*), dstmap(*)                             
      double precision                                                          
     & norm_load, nbeta, dt,                                                    
     & pbar(nodof), ifv(nodof), mdiag(nodof), du(nodof),                        
     & res(nodof)                                                               
c                                                                               
c                    locals                                                     
c                                                                               
      integer :: i                                                              
      double precision                                                          
     &  nfac, sum                                                               
      double precision, parameter :: zero = 0.0d00,                             
     &                               one  = 1.0d00                              
      logical, parameter :: local_debug = .false.                               
                                                                                
      integer :: chknode, chkdof                                                
c                                                                               
c                                                                               
c              set the residual load vector at iter = 0 just before             
c              solving equations.                                               
c                                                                               
c              du will have non-zero terms if using the extrapolated            
c              displacement increment method and/or the user specified          
c              imposed displacements.                                           
c                                                                               
c              pbar then includes the total "eternal" applied load              
c              at n+1 (including estimated inertia).                            
c                                                                               
c              ifv has the "internal" forces for stresses                       
c              at n plus the effects of imposed displacement                    
c              increments/temperatures for step.                                
c                                                                               
c              Section 1.6 of Manual.                                           
c                                                                               
      if( local_debug ) write(iout,9300)                                        
c                                                                               
      nfac = one / (nbeta*dt*dt)                                                
      sum = zero                                                                
c                                                                               
      if( local_debug ) write(iout,9310) dt, nbeta, nfac, nodof                 
c                                                                               
      do i = 1, nodof                                                           
        if( cstmap(i) .eq. 0 ) then                                             
           res(i) = ( pbar(i)- mdiag(i)*du(i)*nfac                              
     &              + contact_force(i) ) - ifv(i)                               
        else                                                                    
           res(i) = zero ! has absolute constraint                              
        end if                                                                  
        sum = sum + res(i) * res(i)                                             
      end do                                                                    
c                                                                               
c              if norm = 0 there is no incremental loading at                   
c              start of step to drive an increment of displacements             
c                                                                               
      norm_load = sqrt( sum )                                                   
c                                                                               
      if( .not. local_debug ) return                                            
c                                                                               
      write(iout,9000) norm_load                                                
      write(iout,9100)                                                          
      do i = 1, nodof                                                           
       write(iout,9200) i, cstmap(i),  res(i), pbar(i),                         
     &               -mdiag(i)*du(i)*nfac,                                      
     &               contact_force(i), -ifv(i)                                  
      end do                                                                    
      write(iout,9320)                                                          
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/,3x,'... norm of res vector (rhs for iter 0): ',e14.6)            
 9100 format(/,6x,"dof",2x,"cstmap",5x, "res",13x,"pbar",                       
     & 9x, "-[M]du*nfact",                                                      
     & 3x,"contact_force",6x,"-ifv" )                                           
 9200 format(3x,2i7,5e16.6)                                                     
 9300 format(/,1x,"... enter upres_iter_0 ...")                                 
 9310 format(/,3x,"... dt, nbeta, nfact, nodof: ",3e14.6,i8)                    
 9320 format(/,1x,"... leave upres_iter_0 ...",//)                              
c                                                                               
      end                                                                       
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine upres                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 05/11/2015 rhd             *          
c     *                                                              *          
c     *     compute structure level, residual force vector for use   *          
c     *     in the next Newton iteration. the next iteration may not *          
c     *     occur if the residuals indicate convergence during the   *          
c     *     convergence tests following this code                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine upres(                                                         
     &    iter, iout, nodof, dt, nbeta, num_term_loads,                         
     &    sum_loads, mgload, mdiag, pbar, du, velocity, accel,                  
     &    ifv, res, cstmap, dstmap, load  )                                     
c                                                                               
      use mod_mpc, only : tied_con_mpcs_constructed, mpcs_exist                 
      use stiffness_data, only : d_lagrange_forces                              
      use contact, only : contact_force                                         
c                                                                               
      implicit none                                                             
c                                                                               
c              parameters                                                       
c                                                                               
      integer :: iter, iout, nodof, num_term_loads, cstmap(nodof),              
     &           dstmap(nodof)                                                  
      double precision ::                                                       
     &  dt, nbeta, sum_loads, mgload, mdiag(nodof), pbar(nodof),                
     &  du(nodof), velocity(nodof), accel(nodof), ifv(nodof),                   
     &  res(nodof), load(nodof)                                                 
c                                                                               
c              locals                                                           
c                                                                               
      integer :: i, chknode, chkdof                                             
      double precision ::                                                       
     &  nfac, zero, one, two, accel_n1, term_load, nfac1, nfac2,                
     &  force_lag, total_external, inertia_contact                              
      logical :: have_mpc_equations                                             
      logical, parameter :: local_debug = .false.                               
      data zero, one, two / 0.0d00, 1.0d00, 2.0d00 /                            
c                                                                               
      mgload = zero                                                             
      nfac   = one / (nbeta*dt*dt)                                              
      nfac1  = (one-two*nbeta)/(two*nbeta)                                      
      nfac2  = one/(nbeta*dt)                                                   
c                                                                               
      num_term_loads = 0                                                        
      sum_loads      = zero                                                     
      have_mpc_equations = tied_con_mpcs_constructed .or. mpcs_exist            
c                                                                               
      if( local_debug ) write(iout,9300) iter, have_mpc_equations               
c                                                                               
c              for iter = 1, 2, 3.. we are computing the                        
c              actual residual load. compute the euclidean                      
c                                                                               
c              residual @ n+1 = applied @ n+1 - internal @ n+1                  
c                                                                               
c              applied = pbar + contact + dlagrange - inertia                   
c                        adjustment                                             
c                                                                               
c                 pbar includes user-applied forces at n+1,                     
c                      no contact, mpc nodal forces at n,                       
c                      inertia forces at n                                      
c                                                                               
c                 Thus, applied here needs the full contact forces,             
c                 the change in mpc nodal forces (sometimes called              
c                 Lagrange multipliers) over the step and the                   
c                 change in inertial forces over the step.                      
c                                                                               
c              internal = simply the assembled int B^T sigma                    
c                                                                               
      if( have_mpc_equations ) then                                             
        do i = 1, nodof                                                         
          res(i) = ( pbar(i) + contact_force(i) + d_lagrange_forces(i)          
     &               - mdiag(i)*du(i)*nfac ) - ifv(i)                           
         end do                                                                 
      else                                                                      
        do i = 1, nodof                                                         
          res(i) = ( pbar(i) + contact_force(i)                                 
     &               - mdiag(i)*du(i)*nfac ) - ifv(i)                           
        end do                                                                  
      end if                                                                    
c                                                                               
      if( local_debug ) then                                                    
        write(iout,9400); write(iout,9405)                                      
        do i = 1, nodof                                                         
          force_lag = zero                                                      
          if( have_mpc_equations ) force_lag = d_lagrange_forces(i)             
          total_external = pbar(i) + contact_force(i) +                         
     &                     force_lag - mdiag(i)*du(i)*nfac                      
          write(iout,9200) i, cstmap(i), res(i), pbar(i),                       
     &                     contact_force(i), force_lag,                         
     &                     mdiag(i)*du(i)*nfac, total_external,                 
     &                     ifv(i)                                               
        end do                                                                  
      end if                                                                    
c                                                                               
c              compute a norm and sum of the current estimate of the            
c              total, applied nodal forces at n+1 for use in                    
c              some convergence tests                                           
c                                                                               
c              treat reaction forces separately. then zero residual             
c              terms for dof w/ absolute constraints.                           
c              ifv is negative of reactions at dof with absolute cons           
c                                                                               
c              Nodal Lagrange forces for MPCs are not included in               
c              the loading measure. They are always a set of                    
c              equilibrating forces.                                            
c                                                                               
      do i = 1, nodof                                                           
        accel_n1 = du(i)*nfac - nfac1*accel(i) - nfac2*velocity(i)              
        inertia_contact = abs( accel_n1) * mdiag(i) +                           
     &                    abs( contact_force(i) )                               
        if( cstmap(i) .eq. 0 ) then ! no abs con                                
           term_load = inertia_contact + abs( load(i) )                         
        else  ! absolute constraint                                             
           term_load = inertia_contact +                                        
     &                 abs( ifv(i) ) !  may want to include + abs( load(i) )    
           res(i)    = zero                                                     
        end if                                                                  
        if( term_load .ne. zero ) then                                          
           mgload = mgload + term_load * term_load                              
           num_term_loads = num_term_loads + 1                                  
           sum_loads = sum_loads + abs(term_load)                               
        end if                                                                  
      end do                                                                    
c                                                                               
      mgload = sqrt( mgload )                                                   
      if( local_debug ) write(iout,9100) mgload                                 
c                                                                               
      if( local_debug ) write(iout,9305)                                        
c                                                                               
      return                                                                    
c                                                                               
 9100 format(5x, "... mgload: ",e14.6)                                          
 9200 format(2x,2i7,7f16.6)                                                     
 9300 format(/,1x,"--- entering upres. iter, have multipoint: ",i3,l2)          
 9305 format(1x,"--- leaving upres ---" )                                       
 9400 format(5x,"... various force vector:")                                    
 9405 format(t5,"sdof", t11,"cstmap", t25,"res", t40, "pbar",                   
     &  t51, "contact force",                                                   
     &  t70, "dlagrange", t82, "inertia adjust",t101,"total ext",t120,          
     &   "ifv" )                                                                
 9990 format('... bugs: ',i6,2f15.6)                                            
c                                                                               
      end                                                                       
                                                                                
                                                                                
