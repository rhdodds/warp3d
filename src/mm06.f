c                                                                               
c           Updated:  17/30/2016 rhd                                            
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06_set_sizes                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 10/1/2015 rhd               *          
c     *                                                              *          
c     *  called by warp3d to obtain history size and                 *          
c     *  other characteristic information about the model            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_set_sizes( info_vector )                                  
      implicit none                                                             
      integer :: info_vector(*)                                                 
c                                                                               
c        set info_data                                                          
c                                                                               
c         1        number of history values per integration                     
c                  point.                                                       
c                                                                               
c         2        number of values in the symmetric part of the                
c                  [D] for each integration point. for solid                    
c                  elements this is 21, for cohesive elements this 6.           
c                                                                               
c         3        = 0, the material model returns "unrotated"                  
c                       Cauchy stresses at n+1                                  
c                  = 1, the material model returns the standard                 
c                       Cauchy stresses at n+1                                  
c                                                                               
c         4        number of history variables per point to be output           
c                  when user requests this type of results                      
c                                                                               
      info_vector(1) = 3                                                        
      info_vector(2) = 21                                                       
      info_vector(3) = 0                                                        
      info_vector(4) = 3                                                        
c                                                                               
      return                                                                    
      end subroutine mm06_set_sizes                                             
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06                              *          
c     *                                                              *          
c     *              Norton power-law creep model                    *          
c     *                                                              *          
c     *           last modified: 3/5/2016 rhd                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine mm06( step, iter, felem, gpn, mxvl, hist_size,                 
     &                 nstrs, nstr, span, kout, dtime, block_props,             
     &                 e_vec, nu_vec, n_power_vec, rtse,                        
     &                 stress_n, stress_np1, dstran,                            
     &                 history_n, history_np1, elem_killed_vec,                 
     &                 do_nonlocal, nonlocal_state, maxnonlocal,                
     &                 cep, compute_creep_strains  )                            
      implicit none                                                             
c                                                                               
c             parameter definitions                                             
c                                                                               
      integer :: step, iter, felem, gpn, mxvl, hist_size,                       
     &           nstrs, nstr, span, kout, maxnonlocal                           
      logical :: elem_killed_vec(mxvl), do_nonlocal,                            
     &           compute_creep_strains                                          
      double precision                                                          
     &  dtime, block_props(mxvl,10), e_vec(mxvl), nu_vec(mxvl),                 
     &  n_power_vec(mxvl), rtse(mxvl,nstr), stress_n(mxvl,nstrs),               
     &  stress_np1(mxvl,nstrs), dstran(mxvl,nstr),                              
     &  history_n(span,hist_size), history_np1(span,hist_size),                 
     &  nonlocal_state(mxvl,maxnonlocal), cep(mxvl,6,6)                         
c                                                                               
c             local variables                                                   
c                                                                               
      integer :: i, j                                                           
      logical :: debug                                                          
      double precision                                                          
     &  zero, half, one, two, three, four, five, six, roothalf,                 
     &  stress(6), statev(10), t_c, ddsdde(6,6),                                
     &  props(10), nrm_dtime, local_dstran(6), local_rtse(6)                    
c                                                                               
      data zero, half, one, two, three, four, five, six, roothalf               
     &  / 0.0d0, 0.5d00, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0,              
     &    0.70710678119d0  /                                                    
c                                                                               
c        this implementation was originally formulated to work with             
c        a normalized stress, modulus and time. That code is                    
c        retained as mm06n and still used as described here.                    
c        Now only time is normalized so that the stress update                  
c        doesn't need B. We can think of stresses as normalized by              
c        a sigma_0 = 1                                                          
c                                                                               
c        The Norton power-law model has no notion of a Sigma_0 -- this          
c        was introduced by the original umat developers for                     
c        purposes of normalization only.                                        
c                                                                               
c        The user defines E, n, nu and B                                        
c        for the creep material.                                                
c        Compute a t_c = 1/B then call the mm06n routine which                  
c        sees properly normalized values.                                       
c                                                                               
c          props(1) - Young's modulus. (often 189 GPa)                          
c          props(2) - Poisson's ratio: nu (often 0.285 in examples)             
c          props(3) - creep exponent: n (often 5.0 in examples)                 
c          props(4) - creep factor B                                            
c                                                                               
c        state variables                                                        
c                                                                               
c          statev(1) - effective creep strain (total)                           
c          statev(2) - effective creep strain rate wrt real time                
c                      (mm06n must have the normalized rate)                    
c          statev(3) - increment of effective creep eps over step               
c                                                                               
      debug = felem .eq. 1  .and. gpn .eq. 8                                    
      debug = .false.                                                           
      if( debug ) write(kout,9010) felem, gpn                                   
      if( debug ) then                                                          
          write(kout,9019)                                                      
          do i = 1, span                                                        
            write(kout,9020) felem+i-1, (dstran(i,j),j=1,6)                     
          end do                                                                
      end if                                                                    
c                                                                               
c             iter = 0 and no global displacement extrapolation.                
c             mm06 just removes estimated creep strain                          
c             increment for load (time) step from dstran and fill               
c             [D] with linear-elastic values.                                   
c                                                                               
      if( compute_creep_strains )  then                                         
        call mm06_step_crp_strains                                              
        if( debug ) then                                                        
           write(kout,9021)                                                     
           do i = 1, span                                                       
               write(kout,9020) felem+i-1, (dstran(i,j),j=1,6)                  
           end do                                                               
        end if                                                                  
        do i = 1, span                                                          
         call cnst6_linear_elastic( e_vec(i), nu_vec(i),                        
     &                              ddsdde, one, two, zero )                    
         cep(i,1:6,1:6) = ddsdde(1:6,1:6)                                       
         stress_np1(i,1:6) = stress_n(i,1:6)                                    
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
      if( step .eq. 1 ) then                                                    
         do i = 1, span                                                         
           history_n(i,1) = zero  ! total creep eff strain                      
           history_n(i,2) = zero  ! current creep rate                          
           history_n(i,3) = zero  ! estimate of delta-creep for step            
         end do                                                                 
      end if                                                                    
c                                                                               
c             loop over all element in block for this integration               
c             point. build simple vectors, constants needed                     
c             by the routine that updates a single point.                       
c                                                                               
      do i = 1, span                                                            
        if( elem_killed_vec(i) ) stress_np1(i,1:6) = zero                       
      end do                                                                    
c                                                                               
      do i = 1, span                                                            
      if( elem_killed_vec(i) ) cycle                                            
      props(1) = e_vec(i)  ! emod                                               
      props(2) = nu_vec(i)                                                      
      props(3) = n_power_vec(i)                                                 
      props(4) = block_props(i,1)                                               
c                                                                               
c             the normalized mm06n code does not know about B.                  
c             we normalize time to accommodate. strain rates                    
c             to-from mm06 must be normalized time                              
c                                                                               
      t_c       = one /  block_props(i,1)  ! 1.0/B                              
      nrm_dtime = dtime / t_c                                                   
      statev(1) = history_n(i,1)                                                
      statev(2) = history_n(i,2) * t_c ! to normalized time                     
      statev(3) = history_n(i,3)                                                
      stress(1:6) = stress_n(i,1:6)                                             
      local_dstran(1:6) = dstran(i,1:6)                                         
c                                                                               
      if( debug ) then                                                          
        write(kout,*) ' ..... before update ......'                             
        write(kout,9200) block_props(i,1), dtime, t_c, nrm_dtime                
        write(kout,9300) (stress(j),j=1,6)                                      
        write(kout,9310) statev(1), statev(2), statev(3)                        
        write(kout,9322) (local_dstran(j),j=1,6)                                
      end if                                                                    
c                                                                               
      call mm06n( stress, local_rtse, statev, local_dstran, nrm_dtime,          
     &            props, felem+i-1, gpn, step, iter, kout )                     
c                                                                               
c             updated stress is in physical units since we used                 
c             an implicit Sigma_0 = 1 for normalization                         
c                                                                               
c             the creep strain rate in statev(2) just set by mm06n              
c             is wrt normalized time as defined for this solid                  
c             creeping material. convert to creep strain rate wrt               
c             to real time                                                      
c                                                                               
      stress_np1(i,1:6) = stress(1:6)                                           
      rtse(i,1:6)       = local_rtse(1:6)                                       
      statev(2)         = statev(2) / t_c                                       
      history_np1(i,1)  = statev(1)                                             
      history_np1(i,2)  = statev(2)                                             
      history_np1(i,3)  = statev(3)                                             
c                                                                               
c            nonlocal state values returned are creep rate wrt real time        
c            and total creep strain                                             
c                                                                               
      if( do_nonlocal ) then                                                    
        nonlocal_state(i,1) = statev(2) ! eff creep eps rate - real time        
        nonlocal_state(i,2) = statev(1) ! eff creep eps                         
        nonlocal_state(i,3) = props(3)  !   n                                   
        nonlocal_state(i,4) = props(4)  !   B                                   
      end if                                                                    
c                                                                               
      call cnst6n( local_rtse, ddsdde, statev, nrm_dtime,                       
     &             props, felem+i-1, gpn, iter, kout )                          
      cep(i,1:6,1:6) = ddsdde(1:6,1:6)                                          
c                                                                               
      call mm06_get_eng_dens( i )                                               
c                                                                               
      if( debug ) then                                                          
        write(kout,*) ' ..... updated values ......'                            
        write(kout,9300) (stress(j),j=1,6)                                      
        write(kout,9310) statev(1), statev(2), statev(3)                        
        write(kout,*) " "                                                       
      end if                                                                    
c                                                                               
      end do  ! over span                                                       
c                                                                               
      if( debug ) write(kout,9012) felem, gpn                                   
      return                                                                    
c                                                                               
 9010 format(/,1x,'.....  entered mm06 creep model, element, gpn: ',            
     &          2i7)                                                            
 9012 format(1x,'.....  leave mm06 creep model, element, gpn: ',                
     &          2i7,/)                                                          
 9019 format(5x,'dstran before stp_crp.... ')                                   
 9020 format(5x,i6,6e14.6)                                                      
 9021 format(5x,'dstran after stp_crp....')                                     
 9022 format(5x,'dstran after  stp_crp : ',6d14.6)                              
 9200 format(5x,'B, dtime, t_c, nrm_dtime: ', 4d14.6)                           
 9300 format(5x,'stresses:     ',6f10.4 )                                       
 9310 format(5x,'statev: ',3d15.8)                                              
 9322 format(5x,'dstran: ',6d14.6)                                              
 9400 format(5x,'stress: ',6f10.5)                                              
c                                                                               
      contains                                                                  
c     =========                                                                 
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06_get_eng_dens                 *          
c     *                                                              *          
c     *        update energy density at gpn for elements in block    *          
c     *                                                              *          
c     *                 last modified: 09/30/2015                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_get_eng_dens( relem )                                     
c                                                                               
c             locally defined variables                                         
c                                                                               
      integer :: relem                                                          
      double precision ::                                                       
     &  dsigma(6), sigma_avg(6), deps_elas(6), emod, nu,                        
     &  shear_mod, loc_dstran(6)                                                
c                                                                               
c             U0(n+1) = U0(n) + dU                                              
c             dU = trans( sigma_{n+1/2} ) .dot.  dstran                         
c             deps_{thermal} already taken out by WARP3D                        
c             => will be under eng_dens in stress output                        
c                                                                               
c             Ue(n+1) = Ue(n) + dUe   (elastic work)                            
c             dUe = trans( sigma_{n+1/2} ) .dot.  deps_elas                     
c             deps_elas = [C] * dsigma                                          
c             => will be under "c3" in stress output                            
c                                                                               
      if( step .eq. 1 ) then                                                    
         stress_n(relem,7) = zero                                               
         stress_n(relem,9) = zero  ! oumm06 moves this at output                
      end if                                                                    
c                                                                               
      loc_dstran(1:6) = dstran(relem,1:6)                                       
      dsigma(1:6) = stress_np1(relem,1:6) - stress_n(relem,1:6)                 
      sigma_avg(1:6) = half * ( stress_np1(relem,1:6) +                         
     &                          stress_n(relem,1:6) )                           
      emod = props(1)                                                           
      nu   = props(2)                                                           
      shear_mod = half * emod / (one + nu )                                     
c                                                                               
      deps_elas(1) =                                                            
     &           ( dsigma(1) - nu*dsigma(2) - nu*dsigma(3) ) / emod             
      deps_elas(2) =                                                            
     &           ( -nu*dsigma(1) + dsigma(2) - nu*dsigma(3) ) / emod            
      deps_elas(3) =                                                            
     &           ( -nu*dsigma(1) - nu*dsigma(2) + dsigma(3) ) / emod            
      deps_elas(4) =  dsigma(4) / shear_mod                                     
      deps_elas(5) =  dsigma(5) / shear_mod                                     
      deps_elas(6) =  dsigma(6) / shear_mod                                     
      stress_np1(relem,7) = stress_n(relem,7) +                                 
     &                      dot_product( sigma_avg, loc_dstran )                
      stress_np1(relem,9) = stress_n(relem,9) +                                 
     &                      dot_product( sigma_avg, deps_elas )                 
c                                                                               
      if( debug ) then                                                          
        write(kout,*) "     update energy density"                              
        write(kout,9010) stress_n(relem,7)                                      
        write(kout,9020) dsigma                                                 
        write(kout,9030) sigma_avg                                              
        write(kout,9040) deps_elas                                              
        write(kout,9050) dot_product( sigma_avg, loc_dstran )                   
        write(kout,9060) dot_product( sigma_avg, deps_elas )                    
       end if                                                                   
c                                                                               
      return                                                                    
c                                                                               
 9010 format(5x,'end dens @n: ',d14.6)                                          
 9020 format(5x,'dsigma:      ',6d14.6)                                         
 9030 format(5x,'sigma_avg:   ',6d14.6)                                         
 9040 format(5x,'deps_elas:   ',6d14.6)                                         
 9050 format(5x,'dU0:         ',d14.6)                                          
 9060 format(5x,'dUelas:      ',d14.6)                                          
c                                                                               
      end subroutine mm06_get_eng_dens                                          
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06_step_crp_strains             *          
c     *                                                              *          
c     *        update dstran for estimate of creep strain increment  *          
c     *        over step                                             *          
c     *                                                              *          
c     *                 last modified: 09/11/2015                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_step_crp_strains                                          
      implicit none                                                             
c                                                                               
c             locally defined variables                                         
c                                                                               
      integer :: i                                                              
      double precision ::                                                       
     &  sig_mean, sig_mises, n_exp, B, factor, t1, t2, t3, t4, t5               
                                                                                
c                                                                               
      do i = 1, span                                                            
c                                                                               
        sig_mean  = ( stress_n(i,1) + stress_n(i,2) +                           
     &                stress_n(i,3) ) / three                                   
        t1 = stress_n(i,1) - stress_n(i,2)                                      
        t2 = stress_n(i,3) - stress_n(i,2)                                      
        t3 = stress_n(i,3) - stress_n(i,1)                                      
        t4 = t1*t1 + t2*t2 + t3*t3                                              
        t5 = stress_n(i,4)**2 + stress_n(i,5)**2 + stress_n(i,6)**2             
        sig_mises = roothalf * sqrt( t4 + six*t5 )                              
        n_exp  = n_power_vec(i)                                                 
        B      = block_props(i,1)                                               
        factor = (three/two) * dtime * B * sig_mises**(n_exp-one)               
        dstran(i,1) = dstran(i,1) -                                             
     &                ( stress_n(i,1) - sig_mean ) * factor                     
        dstran(i,2) = dstran(i,2) -                                             
     &                ( stress_n(i,2) - sig_mean ) * factor                     
        dstran(i,3) = dstran(i,3) -                                             
     &                ( stress_n(i,3) - sig_mean ) * factor                     
        dstran(i,4) = dstran(i,4) - two * stress_n(i,4) * factor                
        dstran(i,5) = dstran(i,5) - two * stress_n(i,5) * factor                
        dstran(i,6) = dstran(i,6) - two * stress_n(i,6) * factor                
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
      end subroutine mm06_step_crp_strains                                      
                                                                                
      end subroutine mm06                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06n                             *          
c     *                                                              *          
c     *           last modified: 9/25/2015 rhd                       *          
c     *           should get inlined                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06n( stress, rtse, statev, dstran, dtime, props,             
     &                  noel, npt, kinc, kiter, kout )                          
c                                                                               
      implicit none                                                             
c                                                                               
c             parameter definitions                                             
c                                                                               
      integer :: noel, npt, kinc, kiter, kout                                   
      double precision ::                                                       
     & stress(6), rtse(6), statev(*), dtime, dstran(6), props(*)                
c                                                                               
c                                                                               
c             key locally defined variables                                     
c                                                                               
c     ssdevnp           : deviatoric stress in present time step t+dt           
c     ssdevnt           : deviatoric stress predictor                           
c     dsndev            : deviatoric strain increment                           
c     smisest           : von mises stress predictor                            
c     smisesnp          : von mises stress in present time step t+dt            
c     sstracen          : stress invariant I1                                   
c     decrnew           : effective strain increment after NR iteration         
c                                                                               
      double precision ::                                                       
     1 re(2), ssdevnt(6), dsndev(6), toler,                                     
     2 ssdevn(6), ssdevnp(6),                                                   
     3 one, two, three, six, zero,                                              
     4 smisest, smisesnp,                                                       
     5 sstracen, decreq1, decreq2, decrnew,                                     
     6 emod, enu, ebulk3, eg2, eg, elam, n,                                     
     7 betastar, f1, re_temp, temp3,                                            
     8 sstracenp, nonlocal_shared(5), tol_null                                  
      logical :: debug                                                          
c                                                                               
      data zero, one, two, three, six, toler, tol_null                          
     & /0.0d0, 1.0d0, 2.0d0, 3.0d0, 6.0d0, 1.d-15, 1.0d-10 /                    
c                                                                               
c        - see time normalization in calling routine                            
c                                                                               
c     props(1) - Young's modulus                                                
c     props(2) - Poisson's ratio: nu (often 0.285 in examples)                  
c     props(3) - creep exponent: n (often 5.0 in examples)                      
c                                                                               
c             state variables                                                   
c                                                                               
c     statev(1) - total effective creep strain                                  
c     statev(2) - effective strain rate in normalized time                      
c     statev(3) - effective creep increment for step                            
c                                                                               
c             material parameters                                               
c                                                                               
      emod    = props(1)                                                        
      enu     = props(2)                                                        
      ebulk3  = emod/(one-two*enu)                                              
      eg2     = emod/(one+enu)                                                  
      eg      = eg2/two                                                         
      elam    = (ebulk3-eg2)/three                                              
      n       = props(3)                                                        
c                                                                               
c                                                                               
c             step (1) -- compute deviatoric                                    
c                         stress at the current time step                       
c                                                                               
      sstracen = stress(1) + stress(2) + stress(3)                              
      ssdevn(1:3) = stress(1:3) - sstracen/three                                
      ssdevn(4:6) = stress(4:6)                                                 
c                                                                               
c             step (2) -- compute deviatoric strain increment                   
c                         and J2 at the current time step                       
c                                                                               
      dsndev(1:3) = dstran(1:3)-(dstran(1)+dstran(2)+dstran(3))/three           
      dsndev(4:6) = dstran(4:6)                                                 
c                                                                               
c             step (3) -- compute predictor stress and mises                    
c                         for predictor state. check for loading with           
c                         null deviator of trial stress. stresses               
c                         then linear elastic. a zero mises                     
c                         stress for trial state will blow-up                   
c                         power-law computations. compare relative              
c                         to emod which is always > 0                           
c                                                                               
      ssdevnt(1:3) = ssdevn(1:3) + eg2*dsndev(1:3)                              
      ssdevnt(4:6) = ssdevn(4:6) + eg*dsndev(4:6)                               
      smisest = (ssdevnt(1)-ssdevnt(2))*(ssdevnt(1)-ssdevnt(2))+                
     1          (ssdevnt(2)-ssdevnt(3))*(ssdevnt(2)-ssdevnt(3))+                
     2          (ssdevnt(3)-ssdevnt(1))*(ssdevnt(3)-ssdevnt(1))                 
      smisest = smisest + six*ssdevnt(4)*ssdevnt(4)                             
      smisest = smisest + six*ssdevnt(5)*ssdevnt(5)                             
      smisest = smisest + six*ssdevnt(6)*ssdevnt(6)                             
      smisest = sqrt(smisest/two)                                               
c                                                                               
c             Update stress, state variables and nonlocal variables.            
c                                                                               
c             If mises for trial elastic stress = 0, the                        
c             creep strain-rate is zero. Updated stresses                       
c             are the trial elastic values. Creep stress update                 
c             cannot have smisest -> 0 (it is in denominator)                   
c                                                                               
      if( smisest .lt. tol_null * emod ) then                                   
        sstracenp = sstracen + ebulk3*(dstran(1)+dstran(2)+dstran(3))           
        stress(1:3) = ssdevnt(1:3) + sstracenp/three                            
        stress(4:6) = ssdevnt(4:6)                                              
        smisesnp = smisest  ! for clarity                                       
        decrnew  = zero     ! for clarity                                       
        statev(1) = statev(1) + decrnew                                         
        statev(2) = smisesnp**n                                                 
        statev(3) = decrnew                                                     
        nonlocal_shared(1) = smisesnp**n                                        
        nonlocal_shared(2) = statev(1)                                          
        rtse(1:6) = ssdevnt(1:6) ! for use in cnst6                             
        return                                                                  
      end if                                                                    
c                                                                               
c             step (4) -- nonlinear scalar consistency equation                 
c                         to compute effective creep strain increment           
c                         (decrnew). Uses Pegasus root finder                   
c                                                                               
      decrnew = zero                                                            
      decreq1 = zero                                                            
      decreq2 = smisest**n*dtime  ! upper estimate                              
      call mm06_derivatives( decreq1, dtime, n, emod, enu,                      
     &                       smisest, f1 )                                      
      re(1) = f1                                                                
      call mm06_derivatives( decreq2, dtime, n, emod, enu,                      
     &                       smisest, f1 )                                      
      re(2) = f1                                                                
      do while( abs(decreq2-decreq1) .gt. toler )                               
         temp3 = (decreq1*re(2)-decreq2*re(1))/(re(2)-re(1))                    
         call mm06_derivatives( temp3, dtime, n, emod, enu,                     
     &                          smisest, f1 )                                   
         re_temp = f1                                                           
         if( re_temp*re(2) .lt. toler ) then                                    
              decreq1 = decreq2                                                 
              call mm06_derivatives( decreq1, dtime, n, emod, enu,              
     &                               smisest, f1 )                              
              re(1) = f1                                                        
         else                                                                   
              re(1) = re(1)*re(2)/(re(2)+re_temp)                               
         end if                                                                 
         decreq2 = temp3                                                        
         re(2) = re_temp                                                        
      end do                                                                    
      decrnew = decreq2                                                         
c                                                                               
c             step (5) -- compute elastic-plastic stress at                     
c                         end of time step.                                     
c                                                                               
      betastar = one - three*emod/two/(one+enu)/smisest*decrnew                 
      smisesnp = betastar*smisest                                               
      ssdevnp(1:3) = betastar*ssdevnt(1:3)                                      
      ssdevnp(4:6) = betastar*ssdevnt(4:6)                                      
      sstracenp = sstracen + ebulk3*(dstran(1)+dstran(2)+dstran(3))             
      stress(1:3) = ssdevnp(1:3) + sstracenp/three                              
      stress(4:6) = ssdevnp(4:6)                                                
      rtse(1:6)   = ssdevnt(1:6) ! for use in [D] computation                   
c                                ! see cnst6                                    
c                                                                               
c             step (7) -- update the state variables                            
c                                                                               
      statev(1) = statev(1) + decrnew                                           
      statev(2) = smisesnp**n ! creep strain rate wrt normalized time           
      statev(3) = decrnew                                                       
c                                                                               
c             step (8) -- nonlocal variables set up by umat and                 
c                         passed to mm04.f                                      
c                                                                               
      nonlocal_shared(1) = smisesnp**n   ! creep strain rate in normalized time 
      nonlocal_shared(2) = statev(1)     ! updated total creep (equiv) strain   
c                                                                               
      return                                                                    
c                                                                               
9000  format(3x,'(',i1,',',i1,'): ', 2e14.6 )                                   
c                                                                               
      end subroutine mm06n                                                      
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine cnst6n                            *          
c     *                                                              *          
c     *            last modified: 5/25/2015 rhd                      *          
c     *            should be inlined                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine cnst6n( ssdevnt, ddsdde, statev, dtime, props,                 
     &                   noel, npt,  kiter, kout )                              
c                                                                               
      implicit none                                                             
c                                                                               
c             parameter definitions                                             
c                                                                               
      integer :: noel, npt, kinc, kiter, kout                                   
      double precision ::                                                       
     & ssdevnt(6), statev(3), dtime, props(3), ddsdde(6,6)                      
c                                                                               
c             locals                                                            
c                                                                               
      integer :: j, k                                                           
      double precision ::                                                       
     1 expon, sigma_e_star, root32, mu, alpha, kmod,                            
     2 one, two, three, zero, aprime, beta1, beta2, beta,                       
     3 decrnew, emod, enu, tol_null                                             
      logical :: debug, null_deviator                                           
                                                                                
      double precision c1, c2, c3, c4                                           
c                                                                               
      data zero, one, two, three, tol_null, root32                              
     & / 0.0d0, 1.0d0, 2.0d0, 3.0d0, 1.0d-15, 1.22474487139d0 /                 
                                                                                
c             material properties from user input                               
c                                                                               
c     props(1) - Young's modulus                                                
c     props(2) - Poisson's ratio: nu (often 0.285 in examples)                  
c     props(3) - creep exponent: n (often 5.0 in examples)                      
c                                                                               
c             state variables                                                   
c                                                                               
c     statev(1) - total effective creep strain                                  
c     statev(2) - effective strain rate in normalized time                      
c     statev(3) - effective creep strain increment over step                    
c                 from last stress update.  consistent                          
c                 with ssdevnt                                                  
c                                                                               
c             material parameters                                               
c                                                                               
      emod    = props(1)                                                        
      enu     = props(2)                                                        
      expon   = props(3)                                                        
      decrnew = statev(3)                                                       
c                                                                               
c             compute the elastic-plastic consistent tangent. couples           
c             normalized stress increments to strain increments                 
c             where time is normalized as tbar                                  
c                                                                               
      mu = emod/two/(one+enu)                                                   
      kmod = emod/three/(one-two*enu)                                           
c                                                                               
      sigma_e_star = ssdevnt(1)**2 + ssdevnt(2)**2 + ssdevnt(3)**2 +            
     &    ( ssdevnt(4)**2 + ssdevnt(5)**2 + ssdevnt(6)**2 ) * two               
      sigma_e_star = root32 * sqrt(sigma_e_star)                                
c                                                                               
      if( sigma_e_star .lt. tol_null*emod ) then                                
        ddsdde(1:6,1:6) = zero                                                  
        c1 = emod/((one+enu)*(one-two*enu))                                     
        c2 = (one-enu)*c1                                                       
        c3 = ((one-two*enu)/two)*c1                                             
        c4 = enu*c1                                                             
        ddsdde(1,1)= c2                                                         
        ddsdde(2,2)= c2                                                         
        ddsdde(3,3)= c2                                                         
        ddsdde(4,4)= c3                                                         
        ddsdde(5,5)= c3                                                         
        ddsdde(6,6)= c3                                                         
        ddsdde(1,2)= c4                                                         
        ddsdde(1,3)= c4                                                         
        ddsdde(2,1)= c4                                                         
        ddsdde(3,1)= c4                                                         
        ddsdde(2,3)= c4                                                         
        ddsdde(3,2)= c4                                                         
        return                                                                  
      end if                                                                    
c                                                                               
      alpha = one - three*mu*decrnew/sigma_e_star                               
      aprime = (one/expon)*(one/dtime)*                                         
     &         (decrnew/dtime)**((one-expon)/expon)                             
      beta1 = three/two/sigma_e_star/sigma_e_star                               
      beta2 = -three*mu +(one-alpha)*(three*mu+aprime)                          
      beta  = beta1 * beta2/(three*mu+aprime)                                   
c                                                                               
      do j = 1, 6                                                               
       do k = 1, 6                                                              
        ddsdde(k,j) = two * mu * beta * ssdevnt(k) * ssdevnt(j)                 
       end do                                                                   
      end do                                                                    
c                                                                               
      ddsdde(1,1) = ddsdde(1,1) + two*mu*alpha                                  
      ddsdde(2,2) = ddsdde(2,2) + two*mu*alpha                                  
      ddsdde(3,3) = ddsdde(3,3) + two*mu*alpha                                  
      ddsdde(4,4) = ddsdde(4,4) + mu*alpha                                      
      ddsdde(5,5) = ddsdde(5,5) + mu*alpha                                      
      ddsdde(6,6) = ddsdde(6,6) + mu*alpha                                      
c                                                                               
      ddsdde(1:3,1:3) = ddsdde(1:3,1:3) +                                       
     &                        (kmod - two*mu*alpha/three)                       
c                                                                               
      return                                                                    
c                                                                               
      end subroutine cnst6n                                                     
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine cnst6_linear_elastic              *          
c     *                                                              *          
c     *       linear-elastic [D] matrix. should be inlined           *          
c     *                                                              *          
c     *                 last modified: 10/18/2015 rhd                *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
      subroutine cnst6_linear_elastic( emod, enu, ddsdde,                       
     &                                 one, two, zero )                         
      implicit none                                                             
      double precision ::                                                       
     &  ddsdde(6,6), emod, enu, one, two, zero, c1, c2, c3, c4                  
c                                                                               
c                       linear isotropic elastic matrix.                        
c                                                                               
      ddsdde(1:6,1:6) = zero                                                    
      c1 = emod/((one+enu)*(one-two*enu))                                       
      c2 = (one-enu)*c1                                                         
      c3 = ((one-two*enu)/two)*c1                                               
      c4 = enu*c1                                                               
      ddsdde(1,1)= c2                                                           
      ddsdde(2,2)= c2                                                           
      ddsdde(3,3)= c2                                                           
      ddsdde(4,4)= c3                                                           
      ddsdde(5,5)= c3                                                           
      ddsdde(6,6)= c3                                                           
      ddsdde(1,2)= c4                                                           
      ddsdde(1,3)= c4                                                           
      ddsdde(2,1)= c4                                                           
      ddsdde(3,1)= c4                                                           
      ddsdde(2,3)= c4                                                           
      ddsdde(3,2)= c4                                                           
c                                                                               
      return                                                                    
      end subroutine cnst6_linear_elastic                                       
                                                                                
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06_derivatives                  *          
c     *                                                              *          
c     *         for power-law (elastic) creep. should be inlined     *          
c     *                                                              *          
c     *                 last modified: 09/21/2015 rhd                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_derivatives( decreq, dtime, n, emod, enu,                 
     1                             smisest, f1 )                                
      implicit none                                                             
      double precision ::                                                       
     & smisest, decreq, dtime, emod, enu, f1, n                                 
c                                                                               
c             locally defined variables                                         
c                                                                               
      double precision ::                                                       
     & one, two, three, temp1                                                   
      data one, two, three                                                      
     & / 1.0d0, 2.0d0, 3.0d0 /                                                  
c                                                                               
      temp1 = one - three*emod*decreq/two/(one+enu)/smisest                     
      f1    = smisest*temp1 - (decreq/dtime)**(one/n)                           
c                                                                               
      return                                                                    
      end subroutine mm06_derivatives                                           
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine mm06_states_values                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 9/24/2015 rhd                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_states_values( itype, elem_states_output,                 
     &                                nrow_states, num_states  )                
      use global_data ! old common.main
c                                                                               
c                       access some global data structures                      
c                                                                               
      use elem_block_data, only: history_blocks, history_blk_list               
      use main_data, only: elems_to_blocks                                      
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: nrow_states, itype, num_states                                 
      double precision :: elem_states_output(nrow_states,*)                     
c                                                                               
c                       locals                                                  
c                                                                               
      double precision,                                                         
     & allocatable :: history_dump(:,:,:), one_elem_states(:)                   
      integer :: relem, elnum, hist_size, blockno                               
      logical :: do_a_block, local_debug                                        
      double precision :: zero                                                  
      data zero / 0.0d00 /                                                      
c                                                                               
c           build creep model states values output.                             
c                                                                               
c              itype > 0 => this is the block number. do all elements           
c                           in the block                                        
c                                                                               
c              itype < 0 => this is an element number. put state                
c                           values into column 1 of results.                    
c                                                                               
      do_a_block = .true.                                                       
      if( itype. gt. 0 ) then                                                   
         do_a_block = .true.                                                    
         blockno = itype                                                        
      else                                                                      
         do_a_block = .false.                                                   
         elnum = -itype                                                         
         blockno = elems_to_blocks(elnum,1)                                     
      end if                                                                    
c                                                                               
      local_debug = .false.                                                     
      felem       = elblks(1,blockno)                                           
      elem_type   = iprops(1,felem)                                             
      mat_type    = iprops(25,felem)                                            
      int_points  = iprops(6,felem)                                             
      span        = elblks(0,blockno)                                           
      hist_size   = history_blk_list(blockno)                                   
                                                                                
      if( local_debug ) write(out,9050) blockno, felem, elem_type,              
     &         mat_type, int_points, span, hist_size                            
c                                                                               
c           temporary block of history so it can be re-organized                
c                                                                               
      allocate( one_elem_states(nrow_states) )                                  
      allocate( history_dump(hist_size,int_points,span) )                       
      history_dump = reshape( history_blocks(blockno)%ptr,                      
     &           (/hist_size,int_points,span/) )                                
c                                                                               
      if( do_a_block ) then                                                     
        do relem = 1, span                                                      
           elnum = felem + relem - 1  ! absolute element number                 
           one_elem_states(1:nrow_states) = zero                                
           call mm06_states_values_a                                            
           elem_states_output(1:nrow_states,relem) =                            
     &                one_elem_states(1:nrow_states)                            
        end do                                                                  
      else                                                                      
        relem = elnum + 1 - felem                                               
        one_elem_states(1:nrow_states) = zero                                   
        call mm06_states_values_a                                               
        elem_states_output(1:nrow_states,1) =                                   
     &                one_elem_states(1:nrow_states)                            
      end if                                                                    
c                                                                               
      deallocate( history_dump, one_elem_states )                               
c                                                                               
      return                                                                    
c                                                                               
 9050 format(10x,"block, felem, etype, mtype:  ",4i7,                           
     &  /,10x,   "int_pts, span, hist_size:    ",3i7 )                          
c                                                                               
      contains                                                                  
c     ========                                                                  
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm06_states_values_a              *          
c     *                                                              *          
c     *                    should get inlined                        *          
c     *                                                              *          
c     *                last modified : 9/24/2015 (rhd)               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_states_values_a                                           
c                                                                               
      implicit none                                                             
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: ipt                                                            
      double precision ::                                                       
     & local_states(3)                                                          
c                                                                               
      local_states(1:3) = zero                                                  
c                                                                               
      do ipt = 1, int_points                                                    
        local_states(1:3) = local_states(1:3)  +                                
     &                      history_dump(1:3,ipt,relem)                         
      end do                                                                    
c                                                                               
      one_elem_states(1:3) =  local_states(1:3) / dble(int_points)              
c                                                                               
      return                                                                    
c                                                                               
      end subroutine mm06_states_values_a                                       
      end subroutine mm06_states_values                                         
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine mm06_states_labels                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 9/24/2015 (rhd)                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm06_states_labels( size_state,                                
     &      num_states, state_labels, state_descriptors, out,                   
     &      comment_lines, max_comment_lines, num_comment_lines )               
c                                                                               
      implicit none                                                             
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: size_state, num_states, out, max_comment_lines,                
     &           num_comment_lines                                              
      character(len=8)  :: state_labels(size_state)                             
      character(len=60) :: state_descriptors(size_state)                        
      character(len=80) :: comment_lines(max_comment_lines)                     
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i                                                              
c                                                                               
      num_states = 3                                                            
      num_comment_lines = 0                                                     
c                                                                               
      state_labels(1) = "crp-eps"                                               
      state_labels(2) = "crp-rate"                                              
      state_labels(3) = "Uelas"                                                 
      state_descriptors(1) = "Eq creep eps"                                     
      state_descriptors(2) = "Eq creep eps dot"                                 
      state_descriptors(3) = "elastic energy density"                           
c                                                                               
      return                                                                    
      end subroutine mm06_states_labels                                         
                                                                                
c *******************************************************************           
c *                                                                 *           
c *        oumm06  material model # 6 -- creep                      *           
c *                                                                 *           
c *           set 3 material model dependent output values          *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine oumm06( gpn, mxvl, span, iout, elestr,                         
     &                   stress, history )                                      
      implicit none                                                             
c                                                                               
c                   parameter declarations                                      
c                   ----------------------                                      
c                                                                               
      integer ::  gpn, mxvl, span, iout                                         
c                                                                               
c                                                                               
      double precision ::                                                       
     & stress(mxvl,*), elestr(mxvl,*), history(mxvl,*)                          
c                                                                               
c               description of parameters                                       
c               -------------------------                                       
c                                                                               
c     gpn               : gauss point number being processed for block          
c     mxvl              : maximum no. elements per block                        
c     span              : number of elements in current block                   
c     iout              : write messages to this device number                  
c     stress            : current stresses for all                              
c                         elements in block for this gauss point                
c (*) elestr            : stresses to be output for elements                    
c     history           : current history values for all                        
c                         elements in block for this gauss point                
c                                                                               
c    (*)  values to be updated by this material model                           
c                                                                               
c                                                                               
c   stress ordering                elestr ordering                              
c     (1) sig-xx                   (1) sig-xx                                   
c     (2) sig-yy                   (2) sig-yy                                   
c     (3) sig-zz                   (3) sig-zz                                   
c     (4) tau-xy                   (4) tau-xy                                   
c     (5) tau-yz                   (5) tau-yz                                   
c     (6) tau-xz                   (6) tau-xz                                   
c     (7) total work density       (7) total work density                       
c                                  (8) mises equiv. stress                      
c                             (*)  (9) mat_val1                                 
c                             (*) (10) mat_val2                                 
c                             (*) (11) mat_val3                                 
c                                                                               
c  NOTE:  do not modify "stress" array or columns 1-8 of the "elestr"           
c         array. only modify columns 9-11 of "elestr". These are the            
c         3 "material model" dependent values output in the stress              
c         values for a gauss point. See Section 2.12 of manual and              
c         description of each material model. The output labels for             
c         columns 9-11 are "c1", "c2", "c3".                                    
c                                                                               
       integer :: i                                                             
c                                                                               
c   mat_va11 = c1 = accumulated, effective creep strain over time               
c   mat_va11 = c2 = current, effective creep strain rate                        
c   mat_va11 = c3 = elastic energy                                              
c                                                                               
c                                                                               
       do i = 1, span                                                           
         elestr(i,7)  = stress(i,7) ! energy density using deps -               
c                                     deps_{thermal}                            
         elestr(i,9)  = history(i,1)                                            
         elestr(i,10) = history(i,2)                                            
         elestr(i,11) = stress(i,9) ! elastic energy density using              
c                                     deps_{elastic}                            
       end do                                                                   
c                                                                               
       return                                                                   
       end subroutine oumm06                                                    
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *         subroutine cnst6  -- no longer used                  *          
c     *                                                              *          
c     *     tangent stiffness for power-law creep model              *          
c     *                                                              *          
c     *           last modified: 10/18/2015 rhd                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine cnst6( span, felem, gpn, iter, kout, mxvl,                     
     &                  nstr, weight, dtime, e_vec, nu_vec,                     
     &                  n_power_vec, rtse, block_props, history_np1,            
     &                  cep, det_jac_block )                                    
      implicit none                                                             
c                                                                               
c             parameter definitions                                             
c                                                                               
      integer :: span, felem, gpn, iter, kout, mxvl, nstr                       
      double precision                                                          
     &  weight, dtime, e_vec(mxvl), nu_vec(mxvl),                               
     &  n_power_vec(mxvl), rtse(mxvl,nstr), block_props(mxvl,*),                
     &  history_np1(span,*), cep(mxvl,6,6), det_jac_block(mxvl)                 
c                                                                               
c             local variables                                                   
c                                                                               
      integer :: i, j, k                                                        
      logical :: debug                                                          
      double precision                                                          
     &  zero, one, two, three, four, five, six, roothalf,                       
     &  statev(10), t_c, stress(6), local_d_mat(6,6),                           
     &  props(10), nrm_dtime, local_rtse(6), B                                  
c                                                                               
      data zero, one / 0.0d0, 1.0d0 /                                           
c                                                                               
      debug = felem .eq. 1  .and. gpn .eq. 8                                    
      debug = .false.                                                           
      if( debug ) write(kout,9010)  felem, gpn                                  
c                                                                               
c             loop over all elements in block for this integration              
c             point.                                                            
c                                                                               
      do i = 1, span                                                            
c                                                                               
      props(1) = e_vec(i)  ! emod                                               
      props(2) = nu_vec(i)                                                      
      props(3) = n_power_vec(i)                                                 
c                                                                               
c             the normalized cnst6n code does not know about B.                 
c             we normalize time to accommodate. strain rates                    
c             to-from cnst6n must be normalized time                            
c                                                                               
      B         = block_props(i,1)                                              
      t_c       = one /  B                                                      
      nrm_dtime = dtime / t_c                                                   
      statev(1) = history_np1(i,1)                                              
      statev(2) = history_np1(i,2) * t_c ! to normalized time                   
      statev(3) = history_np1(i,3)                                              
c                                                                               
      if( debug ) then                                                          
        write(kout,*) ' ..... before [D] computations ......'                   
        write(kout,9200) B                                                      
        write(kout,9210) dtime, t_c                                             
        write(kout,9220) nrm_dtime                                              
        write(kout,9310) statev(1), statev(2), statev(3)                        
      end if                                                                    
c                                                                               
c             deviator of trial elastic stress state                            
c                                                                               
      local_rtse(1:6) = rtse(i,1:6)                                             
c                                                                               
      call cnst6n( local_rtse, local_d_mat, statev, nrm_dtime,                  
     &             props, felem+i-1, gpn, iter, kout )                          
c                                                                               
      if( debug ) then                                                          
       do j = 1, 6                                                              
          write(kout,9400) j, (local_d_mat(j,k),k=1,6)                          
       end do                                                                   
      end if                                                                    
c                                                                               
c             include integration weight and det of [J]                         
c                                                                               
      cep(i,1:6,1:6) = local_d_mat(1:6,1:6) *                                   
     &                 weight * det_jac_block(i)                                
c                                                                               
      end do  ! over span                                                       
c                                                                               
      return                                                                    
c                                                                               
 9010 format(/,1x,'.....  entered cnst6 creep model, element, gpn: ',           
     &          2i7,/)                                                          
 9020 format(5x,i6,6e14.6)                                                      
 9200 format(5x,'B: ',d20.10)                                                   
 9210 format(5x,'dtime, t_c: ',3d15.6)                                          
 9220 format(5x,'nrm_dtime:  ',e14.6)                                           
 9310 format(5x,'statev :    ',3e15.8)                                          
 9400 format(3x,'row',i2,': ',6e14.6 )                                          
c                                                                               
      end subroutine cnst6                                                      
c                                                                               
