c *******************************************************************           
c *                                                                 *           
c *        material model # 2 --  deformation plasticity for power- *           
c *                               law following initial linear      *           
c *                               response with transition region   *           
c *                                                                 *           
c *                               small-strain analysis only        *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine mm02( step, iter, felem, gpn, e, nu, sigyld, exp,              
     &                 stress_n, stress_n1, strain, history,                    
     &                 history1, span, iout, signal_flag )                      
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                   parameter declarations                                      
c                                                                               
      double precision ::                                                       
     & e(*), nu(*), sigyld(*), exp(*), stress_n(mxvl,*),                        
     & stress_n1(mxvl,*), strain(mxvl,*), history(span,*),                      
     & history1(span,*)                                                         
       logical :: signal_flag                                                   
c                                                                               
c                   locally defined                                             
c                                                                               
      double precision ::                                                       
     & k1, k2(mxvl), l2, pi2, twothd, third, root2, rt23, zero, one,            
     & two, half, epsm, e1, e2, e3, epseff(mxvl),                               
     & epsyld(mxvl), epslim(mxvl), g, c, a, b,                                  
     & theta, c1, c2, epsnc(mxvl), signc(mxvl), rnc(mxvl),                      
     & sigeff(mxvl),  eps_elas(mxvl,6), shear_mod(mxvl)                         
      logical :: mm02es, debug, signal, nonlinear_flags(mxvl), test,            
     &           nonlinear_points, linear_flags(mxvl)                           
c                                                                               
      data pi2, twothd, third, root2, rt23                                      
     & / 1.570795d0, 0.66666667d0, 0.333333333d0, 1.4142135623730d0,            
     &   0.816496580927726d0 /                                                  
      data zero, one, two, k1, half                                             
     & / 0.0d0, 1.0d0, 2.0d0, 0.95d0, 0.5d0 /                                   
      data signal / .true. /                                                    
c                                                                               
c              model history data:                                              
c              (1) 0.0 or 1.0 for linear/yielding                               
c              (2) scalar effective stress                                      
c              (3) scalar effective strain                                      
c                                                                               
      debug = .false.                                                           
c                                                                               
      if( step .eq. 1 ) then                                                    
       do i = 1, span                                                           
         history(i,1:5) = zero                                                  
         history1(i,1:5) = zero                                                 
       end do                                                                   
      end if                                                                    
c                                                                               
c             get the effective strain for the total strain. check if           
c             we are still on the linear-reponse below the transition           
c             point. set logical linear vs. nonlinear flag for element.         
c                                                                               
      do i = 1, span                                                            
         epsm = ( strain(i,1) + strain(i,2) + strain(i,3) ) * third             
         e1   = strain(i,1) - epsm                                              
         e2   = strain(i,2) - epsm                                              
         e3   = strain(i,3) - epsm                                              
         epseff(i) = rt23 * sqrt( e1*e1 + e2*e2 + e3*e3 +                       
     &               half*strain(i,4)**2 + half*strain(i,5)**2 +                
     &               half*strain(i,6)**2 )                                      
         epsyld(i) = sigyld(i) / e(i)                                           
         epslim(i) = k1 * epsyld(i) * twothd * (one+nu(i))                      
         nonlinear_flags(i) =  epseff(i) .gt. epslim(i)                         
         linear_flags(i) = .not. nonlinear_flags(i)                             
         history1(i,4) = e(i)                                                   
         history1(i,5) = sigyld(i)                                              
      end do                                                                    
c                                                                               
c             do we have any nonlinear elements at this point                   
c             to process.                                                       
c                                                                               
      nonlinear_points = .false.                                                
      nonlin_point     = 0                                                      
c                                                                               
      do i = 1, span                                                            
         if( nonlinear_flags(i) ) then                                          
            nonlinear_points = .true.                                           
            if( nonlin_point .eq. 0 ) nonlin_point = i                          
         end if                                                                 
      end do                                                                    
c                                                                               
      if( debug ) then                                                          
        write(iout,*) ' '                                                       
        i = max( nonlin_point,1 )                                               
        write(iout,9000) felem+i-1, gpn, e(i), nu(i), sigyld(i),                
     &       exp(i), (strain(i,j),j=1,6),                                       
     &       (stress_n(i,j),j=1,7), history(i,1), history(i,2),                 
     &       history(i,3), epseff(i), epsyld(i), epslim(i)                      
      end if                                                                    
c                                                                               
c             process points still in linear response region                    
c                                                                               
      do i = 1, span                                                            
        if( nonlinear_flags(i) ) cycle                                          
        g = half * e(i) / ( one + nu(i) )                                       
        c = e(i) / ( ( one + nu(i) ) * ( one - two * nu(i) ) )                  
        a = c * ( one - nu(i) )                                                 
        b = c * nu(i)                                                           
        stress_n1(i,1) = a*strain(i,1) + b*(strain(i,2)+strain(i,3))            
        stress_n1(i,2) = b*(strain(i,1)+strain(i,3)) + a*strain(i,2)            
        stress_n1(i,3) = b*(strain(i,1)+strain(i,2)) + a*strain(i,3)            
        stress_n1(i,4) = g*strain(i,4)                                          
        stress_n1(i,5) = g*strain(i,5)                                          
        stress_n1(i,6) = g*strain(i,6)                                          
        stress_n1(i,7) = half * ( stress_n1(i,1)*strain(i,1) +                  
     &                   stress_n1(i,2)*strain(i,2) +                           
     &                   stress_n1(i,3)*strain(i,3) +                           
     &                   stress_n1(i,4)*strain(i,4) +                           
     &                   stress_n1(i,5)*strain(i,5) +                           
     &                   stress_n1(i,6)*strain(i,6) )                           
        history1(i,1) = zero                                                    
        history1(i,2) = zero                                                    
        history1(i,3) = epseff(i)                                               
      end do                                                                    
c                                                                               
      if( .not. nonlinear_points ) go to 1000                                   
c                                                                               
c                                                                               
c             process points beyond linear region. iteratively solve for        
c             effective stress given effective strain. solution for the         
c             effective stress requires a local newtion iteration               
c             done by function mm02es. Function is inlined but loop             
c             not vectorized.                                                   
c                                                                               
      do i = 1, span                                                            
        if( linear_flags(i) ) cycle                                             
        k2(i)  = root2 * ( one-k1 ) / exp(i) + one                              
        theta  = atan( k2(i)**(one-exp(i))/exp(i) )                             
        l2     = -tan(pi2-theta)                                                
        c1     = one / (l2 + one)                                               
        c2     = k2(i)**exp(i)                                                  
        epsnc(i) = ( two*k1 + l2*c2 - k2(i) ) * c1                              
        signc(i) = ( two*l2*k1 - l2*c2 + k2(i) ) * c1                           
        rnc(i)   = root2 * ( k1 - l2*k1 + l2*c2 - k2(i) ) * c1                  
      end do                                                                    
c                                                                               
      do i = 1, span                                                            
       if( linear_flags(i) ) cycle                                              
       test = mm02es( sigyld(i), epsyld(i), nu(i), exp(i),                      
     &                  epseff(i), sigeff(i), k2(i), epsnc(i),                  
     &                  signc(i), rnc(i) )                                      
       if( .not. test ) then                                                    
           write(iout,9250) felem+i-1, gpn                                      
           call die_abort                                                       
           stop                                                                 
       end if                                                                   
      end do                                                                    
c                                                                               
c             compute total stresses for total strain. write a                  
c             yielding message if that option is on.                            
c             get strain energy  density                                        
c                                                                               
      if( signal .and. signal_flag ) then                                       
       do i = 1, span                                                           
        if( nonlinear_flags(i) .and. history(i,1) .eq. zero ) then              
          write(iout,9300) felem+i-1, gpn, sigeff(i)-sigyld(i)*k1               
          history1(i,1) = one                                                   
        end if                                                                  
       end do                                                                   
      end if                                                                    
c                                                                               
      call mm02ns( rnc, signc, k2, k1, exp,                                     
     &             sigyld, epseff, sigeff, nu, e,                               
     &             strain, stress_n1, linear_flags, span )                      
c                                                                               
c             put effective stress and strain energy into stress                
c             vector for point. update history parameters.                      
c                                                                               
      do i = 1, span                                                            
       if( linear_flags(i) ) cycle                                              
       history1(i,1) = one                                                      
       history1(i,2) = sigeff(i)                                                
       history1(i,3) = epseff(i)                                                
      end do                                                                    
      if( debug ) then                                                          
        write(iout,*) ' '                                                       
        i = nonlin_point                                                        
        write(iout,9400) (stress_n1(i,j),j=1,7)                                 
      end if                                                                    
c                                                                               
c             compute plastic work density at gauss pts. and                    
c             the current plastic effective strain.                             
c                                                                               
 1000 continue                                                                  
      do i = 1, span                                                            
         shear_mod(i) = e(i) * half / (one + nu(i))                             
         eps_elas(i,1) = (stress_n1(i,1) - nu(i)*(stress_n1(i,2)+               
     &                   stress_n1(i,3)))/e(i)                                  
         eps_elas(i,2) = (stress_n1(i,2) - nu(i)*(stress_n1(i,1)+               
     &                   stress_n1(i,3)))/e(i)                                  
         eps_elas(i,3) = (stress_n1(i,3) - nu(i)*(stress_n1(i,1)+               
     &                   stress_n1(i,2)))/e(i)                                  
         eps_elas(i,4) = stress_n1(i,4) / shear_mod(i)                          
         eps_elas(i,5) = stress_n1(i,5) / shear_mod(i)                          
         eps_elas(i,6) = stress_n1(i,6) / shear_mod(i)                          
         stress_n1(i,8) =  stress_n1(i,7) -  half *                             
     &    (  eps_elas(i,1) * stress_n1(i,1)                                     
     &     + eps_elas(i,2) * stress_n1(i,2)                                     
     &     + eps_elas(i,3) * stress_n1(i,3)                                     
     &     + eps_elas(i,4) * stress_n1(i,4)                                     
     &     + eps_elas(i,5) * stress_n1(i,5)                                     
     &     + eps_elas(i,6) * stress_n1(i,6) )                                   
         stress_n1(i,9) = zero                                                  
         if( nonlinear_flags(i) )                                               
     &     stress_n1(i,9) = epsyld(i)*(sigeff(i)/sigyld(i))**exp(i)             
     &                       - sigeff(i) / e(i)                                 
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>> debug from mm02. elem, gpn : ',i8,i2,                          
     & /,    '    e, nu, sigyld, exp : ',4f10.2,                                
     & /,    '    strains @ n+1:',                                              
     & /,10x,3e15.6,/,10x,3e15.6,                                               
     & /,    '    stresses @ n :',                                              
     & /,10x,3e15.6,/,10x,4e15.6,                                               
     & /,    '    history : ',f4.1,f10.3,e14.4,                                 
     & /,    '    epseff, epsyld, epslim: ',3e14.6 )                            
 9200 format('    updated elastic stresses: ',                                  
     & /,10x,3e15.6,/,10x,4e15.6 )                                              
 9250 format(//,'>>> Fatal error in mm02',                                      
     &        /,'    mm02es failed to converge for',                            
     &        /,'    element, point: ',i8,i2,                                   
     &        /,'    job terminated....' )                                      
 9300 format(10x,i8,i3,' point yields.  f = ',f12.3 )                           
 9400 format('    updated elastic-plastic stresses: ',                          
     & /,10x,3e15.6,/,10x,4e15.6 )                                              
c                                                                               
      end                                                                       
c ********************************************************************          
c *                                                                  *          
c *   mm02es -- find effective stress for transitional or power-law  *          
c *             portions of uniaxial response curves                 *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      logical function mm02es( sigyld, epsyld, nu, exp, epseff, sigeff,         
     &                         k2, epsnc, signc, rnc )                          
      implicit double precision ( a-h,o-z)                                      
c                                                                               
      double precision ::                                                       
     &   nu, k2                                                                 
      data  maxitr / 30 /, toler / 0.0001d0 /                                   
      data  twothd, third, one, two, half, zero                                 
     &  / 0.666666667d0, 0.333333333d0, 1.0d0, 2.0d0, 0.5d0, 0.0d0 /            
c                                                                               
c             use newton-raphson iterative solution to find the                 
c             effective stress given the effective strain for                   
c             a multi-dimensional deformation model. the uniaxial               
c             curve is a linear region, a circular transition                   
c             region and a power-law region.                                    
c                                                                               
c             sigyld  -- yield stress                                           
c             epsyld  -- yield strain corresponding to sigyld                   
c             nu      -- poisson's ratio                                        
c             exp     -- power-law exponent                                     
c             epseff  -- effective uniaxial strain (given)                      
c             sigeff  -- effective uniaxial stress (computed)                   
c             mm02es  -- .true. if procedure converges                          
c             k2      -- upper point on transition                              
c             epsnc                                                             
c             signc   -- center of circle for transition (given)                
c             rnc     -- radius of transition circle (given)                    
c                                                                               
c                                                                               
      iterno = 1                                                                
      mm02es = .false.                                                          
      expm1  = exp - one                                                        
      trnlim = ( k2**expm1 + twothd * ( nu - half ) ) * k2 * epsyld             
      epsn   = epseff / epsyld                                                  
c                                                                               
c             handle circular transition and power-law in                       
c             separate iteration loops                                          
c                                                                               
      if( epseff .ge. trnlim ) go to 200                                        
c                                                                               
c             perform newton raphson on circular transition part of             
c             response to compute effective stress. starting guess for          
c             normalized effective stress is k2 point on transition curve.      
c             for very low hardening, e.g. n=20, with an effective strain       
c             very close to the upper transition, the iteration process         
c             will blow-up. this is detected by checking the sign of            
c             radical. if a blow-up is detected, simply switch to the power-law 
c             part of curve and convergence will be quick!                      
c                                                                               
      const = ( two * nu - one ) * third                                        
      sign  = k2                                                                
      radical = sqrt( rnc*rnc - (sign - signc)**2 )                             
      resid   = epsn - epsnc + radical - const * sign                           
 100  continue                                                                  
      dfdsn = (signc - sign ) / radical - const                                 
      signew = sign - resid / dfdsn                                             
      radical = rnc*rnc - (signew - signc)**2                                   
      if( radical .lt. zero ) then                                              
         iterno = 1                                                             
         go to 200                                                              
      end if                                                                    
      radical = sqrt( radical )                                                 
      resid   = epsn - epsnc + radical - const * signew                         
      if( abs(signew-sign) .le. toler ) then                                    
        sigeff = signew * sigyld                                                
      else                                                                      
        sign = signew                                                           
        iterno = iterno + 1                                                     
        if ( iterno .gt. maxitr ) return                                        
        go to 100                                                               
      end if                                                                    
      mm02es = .true.                                                           
      return                                                                    
c                                                                               
c             perform newton raphson on power-law part of                       
c             response to compute effective stress. starting guess for          
c             normalized effective stress is n th root of normailzed            
c             effective strain (same approach used for regular deformation      
c             plasticity model).                                                
c                                                                               
 200  continue                                                                  
      sign  = epsn**(one/exp)                                                   
      c     = twothd * ( nu - half )                                            
      resid = sign**exp + c*sign - epsn                                         
 210  continue                                                                  
      dfdsn  = exp*sign**expm1 + c                                              
      signew = sign - resid / dfdsn                                             
      resid  = signew**exp + c*signew - epsn                                    
      if( abs(signew-sign) .le. toler ) then                                    
        sigeff = signew * sigyld                                                
      else                                                                      
        sign = signew                                                           
        iterno = iterno + 1                                                     
        if( iterno .gt. maxitr ) return                                         
        go to 210                                                               
      end if                                                                    
      mm02es = .true.                                                           
      return                                                                    
c                                                                               
      end                                                                       
c ********************************************************************          
c *                                                                  *          
c *   mm02ns -- new stresses and strain energy density               *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine mm02ns( rnc, signc, k2, k1, exp,                               
     &                   sigyld, epseff, sigeff, nu, e,                         
     &                   strain, sig, linear_flags, span )                      
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c             parameter declarations                                            
c                                                                               
      double precision ::                                                       
     & rnc(*), signc(*), k2(*), k1, exp(*), sigyld(*),                          
     & epseff(*), sigeff(*), nu(*), e(*), strain(mxvl,*), sig(mxvl,*)           
      logical :: linear_flags(*)                                                
c                                                                               
c             local declarations                                                
c                                                                               
      double precision ::                                                       
     & one, two, three, half, onefive, six, pt25, twothd,                       
     & third, const1, const2, epsv, sigm, upress(mxvl),                         
     & sigk1, uelask1(mxvl), sncso(mxvl), rncso(mxvl),                          
     & c1(mxvl), rncsq(mxvl), betak1, uk1(mxvl), term1,                         
     & term2, term3, betase, uk2, sigk2, betak2, expp1, c2,                     
     & term4, term5, term6, term7                                               
c                                                                               
      data one, two, three, half, onefive, six, pt25                            
     &   / 1.0d0, 2.0d0, 3.0d0, 0.5d0, 1.5d0, 6.0d0, 0.25d0 /                   
      data twothd, third                                                        
     &  / 0.6666666667d0, 0.33333333333d0 /                                     
c                                                                               
c             compute total stress given total strain and                       
c             effective strain. also the total work density.                    
c                                                                               
c                                                                               
c             compute stress by volumetric and deviatoric                       
c             components.                                                       
c                                                                               
      do i = 1, span                                                            
        if( linear_flags(i) ) cycle                                             
        const1 = twothd * sigeff(i) / epseff(i)                                 
        const2 = e(i) / ( one - two*nu(i) ) - const1                            
        epsv   = ( strain(i,1) + strain(i,2) + strain(i,3) ) * third            
        sig(i,1) = const1 * strain(i,1) + const2 * epsv                         
        sig(i,2) = const1 * strain(i,2) + const2 * epsv                         
        sig(i,3) = const1 * strain(i,3) + const2 * epsv                         
        sig(i,4) = const1 * strain(i,4) * half                                  
        sig(i,5) = const1 * strain(i,5) * half                                  
        sig(i,6) = const1 * strain(i,6) * half                                  
      end do                                                                    
c                                                                               
c             compute exact strain energy density. we have a point on           
c             transition part of curve or power-law part of curve.              
c             energy up to k1*sigyld is same for both cases.                    
c                                                                               
      do i = 1, span                                                            
       if( linear_flags(i) ) cycle                                              
        sigm      = third * ( sig(i,1) + sig(i,2) + sig(i,3) )                  
        upress(i) = ( onefive - three*nu(i)) * sigm * sigm                      
        sigk1     = k1 * sigyld(i)                                              
        uelask1(i)= third * (one+nu(i)) * sigk1 * sigk1                         
        sncso(i)  = signc(i) * sigyld(i)                                        
        rncso(i)  = rnc(i) * sigyld(i)                                          
        c1(i)     = ( two*nu(i) -one ) / six                                    
        rncsq(i)  = rncso(i) * rncso(i)                                         
        betak1    = asin( (sigk1 - sncso(i))/rncso(i) )                         
        term1     = c1(i) * sigk1 * sigk1                                       
        term2     = -sigk1 * sqrt( rncsq(i) - (sigk1-sncso(i))**2 )             
        term3     = rncsq(i) * ( half*betak1 + pt25*sin(two*betak1) )           
        uk1(i) = term1 + term2 + term3                                          
      end do                                                                    
c                                                                               
      do i = 1, span                                                            
       if( linear_flags(i) ) cycle                                              
       if( sigeff(i)/sigyld(i) .lt. k2(i) ) then                                
c                                                                               
c                 Stress point in transition:                                   
c                   1) pressure component (total)                               
c                   2) deviatoric on linear curve up to k1*sigyld               
c                   3) deviatoric transition from k1*sigyld to sigeff.          
c                                                                               
        betase   = asin( (sigeff(i) - sncso(i))/rncso(i) )                      
        term1    = c1(i) * sigeff(i) * sigeff(i)                                
        term2    = -sigeff(i) * sqrt( rncsq(i) -                                
     &              (sigeff(i)-sncso(i))**2 )                                   
        term3    = rncsq(i) * ( half*betase + pt25*sin(two*betase) )            
        uk2      = term1 + term2 + term3                                        
        sig(i,7) = ( upress(i) + uelask1(i) + (uk2-uk1(i)) ) / e(i)             
c                                                                               
       else                                                                     
c                 Stress point in power-law region:                             
c                   1) pressure component (total)                               
c                   2) deviatoric on linear curve up to k1*sigyld               
c                   3) deviatoric transition from k1*sigyld to k2*sigyld        
c                   4) deviatoric power-law from k2*sigyld to sigeff            
c                                                                               
        sigk2    = k2(i) * sigyld(i)                                            
        betak2   = asin( (sigk2 - sncso(i))/rncso(i) )                          
        term1    = c1(i) * sigk2 * sigk2                                        
        term2    = -sigk2 * sqrt( rncsq(i) - (sigk2-sncso(i))**2 )              
        term3    = rncsq(i) * ( half*betak2 + pt25*sin(two*betak2) )            
        uk2      = term1 + term2 + term3                                        
c                                                                               
        expp1    = exp(i) + one                                                 
        c2       = (one/sigyld(i)**(exp(i)-one)) * exp(i) / expp1               
        term4    = c1(i) * sigk2 * sigk2                                        
        term5    = c2 * sigk2**expp1                                            
        term6    = c1(i) * sigeff(i) * sigeff(i)                                
        term7    = c2 * sigeff(i)**expp1                                        
        sig(i,7) = ( upress(i) + uelask1(i) + (uk2-uk1(i)) +                    
     &              (term6+term7) - (term5+term4) ) / e(i)                      
c                                                                               
        end if                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm02_set_sizes                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 12/14/14 rhd                *          
c     *                                                              *          
c     *    called by warp3d for each material model to obtain        *          
c     *    various sizes of data for the model                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm02_set_sizes( info_vector )                                  
      dimension info_vector(*)                                                  
c                                                                               
c        set infor_data                                                         
c                                                                               
c         1        number of history values per integration                     
c                  point. Abaqus calles these "statev". Values                  
c                  double or single precsion based on hardware.                 
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
c         4        number of state variables per point to be output             
c                  when user requests this type of results                      
c                                                                               
      info_vector(1) = 5                                                        
      info_vector(2) = 21                                                       
      info_vector(3) = 0                                                        
      info_vector(4) = 3                                                        
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine mm02_states_values                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/9/15 (rhd)                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm02_states_values( itype, elem_states_output,                 
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
c           build deformation plasticity states values output.                  
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
           call mm02_states_values_a                                            
           elem_states_output(1:nrow_states,relem) =                            
     &                one_elem_states(1:nrow_states)                            
        end do                                                                  
      else                                                                      
        relem = elnum + 1 - felem                                               
        one_elem_states(1:nrow_states) = zero                                   
        call mm02_states_values_a                                               
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
c     *                 subroutine mm02_states_values_a              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 1/9/15 (rhd)               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm02_states_values_a                                           
c                                                                               
      implicit none                                                             
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: ipt                                                            
      double precision ::                                                       
     & nl_flag, sig_eff, eps_eff, ymod, sigyld, epsyld                          
c                                                                               
      nl_flag = zero                                                            
      sig_eff = zero                                                            
      eps_eff = zero                                                            
c                                                                               
      do ipt = 1, int_points                                                    
        ymod    = history_dump(4,ipt,relem)                                     
        sigyld  = history_dump(5,ipt,relem)                                     
        epsyld  = sigyld / ymod                                                 
        nl_flag = nl_flag + history_dump(1,ipt,relem)                           
        sig_eff = sig_eff + history_dump(2,ipt,relem) / sigyld                  
        eps_eff = eps_eff + history_dump(3,ipt,relem) / epsyld                  
                                                                                
      end do                                                                    
c                                                                               
      one_elem_states(1) = nl_flag / dble(int_points)                           
      one_elem_states(2) = sig_eff / dble(int_points)                           
      one_elem_states(3) = eps_eff / dble(int_points)                           
c                                                                               
      return                                                                    
c                                                                               
      end subroutine mm02_states_values_a                                       
      end subroutine mm02_states_values                                         
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm02_states_labels                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/11/2015 (rhd)                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm02_states_labels( size_state,                                
     &      num_states, state_labels, state_descriptors, out,                   
     &      comment_lines, max_comment_lines, num_comment_lines )               
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
      logical, save :: do_print = .false.                                       
c                                                                               
      num_states = 3                                                            
      state_labels(1) = "nl flag"                                               
      state_labels(2) = "seff/s0"                                               
      state_labels(3) = "eff/e0"                                                
      state_descriptors(1) = "=0 all pts linear, =1 all pts nonlinear"          
      state_descriptors(2) = "mises / sig_yld"                                  
      state_descriptors(3) = "eps-eff / eps_yld"                                
c                                                                               
      num_comment_lines = 0                                                     
c                                                                               
      if( do_print ) then                                                       
        do i = 1, 3                                                             
          write(out,9010) i, state_labels(i), state_descriptors(i)              
        end do                                                                  
        do_print = .false.                                                      
      end if                                                                    
c                                                                               
      return                                                                    
 9010 format(2x,i3,2x,a8,2x,a)                                                  
      end                                                                       
c *********************************************************************         
c *                                                                   *         
c *   cnst2 -- tangent modulus matrix for nonlinear elastic model     *         
c *            vectorized version                                     *         
c *                                                                   *         
c *********************************************************************         
c                                                                               
      subroutine cnst2( felem, gpn, e, nu, sigyld, exp, strain_n1,              
     &                  history, cep, span, iout )                              
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                   parameter declarations                                      
c                                                                               
      integer :: felem, gpn, span, iout                                         
      double precision ::                                                       
     & strain_n1(mxvl,*), history(span,*), cep(mxvl,6,6), e(*),                 
     & nu(*), sigyld(*), exp(*)                                                 
c                                                                               
c                   locally defined                                             
c                                                                               
      integer :: nonlin_point, i, j                                             
      double precision ::                                                       
     & k1, k2(mxvl), epsyld(mxvl), epslim(mxvl), sigeff(mxvl),                  
     & epseff(mxvl), l2, c1, c2, c3, c4, theta, pi2,                            
     & root2, zero, one, two, three,                                            
     & epsnc(mxvl), signc(mxvl), rnc(mxvl)                                      
      logical :: debug, nonlinear_points, nonlinear_flags(mxvl)                 
c                                                                               
      data pi2, root2                                                           
     & / 1.570795d00, 1.4142135623730d00  /                                     
      data zero, one, two, k1, three                                            
     & / 0.0d00, 1.0d00, 2.0d00, 0.95d00, 3.0d00 /                              
c                                                                               
c              pull out model properties                                        
c                                                                               
      nonlinear_points = .false.                                                
      nonlin_point     = 0                                                      
      do i = 1, span                                                            
        epsyld(i) = sigyld(i) / e(i)                                            
        k2(i)     = root2 * ( one-k1 ) / exp(i) + one                           
        epslim(i) = k1 * epsyld(i) * (two/three) * (one+nu(i))                  
        sigeff(i) = history(i,2)                                                
        epseff(i) = history(i,3)                                                
      end do                                                                    
      do i = 1, span                                                            
        nonlinear_flags(i) = .false.                                            
        if( epseff(i) .gt. epslim(i) ) then                                     
           nonlinear_points = .true.                                            
           if( nonlin_point .eq. 0 ) nonlin_point = i                           
           nonlinear_flags(i) = .true.                                          
        end if                                                                  
      end do                                                                    
c                                                                               
      debug = .false.                                                           
      if( debug ) then                                                          
        write(iout,*) ' '                                                       
        i = nonlin_point                                                        
        write(iout,9000) felem+i-1, gpn, e(i), nu(i), sigyld(i), exp(i),        
     &                   (strain_n1(i,j),j=1,6), (history(i,j),j=1,3)           
      end if                                                                    
c                                                                               
      do i = 1, span  ! linear elastic points                                   
       if( nonlinear_flags(i) ) cycle                                           
       cep(i,1,4) = zero                                                        
       cep(i,1,5) = zero                                                        
       cep(i,1,6) = zero                                                        
       cep(i,2,4) = zero                                                        
       cep(i,2,5) = zero                                                        
       cep(i,2,6) = zero                                                        
       cep(i,3,4) = zero                                                        
       cep(i,3,5) = zero                                                        
       cep(i,3,6) = zero                                                        
       cep(i,4,1) = zero                                                        
       cep(i,4,2) = zero                                                        
       cep(i,4,3) = zero                                                        
       cep(i,4,5) = zero                                                        
       cep(i,4,6) = zero                                                        
       cep(i,5,1) = zero                                                        
       cep(i,5,2) = zero                                                        
       cep(i,5,3) = zero                                                        
       cep(i,5,4) = zero                                                        
       cep(i,5,6) = zero                                                        
       cep(i,6,1) = zero                                                        
       cep(i,6,2) = zero                                                        
       cep(i,6,3) = zero                                                        
       cep(i,6,4) = zero                                                        
       cep(i,6,5) = zero                                                        
       c1 = (e(i)/((one+nu(i))*(one-two*nu(i))))                                
       c2 = (one-nu(i))*c1                                                      
       c3 = ((one-two*nu(i))/two)*c1                                            
       c4 = nu(i)*c1                                                            
       cep(i,1,1)= c2                                                           
       cep(i,2,2)= c2                                                           
       cep(i,3,3)= c2                                                           
       cep(i,4,4)= c3                                                           
       cep(i,5,5)= c3                                                           
       cep(i,6,6)= c3                                                           
       cep(i,1,2)= c4                                                           
       cep(i,1,3)= c4                                                           
       cep(i,2,1)= c4                                                           
       cep(i,3,1)= c4                                                           
       cep(i,2,3)= c4                                                           
       cep(i,3,2)= c4                                                           
      end do                                                                    
c                                                                               
      if( .not. nonlinear_points ) return                                       
c                                                                               
c             some points are nonlinear                                         
c                                                                               
      do i = 1, span                                                            
        if( .not. nonlinear_flags(i) ) cycle                                    
        theta    = atan( k2(i)**(1.-exp(i))/exp(i) )                            
        l2       = -tan(pi2-theta)                                              
        c1       = one / (l2 + one)                                             
        c2       = k2(i)**exp(i)                                                
        epsnc(i) = ( two*k1 + l2*c2 - k2(i) ) * c1                              
        signc(i) = ( two*l2*k1 - l2*c2 + k2(i) ) * c1                           
        rnc(i)   = root2 * ( k1 - l2*k1 + l2*c2 - k2(i) ) * c1                  
      end do                                                                    
c                                                                               
      call cnst2a( cep, e, nu, sigeff, epseff, strain_n1,                       
     &             sigyld, exp, signc, rnc, k2,                                 
     &             span, nonlinear_flags )                                      
c                                                                               
      if( .not. debug ) return                                                  
      write(iout,9034)                                                          
      i = nonlin_point                                                          
      do j = 1, 6                                                               
        write(iout,9008) cep(i,j,1), cep(i,j,2), cep(i,j,3),                    
     &                     cep(i,j,4), cep(i,j,5), cep(i,j,6)                   
      end do                                                                    
      write(iout,*) ' '                                                         
      return                                                                    
c                                                                               
 9000 format('>> debug from cnst2. elem, gpn : ',i8,i2,                         
     & /,    '    e, nu, sigyld, exp : ',4f10.2,                                
     & /,    '    strains @ n+1 :',                                             
     & /,10x,3e15.6,/,10x,3e15.6,                                               
     & /,    '    history : ',f4.1,f10.3,f10.6 )                                
 9008 format(3x,6e14.4)                                                         
 9020 format(//,5x,'>> update [d] element, point: ',2i6)                        
 9034 format('    elasto-plastic d matrix: ' )                                  
 9038 format('    point is elastic.  epseff, sigeff: ',e15.6,                   
     & f10.3 )                                                                  
c                                                                               
      end                                                                       
c ********************************************************************          
c *                                                                  *          
c *   cnst2a -- tangent modulus matrix                               *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine cnst2a( cep, e, nu, sigeff, epseff, strain, sigyld,            
     &                   exp, signc, rnc, k2, span, nonlinear_flags )           
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
      integer :: span                                                           
      double precision ::                                                       
     & cep(mxvl,6,6), e(*), nu(*), sigeff(*), epseff(*),                        
     & strain(mxvl,*), sigyld(*), exp(*), signc(*),                             
     & rnc(*), k2(*)                                                            
      logical ::  nonlinear_flags(*)                                            
c                                                                               
      integer :: i                                                              
      double precision ::                                                       
     & ev(mxvl), e1(mxvl), e2(mxvl), e3(mxvl),                                  
     & e4(mxvl), e5(mxvl), e6(mxvl), sign(mxvl), radical,                       
     & stiff, const, g(mxvl), expm1, c1,                                        
     & c2(mxvl), c3(mxvl), temp, denom,                                         
     & twothd, fnine, third, half, one, two                                     
c                                                                               
      data  twothd, fnine, third, half, one, two                                
     & / 0.6666666667d00, 0.444444444444444d00, 0.3333333333333d00,             
     &    0.5d00, 1.0d00, 2.0d00 /                                              
c                                                                               
c             compute the tangent modulus matrix                                
c                                                                               
      do i = 1, span                                                            
       if( .not. nonlinear_flags(i) ) cycle                                     
         ev(i) = third * ( strain(i,1) + strain(i,2) + strain(i,3) )            
         e1(i) = strain(i,1) - ev(i)                                            
         e2(i) = strain(i,2) - ev(i)                                            
         e3(i) = strain(i,3) - ev(i)                                            
         e4(i) = strain(i,4) * half                                             
         e5(i) = strain(i,5) * half                                             
         e6(i) = strain(i,6) * half                                             
         sign(i) = sigeff(i) / sigyld(i)                                        
      end do                                                                    
c                                                                               
      do i = 1, span                                                            
        if( .not. nonlinear_flags(i) ) cycle                                    
        if( sign(i) .le. k2(i) ) then                                           
          radical = sqrt( rnc(i)*rnc(i) - (sign(i)-signc(i))**2 )               
          stiff   = (sign(i) - signc(i))/radical -                              
     &                twothd*(half-nu(i))                                       
          const   = one/epseff(i) - e(i)/sigeff(i)/stiff                        
          g(i)    = -fnine * sigeff(i) * const / epseff(i)                      
     &                /epseff(i)                                                
        else                                                                    
          expm1 = exp(i) - one                                                  
          c1    = sign(i)**expm1                                                
          denom = twothd * ( -half + nu(i) ) + exp(i) * c1                      
          temp  = -fnine * sigeff(i) * expm1 * c1 / epseff(i)                   
     &              / epseff(i) /epseff(i)                                      
          g(i)  = temp / denom                                                  
        end if                                                                  
        c2(i) = twothd * sigeff(i) / epseff(i)                                  
        c3(i) = third * ( -c2(i) + e(i)/(one-two*nu(i)) )                       
      end do                                                                    
c                                                                               
c             3-D (row ordering x, y, z, xy, yz, xz)                            
c                                                                               
      do i = 1, span                                                            
       if( .not. nonlinear_flags(i) ) cycle                                     
       cep(i,1,1) = ( c3(i) + g(i)*e1(i)*e1(i) + c2(i) )                        
       cep(i,2,1) = ( c3(i) + g(i)*e2(i)*e1(i) )                                
       cep(i,3,1) = ( c3(i) + g(i)*e3(i)*e1(i) )                                
       cep(i,2,2) = ( c3(i) + g(i)*e2(i)*e2(i) + c2(i) )                        
       cep(i,3,2) = ( c3(i) + g(i)*e3(i)*e2(i) )                                
       cep(i,3,3) = ( c3(i) + g(i)*e3(i)*e3(i) + c2(i) )                        
       cep(i,1,2) = cep(i,2,1)                                                  
       cep(i,1,3) = cep(i,3,1)                                                  
       cep(i,2,3) = cep(i,3,2)                                                  
       cep(i,1,4) = ( g(i)*e1(i)*e4(i) )                                        
       cep(i,1,5) = ( g(i)*e1(i)*e5(i) )                                        
       cep(i,1,6) = ( g(i)*e1(i)*e6(i) )                                        
       cep(i,2,4) = ( g(i)*e2(i)*e4(i) )                                        
       cep(i,2,5) = ( g(i)*e2(i)*e5(i) )                                        
       cep(i,2,6) = ( g(i)*e2(i)*e6(i) )                                        
       cep(i,3,4) = ( g(i)*e3(i)*e4(i) )                                        
       cep(i,3,5) = ( g(i)*e3(i)*e5(i) )                                        
       cep(i,3,6) = ( g(i)*e3(i)*e6(i) )                                        
       cep(i,4,1) = cep(i,1,4)                                                  
       cep(i,5,1) = cep(i,1,5)                                                  
       cep(i,6,1) = cep(i,1,6)                                                  
       cep(i,4,2) = cep(i,2,4)                                                  
       cep(i,5,2) = cep(i,2,5)                                                  
       cep(i,6,2) = cep(i,2,6)                                                  
       cep(i,4,3) = cep(i,3,4)                                                  
       cep(i,5,3) = cep(i,3,5)                                                  
       cep(i,6,3) = cep(i,3,6)                                                  
       cep(i,4,4) = ( g(i)*e4(i)*e4(i) + c2(i)*half )                           
       cep(i,4,5) = ( g(i)*e4(i)*e5(i) )                                        
       cep(i,4,6) = ( g(i)*e4(i)*e6(i) )                                        
       cep(i,5,5) = ( g(i)*e5(i)*e5(i) + c2(i)*half )                           
       cep(i,5,6) = ( g(i)*e5(i)*e6(i) )                                        
       cep(i,6,6) = ( g(i)*e6(i)*e6(i) + c2(i)*half )                           
       cep(i,5,4) = cep(i,4,5)                                                  
       cep(i,6,4) = cep(i,4,6)                                                  
       cep(i,6,5) = cep(i,5,6)                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
