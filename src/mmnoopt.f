c                                                                               
c     Various annoying functions which the Intel compiler breaks with           
c     optimizations on.                                                         
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_solve_bicrystal              *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/11/14                     *          
c     *                                                              *          
c     *    Solve the bicrystal model for a pair                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_solve_bicrystal(props1, np11, n1, props2,                 
     &            np12, n2, sv, lv, tv, cut, iout)                              
      use mm10_defs                                                             
      implicit none                                                             
c                                                                               
      type(crystal_props) :: props1, props2                                     
      type(crystal_state) :: np11, n1, np12, n2                                 
      double precision, dimension(3) :: sv, lv, tv                              
      logical :: cut                                                            
      integer :: iout                                                           
c                                                                               
      double precision, dimension(12) :: x, R, dx, xnew                         
      double precision, dimension(12,12) :: J                                   
      double precision, dimension(3,6) :: Pe, Pc                                
      double precision, dimension(6) :: d_main, d_old, d_rat                    
      double precision, dimension(6,6) :: I                                     
      double precision, dimension(9) :: flatrot                                 
      double precision :: dt, temp                                              
      integer :: step, elem, gp                                                 
      integer :: iter, maxiter, info, k, ls, mls                                
      integer, dimension(12) :: ipiv                                            
      double precision :: nr, inr, tol, atol, c, alpha, red, nlsx, nrs,         
     &      ls1, ls2                                                            
      logical :: aflag, debug, force_atang                                      
c                                                                               
      parameter(maxiter = 250)                                                  
      parameter(tol = 1.0D-6)                                                   
      parameter(atol = 1.0D-10)                                                 
      parameter(debug = .false.)                                                
      parameter(red = 0.5)                                                      
      parameter(c = 1.0e-4)                                                     
      parameter(mls = 10)                                                       
c                                                                               
      force_atang = .true.                                                      
c                                                                               
      flatrot = reshape(np11%R, (/9/))                                          
      dt = np11%tinc                                                            
      temp = np11%temp                                                          
      step = np11%step                                                          
      elem = np11%elem                                                          
      gp = np11%gp                                                              
c                                                                               
      I = 0.0                                                                   
      I(1,1) = 1.0                                                              
      I(2,2) = 1.0                                                              
      I(3,3) = 1.0                                                              
      I(4,4) = 1.0                                                              
      I(5,5) = 1.0                                                              
      I(6,6) = 1.0                                                              
c                                                                               
c           Form our projections                                                
      Pe = 0.0                                                                  
      Pc = 0.0                                                                  
      call mm11_pe(sv, Pe)                                                      
      call mm11_pc(lv, tv, Pc)                                                  
c                                                                               
c           Extract the macroscale deformation from the CP structs,             
c           form the initial guess                                              
      d_main = np11%D                                                           
      d_old = 0.5*(n1%D+n2%D)                                                   
      d_rat = 0.0                                                               
      do k=1,6                                                                  
        if (dabs(d_old(k)) .lt. 1.0e-15) then                                   
              d_rat(k) = 1.0                                                    
        else                                                                    
              d_rat(k) = n1%D(k)/d_old(k)                                       
        end if                                                                  
      end do                                                                    
      x(1:6) = np11%D*d_rat                                                     
      d_rat = 0.0                                                               
      do k=1,6                                                                  
        if (dabs(d_old(k)) .lt. 1.0e-15) then                                   
              d_rat(k) = 1.0                                                    
        else                                                                    
              d_rat(k) = n2%D(k)/d_old(k)                                       
        end if                                                                  
      end do                                                                    
      x(7:12) = np12%D*d_rat                                                    
c           Loop and solve                                                      
c                                                                               
      cut = .false.                                                             
      iter = 0                                                                  
c                                                                               
      if (debug) write(*,*) "Starting bicrystal iterations"                     
                                                                                
c           Just form everything inline...                                      
      call mm10_solve_crystal(props1, np11, n1, cut, iout, force_atang)         
      if (cut) then                                                             
        write(*,*) "CP stress update failed."                                   
        return                                                                  
      end if                                                                    
      call mm10_solve_crystal(props2, np12, n2, cut, iout, force_atang)         
      if (cut) then                                                             
        write(*,*) "CP stress update failed."                                   
        return                                                                  
      end if                                                                    
                                                                                
      R(1:3)  = matmul(Pe, np11%stress - np12%stress)                           
      R(4:6)  = matmul(Pc, np11%D - np12%D)                                     
      R(7:12) = 0.5*(np11%D + np12%D) - d_main                                  
c                                                                               
      nr = dsqrt(dot_product(R,R))                                              
      inr = nr                                                                  
                                                                                
      if (debug) write(*,'("Iter ", i2," norm ",E10.3)') iter, inr              
                                                                                
      do                                                                        
        if ((nr/inr .lt. tol) .or. (nr .lt. atol)) exit                         
c                                                                               
        J(1:3, 1:6)  =  matmul(Pe, np11%tangent)                                
        J(1:3, 7:12) = -matmul(Pe, np12%tangent)                                
        J(4:6, 1:6)  =  matmul(Pc, I)                                           
        J(4:6, 7:12) = -matmul(Pc, I)                                           
        J(7:12,1:6)  = 0.5*I                                                    
        J(7:12,7:12) = 0.5*I                                                    
c                                                                               
        dx = R                                                                  
        call DGESV(12, 1, -J, 12, ipiv, dx, 12, info)                           
c                                                                               
c           Line search                                                         
        alpha = 1.0                                                             
        ls1 = 0.5*dot_product(R,R)                                              
        ls2 = c*dot_product(dx,matmul(transpose(J),R))                          
        ls = 0                                                                  
        do                                                                      
          nlsx = ls1 + ls2*alpha                                                
          xnew = x + alpha*dx                                                   
c                                                                               
          call mm10_setup_np1(flatrot, xnew(1:6), dt, temp, step,               
     &         elem, gp, np11)                                                  
          call mm10_setup_np1(flatrot, xnew(7:12), dt, temp, step,              
     &         elem, gp, np12)                                                  
c           Just form everything inline...                                      
          call mm10_solve_crystal(props1, np11, n1, cut, iout,                  
     &      force_atang)                                                        
          if (cut) then                                                         
            write(*,*) "CP stress update failed."                               
            return                                                              
          end if                                                                
          call mm10_solve_crystal(props2, np12, n2, cut, iout,                  
     &      force_atang)                                                        
          if (cut) then                                                         
            write(*,*) "CP stress update failed."                               
            return                                                              
          end if                                                                
                                                                                
          R(1:3)  = matmul(Pe, np11%stress - np12%stress)                       
          R(4:6)  = matmul(Pc, np11%D - np12%D)                                 
          R(7:12) = 0.5*(np11%D + np12%D) - d_main                              
c                                                                               
          nr = dsqrt(dot_product(R,R))                                          
          nrs = 0.5*dot_product(R,R)                                            
c                                                                               
          if ((nrs .le. nlsx) .or. (ls .gt. mls)) then                          
            x = xnew                                                            
            exit                                                                
          else                                                                  
            alpha = red*alpha                                                   
            ls = ls + 1                                                         
          end if                                                                
                                                                                
        end do                                                                  
c                                                                               
        iter = iter + 1                                                         
        if (debug) write(*,'("Iter ", i2," norm ",E10.3," ls",F10.3)'           
     &      ) iter, nr, alpha                                                   
c                                                                               
        if ((iter .gt. maxiter) .or. isnan(nr)) then                            
            cut = .true.                                                        
            write(*,*) "A bicrystal solution failed"                            
            return                                                              
        end if                                                                  
      end do                                                                    
                                                                                
      end subroutine                                                            
                                                                                
