c                                                                               
c ***************************************************************************** 
c *                                                                           * 
c *   material model #11 -- homogenized grain interface damage model          * 
c *                                                                           * 
c ***************************************************************************** 
c                                                                               
      subroutine mm11(gp, span, hist_sz,                                        
     &            history_n, history_np1,                                       
     &            local_work, uddt, gp_temps,                                   
     &            gp_temp_inc, iout)                                            
            use main_data, only: asymmetric_assembly                            
            use mm10_defs                                                       
            implicit integer (a-z)                                              
      include 'include_sig_up'                                                  
            double precision, intent(in) :: uddt(mxvl,nstr)                     
            double precision, intent(in) :: gp_temps(mxvl),                     
     &                                      gp_temp_inc(mxvl)                   
            integer, intent(in) :: iout                                         
            integer :: span, hist_sz, gp                                        
            double precision :: history_n(span,hist_sz)                         
            double precision :: history_np1(span,hist_sz)                       
c                                                                               
            double precision, dimension(27) :: gradFe                           
            double precision :: dsum, dmax, cutoff                              
            integer :: i,c,cpoff,co_1,cn_1,dci_1, s, p, c1, c2,                 
     &                  co_2, cn_2, e, r_1, r_2                                 
            logical :: debug, cp_cut                                            
            type(crystal_props) :: cc_props_1, cc_props_2                       
            type(crystal_state) :: cc_n_1, cc_np1_1, cc_n_2, cc_np1_2           
c                                                                               
            parameter(cutoff=0.95)                                              
c                                                                               
c            if (.not. asymmetric_assembly) then                                
c              write(*,*) "mm11 requires an asymmetric tangent."                
c              call die_gracefully                                              
c            end if                                                             
c                                                                               
            gradFe = 0.0                                                        
c                                                                               
c                                                                               
c                 Unfortunately we will be repeating a lot of mm10's setup,     
c                 but I really don't want to mess with the local_work           
c                 datastructure more than I have to.                            
c                                                                               
            do i=1,span                                                         
              e = i + local_work%felem                                          
c                                                                               
c                 If the element failed there's no point in doing this          
              if (history_n(i,1+local_work%macro_sz) .gt. cutoff) then          
                 history_np1(i,1+local_work%macro_sz) = 1.0                     
                 goto 1337                                                      
              end if                                                            
                                                                                
c                 The offset to the start of the CP history is:                 
c                 macro_sz + ncrystals + 1                                      
              cp_off=local_work%macro_sz+2*local_work%ncrystals(i)+1            
c                                                                               
c                 Fix the problem with the rotation matrix at small strains     
              if (.not. local_work%geo_non_flg) then                            
                  local_work%rot_blk_n1(i, 1:9, gp) =                           
     &              (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/)             
              end if                                                            
c                                                                               
c                 Loop on crystals                                              
              dsum = 0.0                                                        
              do s=1,local_work%nstacks(i)                                      
              dmax = 0.0                                                        
c                 If any of the crystals in this stack has already hit          
c                 a damage above the cutoff, there's no point in continuing     
c                 with computing the stack                                      
              do p=1,local_work%nper(i)                                         
                c1 = (s-1)*local_work%nper(i)+1                                 
                c2 = c1 + 1                                                     
c                 Where this crystal's damage goes                              
                dci_1 = local_work%macro_sz + 1 + c1                            
                dci_2 = local_work%macro_sz + 1 + c2                            
c                 Where this crystal's RT parameter goes                        
                r_1=local_work%macro_sz+local_work%ncrystals(i)+c1+1            
                r_2=local_work%macro_sz+local_work%ncrystals(i)+c2+1            
                                                                                
                if ((history_n(i,dci_1) .gt. cutoff) .or.                       
     &            (history_n(i,dci_2) .gt. cutoff)) then                        
                  history_np1(i,dci_1) = 1.0                                    
                  history_np1(i,dci_2) = 1.0                                    
                  dmax = 1.0                                                    
                  exit                                                          
                end if                                                          
              end do                                                            
                                                                                
              if (dmax .gt. cutoff) go to 1300                                  
                                                                                
              do p=1,local_work%nper(i)                                         
                c1 = (s-1)*local_work%nper(i)+1                                 
                c2 = c1 + 1                                                     
c                 Where this crystal's damage goes                              
                dci_1 = local_work%macro_sz + 1 + c1                            
                co_1 = 76+max_slip_sys+(c1-1)*(25+max_uhard)+cp_off             
                cn_1 = 76+max_slip_sys+(c1)*(25+max_uhard)-1+cp_off             
                dci_2 = local_work%macro_sz + 1 + c2                            
                co_2 = 76+max_slip_sys+(c2-1)*(25+max_uhard)+cp_off             
                cn_2 = 76+max_slip_sys+(c2)*(25+max_uhard)-1+cp_off             
c                 Where this crystal's RT parameter goes                        
                r_1=local_work%macro_sz+local_work%ncrystals(i)+c1+1            
                r_2=local_work%macro_sz+local_work%ncrystals(i)+c2+1            
c                                                                               
                call mm10_init_cc_props(local_work%c_props(i,c1),               
     &                        local_work%angle_type(i),                         
     &                        local_work%angle_convention(i),cc_props_1)        
                cc_props_1%out = iout                                           
                call mm10_init_cc_props(local_work%c_props(i,c2),               
     &                        local_work%angle_type(i),                         
     &                        local_work%angle_convention(i),cc_props_2)        
                cc_props_2%out = iout                                           
                if (local_work%step .eq. 1) then                                
                  call mm10_init_cc_hist0(cc_props_1,                           
     &                 local_work%c_props(i,c1)%init_angles(1:3),               
     &                 history_n(i,co_1:cn_1))                                  
                  history_n(i,dci_1) = 0.0                                      
                  history_n(i,r_1) =  0.0                                       
                  call mm10_init_cc_hist0(cc_props_2,                           
     &                 local_work%c_props(i,c1)%init_angles(1:3),               
     &                 history_n(i,co_2:cn_2))                                  
                  history_n(i,dci_2) = 0.0                                      
                  history_n(i,r_2) = 0.0                                        
                end if                                                          
                                                                                
                call mm10_copy_cc_hist(history_n(i,co_1:cn_1),                  
     &                   grad_Fe,                                               
     &                   history_n(i,(64+cp_off):(72+cp_off)),                  
     &                   cc_props_1, cc_n_1)                                    
                call mm10_copy_cc_hist(history_n(i,co_2:cn_2),                  
     &                   grad_Fe,                                               
     &                   history_n(i,(64+cp_off):(72+cp_off)),                  
     &                   cc_props_2, cc_n_2)                                    
c                                                                               
c                 Solve the bicrystal stresses for this pair                    
c                                                                               
                call mm10_setup_np1(                                            
     &                  local_work%rot_blk_n1(i,1:9,gp), uddt(i,1:6),           
     &                  local_work%dt, gp_temps(i), local_work%step,            
     &                  i+local_work%felem, gp, cc_np1_1)                       
c                                                                               
                call mm10_setup_np1(                                            
     &                  local_work%rot_blk_n1(i,1:9,gp), uddt(i,1:6),           
     &                  local_work%dt, gp_temps(i), local_work%step,            
     &                  i+local_work%felem, gp, cc_np1_2)                       
c                                                                               
                cp_cut = .false.                                                
c                 Solve the bicrystal pair                                      
                call mm11_solve_bicrystal(cc_props_1,                           
     &            cc_np1_1, cc_n_1,                                             
     &            cc_props_2, cc_np1_2, cc_n_2,                                 
     &            local_work%sv, local_work%lv, local_work%tv,                  
     &            cp_cut, iout)                                                 
                                                                                
c                if ((local_work%step .eq. 100) .or.                            
c     &                  (local_work%step .eq. 1)) then                         
c                  write(*,*) "Crystal 1"                                       
c                  write(*,*) cc_np1_1%stress                                   
c                  write(*,*) cc_np1_1%euler_angles                             
c                  write(*,*) "Crystal 2"                                       
c                  write(*,*) cc_np1_2%stress                                   
c                  write(*,*) cc_np1_2%euler_angles                             
c                end if                                                         
c                                                                               
c                                                                               
                if (cp_cut) then                                                
                  write(*,*) "Bicrystal model failed to converge."              
                  local_work%material_cut_step = .true.                         
                  return                                                        
                end if                                                          
                                                                                
c               Store the CP history for this crystal                           
                call mm10_store_cryhist(cc_props_1, cc_np1_1, cc_n_1,           
     &                  history_np1(i,co_1:cn_1))                               
                call mm10_store_cryhist(cc_props_2, cc_np1_2, cc_n_2,           
     &                  history_np1(i,co_2:cn_2))                               
c                                                                               
c               Calculate and store damage for these crystals                   
                call mm11_calc_cry_damage(cc_props_1, cc_np1_1, cc_n_1,         
     &            local_work%alpha_dmg,                                         
     &            history_n(i,dci_1), history_np1(i,dci_1),                     
     &            history_n(i,r_1), history_np1(i,r_1))                         
                call mm11_calc_cry_damage(cc_props_2, cc_np1_2, cc_n_2,         
     &            local_work%alpha_dmg,                                         
     &            history_n(i,dci_2), history_np1(i,dci_2),                     
     &            history_n(i,r_2), history_np1(i,r_2))                         
c                                                                               
c               Find the max damage in this stack                               
                if (history_np1(i,dci_1) .gt. dmax) then                        
                  dmax = history_np1(i,dci_1)                                   
                end if                                                          
                if (history_np1(i,dci_2) .gt. dmax) then                        
                  dmax = history_np1(i,dci_2)                                   
                end if                                                          
c               If dmax is already 1, we're done                                
                if (dmax .gt. 0.99) then                                        
                  dmax = 1.0                                                    
                  exit                                                          
                end if                                                          
              end do                                                            
c                                                                               
 1300         continue                                                          
              dsum = dsum + dmax                                                
                                                                                
              end do                                                            
c                                                                               
c                       Store R for the next CP run                             
              history_np1(i,(64+cp_off):(72+cp_off)) =                          
     &            local_work%rot_blk_n1(i,1:9,gp)                               
c                                                                               
c                 Store the total damage at state np1                           
              history_np1(i,1+local_work%macro_sz) =                            
     &            dsum / dble(local_work%nstacks(i))                            
c                                                                               
 1337         continue                                                          
c                 Apply the damage to the stress                                
              call mm11_damage_stress(                                          
     &            local_work%urcs_blk_n1(i,1:6,gp),                             
     &            history_n(i,1+local_work%macro_sz),                           
     &            local_work%sv, local_work%lv, local_work%tv)                  
c                                                                               
c                 Store the damage in the user output variables (9 is c1...)    
              local_work%urcs_blk_n1(i,9,gp) =                                  
     &            history_np1(i,1+local_work%macro_sz)                          
            end do                                                              
                                                                                
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_elem_size                    *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/12/14                     *          
c     *                                                              *          
c     *    Calculate the number of crystal pairs per element         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_elem_size(etype, nenodes, node_coords, sv,                
     &            lv, tv, ls, ll, lt, nstacks, nperstack)                       
      implicit none                                                             
c                                                                               
      integer :: etype, nenodes, nstacks, nperstack                             
      double precision, dimension(3) :: sv, lv, tv                              
      double precision :: node_coords(*)                                        
      double precision :: ls, ll, lt                                            
c                                                                               
      integer :: i, j, ns, nl, nt                                               
      double precision, dimension(3) :: size_vec                                
      double precision :: les, lel, let                                         
      double precision :: x1,x2,y1,y2,z1,z2                                     
      logical :: brick, tet, wedge, linear                                      
c                                                                               
      brick = .false.                                                           
      tet = .false.                                                             
      wedge = .false.                                                           
      linear = .false.                                                          
c                                                                               
      brick = etype .ge. 1  .and. etype .le. 5                                  
      tet   = etype .eq. 6 .or. etype .eq. 13                                   
      wedge = etype .eq. 7                                                      
      linear = etype .eq. 2  .or.  etype .eq. 13                                
c                                                                               
      if (.not. brick) then                                                     
        write(*,*) "Cannot determine size for non-brick element"                
        call die_abort                                                          
      end if                                                                    
c                                                                               
c           The brick algorithm is easy -- just figure out the vector           
c           with [max delta_x, max delta_y, max_delta_z], dot it                
c           with the directions to get the element length, and then             
c           divide to get the number                                            
c                                                                               
      if (brick) then                                                           
        size_vec = 0.0                                                          
        do i=1,nenodes                                                          
          x1 = node_coords(i)                                                   
          y1 = node_coords(i+nenodes)                                           
          z1 = node_coords(i+nenodes+nenodes)                                   
          do j=1,nenodes                                                        
            x2 = node_coords(j)                                                 
            y2 = node_coords(j+nenodes)                                         
            z2 = node_coords(j+nenodes+nenodes)                                 
            if (abs(x2-x1) .gt. size_vec(1)) then                               
              size_vec(1) = abs(x2-x1)                                          
            end if                                                              
            if (abs(y2-y1) .gt. size_vec(2)) then                               
              size_vec(2) = abs(y2-y1)                                          
            end if                                                              
            if (abs(z2-z1) .gt. size_vec(3)) then                               
              size_vec(3) = abs(z2-z1)                                          
            end if                                                              
          end do                                                                
        end do                                                                  
        les = dot_product(size_vec, abs(sv))                                    
        lel = dot_product(size_vec, abs(lv))                                    
        let = dot_product(size_vec, abs(tv))                                    
      end if                                                                    
                                                                                
      ns = nint(les/ls)                                                         
      nl = nint(lel/ll)                                                         
      nt = nint(let/lt)                                                         
                                                                                
      if (ns .lt. 1) ns = 1                                                     
      if (nl .lt. 1) nl = 1                                                     
      if (nt .lt. 1) nt = 1                                                     
                                                                                
      if ((ns .lt. 1) .or. (nl .lt. 1) .or. (nt .lt. 1)) then                   
        write(*,*) "Invalid number of crystals"                                 
        call die_abort                                                          
      end if                                                                    
                                                                                
      nstacks = nl*nt                                                           
      nperstack = ns/2                                                          
                                                                                
      if (nperstack .lt. 1) nperstack = 1                                       
                                                                                
      end subroutine                                                            
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_calc_cry_damage              *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/11/14                     *          
c     *                                                              *          
c     *    Calculate the damage for each crystal (state np1)         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_calc_cry_damage(props, np1, n, alpha, dn, dnp1,           
     &      rn, rnp1)                                                           
      use mm10_defs                                                             
      implicit none                                                             
c                                                                               
      type(crystal_props) :: props                                              
      type(crystal_state) :: np1, n                                             
      double precision :: alpha, rn, rnp1, dn, dnp1                             
c                                                                               
      double precision :: rinc, sigma_m, sigma_e                                
      double precision, dimension(6) :: sigma_dev                               
c                                                                               
      sigma_m = np1%stress(1)+np1%stress(2)+np1%stress(3)                       
      sigma_dev = np1%stress                                                    
      sigma_dev(1) = sigma_dev(1) - 1.0/3.0*sigma_m                             
      sigma_dev(2) = sigma_dev(2) - 1.0/3.0*sigma_m                             
      sigma_dev(3) = sigma_dev(3) - 1.0/3.0*sigma_m                             
      sigma_e = dsqrt((dot_product(sigma_dev(1:3), sigma_dev(1:3))+             
     &      2.0*dot_product(sigma_dev(4:6), sigma_dev(4:6)))*3.0/2.0)           
      rinc = dexp(3.0/2.0*sigma_m/sigma_e)                                      
c                                                                               
      rnp1 = rn + rinc*np1%p_strain_inc                                         
c                                                                               
      dnp1 = dn + 1.0/(alpha*dcosh(rnp1/alpha)**2.0)*rinc/alpha*                
     &            np1%p_strain_inc                                              
c                                                                               
c           Could be slightly greater than 1 -- don't let that happen           
      if (dnp1 .gt. 1.0) dnp1 = 1.0                                             
      return                                                                    
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_damage_stress                *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/10/14                     *          
c     *                                                              *          
c     *    apply damage to a stress vector                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_damage_stress(stressv, D, s, l, t)                        
            implicit none                                                       
            double precision, dimension(6), intent(inout) :: stressv            
            double precision, intent(in) :: D                                   
            double precision, dimension(3), intent(in) :: s, l, t               
c                                                                               
            double precision, dimension(6,6) :: pD, I                           
c                                                                               
            I = 0.0                                                             
            I(1,1) = 1.0                                                        
            I(2,2) = 1.0                                                        
            I(3,3) = 1.0                                                        
            I(4,4) = 1.0                                                        
            I(5,5) = 1.0                                                        
            I(6,6) = 1.0                                                        
c                                                                               
            call mm11_pD(s, l, t, pD)                                           
c                                                                               
            stressv = matmul(I-D*pD,stressv)                                    
                                                                                
      end subroutine                                                            
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine cnst11                            *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/11/14                     *          
c     *                                                              *          
c     *    apply damage to a block of tangents                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine cnst11( gp, span, cep_blk, dmg_blk, s_v, l_v, t_v,             
     &                  iout)                                                   
      implicit none                                                             
                                                                                
      integer :: gp, span, iout                                                 
      double precision, dimension(span,6,6) :: cep_blk                          
      double precision, dimension(3) :: s_v, l_v, t_v                           
      double precision, dimension(span) :: dmg_blk                              
                                                                                
      integer :: i                                                              
                                                                                
                                                                                
      do i=1,span                                                               
        call mm11_damage_tangent(cep_blk(i,1:6,1:6), dmg_blk(i),                
     &      s_v, l_v, t_v)                                                      
      end do                                                                    
                                                                                
      end subroutine                                                            
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_damage_tangent               *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/10/14                     *          
c     *                                                              *          
c     *    apply damage to a tangent matrix                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_damage_tangent(cep, D, s, l, t)                           
            use main_data, only: asymmetric_assembly                            
            implicit none                                                       
            double precision, dimension(6,6), intent(inout) :: cep              
            double precision, intent(in) :: D                                   
            double precision, dimension(3), intent(in) :: s, l, t               
c                                                                               
            double precision, dimension(6,6) :: pD, I, fullt                    
c                                                                               
            I = 0.0                                                             
            I(1,1) = 1.0                                                        
            I(2,2) = 1.0                                                        
            I(3,3) = 1.0                                                        
            I(4,4) = 1.0                                                        
            I(5,5) = 1.0                                                        
            I(6,6) = 1.0                                                        
c                                                                               
            call mm11_pD(s, l, t, pD)                                           
c                                                                               
            cep = matmul(I-D*pD,cep)                                            
            if (.not. asymmetric_assembly) then                                 
              cep = 0.5*(cep+transpose(cep))                                    
            end if                                                              
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_pe                           *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/10/14                     *          
c     *                                                              *          
c     *    Return the equilibrium projection                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_pe(n, pe)                                                 
      implicit none                                                             
      double precision, dimension(3,6), intent(out) :: pe                       
      double precision, dimension(3) :: n                                       
      double precision n1, n2, n3, s1, s2, s3, t1, t2, t3                       
c                                                                               
      n1 = n(1)                                                                 
      n2 = n(2)                                                                 
      n3 = n(3)                                                                 
c                                                                               
      pe = 0.0                                                                  
      pe(1,1) = n1                                                              
      pe(1,4) = n2                                                              
      pe(1,6) = n3                                                              
      pe(2,2) = n2                                                              
      pe(2,4) = n1                                                              
      pe(2,5) = n3                                                              
      pe(3,3) = n3                                                              
      pe(3,5) = n2                                                              
      pe(3,6) = n1                                                              
c                                                                               
      return                                                                    
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_pc                           *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/11/14                     *          
c     *                                                              *          
c     *    Return the compatibility projection                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_pc(s, t, pc)                                              
      implicit none                                                             
      double precision, dimension(3,6), intent(out) :: pc                       
      double precision, dimension(3) :: s, t                                    
      double precision s1, s2, s3, t1, t2, t3                                   
c                                                                               
      s1 = s(1)                                                                 
      s2 = s(2)                                                                 
      s3 = s(3)                                                                 
      t1 = t(1)                                                                 
      t2 = t(2)                                                                 
      t3 = t(3)                                                                 
c                                                                               
      pc = 0.0                                                                  
c                                                                               
      pc(1,1) = s1**2.0                                                         
      pc(1,2) = s2**2.0                                                         
      pc(1,3) = s3**2.0                                                         
      pc(1,4) = s1*s2                                                           
      pc(1,5) = s2*s3                                                           
      pc(1,6) = s1*s3                                                           
      pc(2,1) = t1**2.0                                                         
      pc(2,2) = t2**2.0                                                         
      pc(2,3) = t3**2.0                                                         
      pc(2,4) = t1*t2                                                           
      pc(2,5) = t2*t3                                                           
      pc(2,6) = t1*t3                                                           
      pc(3,1) = s1*t1                                                           
      pc(3,2) = s2*t2                                                           
      pc(3,3) = s3*t3                                                           
      pc(3,4) = 0.5*(s2*t1+s1*t2)                                               
      pc(3,5) = 0.5*(s2*t3+s3*t2)                                               
      pc(3,6) = 0.5*(s3*t1+s1*t3)                                               
c                                                                               
      return                                                                    
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_pD                           *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/10/14                     *          
c     *                                                              *          
c     *    Return the stupid huge projection for stress/tangent      *          
c     *    type matrices                                             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_pD(n, s, t, pD)                                           
      implicit none                                                             
      double precision, dimension(6,6), intent(out) :: pD                       
      double precision, dimension(3) :: n, s, t                                 
      double precision n1, n2, n3, s1, s2, s3, t1, t2, t3                       
c                                                                               
      n1 = n(1)                                                                 
      n2 = n(2)                                                                 
      n3 = n(3)                                                                 
      s1 = s(1)                                                                 
      s2 = s(2)                                                                 
      s3 = s(3)                                                                 
      t1 = t(1)                                                                 
      t2 = t(2)                                                                 
      t3 = t(3)                                                                 
c                                                                               
c                                                                               
      pD = 0.0                                                                  
c                                                                               
      pD(1,1) = n1**2.0*(n1**2.0+2*(s1**2.0+t1**2.0))                           
      pD(1,2) = n1*n2*(n1*n2+2*s1*s2+2*t1*t2)                                   
      pD(1,3) = n1*n3*(n1*n3+2*s1*s3+2*t1*t3)                                   
      pD(1,4) = 2*n1*(n1*(s1*s2+t1*t2)+n2*(s1**2.0+t1**2.0)+n2*                 
     &      n1**2.0)                                                            
      pD(1,5) = 2*n1*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3)+n1*n2*n3)               
      pD(1,6) = 2*n1*(n1*(s1*s3+t1*t3)+n3*(s1**2.0+t1**2.0)+n3*                 
     &      n1**2.0)                                                            
      pD(2,1) = n1*n2*(n1*n2+2*s1*s2+2*t1*t2)                                   
      pD(2,2) = n2*(2*n2*(s2**2.0+t2**2.0)+n2**3.0)                             
      pD(2,3) = n2*n3*(n2*n3+2*s2*s3+2*t2*t3)                                   
      pD(2,4) = 2*n2*(n2*(s1*s2+t1*t2)+n1*(n2**2.0+s2**2.0+                     
     &      t2**2.0))                                                           
      pD(2,5) = 2*n2*(n2*(s2*s3+t2*t3)+n3*(s2**2.0+t2**2.0)+                    
     &      n3*n2**2.0)                                                         
      pD(2,6) = 2*n2*(n3*(s1*s2+t1*t2)+n1*(n2*n3+s2*s3+t2*t3))                  
      pD(3,1) = n1*n3*(n1*n3+2*s1*s3+2*t1*t3)                                   
      pD(3,2) = n2*n3*(n2*n3+2*s2*s3+2*t2*t3)                                   
      pD(3,3) = n3**2.0*(n3**2.0+2*(s3**2.0+t3**2.0))                           
      pD(3,4) = 2*n3*(n2*(s1*s3+t1*t3)+n1*(n2*n3+s2*s3+t2*t3))                  
      pD(3,5) = n3*(2*n3*(s2*s3+t2*t3)+2*n2*                                    
     &      (n3**2.0+s3**2.0+t3**2.0))                                          
      pD(3,6) = n3*(2*n3*(s1*s3+t1*t3)+2*n1*                                    
     &      (n3**2.0+s3**2.0+t3**2.0))                                          
      pD(4,1) = n1*(n1*(s1*s2+t1*t2)+n2*(s1**2.0+t1**2.0)+                      
     &      n2*n1**2.0)                                                         
      pD(4,2) = n2*(n2*(s1*s2+t1*t2)+n1*(n2**2.0+s2**2.0+t2**2.0))              
      pD(4,3) = n3*(n2*(s1*s3+t1*t3)+n1*(n2*n3+s2*s3+t2*t3))                    
      pD(4,4) = n1**2.0*(2*n2**2.0+s2**2.0+t2**2.0)+2*n2*n1*                    
     &      (s1*s2+t1*t2)+n2**2.0*(s1**2.0+t1**2.0)                             
      pD(4,5) = n2*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3))+n1*                      
     &      (n2*(s2*s3+t2*t3)+n3*(s2**2.0+t2**2.0)+2*n3*n2**2.0)                
      pD(4,6) = n1**2.0*(2*n2*n3+s2*s3+t2*t3)+n1*(n3*                           
     &      (s1*s2+t1*t2)+n2*(s1*s3+t1*t3))+n2*n3*(s1**2.0+t1**2.0)             
      pD(5,1) = n1*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3)+n1*n2*n3)                 
      pD(5,2) = n2*(n2*(s2*s3+t2*t3)+                                           
     &      n3*(s2**2.0+t2**2.0)+n3*n2**2.0)                                    
      pD(5,3) = n3*(n3*(s2*s3+t2*t3)+n2*(n3**2.0+s3**2.0+t3**2.0))              
      pD(5,4) = n2*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3))+n1*                      
     &      (n2*(s2*s3+t2*t3)+n3*(s2**2.0+t2**2.0)+2*n3*n2**2.0)                
      pD(5,5) = n2**2.0*(2*n3**2.0+s3**2.0+t3**2.0)+2*n3*                       
     &      n2*(s2*s3+t2*t3)+n3**2.0*(s2**2.0+t2**2.0)                          
      pD(5,6) = n3*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3))+                         
     &      n1*(n3*(s2*s3+t2*t3)+n2*(2*n3**2.0+s3**2.0+t3**2.0))                
      pD(6,1) = n1*(n1*(s1*s3+t1*t3)+                                           
     &      n3*(s1**2.0+t1**2.0)+n3*n1**2.0)                                    
      pD(6,2) = n2*(n3*(s1*s2+t1*t2)+n1*(n2*n3+s2*s3+t2*t3))                    
      pD(6,3) = n3*(n3*(s1*s3+t1*t3)+n1*(n3**2.0+s3**2.0+t3**2.0))              
      pD(6,4) = n1**2.0*(2*n2*n3+s2*s3+t2*t3)+                                  
     &      n1*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3))+                             
     &      n2*n3*(s1**2.0+t1**2.0)                                             
      pD(6,5) = n3*(n3*(s1*s2+t1*t2)+n2*(s1*s3+t1*t3))+                         
     &      n1*(n3*(s2*s3+t2*t3)+n2*(2*n3**2.0+s3**2.0+t3**2.0))                
      pD(6,6) = n1**2.0*(2*n3**2.0+s3**2.0+t3**2.0)+                            
     &      2*n3*n1*(s1*s3+t1*t3)+n3**2.0*(s1**2.0+t1**2.0)                     
c                                                                               
      return                                                                    
      end subroutine                                                            
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_set_sizes_special            *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified: 3/10/14                     *          
c     *                                                              *          
c     *    called by warp3d for each material model to obtain        *          
c     *    various sizes of data for the model                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_set_sizes_special( matnum, size_data, local_el )          
      use global_data ! old common.main
      use main_data, only: imatprp                                              
      implicit integer (a-z)                                                    
      dimension size_data(*)                                                    
      integer :: local_el, matnum, ncrystals, nstacks, omatn                    
c                                                                               
      integer, dimension(4) :: fake_sz                                          
c                                                                               
      integer :: old_sz                                                         
c                                                                               
c           We need to store this where the matprops can get it it...           
      old_sz = size_data(1)                                                     
      imatprp(132,matnum) = old_sz                                              
c                                                                               
      fake_sz(1) = 0                                                            
c                                                                               
      ncrystals = imatprp(101,matnum)                                           
c                                                                               
c           New size is going to be:                                            
c                 old_sz                                                        
c               + 1           (damage)                                          
c               + (cp_hist_sz+2*ncrystals) (crystal damage + crystal history)   
c                                                                               
c           Fake out the material model                                         
      omn = iprops(38,local_el)                                                 
      iprops(38,local_el) = matnum                                              
c           Get the CP size                                                     
c      call mm10_set_sizes(fake_sz, local_el)                                   
       write(*,9000)                                                            
       call die_abort                                                           
      imatprp(133,matnum) = fake_sz(1)                                          
c           Sum up                                                              
      size_data(1) = old_sz + 1 + fake_sz(1) + 2*ncrystals                      
c           Reset material data                                                 
      iprops(38,local_el) = omn                                                 
      return                                                                    
 9000 format('>> FATAL ERROR: routine mm11_set_sizes_special',                  
     & /,    '                not implemented',                                 
     & /,    '                job terminated',//)                               
c                                                                               
      end subroutine                                                            
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_set_sizes                    *          
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
      subroutine mm11_set_sizes( info_vector )                                  
      integer, dimension(*) :: info_vector                                      
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
      info_vector(1) = -100    ! set by special version above                   
      info_vector(2) = -100    ! set by special version above                   
      info_vector(3) = -100    ! set by special version above                   
      info_vector(4) = 0                                                        
c                                                                               
      return                                                                    
      end                                                                       
c            dummy routines for model not yet supporting                        
c            states output                                                      
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine mm11_states_values                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/3/2015 (rhd))                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_states_values( itype, elem_states_output,                 
     &                                 nrow_states, num_states  )               
      use global_data ! old common.main
c                                                                               
c                       access some global data structures                      
c                                                                               
      use elem_block_data, only: history_blocks, history_blk_list               
      use main_data, only: elems_to_blocks, cohesive_ele_types                  
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: nrow_states, itype, num_states                                 
      double precision :: elem_states_output(nrow_states,*)                     
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm11_states_labels                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/11/2015 (rhd)                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm11_states_labels( size_state,                                
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
c                                                                               
      num_states = 0                                                            
      num_comment_lines = 0                                                     
      state_labels(1) = "..."                                                   
      state_descriptors(1) = "...."                                             
c                                                                               
      return                                                                    
      end                                                                       
