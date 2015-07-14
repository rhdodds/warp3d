c     ****************************************************************
c     *                                                              *
c     *                      subroutine incrystal                    *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 7/14/2015 rhd              *
C     *                                                              *
C     *     this subroutine supervises and conducts the input of the *
c     *     properties of crystals into the crystal library          *
c     *                                                              *
c     ****************************************************************
c
      subroutine incrystal( sbflg1, sblfg2, chkparam, out )
            use crystal_data
            implicit integer (a-z)
            logical :: sbflg1, sbflg2
c
c           Local variables
            integer :: cnum
            integer :: dumi, nc
            real :: dumr
            double precision :: dumd
            character :: dums*10, lab*15
            logical :: reading
c
c           Check to make sure we have a valid crystal number
            if (.not. numi(cnum)) then
                  call errmsg(350,dumi,dums,dumr,dumd)
                  return
            end if
            if ((cnum .gt. max_crystals) .or. (cnum .lt. 1)) then
                  call errmsg(351,cnum,dums,dumr,dumd)
                  return
            end if
            call readsc()
c           Next line should start with properties
            if ( .not. matchs('properties',10) ) then
                  call errmsg(352,dumi,dums,dumr,dumd)
                  return
            end if
c           Initialize crystal
            call initialize_new_crystal(cnum, out)
c           Read in properties
            reading = .true.
            do while (reading)
                  if ( matchs_exact('slip_type')) then
                        if (.not. label(dumi)) then
                              call errmsg(353,dumi,dums,dumr,dumd)
                        else
                              lab = ' '
                              call entits(lab,nc)
                        end if
                        if (lab(1:nc) .eq. 'fcc') then
                              c_array(cnum)%slip_type = 1
                        elseif (lab(1:nc) .eq. 'bcc') then
                              c_array(cnum)%slip_type = 2
                        elseif (lab(1:nc) .eq. 'single') then
                              c_array(cnum)%slip_type = 3
                        else
                              call errmsg(353,dumi,dums,dumr,dumd)
                        end if
                  elseif ( matchs_exact('elastic_type')) then
                        if ( .not. label(dumi)) then
                              call errmsg(354,dumi,dums,dumr,dumd)
                        else
                              lab = ' '
                              call entits(lab,nc)
                        end if
                        if (lab(1:nc) .eq. 'isotropic') then
                              c_array(cnum)%elastic_type = 1
                        elseif (lab(1:nc) .eq. 'cubic') then
                              c_array(cnum)%elastic_type = 2
                        else
                              call errmsg(354,dumi,dums,dumr,dumd)
                        end if
                  elseif ( matchs_exact('e')) then
                        if (.not. numd(c_array(cnum)%e)) then
                              call errmsg(5,dumi,'e',dumr,dumd)
                        end if
                  elseif ( matchs_exact('nu')) then
                        if (.not. numd(c_array(cnum)%nu)) then
                              call errmsg(5,dumi,'nu',dumr,dumd)
                        end if
                  elseif ( matchs_exact('mu')) then
                        if (.not. numd(c_array(cnum)%mu)) then
                              call errmsg(5,dumi,'mu',dumr,dumd)
                        end if
                  elseif ( matchs_exact('harden_n')) then
                     if (.not. numd(c_array(cnum)%harden_n)) then
                           call errmsg(5,dumi,'harden_n',dumr,dumd)
                     end if
                  elseif ( matchs_exact('tau_a')) then
                     if (.not. numd(c_array(cnum)%tau_a)) then
                           call errmsg(5,dumi,'tau_a',dumr,dumd)
                     end if
                  elseif ( matchs_exact('tau_hat_y')) then
                     if (.not. numd(c_array(cnum)%tau_hat_y)) then
                           call errmsg(5,dumi,'tau_hat_y',dumr,dumd)
                     end if
                  elseif ( matchs_exact('g_0_y')) then
                     if (.not. numd(c_array(cnum)%g_o_y)) then
                           call errmsg(5,dumi,'g_0_y',dumr,dumd)
                     end if
                  elseif ( matchs_exact('tau_hat_v')) then
                     if (.not. numd(c_array(cnum)%tau_hat_v)) then
                           call errmsg(5,dumi,'tau_hat_v',dumr,dumd)
                     end if
                  elseif ( matchs_exact('g_0_v')) then
                     if (.not. numd(c_array(cnum)%g_o_v)) then
                           call errmsg(5,dumi,'g_0_v',dumr,dumd)
                     end if
                  elseif ( matchs_exact('b')) then
                     if (.not. numd(c_array(cnum)%b)) then
                           call errmsg(5,dumi,'b',dumr,dumd)
                     end if
                  elseif ( matchs_exact('p_v')) then
                     if (.not. numd(c_array(cnum)%p_v)) then
                           call errmsg(5,dumi,'p_v',dumr,dumd)
                     end if
                  elseif ( matchs_exact('q_v')) then
                     if (.not. numd(c_array(cnum)%q_v)) then
                           call errmsg(5,dumi,'q_v',dumr,dumd)
                     end if
                  elseif ( matchs_exact('p_y')) then
                     if (.not. numd(c_array(cnum)%p_y)) then
                           call errmsg(5,dumi,'p_y',dumr,dumd)
                     end if
                  elseif ( matchs_exact('q_y')) then
                     if (.not. numd(c_array(cnum)%q_y)) then
                           call errmsg(5,dumi,'q_y',dumr,dumd)
                     end if
                  elseif ( matchs_exact('boltz')) then
                     if (.not. numd(c_array(cnum)%boltz)) then
                           call errmsg(5,dumi,'boltz',dumr,dumd)
                     end if
                  elseif ( matchs_exact('eps_dot_0_v')) then
                     if (.not. numd(c_array(cnum)%eps_dot_o_v)) then
                           call errmsg(5,dumi,'eps_dot_0_v',dumr,dumd)
                     end if
                  elseif ( matchs_exact('eps_dot_0_y')) then
                     if (.not. numd(c_array(cnum)%eps_dot_o_y)) then
                           call errmsg(5,dumi,'eps_dot_0_y',dumr,dumd)
                     end if
                  elseif ( matchs_exact('T_0')) then
                     if (.not. numd(c_array(cnum)%t_o)) then
                           call errmsg(5,dumi,'T_0',dumr,dumd)
                     end if
                  elseif ( matchs_exact('mu_0')) then
                     if (.not. numd(c_array(cnum)%mu_o)) then
                           call errmsg(5,dumi,'mu_0',dumr,dumd)
                     end if
                  elseif ( matchs_exact('D_0')) then
                      if (.not. numd(c_array(cnum)%D_o)) then
                           call errmsg(5,dumi,'D_0',dumr,dumd)
                      end if
                  elseif ( matchs_exact('theta_0')) then
                     if (.not. numd(c_array(cnum)%theta_o)) then
                           call errmsg(5,dumi,'theta_0',dumr,dumd)
                     end if
                  elseif ( matchs_exact('tau_y')) then
                     if (.not. numd(c_array(cnum)%tau_y)) then
                           call errmsg(5,dumi,'tau_y',dumr,dumd)
                     end if
                  elseif ( matchs_exact('tau_v')) then
                     if (.not. numd(c_array(cnum)%tau_v)) then
                           call errmsg(5,dumi,'tau_v',dumr,dumd)
                     end if
                  elseif ( matchs_exact('k_0')) then
                     if (.not. numd(c_array(cnum)%k_o)) then
                           call errmsg(5,dumi,'k_0',dumr,dumd)
                     end if
                  elseif ( matchs_exact('voche_m')) then
                     if (.not. numd(c_array(cnum)%voche_m)) then
                           call errmsg(5,dumi,'voche_m',dumr,dumd)
                     end if
                  elseif ( matchs_exact('voce_m')) then
                     if (.not. numd(c_array(cnum)%voche_m)) then
                           call errmsg(5,dumi,'voche_m',dumr,dumd)
                     end if
                  elseif ( matchs_exact('u_1')) then
                     if (.not. numd(c_array(cnum)%u1)) then
                           call errmsg(5,dumi,'u_1',dumr,dumd)
                     end if
                  elseif ( matchs_exact('u_2')) then
                     if (.not. numd(c_array(cnum)%u2)) then
                           call errmsg(5,dumi,'u_2',dumr,dumd)
                     end if
                  elseif ( matchs_exact('u_3')) then
                     if (.not. numd(c_array(cnum)%u3)) then
                           call errmsg(5,dumi,'u_3',dumr,dumd)
                     end if
                  elseif ( matchs_exact('u_4')) then
                     if (.not. numd(c_array(cnum)%u4)) then
                           call errmsg(5,dumi,'u_4',dumr,dumd)
                     end if
                  elseif ( matchs_exact('u_5')) then
                     if (.not. numd(c_array(cnum)%u5)) then
                           call errmsg(5,dumi,'u_5',dumr,dumd)
                     end if
                  elseif ( matchs_exact('u_6')) then
                     if (.not. numd(c_array(cnum)%u6)) then
                           call errmsg(5,dumi,'u_6',dumr,dumd)
                     end if
                  elseif ( matchs_exact('hardening')) then
                     if (.not. label(dumi)) then
                           call errmsg(5,dumi,'hardening',dumr,dumd)
                     else
                           lab = ' '
                           call entits(lab,nc)
                           if (lab(1:nc) .eq. 'voche') then
                              c_array(cnum)%h_type = 1
                           elseif (lab(1:nc) .eq. 'voce') then
                              c_array(cnum)%h_type = 1
                           elseif (lab(1:nc) .eq. 'mts') then
                              c_array(cnum)%h_type = 2
                           elseif (lab(1:nc) .eq. 'user') then
                              c_array(cnum)%h_type = 3
                           else
                              call errmsg(364,dumi,dums,dumr,dumd)
                           end if
                     end if
                  elseif ( endcrd(dum) ) then
                        reading = .false.
                        cycle
                  elseif ( matchs(',',1) ) then 
                        call readsc()
                  else
                        call entits(lab,nc)
                        call errmsg(355,dumi,lab(1:nc),dumr,dumd)
                        
                        call scan()
                        cycle
                  end if
            end do
c
            call finalize_new_crystal(cnum, out)
c
c           Let further input systems know we have a valid crystal
            defined_crystal = .true.
c            call print_crystal(cnum)
            chkparam = cnum
            sbflg1 = .true.
            sbflg2 = .true.
            return
      end subroutine
