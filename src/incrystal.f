c     ****************************************************************
c     *                                                              *
c     *                      subroutine incrystal                    *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 12/20/2016 rhd             *
C     *                                                              *
C     *     input properties of crystals into the crystal library    *
c     *                                                              *
c     ****************************************************************
c
      subroutine incrystal( sbflg1, sbflg2, chkparam, out )
      use crystal_data
      implicit none 
c
c
c              parameters
c      
      logical :: sbflg1, sbflg2
      integer :: out, chkparam
c
c              locals

      integer :: cnum, dumi, nc
      real :: dumr
      double precision :: dumd
      character :: dums*10, lab*15
      logical :: reading
      logical, external :: numi, matchs, label, matchs_exact, numd,
     &                     endcrd 
c
c               make sure we have a valid crystal number
c
      if( .not. numi(cnum) ) then
            call errmsg(350,dumi,dums,dumr,dumd)
            return
      end if
      if( (cnum .gt. max_crystals) .or. (cnum .lt. 1) ) then
            call errmsg(351,cnum,dums,dumr,dumd)
            return
      end if
      call readsc()
c
c               next line should start with properties
c
      if(  .not. matchs('properties',10) ) then
            call errmsg(352,dumi,dums,dumr,dumd)
            return
      end if
c      
c                initialize crystal
c
      call initialize_new_crystal( cnum, out )
      reading = .true.
      do while( reading )
       if(  matchs_exact('slip_type') ) then
             if( .not. label(dumi) ) then
                   call errmsg(353,dumi,dums,dumr,dumd)
             else
                   lab = ' '
                   call entits(lab,nc)
             end if
             if( lab(1:nc) .eq. 'fcc') then
                   c_array(cnum)%slip_type = 1
             elseif( lab(1:nc) .eq. 'bcc') then
                   c_array(cnum)%slip_type = 2
             elseif( lab(1:nc) .eq. 'single') then
                   c_array(cnum)%slip_type = 3
             elseif( lab(1:nc) .eq. 'roters') then
                   c_array(cnum)%slip_type = 6
             elseif( lab(1:nc) .eq. 'bcc12') then
                   c_array(cnum)%slip_type = 7
             elseif( lab(1:nc) .eq. 'bcc48') then
                   c_array(cnum)%slip_type = 8
             elseif( lab(1:nc) .eq. 'hcp6') then
                   c_array(cnum)%slip_type = 9
             elseif( lab(1:nc) .eq. 'hcp18') then
                   c_array(cnum)%slip_type = 10
             else
                   call errmsg(353,dumi,dums,dumr,dumd)
             end if
       elseif(  matchs_exact('elastic_type') ) then
             if(  .not. label(dumi) ) then
                   call errmsg(354,dumi,dums,dumr,dumd)
             else
                   lab = ' '
                   call entits(lab,nc)
             end if
             if( lab(1:nc) .eq. 'isotropic') then
                   c_array(cnum)%elastic_type = 1
             elseif( lab(1:nc) .eq. 'cubic') then
                   c_array(cnum)%elastic_type = 2
             elseif( lab(1:nc) .eq. 'ti6242') then
                   c_array(cnum)%elastic_type = 3
             else
                   call errmsg(354,dumi,dums,dumr,dumd)
             end if
       elseif(  matchs_exact('alter_mode') ) then
             if( .not. label(dumi) ) then
                   call errmsg(364,dumi,dums,dumr,dumd)
             else
                   lab = ' '
                   call entits(lab,nc)
             end if
             if( lab(1:nc) .eq. 'true') then
                   c_array(cnum)%alter_mode = .true.
             elseif( lab(1:nc) .eq. 'false') then
                   c_array(cnum)%alter_mode = .false.
             else
                   call errmsg(364,dumi,dums,dumr,dumd)
             end if
       elseif(  matchs_exact('e') ) then
             if( .not. numd(c_array(cnum)%e) ) then
                   call errmsg(5,dumi,'e',dumr,dumd)
             end if
       elseif(  matchs_exact('nu') ) then
             if( .not. numd(c_array(cnum)%nu) ) then
                   call errmsg(5,dumi,'nu',dumr,dumd)
             end if
       elseif(  matchs_exact('mu') ) then
             if( .not. numd(c_array(cnum)%mu) ) then
                   call errmsg(5,dumi,'mu',dumr,dumd)
             end if
       elseif(  matchs_exact('harden_n') ) then
          if( .not. numd(c_array(cnum)%harden_n) ) then
                call errmsg(5,dumi,'harden_n',dumr,dumd)
          end if
       elseif(  matchs_exact('tau_a') ) then
          if( .not. numd(c_array(cnum)%tau_a) ) then
                call errmsg(5,dumi,'tau_a',dumr,dumd)
          end if
       elseif(  matchs_exact('tau_hat_y') ) then
          if( .not. numd(c_array(cnum)%tau_hat_y) ) then
                call errmsg(5,dumi,'tau_hat_y',dumr,dumd)
          end if
       elseif(  matchs_exact('g_0_y') ) then
          if( .not. numd(c_array(cnum)%g_o_y) ) then
                call errmsg(5,dumi,'g_0_y',dumr,dumd)
          end if
       elseif(  matchs_exact('tau_hat_v') ) then
          if( .not. numd(c_array(cnum)%tau_hat_v) ) then
                call errmsg(5,dumi,'tau_hat_v',dumr,dumd)
          end if
       elseif(  matchs_exact('g_0_v') ) then
          if( .not. numd(c_array(cnum)%g_o_v) ) then
                call errmsg(5,dumi,'g_0_v',dumr,dumd)
          end if
       elseif(  matchs_exact('b') ) then
          if( .not. numd(c_array(cnum)%b) ) then
                call errmsg(5,dumi,'b',dumr,dumd)
          end if
       elseif(  matchs_exact('p_v') ) then
          if( .not. numd(c_array(cnum)%p_v) ) then
                call errmsg(5,dumi,'p_v',dumr,dumd)
          end if
       elseif(  matchs_exact('q_v') ) then
          if( .not. numd(c_array(cnum)%q_v) ) then
                call errmsg(5,dumi,'q_v',dumr,dumd)
          end if
       elseif(  matchs_exact('p_y') ) then
          if( .not. numd(c_array(cnum)%p_y) ) then
                call errmsg(5,dumi,'p_y',dumr,dumd)
          end if
       elseif(  matchs_exact('q_y') ) then
          if( .not. numd(c_array(cnum)%q_y) ) then
                call errmsg(5,dumi,'q_y',dumr,dumd)
          end if
       elseif(  matchs_exact('boltz') ) then
          if( .not. numd(c_array(cnum)%boltz) ) then
                call errmsg(5,dumi,'boltz',dumr,dumd)
          end if
       elseif(  matchs_exact('eps_dot_0_v') ) then
          if( .not. numd(c_array(cnum)%eps_dot_o_v) ) then
                call errmsg(5,dumi,'eps_dot_0_v',dumr,dumd)
          end if
       elseif(  matchs_exact('eps_dot_0_y') ) then
          if( .not. numd(c_array(cnum)%eps_dot_o_y) ) then
                call errmsg(5,dumi,'eps_dot_0_y',dumr,dumd)
          end if
       elseif(  matchs_exact('T_0') ) then
          if( .not. numd(c_array(cnum)%t_o) ) then
                call errmsg(5,dumi,'T_0',dumr,dumd)
          end if
       elseif(  matchs_exact('mu_0') ) then
          if( .not. numd(c_array(cnum)%mu_o) ) then
                call errmsg(5,dumi,'mu_0',dumr,dumd)
          end if
       elseif(  matchs_exact('D_0') ) then
           if( .not. numd(c_array(cnum)%D_o) ) then
                call errmsg(5,dumi,'D_0',dumr,dumd)
           end if
       elseif(  matchs_exact('theta_0') ) then
          if( .not. numd(c_array(cnum)%theta_o) ) then
                call errmsg(5,dumi,'theta_0',dumr,dumd)
          end if
       elseif(  matchs_exact('tau_y') ) then
          if( .not. numd(c_array(cnum)%tau_y) ) then
                call errmsg(5,dumi,'tau_y',dumr,dumd)
          end if
       elseif(  matchs_exact('tau_v') ) then
          if( .not. numd(c_array(cnum)%tau_v) ) then
                call errmsg(5,dumi,'tau_v',dumr,dumd)
          end if
       elseif(  matchs_exact('k_0') ) then
          if( .not. numd(c_array(cnum)%k_o) ) then
                call errmsg(5,dumi,'k_0',dumr,dumd)
          end if
       elseif(  matchs_exact('voche_m') ) then
          if( .not. numd(c_array(cnum)%voche_m) ) then
                call errmsg(5,dumi,'voche_m',dumr,dumd)
          end if
       elseif(  matchs_exact('voce_m') ) then
          if( .not. numd(c_array(cnum)%voche_m) ) then
                call errmsg(5,dumi,'voche_m',dumr,dumd)
          end if
       elseif(  matchs_exact('iD_v') ) then
          if( .not. numd(c_array(cnum)%iD_v) ) then
                call errmsg(5,dumi,'iD_v',dumr,dumd)
          end if
       elseif(  matchs_exact('u_1') ) then
          if( .not. numd(c_array(cnum)%u1) ) then
                call errmsg(5,dumi,'u_1',dumr,dumd)
          end if
       elseif(  matchs_exact('u_2') ) then
          if( .not. numd(c_array(cnum)%u2) ) then
                call errmsg(5,dumi,'u_2',dumr,dumd)
          end if
       elseif(  matchs_exact('u_3') ) then
          if( .not. numd(c_array(cnum)%u3) ) then
                call errmsg(5,dumi,'u_3',dumr,dumd)
          end if
       elseif(  matchs_exact('u_4') ) then
          if( .not. numd(c_array(cnum)%u4) ) then
                call errmsg(5,dumi,'u_4',dumr,dumd)
          end if
       elseif(  matchs_exact('u_5') ) then
          if( .not. numd(c_array(cnum)%u5) ) then
                call errmsg(5,dumi,'u_5',dumr,dumd)
          end if
       elseif(  matchs_exact('u_6') ) then
          if( .not. numd(c_array(cnum)%u6) ) then
                call errmsg(5,dumi,'u_6',dumr,dumd)
          end if
       elseif(  matchs_exact('u_7') ) then
          if( .not. numd(c_array(cnum)%u7) ) then
                call errmsg(5,dumi,'u_7',dumr,dumd)
          end if
       elseif(  matchs_exact('u_8') ) then
          if( .not. numd(c_array(cnum)%u8) ) then
                call errmsg(5,dumi,'u_8',dumr,dumd)
          end if
       elseif(  matchs_exact('u_9') ) then
          if( .not. numd(c_array(cnum)%u9) ) then
                call errmsg(5,dumi,'u_9',dumr,dumd)
          end if
       elseif(  matchs_exact('u_10') ) then
          if( .not. numd(c_array(cnum)%u10) ) then
                call errmsg(5,dumi,'u_10',dumr,dumd)
          end if
       elseif(  matchs_exact('hardening') ) then
          if( .not. label(dumi) ) then
                call errmsg(5,dumi,'hardening',dumr,dumd)
          else
                lab = ' '
                call entits(lab,nc)
                if( lab(1:nc) .eq. 'voche') then
                   c_array(cnum)%h_type = 1
                elseif( lab(1:nc) .eq. 'voce') then
                   c_array(cnum)%h_type = 1
                elseif( lab(1:nc) .eq. 'mts') then
                   c_array(cnum)%h_type = 2
                elseif( lab(1:nc) .eq. 'user') then
                   c_array(cnum)%h_type = 3
                elseif( lab(1:nc) .eq. 'ornl') then
                   c_array(cnum)%h_type = 4
                elseif( lab(1:nc) .eq. 'roters') then
                   c_array(cnum)%h_type = 7
                elseif( lab(1:nc) .eq. 'djgm') then
                   c_array(cnum)%h_type = 9
                else
                   call errmsg(364,dumi,dums,dumr,dumd)
                end if
          end if
       elseif(  matchs_exact('cp_001') ) then
          if( .not. numd(c_array(cnum)%cp_001) ) then
                call errmsg(5,dumi,'cp_001',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_002') ) then
          if( .not. numd(c_array(cnum)%cp_002) ) then
                call errmsg(5,dumi,'cp_002',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_003') ) then
          if( .not. numd(c_array(cnum)%cp_003) ) then
                call errmsg(5,dumi,'cp_003',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_004') ) then
          if( .not. numd(c_array(cnum)%cp_004) ) then
                call errmsg(5,dumi,'cp_004',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_005') ) then
          if( .not. numd(c_array(cnum)%cp_005) ) then
                call errmsg(5,dumi,'cp_005',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_006') ) then
          if( .not. numd(c_array(cnum)%cp_006) ) then
                call errmsg(5,dumi,'cp_006',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_007') ) then
          if( .not. numd(c_array(cnum)%cp_007) ) then
                call errmsg(5,dumi,'cp_007',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_008') ) then
          if( .not. numd(c_array(cnum)%cp_008) ) then
                call errmsg(5,dumi,'cp_008',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_009') ) then
          if( .not. numd(c_array(cnum)%cp_009) ) then
                call errmsg(5,dumi,'cp_009',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_010') ) then
          if( .not. numd(c_array(cnum)%cp_010) ) then
                call errmsg(5,dumi,'cp_010',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_011') ) then
          if( .not. numd(c_array(cnum)%cp_011) ) then
                call errmsg(5,dumi,'cp_011',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_012') ) then
          if( .not. numd(c_array(cnum)%cp_012) ) then
                call errmsg(5,dumi,'cp_012',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_013') ) then
          if( .not. numd(c_array(cnum)%cp_013) ) then
                call errmsg(5,dumi,'cp_013',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_014') ) then
          if( .not. numd(c_array(cnum)%cp_014) ) then
                call errmsg(5,dumi,'cp_014',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_015') ) then
          if( .not. numd(c_array(cnum)%cp_015) ) then
                call errmsg(5,dumi,'cp_015',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_016') ) then
          if( .not. numd(c_array(cnum)%cp_016) ) then
                call errmsg(5,dumi,'cp_016',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_017') ) then
          if( .not. numd(c_array(cnum)%cp_017) ) then
                call errmsg(5,dumi,'cp_017',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_018') ) then
          if( .not. numd(c_array(cnum)%cp_018) ) then
                call errmsg(5,dumi,'cp_018',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_019') ) then
          if( .not. numd(c_array(cnum)%cp_019) ) then
                call errmsg(5,dumi,'cp_019',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_020') ) then
          if( .not. numd(c_array(cnum)%cp_020) ) then
                call errmsg(5,dumi,'cp_020',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_021') ) then
          if( .not. numd(c_array(cnum)%cp_021) ) then
                call errmsg(5,dumi,'cp_021',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_022') ) then
          if( .not. numd(c_array(cnum)%cp_022) ) then
                call errmsg(5,dumi,'cp_022',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_023') ) then
          if( .not. numd(c_array(cnum)%cp_023) ) then
                call errmsg(5,dumi,'cp_023',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_024') ) then
          if( .not. numd(c_array(cnum)%cp_024) ) then
                call errmsg(5,dumi,'cp_024',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_025') ) then
          if( .not. numd(c_array(cnum)%cp_025) ) then
                call errmsg(5,dumi,'cp_025',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_026') ) then
          if( .not. numd(c_array(cnum)%cp_026) ) then
                call errmsg(5,dumi,'cp_026',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_027') ) then
          if( .not. numd(c_array(cnum)%cp_027) ) then
                call errmsg(5,dumi,'cp_027',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_028') ) then
          if( .not. numd(c_array(cnum)%cp_028) ) then
                call errmsg(5,dumi,'cp_028',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_029') ) then
          if( .not. numd(c_array(cnum)%cp_029) ) then
                call errmsg(5,dumi,'cp_029',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_030') ) then
          if( .not. numd(c_array(cnum)%cp_030) ) then
                call errmsg(5,dumi,'cp_030',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_031') ) then
          if( .not. numd(c_array(cnum)%cp_031) ) then
                call errmsg(5,dumi,'cp_031',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_032') ) then
          if( .not. numd(c_array(cnum)%cp_032) ) then
                call errmsg(5,dumi,'cp_032',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_033') ) then
          if( .not. numd(c_array(cnum)%cp_033) ) then
                call errmsg(5,dumi,'cp_033',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_034') ) then
          if( .not. numd(c_array(cnum)%cp_034) ) then
                call errmsg(5,dumi,'cp_034',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_035') ) then
          if( .not. numd(c_array(cnum)%cp_035) ) then
                call errmsg(5,dumi,'cp_035',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_036') ) then
          if( .not. numd(c_array(cnum)%cp_036) ) then
                call errmsg(5,dumi,'cp_036',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_037') ) then
          if( .not. numd(c_array(cnum)%cp_037) ) then
                call errmsg(5,dumi,'cp_037',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_038') ) then
          if( .not. numd(c_array(cnum)%cp_038) ) then
                call errmsg(5,dumi,'cp_038',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_039') ) then
          if( .not. numd(c_array(cnum)%cp_039) ) then
                call errmsg(5,dumi,'cp_039',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_040') ) then
          if( .not. numd(c_array(cnum)%cp_040) ) then
                call errmsg(5,dumi,'cp_040',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_041') ) then
          if( .not. numd(c_array(cnum)%cp_041) ) then
                call errmsg(5,dumi,'cp_041',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_042') ) then
          if( .not. numd(c_array(cnum)%cp_042) ) then
                call errmsg(5,dumi,'cp_042',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_043') ) then
          if( .not. numd(c_array(cnum)%cp_043) ) then
                call errmsg(5,dumi,'cp_043',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_044') ) then
          if( .not. numd(c_array(cnum)%cp_044) ) then
                call errmsg(5,dumi,'cp_044',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_045') ) then
          if( .not. numd(c_array(cnum)%cp_045) ) then
                call errmsg(5,dumi,'cp_045',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_046') ) then
          if( .not. numd(c_array(cnum)%cp_046) ) then
                call errmsg(5,dumi,'cp_046',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_047') ) then
          if( .not. numd(c_array(cnum)%cp_047) ) then
                call errmsg(5,dumi,'cp_047',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_048') ) then
          if( .not. numd(c_array(cnum)%cp_048) ) then
                call errmsg(5,dumi,'cp_048',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_049') ) then
          if( .not. numd(c_array(cnum)%cp_049) ) then
                call errmsg(5,dumi,'cp_049',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_050') ) then
          if( .not. numd(c_array(cnum)%cp_050) ) then
                call errmsg(5,dumi,'cp_050',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_051') ) then
          if( .not. numd(c_array(cnum)%cp_051) ) then
                call errmsg(5,dumi,'cp_051',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_052') ) then
          if( .not. numd(c_array(cnum)%cp_052) ) then
                call errmsg(5,dumi,'cp_052',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_053') ) then
          if( .not. numd(c_array(cnum)%cp_053) ) then
                call errmsg(5,dumi,'cp_053',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_054') ) then
          if( .not. numd(c_array(cnum)%cp_054) ) then
                call errmsg(5,dumi,'cp_054',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_055') ) then
          if( .not. numd(c_array(cnum)%cp_055) ) then
                call errmsg(5,dumi,'cp_055',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_056') ) then
          if( .not. numd(c_array(cnum)%cp_056) ) then
                call errmsg(5,dumi,'cp_056',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_057') ) then
          if( .not. numd(c_array(cnum)%cp_057) ) then
                call errmsg(5,dumi,'cp_057',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_058') ) then
          if( .not. numd(c_array(cnum)%cp_058) ) then
                call errmsg(5,dumi,'cp_058',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_059') ) then
          if( .not. numd(c_array(cnum)%cp_059) ) then
                call errmsg(5,dumi,'cp_059',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_060') ) then
          if( .not. numd(c_array(cnum)%cp_060) ) then
                call errmsg(5,dumi,'cp_060',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_061') ) then
          if( .not. numd(c_array(cnum)%cp_061) ) then
                call errmsg(5,dumi,'cp_061',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_062') ) then
          if( .not. numd(c_array(cnum)%cp_062) ) then
                call errmsg(5,dumi,'cp_062',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_063') ) then
          if( .not. numd(c_array(cnum)%cp_063) ) then
                call errmsg(5,dumi,'cp_063',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_064') ) then
          if( .not. numd(c_array(cnum)%cp_064) ) then
                call errmsg(5,dumi,'cp_064',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_065') ) then
          if( .not. numd(c_array(cnum)%cp_065) ) then
                call errmsg(5,dumi,'cp_065',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_066') ) then
          if( .not. numd(c_array(cnum)%cp_066) ) then
                call errmsg(5,dumi,'cp_066',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_067') ) then
          if( .not. numd(c_array(cnum)%cp_067) ) then
                call errmsg(5,dumi,'cp_067',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_068') ) then
          if( .not. numd(c_array(cnum)%cp_068) ) then
                call errmsg(5,dumi,'cp_068',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_069') ) then
          if( .not. numd(c_array(cnum)%cp_069) ) then
                call errmsg(5,dumi,'cp_069',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_070') ) then
          if( .not. numd(c_array(cnum)%cp_070) ) then
                call errmsg(5,dumi,'cp_070',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_071') ) then
          if( .not. numd(c_array(cnum)%cp_071) ) then
                call errmsg(5,dumi,'cp_071',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_072') ) then
          if( .not. numd(c_array(cnum)%cp_072) ) then
                call errmsg(5,dumi,'cp_072',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_073') ) then
          if( .not. numd(c_array(cnum)%cp_073) ) then
                call errmsg(5,dumi,'cp_073',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_074') ) then
          if( .not. numd(c_array(cnum)%cp_074) ) then
                call errmsg(5,dumi,'cp_074',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_075') ) then
          if( .not. numd(c_array(cnum)%cp_075) ) then
                call errmsg(5,dumi,'cp_075',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_076') ) then
          if( .not. numd(c_array(cnum)%cp_076) ) then
                call errmsg(5,dumi,'cp_076',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_077') ) then
          if( .not. numd(c_array(cnum)%cp_077) ) then
                call errmsg(5,dumi,'cp_077',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_078') ) then
          if( .not. numd(c_array(cnum)%cp_078) ) then
                call errmsg(5,dumi,'cp_078',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_079') ) then
          if( .not. numd(c_array(cnum)%cp_079) ) then
                call errmsg(5,dumi,'cp_079',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_080') ) then
          if( .not. numd(c_array(cnum)%cp_080) ) then
                call errmsg(5,dumi,'cp_080',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_081') ) then
          if( .not. numd(c_array(cnum)%cp_081) ) then
                call errmsg(5,dumi,'cp_081',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_082') ) then
          if( .not. numd(c_array(cnum)%cp_082) ) then
                call errmsg(5,dumi,'cp_082',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_083') ) then
          if( .not. numd(c_array(cnum)%cp_083) ) then
                call errmsg(5,dumi,'cp_083',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_084') ) then
          if( .not. numd(c_array(cnum)%cp_084) ) then
                call errmsg(5,dumi,'cp_084',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_085') ) then
          if( .not. numd(c_array(cnum)%cp_085) ) then
                call errmsg(5,dumi,'cp_085',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_086') ) then
          if( .not. numd(c_array(cnum)%cp_086) ) then
                call errmsg(5,dumi,'cp_086',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_087') ) then
          if( .not. numd(c_array(cnum)%cp_087) ) then
                call errmsg(5,dumi,'cp_087',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_088') ) then
          if( .not. numd(c_array(cnum)%cp_088) ) then
                call errmsg(5,dumi,'cp_088',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_089') ) then
          if( .not. numd(c_array(cnum)%cp_089) ) then
                call errmsg(5,dumi,'cp_089',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_090') ) then
          if( .not. numd(c_array(cnum)%cp_090) ) then
                call errmsg(5,dumi,'cp_090',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_091') ) then
          if( .not. numd(c_array(cnum)%cp_091) ) then
                call errmsg(5,dumi,'cp_091',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_092') ) then
          if( .not. numd(c_array(cnum)%cp_092) ) then
                call errmsg(5,dumi,'cp_092',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_093') ) then
          if( .not. numd(c_array(cnum)%cp_093) ) then
                call errmsg(5,dumi,'cp_093',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_094') ) then
          if( .not. numd(c_array(cnum)%cp_094) ) then
                call errmsg(5,dumi,'cp_094',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_095') ) then
          if( .not. numd(c_array(cnum)%cp_095) ) then
                call errmsg(5,dumi,'cp_095',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_096') ) then
          if( .not. numd(c_array(cnum)%cp_096) ) then
                call errmsg(5,dumi,'cp_096',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_097') ) then
          if( .not. numd(c_array(cnum)%cp_097) ) then
                call errmsg(5,dumi,'cp_097',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_098') ) then
          if( .not. numd(c_array(cnum)%cp_098) ) then
                call errmsg(5,dumi,'cp_098',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_099') ) then
          if( .not. numd(c_array(cnum)%cp_099) ) then
                call errmsg(5,dumi,'cp_099',dumr,dumd)
          end if
       elseif(  matchs_exact('cp_100') ) then
          if( .not. numd(c_array(cnum)%cp_100) ) then
                call errmsg(5,dumi,'cp_100',dumr,dumd)
          end if
       elseif(  matchs_exact('atol') ) then
          if( .not. numd(c_array(cnum)%atol) ) then
                call errmsg(5,dumi,'atol',dumr,dumd)
          end if
       elseif(  matchs_exact('atol1') ) then
          if( .not. numd(c_array(cnum)%atol1) ) then
                call errmsg(5,dumi,'atol1',dumr,dumd)
          end if
       elseif(  matchs_exact('rtol') ) then
          if( .not. numd(c_array(cnum)%rtol) ) then
                call errmsg(5,dumi,'rtol',dumr,dumd)
          end if
       elseif(  matchs_exact('rtol1') ) then
          if( .not. numd(c_array(cnum)%rtol1) ) then
                call errmsg(5,dumi,'rtol1',dumr,dumd)
          end if
       elseif(  matchs_exact('xtol') ) then
          if( .not. numd(c_array(cnum)%xtol) ) then
                call errmsg(5,dumi,'xtol',dumr,dumd)
          end if
       elseif(  matchs_exact('xtol1') ) then
          if( .not. numd(c_array(cnum)%xtol1) ) then
                call errmsg(5,dumi,'xtol1',dumr,dumd)
          end if
       elseif(  matchs_exact('miter') ) then
          if( .not. numi(c_array(cnum)%miter) ) then
                call errmsg(5,dumi,'miter',dumr,dumd)
          end if
       elseif(  matchs_exact('gpp') ) then
          if( .not. numi(c_array(cnum)%gpp) ) then
                call errmsg(5,dumi,'gpp',dumr,dumd)
          end if
       elseif(  matchs_exact('method') ) then
          if( .not. numi(c_array(cnum)%method) ) then
                call errmsg(5,dumi,'method',dumr,dumd)
          end if
       elseif(  matchs_exact('dstep') ) then
          if( .not. numi(c_array(cnum)%st_it(1) )) then
                call errmsg(5,dumi,'dstep',dumr,dumd)
          end if
       elseif(  matchs_exact('diter') ) then
          if( .not. numi(c_array(cnum)%st_it(2) )) then
                call errmsg(5,dumi,'diter',dumr,dumd)
          end if
       elseif(  matchs_exact('delem') ) then
          if( .not. numi(c_array(cnum)%st_it(3) )) then
                call errmsg(5,dumi,'delem',dumr,dumd)
          end if
       elseif(  matchs_exact('tang_calc') ) then
          if( .not. numi(c_array(cnum)%tang_calc) ) then
                call errmsg(5,dumi,'tang_calc',dumr,dumd)
          end if
       elseif(  matchs_exact('method') ) then
          if( .not. numi(c_array(cnum)%method) ) then
                call errmsg(5,dumi,'method',dumr,dumd)
          end if
       elseif(  matchs_exact('solver') ) then
          if( .not. label(dumi) ) then
                call errmsg(5,dumi,'solver',dumr,dumd)
          else
                lab = ' '
                call entits(lab,nc)
                if( lab(1:nc) .eq. 'nr') then
                   c_array(cnum)%solver = .true.
                elseif( lab(1:nc) .eq. 'tr') then
                   c_array(cnum)%solver = .false.
                else
                   call errmsg(364,dumi,dums,dumr,dumd)
                end if
          end if
       elseif(  matchs_exact('strategy') ) then
          if( .not. label(dumi) ) then
                call errmsg(5,dumi,'strategy',dumr,dumd)
          else
                lab = ' '
                call entits(lab,nc)
                if( lab(1:nc) .eq. 'geom') then
                   c_array(cnum)%strategy = .true.
                elseif( lab(1:nc) .eq. 'cubic') then
                   c_array(cnum)%strategy = .false.
                else
                   call errmsg(364,dumi,dums,dumr,dumd)
                end if
          end if
       elseif(  matchs_exact('gpall') ) then
          if( .not. label(dumi) ) then
                call errmsg(5,dumi,'gpall',dumr,dumd)
          else
                lab = ' '
                call entits(lab,nc)
                if( lab(1:nc) .eq. 'on') then
                   c_array(cnum)%gpall = .true.
                elseif( lab(1:nc) .eq. 'off') then
                   c_array(cnum)%gpall = .false.
                else
                   call errmsg(364,dumi,dums,dumr,dumd)
                end if
          end if
       elseif(  endcrd(dumi) ) then
             reading = .false.
             cycle
       elseif(  matchs(',',1) ) then 
             call readsc()
       else
             call entits(lab,nc)
             call errmsg(355,dumi,lab(1:nc),dumr,dumd)
             
             call scan()
             cycle
       end if
      end do  !  big do while
c
      call finalize_new_crystal(cnum, out)
c
c              let further input systems know we have a valid crystal
c
      defined_crystal = .true.
c      call print_crystal(cnum)
      chkparam = cnum
      sbflg1 = .true.
      sbflg2 = .true.
      return
c      
      end 
