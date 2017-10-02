c
c            mm10_output  moved to mm10_a.f to allow expanded use
c            of inlining.
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_history_locs             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/30/2017 rhd               *
c     *                                                              *
c     *  set locations (indexes) of various data within the          *
c     *  history vectory for a single integration point.             *
c     *  makes changing history much cleaner                         *
c     *                                                              *
c     *  WARP3D insures that this routine is called before any mm10  *
c     *  routine is used at runtime                                  *
c     *                                                              *
c     *  usually called once per analysis                            *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_history_locs
      use mm10_defs, only : indexes_common, index_crys_hist,
     & num_common_indexes, num_crystal_terms, length_crys_hist,
     & one_crystal_hist_size, common_hist_size, length_comm_hist
      use crystal_data, only : c_array, crystal_input,
     &            data_offset
      use global_data  ! old common.main
      use main_data
      implicit none
c
      integer :: i, j, k, ncrystals, cnum, osn, num, ecount,
     &  total_hist_size, crystal,
     &  loc_start, loc_cauchy, loc_euler, loc_pls_R, loc_uddt,
     &  loc_els_eps, loc_cur_slip_incr, loc_tau_tilde, loc_user_hist,
     &  loc_tt_rate, loc_ep, loc_ed, indev, outdev, idummy, jdummy,
     &  nslip, num_hard, cur_slip, cur_hard
      logical :: local_debug, use_max
c
c
c              the top part of the history vectory contains data
c              not dependent on the number of crystals defined
c              at an integration point.
c
c              These data are followed by a set of data for each crystal.
c              The number, type, layout of data for each crystalis the
c              same.
c
c              1. Data not dependent of crystals
c
c               1 -- [D] 6 x 6 non-symmetric
c               2 -- grad_Fe(3,3,3)
c               3 -- R(3,3) from F = R U
c               4 -- work density, plastic work density, equiv eps-pls
c               5 -- slip history: sum over time of (signed) plastic
c                    slip on each slip system. Makes sense for material
c                    w/ 1 crystal - slip on on a system is summed for
c                    each crystal then averaged w/o accounting for
c                    crystal orientations
c
c               for current max slip systems = 48, this top part
c               of the history vector uses: 123 dP words
c
c
c               allocate variable sized arrays to be less than max
c               note: this code is from avg_cry_elast_props
c               in mod_crystals.f
c
      local_debug = .false.
      call iodevn( indev, outdev, idummy, jdummy )
c
      num_common_indexes = 5
      num_crystal_terms  = 11
c
c                run through materials employed in model.
c                find the CP materials, determine number of variables
c
      cur_slip = 0
      cur_hard = 0
c
      do i = 1, nummat
       if( matprp(9,i) /= 10 ) cycle
       matprp(1,i) = 0.0  !  matprp is single precision
       matprp(2,i) = 0.0
       ecount = 0
       do j = 1, noelem  ! all in model
         if( iprops(38,j)  /=  i) cycle   ! element has this material?
           ecount = ecount + 1
            ncrystals = imatprp(101, i)
            do k = 1, ncrystals
              if( imatprp(104,i) == 1) then ! crystal number
                cnum = imatprp(105,i)
              elseif( imatprp(104,i) == 2 ) then
                osn  = data_offset(j)
                cnum = crystal_input(osn,k)
c                  couldn't do this earlier, so check here
                if( (cnum > max_crystals) .or. (cnum < 0) ) then
                  write (outdev,9505) cnum
                  call die_gracefully
                 end if
              else ! something bad wrong
                write(outdev,9502)
                call die_gracefully
              end if
              if( cur_slip < c_array(cnum)%nslip )
     &             cur_slip = c_array(cnum)%nslip
              if( cur_hard < c_array(cnum)%num_hard )
     &             cur_hard = c_array(cnum)%num_hard
            end do ! over k, crystals
        end do ! over j ( noelem)
      end do  ! over i, materials
c
c              consistency checks
c
      if( cur_hard > max_uhard ) then
        write(outdev,9200) 3
        call die_gracefully
      end if
      if( cur_hard == 0 ) then
        write(outdev,9200) 5
        call die_gracefully
      end if
      if( cur_slip > max_slip_sys ) then
        write(outdev,9200) 4
        call die_gracefully
      end if
c
c              use max value if either one is close
c
      if( (cur_hard == max_uhard) .or.
     &    (cur_slip == max_slip_sys) ) then
        use_max = .true.
      else
        use_max = .false.
      end if
      num_hard = cur_hard
      nslip    = cur_slip
c
      if( .not. allocated( indexes_common ) )
     &   allocate( indexes_common(num_common_indexes,2) )
      if( .not. allocated( index_crys_hist ) )
     &   allocate( index_crys_hist(max_crystals,num_crystal_terms,2) )
      if( .not. allocated( length_comm_hist ) )
     &   allocate( length_comm_hist(num_common_indexes) )
      if( .not. allocated( length_crys_hist ) )
     &   allocate( length_crys_hist(num_crystal_terms) )
c
c              length of common history. length of common history
c
      length_comm_hist(1) = 36
      length_comm_hist(2) = 27
      length_comm_hist(3) = 9
      length_comm_hist(4) = 3
      if( use_max ) then
        length_comm_hist(5) = max_slip_sys
      else
        length_comm_hist(5) = nslip
      end if
      common_hist_size = 0
      do i = 1, num_common_indexes
        common_hist_size = common_hist_size
     &                           + length_comm_hist(i)
      end do
c
c              start, last index for quantities
c
      indexes_common(1,1) = 1
      indexes_common(2,1) = indexes_common(1,1) + length_comm_hist(1)
      indexes_common(3,1) = indexes_common(2,1) + length_comm_hist(2)
      indexes_common(4,1) = indexes_common(3,1) + length_comm_hist(3)
      indexes_common(5,1) = indexes_common(4,1) + length_comm_hist(4)
c
      indexes_common(1,2) = length_comm_hist(1)
      indexes_common(2,2) = indexes_common(1,2) + length_comm_hist(2)
      indexes_common(3,2) = indexes_common(2,2) + length_comm_hist(3)
      indexes_common(4,2) = indexes_common(3,2) + length_comm_hist(4)
      indexes_common(5,2) = indexes_common(4,2) + length_comm_hist(5)
c
c              2. Data for each crystal specified at the integraion
c                 point
c
c                 index_crys_hist(max_crystals,num_crys_terms,2)
c
c                 use: starting location in history vector for
c                      crystal 42, elastic strains(term=5) is

c                      estrain_start = index_crys_hist(42,5,1)
c                      estrain_end   = index_crys_hist(42,5,2)
c
c                |<- crystal term #
c             %  1 -- unrotated Cauchy stress ( 6 x 1 )
c             %  2 -- updated (current) Euler angles ( 3 x 1 )
c             %  3 -- plastic rotation tensor ( 3 x 3 )
c             %  4 -- uddt. unrotation deformation tensor ( 6 x 1 )
c             %  5 -- elastic lattice strain tensor ( 6 x 1 )
c                6 -- current slip increment for each slip system
c                     ( max_slip_sys x 1 )
c                7 -- tau_tilde ( max_uhard x 1 )
c                8 -- more user history for the crystal
c                     ( max_uhard x 1 )
c                9 -- tt_rate. ( max_uhard x 1 ). rate of change of
c                     hardening variables. used for prediction/init of
c                     local NR
c
c               Current size per crystal:
c                   30 (=% terms) + max_slip_sys + 2 * max_uhard + max_uhard
c
c               more user history currently ( 14 x 1 )
c                  1 -- time OR tau_y  ( 1 x 1 ) [tau_y is for MTS model]
c                  2 -- mu_harden (1 x 1 ) [mu_harden is for MTS model]
c                  3 -- unused (3 x 1 )
c                  6 -- maxslip / np1%tinc ( 1 x 1 ) = maximum sliprate
c                                                     over all systems
c                  7 -- dble(sysID) ( 1 x 1 ) = which system has the
c                                               maximum rate
c                  8 -- dble(numAct) ( 1 x 1 ) = number of systems
c                                          with rate >= 0.1*max_rate
c                  9 -- unused (2 x 1 )
c                 - Variables for equivalent power-law creep -
c                 11 -- ec_dot = np1%p_strain_inc/np1%tinc ( 1 x 1 )
c                 12 -- n_eff ( 1 x 1 )
c                 13 -- s_trace (1 x 1)
c                 14 -- B_eff ( 1 x 1 )
c                 15 -- ed_dot = diffusion rate ( 1 x 1 )
c                 16-21 -- cp strain rate tensor 6x1
c                 22-27 -- diff strain rate tensor 6x1
c
c
      length_crys_hist(1) = 6
      length_crys_hist(2) = 3
      length_crys_hist(3) = 9
      length_crys_hist(4) = 6
      length_crys_hist(5) = 6
      if( use_max ) then
        length_crys_hist(6) = max_slip_sys
        length_crys_hist(7) = max_uhard
        length_crys_hist(8) = max_uhard
        length_crys_hist(9) = max_uhard
        length_crys_hist(10) = 6
        length_crys_hist(11) = 6
      else
        length_crys_hist(6) = nslip
        length_crys_hist(7) = num_hard
        length_crys_hist(8) = 15
        length_crys_hist(9) = num_hard
        length_crys_hist(10) = 6
        length_crys_hist(11) = 6
      end if
c
c              length of one crystal history
c
      one_crystal_hist_size = 0
      do i = 1, num_crystal_terms
        one_crystal_hist_size = one_crystal_hist_size
     &                           + length_crys_hist(i)
      end do
      total_hist_size   = indexes_common(num_common_indexes,2) +
     &                    ( max_crystals *  one_crystal_hist_size )
c
      if( local_debug ) then
       write(outdev,*) "... inside mm10_set_history_locs ..."
       write(outdev,9100) max_crystals, max_slip_sys, max_uhard,
     &                    num_common_indexes, num_crystal_terms,
     &                    one_crystal_hist_size, total_hist_size
       write(outdev,*) "      .... indexes_common ...."
       write(outdev,9110) (i, indexes_common(i,1),indexes_common(i,2),
     &                     i = 1, num_common_indexes )
      end if
c
      do crystal = 1, max_crystals
c
        loc_start         = ( indexes_common(num_common_indexes,2) + 1 )
     &                      + (crystal-1)*one_crystal_hist_size
        loc_cauchy        = loc_start
        loc_euler         = loc_cauchy + length_crys_hist(1)
        loc_pls_R         = loc_euler + length_crys_hist(2)
        loc_uddt          = loc_pls_R + length_crys_hist(3)
        loc_els_eps       = loc_uddt + length_crys_hist(4)
        loc_cur_slip_incr = loc_els_eps + length_crys_hist(5)
        loc_tau_tilde     = loc_cur_slip_incr  + length_crys_hist(6)
        loc_user_hist     = loc_tau_tilde      + length_crys_hist(7)
        loc_tt_rate       = loc_user_hist      + length_crys_hist(8)
        loc_ep            = loc_tt_rate        + length_crys_hist(9)
        loc_ed            = loc_ep             + length_crys_hist(10)
c
        index_crys_hist(crystal,1,1) = loc_cauchy
        index_crys_hist(crystal,1,2) = loc_cauchy +
     &                                 length_crys_hist(1)-1
c
        index_crys_hist(crystal,2,1) = loc_euler
        index_crys_hist(crystal,2,2) = loc_euler + length_crys_hist(2)-1
c
        index_crys_hist(crystal,3,1) = loc_pls_R
        index_crys_hist(crystal,3,2) = loc_pls_R + length_crys_hist(3)-1
c
        index_crys_hist(crystal,4,1) = loc_uddt
        index_crys_hist(crystal,4,2) = loc_uddt + length_crys_hist(4)-1
c
        index_crys_hist(crystal,5,1) = loc_els_eps
        index_crys_hist(crystal,5,2) = loc_els_eps +
     &                                 length_crys_hist(5)-1
c
        index_crys_hist(crystal,6,1) = loc_cur_slip_incr
        index_crys_hist(crystal,6,2) = loc_cur_slip_incr +
     &                                 length_crys_hist(6)-1
c
        index_crys_hist(crystal,7,1) = loc_tau_tilde
        index_crys_hist(crystal,7,2) = loc_tau_tilde +
     &                                 length_crys_hist(7)-1
c
        index_crys_hist(crystal,8,1) = loc_user_hist
        index_crys_hist(crystal,8,2) = loc_user_hist +
     &                                 length_crys_hist(8)-1
c
        index_crys_hist(crystal,9,1) = loc_tt_rate
        index_crys_hist(crystal,9,2) = loc_tt_rate   +
     &                                 length_crys_hist(9)-1
c
        index_crys_hist(crystal,10,1) = loc_ep
        index_crys_hist(crystal,10,2) = loc_ep   +
     &                                 length_crys_hist(10)-1
c
        index_crys_hist(crystal,11,1) = loc_ed
        index_crys_hist(crystal,11,2) = loc_ed   +
     &                                 length_crys_hist(11)-1
c
c                consistency check
c
        if( local_debug ) then
          write(outdev,*) "      .... index_crys_hist: crystal: ",
     &                    crystal
          write(outdev,9120)  ( i, index_crys_hist(crystal,i,1),
     &                          index_crys_hist(crystal,i,2),
     &                          i = 1, num_crystal_terms )
        end if
c
        if( index_crys_hist(crystal,num_crystal_terms,2) ==
     &    indexes_common(num_common_indexes,2) +
     &    crystal * one_crystal_hist_size  ) cycle
        write(outdev,9200) 1
        call die_gracefully
c
      end do  ! over crystals
c
c                final consistency check
c
      if( index_crys_hist(max_crystals,num_crystal_terms,2)
     &      /= total_hist_size ) then
        write(outdev,9200) 2
        call die_gracefully
      end if
c
      return
c
 9100 format(10x,"max_crys, max_slip, max_uhard:      ",3i6,
     & /,    10x,"num_common_indexes, num_crys_terms: ",2i6,
     & /,    10x,"one_crys_his_size, total_hist_size: ", i6, i10 )
 9110 format(10x,i2,2i6)
 9120 format( 10x,i2,2i10)
 9200 format(/,1x">>>>> FATAL ERROR:",
     &       /,1x"      routine mm10_set_history_locs @",i2,
     &       /,1x"      aborting job..." )

 9502 format(/,1x,
     & '>>>> FATAL ERROR: unexpected input type in',
     &    ' routine  mm10_set_history_locs',
     &   /,    ' Aborting.',/)
 9505 format(/,1x,
     & '>>>> FATAL ERROR: invalid crystal number:',i10,' detected in',
     &    ' routine  mm10_set_history_locs',
     &   /,    ' Aborting.',/)
c
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_state_sizes              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 6/15/2016 tjt               *
c     *                                                              *
c     *    called by oustates to get # states output values per      *
c     *    element to reflect various options                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_state_sizes( info_vector,
     &                           ncrystals, user_output_opt1,
     &                           user_output_opt2, iout  )
      use mm10_defs, only : length_crys_hist, length_comm_hist
      implicit none
      include 'param_def'
      integer :: ncrystals, iout, user_output_opt1, user_output_opt2,
     &           info_vector(*)
c
      logical, parameter :: here_debug = .false.
      integer :: num_output_values, n
c
c        set info data for state variable output only
c
c
c         1        number of state variables per integration
c                  point to be output based on user
c                  options
c
c         user_output_opt1, 2 are values set by user on the
c
c             output... states type <int> <int> command
c
c         They indicate which grouping of states values to output
c         and for which crystal: all or just 1
c
c         Default values = 6 mean all values for all crystals of
c         all CP materials in FE model (set by oudrive))
c
c
      if( here_debug ) then
        write(iout,9000) ncrystals, user_output_opt1, user_output_opt2
      end if
c
      if( user_output_opt1 <= 0 .or. user_output_opt1 > 7 ) then
        write(iout,9010) user_output_opt1
        user_output_opt1 = 1
      end if
c
      info_vector(1) = -100000000  ! insanity prevention
c
      select case( user_output_opt1 )
       case( 1 )
        n = ncrystals
        if( user_output_opt2 /= 0 ) n = 1
        info_vector(1) = n * 3
       case( 2 )
        n = ncrystals
        if( user_output_opt2 /= 0 ) n = 1
        info_vector(1) = n * 7
       case( 3 )
        n = ncrystals
        if( user_output_opt2 /= 0 ) n = 1
        info_vector(1) = n * 12 + length_comm_hist(5) +1
       case( 4 )
        n = ncrystals
        if( user_output_opt2 /= 0 ) n = 1
        info_vector(1) = 9 + n * 12
       case( 5 )
        info_vector(1) = 18 + 3 + length_comm_hist(5)
       case( 6 )
        n = ncrystals
        if( user_output_opt2 /= 0 ) n = 1
        info_vector(1) = n * ( 3 + 4 + length_crys_hist(7)  )
       case( 7 )
        n = ncrystals
        if( user_output_opt2 /= 0 ) n = 1
        info_vector(1) = 18 + 3 + length_comm_hist(5)
     &  + n * ( 6 + 3 + 9 + 6 + length_crys_hist(6)
     &        + length_crys_hist(7) + 5 + 5
     &        + length_crys_hist(10) + length_crys_hist(11))
       case default
        write(iout,9100) 1
        call die_abort
      end select
      return
c
 9000 format(5x,".... inside mm10_set_state_sizes ....",
     & /,15x,"ncrystals, uopt1,2: ",3i6)
 9010 format(/1x,'>>>>> Warning: invalid states output type: ',i5,
     & /,     1x,'               value of 1 used.',//)
 9100 format('>> FATAL ERROR: routine  mm10_set_state_sizes @ ',i2,
     & /,    '                job terminated',//)
c
      end
