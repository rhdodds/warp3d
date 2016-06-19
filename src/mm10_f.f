c
c             routines to support various states output options
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_history_locs             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 6/14/2016 tjt               *
c     *                                                              *
c     *  set locations (indexes) of various data within the          *
c     *  history vectory for a single integration point.             *
c     *  makes changing history much cleaner                         *
c     *                                                              *
c     *  WARP3D insures that this routine is called before any mm10  *
c     *  routine is used at runtime                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_history_locs
      use mm10_defs, only : indexes_common, index_crys_hist,
     & num_common_indexes, num_crystal_terms, length_crys_hist,
     & one_crystal_hist_size, common_hist_size, length_comm_hist
      use crystal_data, only : c_array, crystal_input,
     &            data_offset
      use main_data
      implicit none
c   $add param_def
$add common.main
c
      integer :: i, j, k, ncrystals, cnum, osn, num, ecount, 
     &  total_hist_size, crystal, elnum,
     &  loc_start, loc_cauchy, loc_euler, loc_pls_R, loc_uddt,
     &  loc_els_eps, loc_cur_slip_incr, loc_tau_tilde, loc_user_hist,
     &  loc_tt_rate, indev, outdev, idummy, jdummy,
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
      local_debug = .false.
      if( local_debug ) call iodevn( indev, outdev, idummy, jdummy )
c
      num_common_indexes = 5
      num_crystal_terms  = 9
c
c     Allocate variable sized arrays to be less than max
c      Note: this code is from avg_cry_elast_props in mod_crystals.f
c
c     Run through our materials, find the CP materials, determine number of variables
      cur_slip = 0
      cur_hard = 0
      do i=1,nummat
      if (matprp(9,i) .eq. 10) then
        matprp(1,i) = 0.0
        matprp(2,i) = 0.0
        ecount = 0
        do j = 1, noelem
          if (iprops(38,j)  .eq. i) then
            ecount = ecount + 1
            ncrystals = imatprp(101, i)
            do k = 1, ncrystals
c             Get the local crystal number
              if (imatprp(104,i) .eq. 1) then
                cnum = imatprp(105,i)
              elseif (imatprp(104,i) .eq. 2) then
                osn = data_offset(elnum)
                cnum = crystal_input(osn,k)
c               Couldn't do this earlier, so check here
                if ((cnum .gt. max_crystals) .or. 
     &                (cnum .lt. 0)) then
                  write (out,'("Crystal ", i3, " not valid")')
     &                 cnum
                  call die_gracefully
                 end if
              else
                write(out,9502) 
                call die_gracefully
              end if
c              
c             change maximums
c   
              if( cur_slip .lt. c_array(cnum)%nslip ) then
                cur_slip = c_array(cnum)%nslip
              endif
              if( cur_hard .lt. c_array(cnum)%num_hard ) then
                cur_hard = c_array(cnum)%num_hard
              endif

            end do

          end if
        end do

      end if
      end do
c
c     checks
      if( cur_hard .gt. max_uhard ) then
        write(outdev,9200) 3
        call die_gracefully
      end if
      if( cur_slip .gt. max_slip_sys ) then
        write(outdev,9200) 4
        call die_gracefully
      end if
c     use max value if either one is close
      if( (cur_hard .eq. max_uhard) .or.
     &    (cur_slip .eq. max_slip_sys) ) then
        use_max = .true.
      else
        use_max = .false. 
      end if
      num_hard = cur_hard
      nslip = cur_slip
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
c      Length of common history
      length_comm_hist(1) = 36
      length_comm_hist(2) = 27
      length_comm_hist(3) = 9
      length_comm_hist(4) = 3
      if( use_max ) then
        length_comm_hist(5) = max_slip_sys
      else
        length_comm_hist(5) = nslip
      endif
c      Sum up length of common history
      common_hist_size = 0
      do i = 1,num_common_indexes
        common_hist_size = common_hist_size
     &                           + length_comm_hist(i)
      end do
c
c       * start index *
      indexes_common(1,1) = 1
      indexes_common(2,1) = indexes_common(1,1) + length_comm_hist(1)
      indexes_common(3,1) = indexes_common(2,1) + length_comm_hist(2)
      indexes_common(4,1) = indexes_common(3,1) + length_comm_hist(3)
      indexes_common(5,1) = indexes_common(4,1) + length_comm_hist(4)
c       * last index *
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
      else
        length_crys_hist(6) = nslip
        length_crys_hist(7) = num_hard
        length_crys_hist(8) = 14
        length_crys_hist(9) = num_hard
      endif
c      Sum up length of one crystal history
      one_crystal_hist_size = 0
      do i = 1,num_crystal_terms
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
     & '>>>> System error: unexpected input type in avg_elast_props!',
     &       ' Aborting.'/)
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
$add param_def
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
        info_vector(1) = n * 12
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
     &        + length_crys_hist(7) + 5 + 4 )
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

c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/3/2016 (rhd)                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels( max_num_crystals, user_output_opt1,
     &      user_output_opt2, max_cp_states_values, size_st_labels,
     &      state_labels, state_descriptors, iout, comment_lines,
     &      max_comment_lines, num_comment_lines,
     &      nterms_crystal_list, crystal_list )
      implicit none
$add param_def
c
c                       parameters
c
      integer ::  max_num_crystals, user_output_opt1,
     &            user_output_opt2, max_cp_states_values,
     &            size_st_labels, iout, max_comment_lines,
     &            num_comment_lines, nterms_crystal_list,
     &            crystal_list(*)
      character(len=8)  :: state_labels(size_st_labels)
      character(len=60) :: state_descriptors(size_st_labels)
      character(len=80) :: comment_lines(max_comment_lines)
c
c                       locals
c
      integer :: i
      logical, save :: do_print = .false.
      logical, parameter :: here_debug = .false.
      character(len=24)  :: sdate_time_tmp
c
c                       the CP states output has options for multiple
c                       levels of detail.
c
      if( here_debug ) write(iout,9000) max_num_crystals,
     &      user_output_opt1,
     &      user_output_opt2,  max_cp_states_values
c
      state_labels(1:size_st_labels) = " "
      state_descriptors = " "
      comment_lines     = " "
c
      num_comment_lines = 3
      call fdate( sdate_time_tmp )
      write(comment_lines(1),9040) sdate_time_tmp
      write(comment_lines(2),9050) user_output_opt1, user_output_opt2
      write(comment_lines(3),9060) max_num_crystals
c
      select case( user_output_opt1 )
      case( 1 )
        call mm10_states_labels_type_1
      case( 2 )
        call mm10_states_labels_type_2
      case( 3 )
        call mm10_states_labels_type_3
      case( 4 )
        call mm10_states_labels_type_4
      case( 5 )
        call mm10_states_labels_type_5
      case( 6 )
        call mm10_states_labels_type_6
      case( 7 )
        call mm10_states_labels_type_7
      end select   ! no default needed. value chkd earlier
c
      return

 9000 format(5x,'.... inside mm10_states_labels ....',
     & /,10x,'max_num_crystals, uopt1, uopt2, max_cp_states_values: ',
     &  4i5 )
 9010 format('>> FATAL ERROR: routine mm10_states_labels @ ',i2,
     & /,    '                job aborted',//)
 9040 format(a24)
 9050 format('==> User states type: ',i2,' option: ',i2)
 9060 format('==> Max # crystals used by any CP material: ',i5 )
c
      contains
c     ========

c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_6             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/14/2016  (tjt)               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_6
      use mm10_defs, only : length_crys_hist
      implicit none
c
      integer :: num_states_here, i, s, nc, cry_id, hard_no, slip_sys
      character(len=4) :: crystal_id
      character(len=2) :: hard_id, slip_sys_id
c
c              sanity check
c
      num_states_here =  nterms_crystal_list *
     &  ( 3 + 4 + length_crys_hist(7)  )
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c
c         **   1. for each crystal:
c                  a. updated Euler angles (3x1)
c                  b. user hardening variables (each 1 value))
c                       equivalent creep strain rate
c                       effective Norton exponent n
c                       effective_Norton_stress
c                       effective Norton B value
c                  c. tau_tilde (max_uhard x 1)
c
c         **   values are averages computed over all integration points
c 
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model
c
      s = 0
      do nc = 1, nterms_crystal_list
c
         cry_id = crystal_list(nc)
         write(crystal_id,fmt="(i4.4)") cry_id
c
         state_labels(s+1) = "EA1-" // crystal_id
         state_labels(s+2) = "EA2-" // crystal_id
         state_labels(s+3) = "EA3-" // crystal_id
c
         state_descriptors(s+1) =
     &    "Eul angle- " // "cry # " // crystal_id
         state_descriptors(s+2) = " "
         state_descriptors(s+3) = " "
         s = s + 3
c
         state_labels(s+1) = "crt-" // crystal_id 
         state_labels(s+2) = "nef-" // crystal_id 
         state_labels(s+3) = "sig-" // crystal_id
         state_labels(s+4) = "Bef-" // crystal_id 
         state_descriptors(s+1) = "creep rate cry # " // crystal_id
         state_descriptors(s+2) = "eff_Norton_n"
         state_descriptors(s+3) = "eff Norton_stress"
         state_descriptors(s+4) = "Norton_B_eff"
         s = s + 4         
c
         do hard_no = 1, length_crys_hist(7)
           write(hard_id,fmt="(i2.2)") hard_no
           state_labels(s+hard_no) = "h" // hard_id // "-" 
     &      // crystal_id
         end do
         state_descriptors(s+1) = "hard val cry # " 
     &    // crystal_id         
         s = s + length_crys_hist(7)
c
      end do
c
      if( s .ne. num_states_here ) then
         write(iout,9030) 3
         call die_gracefully
      end if
c
      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i5,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_6 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_6
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_7             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/16/2016  (tjt)               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_7
      use mm10_defs, only : length_crys_hist, length_comm_hist
      implicit none
c
      integer :: num_states_here, i, s, nc, cry_id, hard_no, slip_sys
      character(len=4) :: crystal_id
      character(len=2) :: hard_id, slip_sys_id
c
c              sanity check
c
      num_states_here = 21 + length_comm_hist(5) 
     &  + nterms_crystal_list *
     &  ( 6 + 3 + 9 + 6 + length_crys_hist(6)
     &      + length_crys_hist(7) + 5 + 4 )
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c
c              1. 21+nslip common values not for a specific crystal
c                  a. -curl(Fe^-1) psuedo Nye tensor (3x3)
c                       gradFe (3x3x3) averaged over all element
c                       integration points then 3x3 Nye computed
c                  b. R(3,3) from F = R U
c                  c. work density (1x1)
c                  d. plastic work density (1x1
c                  e. equiv eps-pls (1x1
c                  f. slip history: sum over time of (signed) plastic
c                     slip on each slip system. Makes sense for material
c                     w/ 1 crystal - slip on on a system is summed for
c                     each crystal then averaged w/o accounting for
c                     crystal orientations (max_slip_sys x1)

c         **   2. for each crystal:
c                  a. unrotated cauchy stress (6x1)
c                  b. updated Euler angles (3x1)
c                  c. Rp (3x3)
c                  d. elastic lattice strain (6x1)
c                  e. current increment each slip sys (max_slip_sys x 1)
c                  f. tau_tilde (max_uhard x 1)
c                  g. user hardening variables (each 1 value))
c                      (1) time OR tau_y [tau_y for MTS model]
c                      (2) mu_harden [mu_harden for MTS model]
c                      (3) max slip rate over all systems
c                      (4) system w/ max rate
c                      (5) # systems w/ rate >= 0.1*max_rate
c                      (6) equivalent creep strain rate
c                      (7) effective Norton exponent n
c                      (8) effective_Norton_stress
c                      (9) effective Norton B value
c
c         **   values are averages computed over all integration points
c 
      s = 0
      state_labels(s+1) = "nye-11"    
      state_labels(s+2) = "nye-21"    
      state_labels(s+3) = "nye-31"    
      state_labels(s+4) = "nye-12"    
      state_labels(s+5) = "nye-22"    
      state_labels(s+6) = "nye-32"    
      state_labels(s+7) = "nye-13"    
      state_labels(s+8) = "nye-23"    
      state_labels(s+9) = "nye-33" 
      state_descriptors(s+1) = "pseudo Nye tensor"         
      s = s + 9
c
      state_labels(s+1) = "R-11"
      state_labels(s+2) = "R-21"
      state_labels(s+3) = "R-31"
      state_labels(s+4) = "R-12"
      state_labels(s+5) = "R-22"
      state_labels(s+6) = "R-32"
      state_labels(s+7) = "R-13"
      state_labels(s+8) = "R-23"
      state_labels(s+9) = "R-33"
      state_descriptors(s+1) = "R of RU rotation"
      s = s + 9
c
      state_labels(s+1) = "U-total"
      state_labels(s+2) = "U-plast"
      state_labels(s+3) = "eq-epspl"
      state_descriptors(s+1) = "total work density"
      state_descriptors(s+2) = "plastic work density"
      state_descriptors(s+3) = "equiv plastic eps"
      s = s + 3
c
      do i = 1, length_comm_hist(5) 
         write(state_labels(s+i), 9000) i
      end do
      state_descriptors(s+1) = "integrated slp  ea sys"      
      s = s + length_comm_hist(5) 
c
      if( s .ne. (21 + length_comm_hist(5)) ) then
        write(iout,9030) 2
        call die_gracefully
      end if
c
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model
c
      do nc = 1, nterms_crystal_list
c
         cry_id = crystal_list(nc)
         write(crystal_id,fmt="(i4.4)") cry_id
c
         state_labels(s+1) = "ur1-" // crystal_id
         state_labels(s+2) = "ur2-" // crystal_id
         state_labels(s+3) = "ur3-" // crystal_id
         state_labels(s+4) = "ur4-" // crystal_id
         state_labels(s+5) = "ur5-" // crystal_id
         state_labels(s+6) = "ur6-" // crystal_id
         state_descriptors(s+1) = "urot c sig " // "cry # " 
     &    // crystal_id
         state_descriptors(s+2:s+6) = " "
         s = s + 6
c
         state_labels(s+1) = "EA1-" // crystal_id
         state_labels(s+2) = "EA2-" // crystal_id
         state_labels(s+3) = "EA3-" // crystal_id
c
         state_descriptors(s+1) =
     &    "Eul angle- " // "cry # " // crystal_id
         state_descriptors(s+2) = " "
         state_descriptors(s+3) = " "
         s = s + 3
c
         state_labels(s+1) = "p11-" // crystal_id
         state_labels(s+2) = "p21-" // crystal_id
         state_labels(s+3) = "p31-" // crystal_id
         state_labels(s+4) = "p12-" // crystal_id
         state_labels(s+5) = "p22-" // crystal_id
         state_labels(s+6) = "p32-" // crystal_id
         state_labels(s+7) = "p13-" // crystal_id
         state_labels(s+8) = "p23-" // crystal_id
         state_labels(s+9) = "p33-" // crystal_id
         state_descriptors(s+1) = "pls rot cry # " // crystal_id
         s = s + 9
c
         state_labels(s+1) = "e11-" // crystal_id
         state_labels(s+2) = "e22-" // crystal_id
         state_labels(s+3) = "e33-" // crystal_id
         state_labels(s+4) = "e13-" // crystal_id
         state_labels(s+5) = "e23-" // crystal_id
         state_labels(s+6) = "e12-" // crystal_id
         state_descriptors(s+1) = "latt. eps cry # " // crystal_id
         state_descriptors(s+2:s+6) = "  "
         s = s + 6
c
         do slip_sys = 1, length_crys_hist(6)
           write(slip_sys_id,fmt="(i2.2)") slip_sys
           state_labels(s+slip_sys) = "s" // slip_sys_id // "-"
     &            // crystal_id
         end do
         state_descriptors(s+1) = "slip incr cry # " 
     &            // crystal_id
         s = s + length_crys_hist(6)
c
         do hard_no = 1, length_crys_hist(7)
           write(hard_id,fmt="(i2.2)") hard_no
           state_labels(s+hard_no) = "h" // hard_id // "-" 
     &      // crystal_id
         end do
         state_descriptors(s+1) = "hard val cry # " 
     &    // crystal_id         
         s = s + length_crys_hist(7)
c
         state_labels(s+1) = "tty-" // crystal_id 
         state_descriptors(s+1) = "time OR tau_y"
         state_labels(s+2) = "muh-"  // crystal_id 
         state_descriptors(s+2) = "MTS model mu-hard"
         state_labels(s+3) = "msr-" // crystal_id
         state_descriptors(s+3) = "max sl rate all sys"
         state_labels(s+4) = "mID-" // crystal_id
         state_descriptors(s+4) = "sys # w/ max rate"
         state_labels(s+5) = "#ac-" // crystal_id
         state_descriptors(s+5) = "no.sys >0.1*max_rate"
         s = s + 5
c
         state_labels(s+1) = "crt-" // crystal_id 
         state_labels(s+2) = "nef-" // crystal_id 
         state_labels(s+3) = "sig-" // crystal_id
         state_labels(s+4) = "Bef-" // crystal_id 
         state_descriptors(s+1) = "creep rate cry # " // crystal_id
         state_descriptors(s+2) = "eff_Norton_n"
         state_descriptors(s+3) = "eff Norton_stress"
         state_descriptors(s+4) = "Norton_B_eff"
         s = s + 4
c
      end do
c
      if( s .ne. num_states_here ) then
         write(iout,9030) 3
         call die_gracefully
      end if
c
      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i5,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_7 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_7

c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_1             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/9/2016  (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_1
      implicit none
c
      integer :: s, nc, cry_id, num_states_here
      character(len=4) :: crystal_id
c
c              sanity check
c
      num_states_here = 3 * nterms_crystal_list
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c
c              3 updated euler angles for each crystal
c
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model
c
      s = 0
      do nc = 1, nterms_crystal_list
c
         cry_id = crystal_list(nc)
         write(crystal_id,fmt="(i4.4)") cry_id
         state_labels(s+1) = "EA1-" // crystal_id
         state_labels(s+2) = "EA2-" // crystal_id
         state_labels(s+3) = "EA3-" // crystal_id
c
         state_descriptors(s+1) =
     &    "Eu angle- " // "cry # " // crystal_id
         state_descriptors(s+2) = " "
         state_descriptors(s+3) = " "
         s = s + 3
c
      end do
c
      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i5,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_1 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_1
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_2             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/10/2016  (rhd)               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_2 
      implicit none
c
      integer :: s, nc, cry_id, num_states_here
      character(len=4) :: crystal_id
c
c              sanity check
c
      num_states_here = 7 * nterms_crystal_list
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c
c              states has for each crystal
c               (1)  updated euler angles 3x1
c               (2)  scalar creep values (4x1)
c
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model
c
      s = 0
      
      do nc = 1, nterms_crystal_list
c
         cry_id = crystal_list(nc)
         write(crystal_id,fmt="(i4.4)") cry_id
c
         state_labels(s+1) = "EA1-" // crystal_id
         state_labels(s+2) = "EA2-" // crystal_id
         state_labels(s+3) = "EA3-" // crystal_id
c
         state_descriptors(s+1) =
     &    "Eul angle- " // "cry # " // crystal_id
         state_descriptors(s+2) = " "
         state_descriptors(s+3) = " "
         s = s + 3
c
         state_labels(s+1) = "crt-" // crystal_id 
         state_labels(s+2) = "nef-" // crystal_id 
         state_labels(s+3) = "sig-" // crystal_id
         state_labels(s+4) = "Bef-" // crystal_id 
         state_descriptors(s+1) = "creep rate cry # " // crystal_id
         state_descriptors(s+2) = "eff_Norton_n"
         state_descriptors(s+3) = "eff Norton_stress"
         state_descriptors(s+4) = "Norton_B_eff"
         s = s + 4
c
      end do
c
      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i5,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_2 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_2
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_3             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/9/2016 (rhd)                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_3
      implicit none
c
      integer :: s, nc, cry_id, num_states_here
      character(len=4) :: crystal_id
c
c              sanity check
c
      num_states_here = 12 * nterms_crystal_list
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c              states has for each crystal
c               (1)  updated euler angles 3x1
c               (2)  lattice strain 6x1
c               (3)  max-slip rate over all slip systems 1x1
c               (4)  id for slip system with max rate 1x1
c               (5)  number of active slip systems  1x1
c
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model
c
      s = 0
      do nc = 1, nterms_crystal_list
c
         cry_id = crystal_list(nc)
         write(crystal_id,fmt="(i4.4)") cry_id
c
         state_labels(s+1) = "EA1-" // crystal_id
         state_labels(s+2) = "EA2-" // crystal_id
         state_labels(s+3) = "EA3-" // crystal_id
c
         state_descriptors(s+1) =
     &    "Eul angle- " // "cry # " // crystal_id
         state_descriptors(s+2) = " "
         state_descriptors(s+3) = " "
         s = s + 3
c
         state_labels(s+1) = "e11-" // crystal_id
         state_labels(s+2) = "e22-" // crystal_id
         state_labels(s+3) = "e33-" // crystal_id
         state_labels(s+4) = "e13-" // crystal_id
         state_labels(s+5) = "e23-" // crystal_id
         state_labels(s+6) = "e12-" // crystal_id
         state_descriptors(s+1) = "latt. eps cry # " // crystal_id
         state_descriptors(s+2:s+6) = "  "
         s = s + 6
c
         state_labels(s+1) = "msr-" // crystal_id
         state_descriptors(s+1) = "max sl rate all sys"
         state_labels(s+2) = "mID-"  // crystal_id
         state_descriptors(s+2) = "sys # w/ max rate"
         state_labels(s+3) = "#ac-" // crystal_id
         state_descriptors(s+3) = "no.sys >0.1*max_rate"
         s = s + 3
c
      end do
c
      if( s .ne. num_states_here ) then
         write(iout,9030) 3
         call die_gracefully
      end if

      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i5,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_3 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_3
c
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_4             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/9/2016 rhd)                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_4
      implicit none
c
      integer :: s, nc, cry_id, num_states_here
      character(len=4) :: crystal_id
c
c              sanity check
c
      num_states_here = 9 + 12 * nterms_crystal_list
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c
c           (1) nye tensor computed from averaged gradient
c                over int_points gradient for each point at
c                top of history not connected to a crystal
c                3x3

c           (2) states has for each crystal
c                (1)  updated euler angles 3x1
c                (2)  lattice strain 6x1
c                (3)  max-slip rate over all slip systems 1x1
c                (4)  id for slip system with max rate
c                (5)  number of active slip systems
c                      /=0, labels for just this crystal number
c
      state_labels(1) = "nye-11" 
      state_labels(2) = "nye-21" 
      state_labels(3) = "nye-31" 
      state_labels(4) = "nye-12" 
      state_labels(5) = "nye-22" 
      state_labels(6) = "nye-32" 
      state_labels(7) = "nye-13" 
      state_labels(8) = "nye-23" 
      state_labels(9) = "nye-33" 
      state_descriptors(1) = "pseudo Nye tensor"
c
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model

      s = 9
      do nc = 1, nterms_crystal_list
c
         cry_id = crystal_list(nc)
         write(crystal_id,fmt="(i4.4)") cry_id
c
         state_labels(s+1) = "EA1-" // crystal_id
         state_labels(s+2) = "EA2-" // crystal_id
         state_labels(s+3) = "EA3-" // crystal_id
c
         state_descriptors(s+1) =
     &    "Eul angle- " // "cry # " // crystal_id
         state_descriptors(s+2) = " "
         state_descriptors(s+3) = " "
         s = s + 3
c
         state_labels(s+1) = "e11-" // crystal_id
         state_labels(s+2) = "e22-" // crystal_id
         state_labels(s+3) = "e33-" // crystal_id
         state_labels(s+4) = "e13-" // crystal_id
         state_labels(s+5) = "e23-" // crystal_id
         state_labels(s+6) = "e12-" // crystal_id
         state_descriptors(s+1) = "latt. eps cry # " // crystal_id
         state_descriptors(s+2:s+6) = "  "
         s = s + 6
c
         state_labels(s+1) = "msr-" // crystal_id
         state_descriptors(s+1) = "max sl rate all sys"
         state_labels(s+2) = "mID-"  // crystal_id
         state_descriptors(s+2) = "sys # w/ max rate"
         state_labels(s+3) = "#ac-" // crystal_id
         state_descriptors(s+3) = "no.sys >0.1*max_rate"
         s = s + 3         
c
      end do
c
      if( s .ne. num_states_here ) then ! sanity check
         write(iout,9030) 3
         call die_gracefully
      end if
c
      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i5,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_4 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_4
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels_type_5             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/16/2016 tjt                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels_type_5
      use mm10_defs, only : length_comm_hist
      implicit none
c
      integer :: num_states_here, i, s
c
      num_states_here = 21 + length_comm_hist(5)  ! only selected common values 
c                             independent of number of crystals
c
c              sanity check
c
      if( num_states_here .ne. max_cp_states_values ) then
         write(iout,9030) 1
         call die_gracefully
      end if
c
c
c              1. 21 + nslip common values not for a specific crystal
c                  a. -curl(Fe^-1) psuedo Nye tensor (3x3)
c                       gradFe (3x3x3) averaged over all element
c                       integration points then 3x3 Nye computed
c                  b. R(3,3) from F = R U
c                  c. work density (1x1)
c                  d. plastic work density (1x1
c                  e. equiv eps-pls (1x1
c                  f. slip history: sum over time of (signed) plastic
c                     slip on each slip system. Makes sense for material
c                     w/ 1 crystal - slip on on a system is summed for
c                     each crystal then averaged w/o accounting for
c                     crystal orientations (max_slip_sys x1)
c
c         **   values are averages computed over all integration points
c
      s = 0
      state_labels(s+1) = "nye-11"    
      state_labels(s+2) = "nye-21"    
      state_labels(s+3) = "nye-31"    
      state_labels(s+4) = "nye-12"    
      state_labels(s+5) = "nye-22"    
      state_labels(s+6) = "nye-32"    
      state_labels(s+7) = "nye-13"    
      state_labels(s+8) = "nye-23"    
      state_labels(s+9) = "nye-33" 
      state_descriptors(s+1) = "pseudo Nye tensor"         
      s = s + 9
c
      state_labels(s+1) = "R-11"
      state_labels(s+2) = "R-21"
      state_labels(s+3) = "R-31"
      state_labels(s+4) = "R-12"
      state_labels(s+5) = "R-22"
      state_labels(s+6) = "R-32"
      state_labels(s+7) = "R-13"
      state_labels(s+8) = "R-23"
      state_labels(s+9) = "R-33"
      state_descriptors(s+1) = "R of RU rotation"
      s = s + 9
c
      state_labels(s+1) = "U-total"
      state_labels(s+2) = "U-plast"
      state_labels(s+3) = "eq-epspl"
      state_descriptors(s+1) = "total work density"
      state_descriptors(s+2) = "plastic work density"
      state_descriptors(s+3) = "equiv plastic eps"
      s = s + 3
c
      do i = 1, length_comm_hist(5) 
         write(state_labels(s+i), 9000) i
         state_descriptors(s+i) = "integrated slip"
      end do
      s = s + length_comm_hist(5) 
c
      if( s .ne. (21 + length_comm_hist(5) ) ) then  ! sanity check
        write(iout,9030) 2
        call die_gracefully
      end if
d
      if( do_print ) then
        do i = 1, num_states_here
          write(iout,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
c
 9000 format("slip-",i2.2)
 9010 format(2x,i3,2x,a8,2x,a)
 9020 format("hard.-",i2.2)
 9030 format('>> FATAL ERROR: routine mm10_states_labels_type_5 @ ',i2,
     & /,    '                job aborted',//)
c
      end subroutine mm10_states_labels_type_5
      end subroutine mm10_states_labels

c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 6/15/2016  tjt                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values( itype, elem_states_output,
     &                               num_states, iout,
     &                               iprops, felem, hist_size,
     &                               int_points, span, history_blk,
     &                               ncrystals_blk, max_num_crystals,
     &                               user_output_opt1, user_output_opt2,
     &                               nterms_crystal_list, crystal_list )
c
      implicit none
c
c                       parameters
c
      integer :: itype, num_states, felem, iprops(*),
     &           hist_size, int_points, span, iout, ncrystals_blk,
     &           max_num_crystals, user_output_opt1, user_output_opt2,
     &           nterms_crystal_list, crystal_list(*)
#dbl      double precision ::
#sgl      real  ::
     & elem_states_output(num_states,span),
     & history_blk(hist_size,int_points,span)
c
c                       locals
c
      integer :: relem, elnum, elem_type, mat_type, nc, cry_id,
     &           num_crystals_to_process
      integer :: crystal_list_for_blk(max_num_crystals) ! automatic
      logical :: do_a_block, local_debug
      double precision :: zero
      data zero / 0.0d00 /
c
@!DIR$ ASSUME_ALIGNED elem_states_output:64,history_blk:64
c
c           build CP states values output.
c
c              itype > 0 => this is the block number. do all elements
c                           in the block
c
c              itype < 0 => this is an element number. put state
c                           values into column 1 of results.
c
      do_a_block = .true.
      if( itype < 0 ) then
         do_a_block = .false.
         elnum = -itype
      end if
c
      local_debug = .false.
      elem_type   = iprops(1)
      mat_type    = iprops(25)
c
      if( local_debug ) write(iout,9050) felem, elem_type,
     &         mat_type, int_points, span, hist_size,
     &         ncrystals_blk, max_num_crystals, user_output_opt1,
     &         user_output_opt2

c              sanity checks. crystal_list has crystal #s for states
c              output (increasing order).
c              ncrystals_blk = # of crystals for elements in this blk
c              it is thus possible for a crystal # in the list to
c              be more than ncrystals_blk.
c              build a crystal list just for this block with allowable
c              crystals so lower level can be simpler.
c              state entries = zero for crystals not present for CP
c              material of this block
c
      num_crystals_to_process = 0
      do nc = 1, nterms_crystal_list
        cry_id = crystal_list(nc)
        if( cry_id > ncrystals_blk ) exit
        num_crystals_to_process =  num_crystals_to_process + 1
        crystal_list_for_blk(nc) = cry_id
      end do
      if( num_crystals_to_process == 0 ) then ! something wrong
        write(iout,9000) 1
        call die_gracefully
      end if
c
      if( do_a_block ) then
        do relem = 1, span
           elnum = felem + relem - 1  ! absolute element number
           select case( user_output_opt1 )
           case(1)
             call mm10_states_values_1( elem_states_output(1,relem),
     &                                   history_blk )
           case(2)
             call mm10_states_values_2( elem_states_output(1,relem),
     &                                   history_blk )
           case(3)
             call mm10_states_values_3( elem_states_output(1,relem),
     &                                   history_blk )
           case(4)
             call mm10_states_values_4( elem_states_output(1,relem),
     &                                   history_blk )
           case(5)
             call mm10_states_values_5( elem_states_output(1,relem),
     &                                   history_blk )
           case(6)
             call mm10_states_values_6( elem_states_output(1,relem),
     &                                   history_blk )
           case(7)
             call mm10_states_values_7( elem_states_output(1,relem),
     &                                   history_blk )
           end select
        end do
      else
           relem = elnum + 1 - felem
           select case( user_output_opt1 )
           case(1)
             call mm10_states_values_1( elem_states_output(1,1),
     &                                   history_blk )
           case(2)
             call mm10_states_values_2( elem_states_output(1,1),
     &                                   history_blk )
           case(3)
             call mm10_states_values_3( elem_states_output(1,1),
     &                                   history_blk )
           case(4)
             call mm10_states_values_4( elem_states_output(1,1),
     &                                   history_blk )
           case(5)
             call mm10_states_values_5( elem_states_output(1,1),
     &                                   history_blk )
           case(6)
             call mm10_states_values_6( elem_states_output(1,1),
     &                                   history_blk )
           case(7)
             call mm10_states_values_7( elem_states_output(1,1),
     &                                   history_blk )
           end select
      end if
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values @: ',i2,
     &   /,  '                job termindated.'// )

 9050 format(10x,"felem, etype, mtype:  ",3i7,
     &  /,10x,   "int_pts, span, hist_size:    ",3i7,
     &  /,10x,   "ncry_blk, mxcry, uopt1, uopt2: ",4i5)
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_6                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/11/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_6( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist,
     &   length_crys_hist
c
c              parameters. history_blk passed to declare alignment
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c              locals
c
      integer :: num_states_here, s, e, sh, eh, cry_id, nc
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64
@!DIR$ ASSUME_ALIGNED hisblk:64
c
c              locations of data in the history vector for each
c              integration point obtained for locations array. see
c              routine mm10_set_history_locs and how to
c              use the index arrays.
c
c              note: one_elem_states zeroed before entry.
c
c              sanity check
c
      num_states_here = num_crystals_to_process *
     &  ( 3 + 4 + length_crys_hist(7) )
c
      if( num_states_here > num_states ) then
         write(iout,9000) 1
         write(iout,*) '.... num_states: ', num_states
         call die_gracefully
      end if
c
c             see mm10_states_labels_type_7 for detailed summary &
c             labels for states values loaded here into the output
c             vector
c
c         **   1. for each crystal:
c                  a. updated Euler angles (3x1)
c                  b. user hardening variables (each 1 value))
c                      (6) equivalent creep strain rate
c                      (7) effective Norton exponent n
c                      (8) effective_Norton_stress
c                      (9) effective Norton B value
c                  c. tau_tilde (max_uhard x 1)
c
c         **   values are averages computed over all integration points
c 
c              nterms_crystal_list will = 1 or the maximum number of
c              crystals used in all the CP materials of the FE  model

c           1. loop over all required crystals in the block list
c
      s = 1
      do nc = 1, num_crystals_to_process
c
c           number of crystal being processed
c
        cry_id = crystal_list_for_blk(nc)
        if( cry_id .gt. ncrystals_blk ) cycle
c
c           1a. updated Euler angles (3x1)
c
        e = s + 2
        sh = index_crys_hist(cry_id,2,1)
        eh = index_crys_hist(cry_id,2,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c        
c           1b. user hardening variables (each 1 value))
c                      (6) equivalent creep strain rate
c                      (7) effective Norton exponent n
c                      (8) effective_Norton_stress
c                      (9) effective Norton B value
c
        e = s + 3
        sh = index_crys_hist(cry_id,8,1) + 11 - 1
        eh = sh + 3
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 4        
c
c           1c. slip system hardening variable (max_uhard x 1)
c               [ tau tilde ]
c
        e = s + (length_crys_hist(7) - 1)
        sh = index_crys_hist(cry_id,7,1)
        eh = index_crys_hist(cry_id,7,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + length_crys_hist(7)
c
      end do ! over crystals
c
      if( s-1 .ne. num_states_here ) then
         write(iout,9000) 3
         call die_gracefully
      end if
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_6 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_6
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_7                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/16/2016 tjt              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_7( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist,
     &    length_crys_hist, length_comm_hist
c
c              parameters. history_blk passed to declare alignment
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c              locals
c
      integer :: num_states_here, s, e, sh, eh, cry_id, nc
#dbl      double precision :: nyevec(9),  gradvec(27)
#sgl      real :: nyevec(9),  gradvec(27)
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64
@!DIR$ ASSUME_ALIGNED gradvec:64, nyevec:64, hisblk:64
c
c              locations of data in the history vector for each
c              integration point obtained for locations array. see
c              routine mm10_set_history_locs and how to
c              use the index arrays.
c
c              note: one_elem_states zeroed before entry.
c
c              sanity check
c
      num_states_here = 21 + length_comm_hist(5) 
     &  + num_crystals_to_process *
     &  ( 6 + 3 + 9 + 6 + length_crys_hist(6)
     &      + length_crys_hist(7) + 5 + 4 )
c
      if( num_states_here > num_states ) then
         write(iout,9000) 1
         write(iout,*) '.... num_states: ', num_states
         call die_gracefully
      end if
c
c             see mm10_states_labels_type_6 for detailed summary &
c             labels for states values loaded here into the output
c             vector
c
c              1. 21 + nslip common values not for a specific crystal
c                  a. -curl(Fe^-1) psuedo Nye tensor (3x3)
c                       gradFe (3x3x3) averaged over all element
c                       integration points then 3x3 Nye computed
c                  b. R(3,3) from F = R U
c                  c. work density (1x1)
c                  d. plastic work density (1x1
c                  e. equiv eps-pls (1x1
c                  f. slip history: sum over time of (signed) plastic
c                     slip on each slip system. Makes sense for material
c                     w/ 1 crystal - slip on on a system is summed for
c                     each crystal then averaged w/o accounting for
c                     crystal orientations (max_slip_sys x1)

c         **   2. for each crystal:
c                  a. unrotated cauchy stress (6x1)
c                  b. updated Euler angles (3x1)
c                  c. Rp (3x3)
c                  d. elastic lattice strain (6x1)
c                  e. current increment each slip sys (max_slip_sys x 1)
c                  f. tau_tilde (max_uhard x 1)
c                  g. user hardening variables (each 1 value))
c                      (1) time OR tau_y [tau_y for MTS model]
c                      (2) mu_harden [mu_harden for MTS model]
c                      (3) max slip rate over all systems
c                      (4) system w/ max rate
c                      (5) # systems w/ rate >= 0.1*max_rate
c                      (6) equivalent creep strain rate
c                      (7) effective Norton exponent n
c                      (8) effective_Norton_stress
c                      (9) effective Norton B value
c
c         **   values are averages computed over all integration points
c
c
c             1. common values ( 69 total ))
c
c             1a. Nye tensor, averaged over int. pts.
c                 gradFe = reshape(
c                     sum(hisblk(s:e,1:int_points,relem),2 ) /
c                     dble(int_points), (/3,3,3/) )
c
      s  = 1
      e  = s + 8
      sh = indexes_common(2,1)
      eh = indexes_common(2,2)
      gradvec = zero
      call mm10_avg_states( gradvec, 1, 27, sh, eh, hisblk(1,1,relem),
     &                      hist_size, int_points )
      call mm10_states_nye( gradvec, nyevec )
      one_elem_states(s:e) = nyevec
      s = s + 9
c
c             1b. R of F = RU averaged over element int points
c
      e = s + 8
      sh = indexes_common(3,1)
      eh = indexes_common(3,2)
      call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
      s = s + 9
c
c             1c,d,e. work density (1x1), plastic work density (1x1
c                     equiv eps-pls (1x1)
c
      e = s + 2
      sh = indexes_common(4,1)
      eh = indexes_common(4,2)
      call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
      s = s + 3
c
c             1f. current (total) slip on each system. averaged over
c                 integration points
c
      sh  = indexes_common(5,1)
      eh  = indexes_common(5,2)
      e   = s + (length_comm_hist(5) - 1)
      call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
      s = s + length_comm_hist(5)
c
      if( s .ne. (21 + length_comm_hist(5) + 1) ) then
        write(iout,9000) 2
        call die_gracefully
      end if
c
c           2. loop over all required crystals in the block list
c
      do nc = 1, num_crystals_to_process
c
c           number of crystal being processed
c
        cry_id = crystal_list_for_blk(nc)
        if( cry_id .gt. ncrystals_blk ) cycle
c
c           2a. unrotated cauchy stresses (6x1)
c
        e = s + 5
        sh = index_crys_hist(cry_id,1,1)
        eh = index_crys_hist(cry_id,1,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 6
c
c           2b. updated Euler angles (3x1)
c
        e = s + 2
        sh = index_crys_hist(cry_id,2,1)
        eh = index_crys_hist(cry_id,2,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
c           2c. Rp (3x3)
c
        e = s + 8
        sh = index_crys_hist(cry_id,3,1)
        eh = index_crys_hist(cry_id,3,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 9
c
c           2d. elastic lattice strain (6x1)
c
        e = s + 5
        sh = index_crys_hist(cry_id,5,1)
        eh = index_crys_hist(cry_id,5,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 6
c
c           2e. current increment of slip each system (length_crys_hist(6) x 1)
c
        e = s + (length_crys_hist(6) -1)
        sh = index_crys_hist(cry_id,6,1)
        eh = index_crys_hist(cry_id,6,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + length_crys_hist(6)
c
c           2f. slip system hardening variable (max_uhard x 1)
c               [ tau tilde ]
c
        e = s + (length_crys_hist(7) - 1)
        sh = index_crys_hist(cry_id,7,1)
        eh = index_crys_hist(cry_id,7,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + length_crys_hist(7)
c
c           2g. user hardening variables (each 1 value))
c                      (1) time OR tau_y [tau_y for MTS model]
c                      (2) mu_harden [mu_harden for MTS model]
c                      (3) max slip rate over all systems
c                      (4) system w/ max rate
c                      (5) # systems w/ rate >= 0.1*max_rate
c                      (6) equivalent creep strain rate
c                      (7) effective Norton exponent n
c                      (8) effective_Norton_stress
c                      (9) effective Norton B value
c
        e = s + 1
        sh = index_crys_hist(cry_id,8,1) + 0
        eh = sh + 1
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 2
c
        e = s + 2
        sh = index_crys_hist(cry_id,8,1) + 6 - 1
        eh = sh + 2
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
        e = s + 3
        sh = index_crys_hist(cry_id,8,1) + 11 - 1
        eh = sh + 3
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 4
c
      end do ! over crystals
c
      if( s-1 .ne. num_states_here ) then
         write(iout,9000) 3
         call die_gracefully
      end if
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_7 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_7
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_1                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/5/2016  rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_1( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist
c
c               parameters. history_blk passed so alignment
c               can be declared
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c               locals
c
      integer :: s, e, sh, eh, num_states_here, nc, cry_id
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64, hisblk:64
c
c              locations of data in the history vector for each
c              integration point obtained for locations array. see
c              routine mm10_set_history_locs and how to
c              use the index arrays.
c
c              note: one_elem_states zeroed before entry.
c
c              for each required crystal:
c               - each of 3 Euler angles avgeraed over int points
c               - 4 x1 creep parameters averaged over int points
c
      num_states_here = num_crystals_to_process * 3
      if( num_states_here > num_states ) then ! sanity check
         write(iout,9000) 1
         write(iout,*) '.... num_states: ', num_states
         call die_gracefully
      end if
c
c           1. loop over all required crystals in the global list.
c
      s = 1
      do nc = 1, num_crystals_to_process
c
c           2a. number of crystal being processed
c
        cry_id = crystal_list(nc)
        if( cry_id .gt. ncrystals_blk ) cycle
c
c           2b. 3 Euler angles
c
        e = s + 2
        sh = index_crys_hist(cry_id,2,1)
        eh = index_crys_hist(cry_id,2,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
      end do ! over crystals
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_1 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_1
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_2                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/3/2016  rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_2( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist
c
c               parameters. history_blk passed so alignment
c               can be declared
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c               locals
c
      integer :: s, e, sh, eh, num_states_here, nc, cry_id
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64, hisblk:64
c
c              locations of data in the history vector for each
c              integration point obtained for locations array. see
c              routine mm10_set_history_locs and how to
c              use the index arrays.
c
c              note: one_elem_states zeroed before entry.
c
c              for each required crystal:
c               - each of 3 Euler angles avgeraed over int points
c               - 4 x 1 creep parameters averaged over int points
c
      num_states_here = num_crystals_to_process * 7
      if( num_states_here > num_states ) then ! sanity check
         write(iout,9000) 1
         write(iout,*) '.... num_states: ', num_states
         call die_gracefully
      end if
c
c           1. loop over all required crystals in the global list.
c
      s = 1
      do nc = 1, num_crystals_to_process
c
c           2a. number of crystal being processed
c
        cry_id = crystal_list(nc)
        if( cry_id .gt. ncrystals_blk ) cycle
c
c           2b. 3 Euler angles
c
        e = s + 2
        sh = index_crys_hist(cry_id,2,1)
        eh = index_crys_hist(cry_id,2,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
c           2c. 4 x 1 average creep values over int points
c
        e = s + 3
        sh = index_crys_hist(cry_id,8,1) + 10
        eh = sh + 3
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 4
c
      end do ! over crystals
c
      if( s-1 .ne. num_states_here ) then
         write(iout,9000) 3
         call die_gracefully
      end if
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_2 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_2
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_3                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/5/2016  rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_3( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist
c
c               parameters. history_blk passed so alignment
c               can be declared
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c               locals
c
      integer :: s, e, sh, eh, num_states_here, nc, cry_id
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64, hisblk:64
c
c              locations of data in the history vector for each
c              integration point obtained for locations array. see
c              routine mm10_set_history_locs and how to
c              use the index arrays.

c              for history contents per crystal, routine
c              mm10_set_history_locs and how to use the index arrays.
c
c              note: one_elem_states zeroed before entry.
c
c              for each required crystal:
c               - each of 3 Euler angles averaged over int points
c               - 6 x 1 lattice strains averaged over int points
c               - 3 x 1 values for max slip rates (see notes below)
c
      num_states_here = num_crystals_to_process * 12
      if( num_states_here > num_states ) then ! sanity check
         write(iout,9000) 1
         write(iout,*) '.... num_states: ', num_states
         call die_gracefully
      end if
c
c           1. loop over all required crystals in the global list.
c
      s = 1
      do nc = 1, num_crystals_to_process
c
c           2a. number of crystal being processed
c
        cry_id = crystal_list(nc)
        if( cry_id .gt. ncrystals_blk ) cycle
c
c           2b. 3 updated euler analges
c
        e = s + 2
        sh = index_crys_hist(cry_id,2,1)
        eh = index_crys_hist(cry_id,2,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
c           2c. average lattice strain over int points
c
        e = s + 5
        sh = index_crys_hist(cry_id,5,1)
        eh = index_crys_hist(cry_id,5,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 6
c
c           2d. average of (a) maximum sliprate over all systems,
c               (b) which slip system has the maximum rate
c               (c) nmber of slip systems  with rate >=
c                   0.1*max_rate lattice strain over int points
c               => meaning of (2) & (3) averaged over over int points is
c                  questionable
c
        e = s + 2
        sh = index_crys_hist(cry_id,8,1) + 5
        eh = sh + 2
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
      end do ! over crystals
c
      if( s-1 .ne. num_states_here ) then
         write(iout,9000) 3
         call die_gracefully
      end if
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_3 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_3
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_4                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/5/2016  rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_4( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist
c
c               parameters. history_blk passed so alignment
c               can be declared
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c               locals
c
      integer :: s, e, sh, eh, num_states_here, nc, cry_id
#dbl      double precision :: nyevec(9),  gradvec(27)
#sgl      real :: nyevec(9),  gradvec(27)
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64
@!DIR$ ASSUME_ALIGNED gradvec:64, nyevec:64, hisblk:64
c
c              locations of data in the history vector for each
c              integration point obtained for locations array. see
c              routine mm10_set_history_locs and how to
c              use the index arrays.

c              note: one_elem_states zeroed before entry.
c
c              sanity check
c
      num_states_here = 9 + num_crystals_to_process * 12
      if( num_states_here > num_states ) then ! sanity check
         write(iout,9000) 1
         write(iout,*) '.... num_states: ', num_states
         call die_gracefully
      end if
c
c              states vector constructed here
c                1 -> 9 nye tensor
c               for each crystal:
c               - each of 3 Euler angles avgeraed over int points
c               - 6 x 1 lattice strains averaged over int points
c               - 3 x 1 values for max slip rates (see notes below)
c
c           1. Nye tensor, averaged over int. pts.
c              gradFe = reshape(
c                     sum(hisblk(s:e,1:int_points,relem),2 ) /
c                     dble(int_points), (/3,3,3/) )

      sh = indexes_common(2,1)
      eh = indexes_common(2,2)
      gradvec = zero
      call mm10_avg_states( gradvec, 1, 27, sh, eh, hisblk(1,1,relem),
     &                      hist_size, int_points )
      call mm10_states_nye( gradvec, nyevec )
      one_elem_states(1:9) = nyevec
c

c           2. loop over all required crystals in the global list.
c              the CP material for elements in this block has
c              ncrystals_blk which may be < than the material in
c              FE model with largest number of crystals
c
      s = 10
      do nc = 1, num_crystals_to_process
c
c           2a. number of crystal being processed
c
        cry_id = crystal_list(nc)
        if( cry_id .gt. ncrystals_blk ) cycle
c
c           2b. 3 updated Euler angles
c
        e = s + 2
        sh = index_crys_hist(cry_id,2,1)
        eh = index_crys_hist(cry_id,2,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
c           2c. average lattice strain over int points
c
        e = s + 5
        sh = index_crys_hist(cry_id,5,1)
        eh = index_crys_hist(cry_id,5,2)
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 6
c
c           2d. average of (1) maximum sliprate over all systems,
c               (2) which slip system has the maximum rate
c               (3) nmber of slip systems  with rate >=
c                   0.1*max_rate lattice strain over int points
c               => meaning of (2) & (3) averaged over over int points is
c                  questionable
c
        e = s + 2
        sh = index_crys_hist(cry_id,8,1) + 5
        eh = sh + 2
        call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
        s = s + 3
c
      end do ! over crystals
c
      if( s-1 .ne. num_states_here ) then
         write(iout,9000) 3
         write(iout,*) '... s:', s
         call die_gracefully
      end if
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_4 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_4
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_5                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/16/2016 tjt              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_5( one_elem_states, hisblk )
c
      use mm10_defs, only : indexes_common, index_crys_hist,
     &                      length_comm_hist
c
c                       parameters. history_blk passed so alignment
c                       can be declared
c
      implicit none
$add param_def
#dbl      double precision ::
#sgl      real ::
     & one_elem_states(num_states), hisblk(hist_size,int_points,span)
c
c                       locals
c
      integer :: num_states_here, s, e, sh, eh
#dbl      double precision :: nyevec(9),  gradvec(27)
#sgl      real :: nyevec(9),  gradvec(27)
c
@!DIR$ ASSUME_ALIGNED one_elem_states:64
@!DIR$ ASSUME_ALIGNED gradvec:64, nyevec:64, hisblk:64

c             see mm10_states_labels_type_6 for detailed summary &
c             labels for states values loaded here into the output
c             vector
c
c             21+nslip common values not for a specific crystal
c                  a. -curl(Fe^-1) psuedo Nye tensor (3x3)
c                       gradFe (3x3x3) averaged over all element
c                       integration points then 3x3 Nye computed
c                  b. R(3,3) from F = R U
c                  c. work density (1x1)
c                  d. plastic work density (1x1
c                  e. equiv eps-pls (1x1
c                  f. slip history: sum over time of (signed) plastic
c                     slip on each slip system. Makes sense for material
c                     w/ 1 crystal - slip on on a system is summed for
c                     each crystal then averaged w/o accounting for
c                     crystal orientations (max_slip_sys x1)
c        **   values are averages computed over all integration points
c             note: one_elem_states zeroed before entry.
c
      num_states_here = 21 + length_comm_hist(5)
c
      if( num_states_here .ne. num_states ) then
         write(iout,9000) 1
         write(iout,*) ' ... num_states: ', num_states
         call die_gracefully
      end if
c
c             1a. Nye tensor, averaged over int. pts.
c                 gradFe = reshape(
c                     sum(hisblk(s:e,1:int_points,relem),2 ) /
c                     dble(int_points), (/3,3,3/) )
c
      s  = 1
      e  = s + 8
      sh = indexes_common(2,1)
      eh = indexes_common(2,2)
      gradvec = zero
      call mm10_avg_states( gradvec, 1, 27, sh, eh, hisblk(1,1,relem),
     &                      hist_size, int_points )
      call mm10_states_nye( gradvec, nyevec )
      one_elem_states(s:e) = nyevec(1:9)
      s = s + 9
c
c             1b. R of F = RU averaged over element int points
c
      e = s + 8
      sh = indexes_common(3,1)
      eh = indexes_common(3,2)
      call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
      s = s + 9
c
c             1c,d,e. work density (1x1), plastic work density (1x1
c                     equiv eps-pls (1x1)
c
      e = s + 2
      sh = indexes_common(4,1)
      eh = indexes_common(4,2)
      call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
      s = s + 3
c
c             1f. current (total) slip on each system. averaged over
c                 integration points
c
      sh  = indexes_common(5,1)
      eh  = indexes_common(5,2)
      e   = s + (length_comm_hist(5) - 1)
      call mm10_avg_states( one_elem_states, s, e, sh, eh,
     &                  hisblk(1,1,relem), hist_size, int_points )
      s = s + length_comm_hist(5)
c
      if( s .ne. (21 + length_comm_hist(5) + 1) ) then
        write(iout,9000) 2
        call die_gracefully
      end if
      write(iout,*) '... inside mm10_states_values_6 @ 1 '
c
      return
c
 9000 format('>> FATAL ERROR: mm10_states_values_5 @: ',i2,
     &   /,  '                job termindated.'// )
c
      end subroutine mm10_states_values_5
      end subroutine mm10_states_values
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_states_nye                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 6/5/2016 rhd                *
c     *                                                              *
c     *    service routine to get nye tensor for states output       *                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_nye( gradFe, nye )
      implicit none
c
      double precision :: gradFe(3,3,3), nye(3,3)
c
      double precision, parameter :: zero = 0.0d00
      integer :: d, z, b, w, kk
c
@!DIR$ ASSUME_ALIGNED gradFe:64, nye:64
c
      nye = zero
c
      do d = 1,3
          do z = 1,3
            do b = 1,3
@!DIR$ IVDEP
                do w = 1,3
                  kk = (z-b)*(b-w)*(w-z)/2
                  nye(d,z) = nye(d,z) - dble(kk)*gradFe(d,b,w)
                end do
             end do
          end do
      end do
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_avg_states                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 6/2/2016 rhd                *
c     *                                                              *
c     *    service routine to help putting values in the states      *
c     *    output vector                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_avg_states( states_vec, s, e, sh, eh, hisblk,
     &                            hist_size, int_points )
      implicit none

      integer :: s, e, sh, eh,  hist_size, int_points
#dbl      double precision ::
#sgl      real  ::
     & states_vec(*), hisblk(hist_size,int_points)

      integer :: igp, len1, len2
#dbl      double precision :: one, dint_pts
#sgl      real  :: one, dint_pts
c
      data one / 1.0d00 /
@!DIR$ ASSUME_ALIGNED states_vec:64, hisblk:64
c
      len1 = e - s + 1
      len2 = eh - sh + 1
      if( len1 .ne. len2 ) then  ! sanity check
        write(*,9100) len1, len2
        call die_abort
      endif
c
      dint_pts = one / dble( int_points )
c
@!DIR$ IVDEP
      do igp = 1, int_points
         states_vec(s:e) = states_vec(s:e) + hisblk(sh:eh,igp)
      end do
      states_vec(s:e) = states_vec(s:e) * dint_pts
c
      return
 9100 format('>>> FATAL ERROR: routine. mm10_avg_states.',
     & /,    '                 len1, len2: ', 2i6)
      end
