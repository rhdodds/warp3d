c     ****************************************************************
c     *                                                              *
c     *                      subroutine adapt_check                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/11/2018 RHD              *
c     *                                                              *
c     *     executes all logic to drive the adaptive solution        *
c     *     process for the static analysis of a load step           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine adapt_check( scaling_adapt, process_type, step,
     &                        iout )
      use adaptive_steps
      implicit none
c
c              parameters
c
      integer :: process_type, step, iout
      double precision :: scaling_adapt
c
c              locals
c
      integer :: num_subs, i, j
      logical :: local_debug
      double precision ::
     & start_scale, end_scale, incr_scale, start, zero, one
c
      data zero, one / 0.0d00, 1.0d00 /
      data local_debug / .false. /
c
      select case( process_type )
c
c                         process type - 1
c                         ----------------
c
c            iterations for an increment failed to converge or
c            material models requested immediate step size cut.
c            subdivide increment and start over for increment.
c            we may be at the limit of subdivision now.
c
      case( 1 )
      adapt_result = 1
      num_subs     = adapt_divisions
      if( adapt_disp_fact .le. adapt_min_fact ) then
           write(iout,*) ' '
           write(iout,9420) step, adaptive_stack(1,adapt_level),
     &                      adaptive_stack(2,adapt_level)
           write(iout,*) ' '
           adapt_result = 3
           call adapt_check_a
           return
      end if
c
      start_scale  = adaptive_stack(1,adapt_level)
      end_scale    = adaptive_stack(2,adapt_level)
      incr_scale   = (end_scale - start_scale) / dble(num_subs)
c
      if( adapt_level+num_subs .ge. adapt_cols ) then
        write(iout,9010) adapt_level+num_subs
        call die_abort
      end if
      do i = adapt_level+1, adapt_level+num_subs, 1
        adaptive_stack(3,i) = incr_scale
      end do
c
      start = start_scale
      do i = adapt_level+num_subs, adapt_level+1, -1
        adaptive_stack(1,i) = start
        adaptive_stack(2,i) = start + incr_scale
        start               = start + incr_scale
      end do
c
      if( local_debug ) write (iout,9500)
     &    adapt_level+num_subs, adaptive_stack(4,adapt_level+num_subs),
     &    adapt_level, adaptive_stack(4,adapt_level),
     &    adaptive_stack(4,adapt_level+num_subs)
c
c        the adaptive_stack contains a scaling factor for extrapolating
c        solution parameters from one subincremented step to the next.
c        this section of code calculates the beginning and ending scaling
c        factors for this adaptive cut and then fill in the remaining
c        values inbetween with 1.0
c
      adaptive_stack(4,adapt_level+num_subs) =
     &          adaptive_stack(4,adapt_level) / num_subs
      adaptive_stack(4,adapt_level) = adaptive_stack(4,adapt_level) /
     &          adaptive_stack(4,adapt_level+num_subs)
      adaptive_stack(4,adapt_level+1:adapt_level+num_subs-1) = one
c
c
c
      adapt_level       = adapt_level + num_subs
      scaling_adapt     = adaptive_stack(4,adapt_level)
      adapt_disp_fact   = adaptive_stack(3,adapt_level)
      adapt_temper_fact = adaptive_stack(3,adapt_level)
      adapt_load_fact   = adaptive_stack(2,adapt_level)
      adapt_result      = 2
      write(iout,*) ' '
      write(iout,9300) step, start_scale, end_scale,
     &                 incr_scale, adaptive_stack(1,adapt_level),
     &                 adaptive_stack(2,adapt_level)
      write(iout,*) ' '
      call adapt_check_a
c
c                         process type - 2
c                         ----------------
c
c            iterations for an increment converged.
c            start processing previous increment in stack
c            unless we are done with step.
c
c            when the end_scale = the previous adaptive stack value
c            then we are at the end of a subincrementation and we can
c            progress to the next level of subincrementation.  At this
c            point, we need to check if this first energy in the stack.
c            if it is, then we can pass the current value on the stack
c            to the scaling_adapt.  If not, we need to multiply the
c            stack value by the number of subincrements.
c
c
      case( 2 )
      adapt_result = 1
      num_subs     = adapt_divisions
      if( adapt_level .le. 1  ) then
          call adapt_check_a
          return
      end if
      end_scale = adaptive_stack(2,adapt_level)
      if( end_scale .eq. adaptive_stack(2,adapt_level-1) ) then
          adapt_level  = adapt_level - 2
          if( adapt_level .eq. 0 ) then
             scaling_adapt = adaptive_stack(4,1)
             call adapt_check_a
             return
          end if
          adaptive_stack(4,adapt_level) = num_subs *
     &               adaptive_stack(4,adapt_level)
          if( local_debug )
     &        write (iout,9510) adapt_level, num_subs,
     &        adaptive_stack(4,adapt_level), num_subs
          scaling_adapt = adaptive_stack(4,adapt_level)
      else
          adapt_level  = adapt_level - 1
      end if
      if( adapt_level .le. 1 ) then
         call adapt_check_a
         return
      end if
c
c                  more work to do.
c
      scaling_adapt     = adaptive_stack(4,adapt_level)
      adapt_disp_fact   = adaptive_stack(3,adapt_level)
      adapt_temper_fact = adaptive_stack(3,adapt_level)
      adapt_load_fact   = adaptive_stack(2,adapt_level)
      adapt_result      = 2
      write(iout,9200) step, adaptive_stack(1,adapt_level),
     &                 adaptive_stack(2,adapt_level)
      write(iout,*) ' '
      call adapt_check_a
c
c                         process type - 3
c                         ----------------
c
c            initialize adaptive stack for step
c
      case( 3 )
      adaptive_stack(1:adapt_rows,1:adapt_cols) = zero
      adaptive_stack(1,1) = zero
      adaptive_stack(2,1) = one
      adaptive_stack(3,1) = one
      adaptive_stack(4,1) = one
      adapt_level         = 1
      adapt_disp_fact     = one
      adapt_temper_fact   = one
      adapt_load_fact     = one
      adapt_divisions     = 4
      adapt_min_fact      = one / dble(adapt_divisions**2)
      scaling_adapt       = one
c
      case default
        write(iout,9000)
        call die_abort
      end select
c
      return
c
 9000 format( 1x,'>>>>> FATAL Error: adapt_check. bad process_type',
     &   /,   1x,'                   job aborted.',/)
 9010 format( 1x,'>>>>> FATAL Error: adapt_check. table overflow: ',i5,
     &   /,   1x,'                   job aborted.',/)
 9200 format(/,4x,'adaptive solution for load step: ',i7,
     & /,    7x,'now advancing step solution from: ',f7.4,' to: ',f7.4)
 9300 format(/,6x,'adaptive solution driver for step: ',i7,
     & /, 9x,'subdividing current increment from: ',f7.4,' to: ',f7.4,
     & /, 9x,'into smaller increments of:         ',f7.4,
     & /, 9x,'now advancing step solution from:   ',f7.4,' to: ',f7.4)
 9420 format(/,6x,'adaptive solution for load step: ',i7,
     & /, 9x,'failed in advancing increment from: ',f7.4,' to: ',f7.4,
     & /, 9x,'subdivision limit has been reached.',
     & /, 9x,'the analysis will be terminated; restart from last ',
     & /, 9x,'converged database file.')
 9500 format ( ' Divide adaptive scaling factor (', i3, ') = ',
     &           e16.6, ' by ', i3  /
     & ' Divide adaptive scaling factor (', i3, ') = ',
     &           e16.6, ' by ', e16.6  )
 9510 format ( 3x, ' Multiply the adaptive scale factor in row ',
     &               i3 / ' by to obtain the following ',
     &             'adaptive_stack = ',e16.6 )

      contains
c     ========

      subroutine adapt_check_a
c
      if ( .not. local_debug ) return
      write(iout,*) '>> returning from adapt_check:'
      write(iout,*) '     process type: ', process_type
      write(iout,*) '     adapt_result: ', adapt_result
      write(iout,9000) adapt_level, adapt_disp_fact, adapt_load_fact
      write(iout,9100) (i, (adaptive_stack(j,i),j=1,4),i=1,20)
      return
c
 9000 format(' new level, disp_fact, load_fact: ',i3,3f10.4)
 9100 format(1x,i4,4f12.4)
c
      end subroutine adapt_check_a
      end subroutine adapt_check

