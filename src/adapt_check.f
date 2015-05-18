c     ****************************************************************
c     *                                                              *
c     *                      subroutine adapt_check                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/08/94                   *
c     *                   last modified : 02/15/94                   *
c     *                                   09/27/94 kck               *
c     *                                   04/07/95 kck               *
c     *                                                              *
c     *     executes all logic to drive the adaptive solution        *
c     *     process for the static analysis of a load step           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine adapt_check( scaling_adapt, process_type, step )
      use adaptive_steps
      implicit integer (a-z)
c
      logical local_debug
#dbl      double precision
#sgl      real
     & min_fact, start_scale, end_scale, incr_scale, start, zero, one,
     & scaling_adapt
c     
      data min_fact, zero, one / 0.0625, 0.0, 1.0 /
      data local_debug / .false. /
c
      call iodevn( idummy, iout, dummy, 1 )
c
      go to ( 100, 200, 300 ), process_type
c
c                         process type - 1
c                         ----------------
c
c            iterations for an increment failed to converge or
c            material models requested immediate step size cut.
c            subdivide increment and start over for increment.
c            we may be at the limit of subdivision now.
c
 100  continue
      adapt_result = 1
      num_subs     = adapt_divisions
      if ( adapt_disp_fact .le. min_fact ) then
           write(iout,*) ' '
           write(iout,9420) step, adaptive_stack(1,adapt_level),
     &                      adaptive_stack(2,adapt_level)
           write(iout,*) ' '
           adapt_result = 3
           go to 1000
      end if
c
      start_scale  = adaptive_stack(1,adapt_level)
      end_scale    = adaptive_stack(2,adapt_level)
#sgl      incr_scale   = (end_scale - start_scale) / real(num_subs)
#dbl      incr_scale   = (end_scale - start_scale) / dble(num_subs)
c
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
      if ( local_debug ) write (iout,9500)
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
      go to 1000
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
 200  continue
      adapt_result = 1
      num_subs     = adapt_divisions
      if ( adapt_level .le. 1  ) go to 1000
      end_scale = adaptive_stack(2,adapt_level)
      if ( end_scale .eq. adaptive_stack(2,adapt_level-1) ) then
          adapt_level  = adapt_level - 2
          if ( adapt_level .eq. 0 ) then
             scaling_adapt = adaptive_stack(4,1)
             go to 1000
          end if
          adaptive_stack(4,adapt_level) = num_subs *
     &               adaptive_stack(4,adapt_level)
          if (local_debug )
     &        write (iout,9510) adapt_level, num_subs,
     &        adaptive_stack(4,adapt_level), num_subs
          scaling_adapt = adaptive_stack(4,adapt_level)
      else
          adapt_level  = adapt_level - 1
      end if
      if ( adapt_level .le. 1 ) go to 1000
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
      go to 1000
c
c                         process type - 3
c                         ----------------
c
c            initialize adaptive stack for step
c
 300  continue
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
      scaling_adapt       = one
      return
c
c
 1000 continue
      if ( .not. local_debug ) return
      write(iout,*) '>> returning from adapt_check:'
      write(iout,*) '     process type: ', process_type
      write(iout,*) '     adapt_result: ', adapt_result
      write(iout,9000) adapt_level, adapt_disp_fact, adapt_load_fact
      write(iout,9100) (i, (adaptive_stack(j,i),j=1,4),i=1,20)
c
      return
c
 9000 format(' new level, disp_fact, load_fact: ',i3,3f10.4)
 9100 format(1x,i4,4f12.4)
 9200 format(/,4x,'adaptive solution for load step: ',i5,
     & /,    7x,'now advancing step solution from: ',f7.4,' to: ',f7.4)
 9300 format(/,6x,'adaptive solution driver for step: ',i5,
     & /, 9x,'subdividing current increment from: ',f7.4,' to: ',f7.4,
     & /, 9x,'into smaller increments of:         ',f7.4,
     & /, 9x,'now advancing step solution from:   ',f7.4,' to: ',f7.4)
 9420 format(/,6x,'adaptive solution for load step: ',i5,
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


c
      end


