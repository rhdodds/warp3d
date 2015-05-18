c     ****************************************************************
c     *                                                              *
c     *                      subroutine stifup                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 6/25/12                    *
c     *                                                              *
c     *     decide and invoke the process to compute element         *
c     *     element stiffness matrices.                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine stifup( step, iter, force_lin_k_iter_1 )
      implicit integer (a-z)
$add common.main
c
      logical force_lin_k_iter_1
      newstf = .true.
c
c               stepfor step 1, iter 0, 1 the linear [k]s have already
c               been computed by system initialization.
c
      if ( step .eq. 1 .and. iter .le. 1 ) return
c
c               for iteration 0 of a load step we have several cases to
c               consider:
c                a) if the user wants the linear [k]
c                   just do it and leave.
c                b) the force_k_update = 2 case means the adaptive
c                   algorithm forced a restart of the step with
c                   increased number of subincrements. we can force
c                   a tangent or linear stiffness update based on the solution
c                   at start of step since all state variables
c                   are now saved and will be reset.
c                   Experimenting with best updating algorithm to handle
c                   first adaptive increment following non-convergence.
c                c) force_k_update = 1 means we just did a restart
c                d) if not a or b, we just get the tangent stiffness
c
      local_step = step   !  protects system level values
      local_iter = iter
c
      if( local_iter .eq. 1 ) then
         write(out,9200)
         call die_abort
      end if
c
      if( local_iter .gt. 1 ) then
        if ( show_details ) write(out,9110) local_step, local_iter
        call tanstf( .false., local_step, local_iter )
        force_k_update = 0
        return
      end if
c
c               cases for iteration 0 of a load step
c
      if( lnkit1 .or. force_lin_k_iter_1 ) then
         if( show_details ) write(out,9100) local_step
         call lnstff( local_step, local_iter )
         force_k_update = 0
         return
      end if
c
      if( force_k_update .eq. 2 ) then
         if( show_details ) write(out,9120) local_step
         call tanstf( .false., local_step, local_iter )
         force_k_update = 0
         return
      end if
c
      if( growth_k_flag ) then
         if( show_details ) write(out,9100) local_step
         call lnstff( local_step, local_iter )
         growth_k_flag = .false.
         return
      end if
c
      if ( show_details ) write(out,9120) local_step
      call tanstf( .false., local_step, local_iter )
      force_k_update = 0
      return
c
 9100 format(7x,
     & '>> computing linear stiffness start of step        ',i5)
 9110 format(7x,
     & '>> computing tangent stiffness step, iteration:    ',i5,i3)
 9120 format(7x,
     & '>> computing tangent stiffness start of step       ',i5)
 9200 format(
     & '>>> FATAL ERROR: stifup. Please contact WARP3D developers...',/,
     & '                 Job terminated.',/)
c
      end






