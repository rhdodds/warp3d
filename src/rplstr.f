c     ****************************************************************
c     *                                                              *
c     *                      subroutine rplstr                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/18/04 rhd               *
c     *                                                              *
c     *     this subroutine stores globally the recovered material   *
c     *     and stress states.                                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rplstr( span, felem, ngp, mat_type, iter, geonl,
     &                   local_work, blk )  
      use elem_block_data, only : history1_blocks, rot_n1_blocks,
     &                            rts_blocks, eps_n1_blocks,
     &                            urcs_n1_blocks, history_blk_list

      use segmental_curves, only : max_seg_points, max_seg_curves 
      implicit integer (a-z)
$add common.main
$add include_sig_up
      logical process_rts, process_hist

c
c                   parameter declarations
c
c      
      logical geonl
c
c                                         
c             replace stress/strain data:
c
c               1) unrotated cauchy stresses at n+1 (all models)
c
c               2) deviatoric components of trial elastic stress
c                  increment over step (if model uses them)
c
c               3) ddtse comes back from element strain routine with
c                  total strains at n+1 (all models)
c 
c               4) [R,n+1] for geonl models
c
c               5) history data updated to n+1
c
c             iter = 0 indicates the "pre" step computations used
c             by the load step processor to get estimated stresses
c             for non-zero imposed displacements and/or the
c             extrapolated displacement increment for step.
c             we do no want to update the material state variables
c             especially (2) above since they most likely will be
c             used immediately after this for stiffness
c             computation.
c
c
      process_rts  = mat_type .ne. 2  ! note the .ne. 
c
      process_hist = mat_type .eq. 1  .or.  mat_type .eq. 2 .or.
     &               mat_type .eq. 3  .or.  mat_type .eq. 4 .or.
     &               mat_type .eq. 5  .or.  mat_type .eq. 6 .or.
     &               mat_type .eq. 7  .or.  mat_type .eq. 8 .or.
     &               mat_type .eq. 10
      
c
      call scstr( local_work%urcs_blk_n1, urcs_n1_blocks(blk)%ptr(1),
     &            ngp, nstrs,  span)
c
      call scstr( local_work%ddtse, eps_n1_blocks(blk)%ptr(1),
     &            ngp, nstr, span )
c
      if ( iter .gt. 0 ) then
         if ( process_rts ) then
           call scstr( local_work%rtse, rts_blocks(blk)%ptr(1), ngp,
     &                 nstr, span )
         end if
         if ( geonl ) then
           call scstr( local_work%rot_blk_n1, rot_n1_blocks(blk)%ptr(1),
     &                 ngp, 9, span )
         end if
      else
         if ( geonl ) then
           call gastr( local_work%rot_blk_n1, rot_n1_blocks(blk)%ptr(1),
     &                 ngp, 9, span )
          end if
      end if
c
      if ( process_hist .and. iter .gt. 0 ) then
        hist_size = history_blk_list(blk)
        call scstr_history( local_work%elem_hist1(1,1,1),
     &                      history1_blocks(blk)%ptr(1), ngp, 
     &                      hist_size, span )
      end if
c
      return
      end
