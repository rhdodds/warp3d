c     ****************************************************************
c     *                                                              *
c     *                      subroutine rplstr                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/17/2017 rhd             *
c     *                                                              *
c     *     stores globally the recovered material                   *
c     *     and stress states.                                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rplstr( span, felem, ngp, mat_type, iter, geonl,
     &                   local_work, blk, block_stress_sums,
     &                   block_strain_sums, local_debug_sums )
      use global_data ! old common.main
      use elem_block_data, only : history1_blocks, rot_n1_blocks,
     &                            eps_n1_blocks,
     &                            urcs_n1_blocks, history_blk_list

      use segmental_curves, only : max_seg_points, max_seg_curves
      implicit none
c
      include 'include_sig_up'
c
c                   parameter declarations
c
      integer :: span, felem, ngp, iter, blk
      logical :: geonl, local_debug_sums
      double precision :: block_stress_sums(1), block_strain_sums(1)
c
c                   local declarations
c
      integer :: mat_type, hist_size
      logical :: process_hist, save_history_1, save_history_2
c
c
c             replace stress/strain data:
c
c               1) unrotated cauchy stresses at n+1 (all models)
c
c               2) ddtse comes back from element strain routine with
c                  total strains at n+1 (all models)
c
c               3) [R,n+1] for geonl models
c
c               4) history data updated to n+1
c
c             iter = 0 indicates the "pre" step computations used
c             by the load step processor to get estimated stresses
c             for non-zero imposed displacements and/or the
c             extrapolated displacement increment for step.
c             we do no want to update the material state variables
c             especially (2) above since they most likely will be
c             used immediately after this for stiffness
c             computation. Except for the CP model since it hides the
c             [D] matrix in there.
c
c
      process_hist = mat_type .eq. 1  .or.  mat_type .eq. 2 .or.
     &               mat_type .eq. 3  .or.  mat_type .eq. 4 .or.
     &               mat_type .eq. 5  .or.  mat_type .eq. 6 .or.
     &               mat_type .eq. 7  .or.  mat_type .eq. 8 .or.
     &               mat_type .eq. 10
c
      call rp_scstr( local_work%urcs_blk_n1,
     &            urcs_n1_blocks(blk)%ptr(1), ngp, nstrs,  span)
      if( local_debug_sums )
     &   block_stress_sums(blk) =
     &          norm2( urcs_n1_blocks(blk)%ptr(1:ngp*nstrs*span) )
c
      call rp_scstr( local_work%ddtse, eps_n1_blocks(blk)%ptr(1),
     &               ngp, nstr, span )
      if( local_debug_sums )
     &   block_strain_sums(blk) =
     &          norm2( eps_n1_blocks(blk)%ptr(1:ngp*nstr*span) )
c
      if( iter .gt. 0 .and. geonl )
     &      call rp_scstr( local_work%rot_blk_n1,
     &               rot_n1_blocks(blk)%ptr(1), ngp, 9, span )
c
      save_history_1 = process_hist .and. iter > 0
      save_history_2 = mat_type .eq. 10 ! [D] is hidden in history ...
      if( save_history_1  .or.  save_history_2 ) then
        hist_size = history_blk_list(blk)
        call rp_scstr_history( local_work%elem_hist1(1,1,1),
     &                      history1_blocks(blk)%ptr(1), ngp,
     &                      hist_size, span )
      end if
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rp_gastr                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/30/2016 rhd              *
c     *                                                              *
c     *     gathers element stresses from the global                 *
c     *     stress data structure to a block of similar,             *
c     *     elements for all gauss points.                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rp_gastr( ml, mg, ngp, nprm, span )
      implicit none
      include 'param_def'
c
c               parameter declarations
c
      integer :: ngp, nprm, span
      double precision :: ml(mxvl,nprm,*), mg(nprm,ngp,*)
c
c               local declarations
c
      integer :: i, j, k
c
c               on current hardware, the manual unroll runs
c               slower
c
c      if( ngp .ne. 8 ) then
!DIR$ LOOP COUNT MAX=27
        do k = 1, ngp
         do  j = 1, nprm
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
            do i = 1, span
               ml(i,j,k) = mg(j,k,i)
            end do
         end do
        end do
        return
c      end if
c
c                number of gauss points = 8, unroll.
c
c      do j = 1, nprm
c@!DIR$ IVDEP
c        do i = 1, span
c            ml(i,j,1) = mg(j,1,i)
c            ml(i,j,2) = mg(j,2,i)
c            ml(i,j,3) = mg(j,3,i)
c            ml(i,j,4) = mg(j,4,i)
c            ml(i,j,5) = mg(j,5,i)
c            ml(i,j,6) = mg(j,6,i)
c            ml(i,j,7) = mg(j,7,i)
c            ml(i,j,8) = mg(j,8,i)
c        end do
c      end do
c
c      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine rp_scstr                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/30/2016 rhd              *
c     *                                                              *
c     *     scatters element stresses to the global                  *
c     *     stress data structure from a block of similar            *
c     *     elements for all gauss points.                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rp_scstr( ml, mg, ngp, nprm, span )
      implicit none
      include 'param_def'
c
c               parameter declarations
c
      integer :: ngp, nprm, span
      double precision ::
     &     ml(mxvl,nprm,*),mg(nprm,ngp,*)
c
c               local declarations
c
      integer :: i, j, k
c
c               on current hardware, the manual unroll runs
c               slower
c
c
c      if( ngp .ne. 8 ) then
!DIR$ LOOP COUNT MAX=27
        do k = 1, ngp
           do j = 1, nprm
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              do i = 1, span
                 mg(j,k,i) = ml(i,j,k)
              end do
           end do
        end do
        return
c      end if
c
c                       number of gauss points = 8
c
c      do j = 1, nprm
c@!DIR$ IVDEP
c        do i = 1, span
c          mg(j,1,i) = ml(i,j,1)
c          mg(j,2,i) = ml(i,j,2)
c          mg(j,3,i) = ml(i,j,3)
c          mg(j,4,i) = ml(i,j,4)
c          mg(j,5,i) = ml(i,j,5)
c          mg(j,6,i) = ml(i,j,6)
c          mg(j,7,i) = ml(i,j,7)
c          mg(j,8,i) = ml(i,j,8)
c        end do
c      end do
c
c      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rp_scstr_history             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/28/2016 rhd             *
c     *                                                              *
c     *     scatters element history to the global                   *
c     *     history data structure from a block of similar           *
c     *     elements for all gauss points.                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rp_scstr_history( local_history, global_history,
     &                          ngp, hist_size, span )
      implicit none
c
c               parameter declarations
c
      integer :: ngp, hist_size, span
      double precision
     &      local_history(span,hist_size,ngp),
     &      global_history(hist_size,ngp,span)
c
c               local declarations
c
      integer :: i, j, k
c
c               on current hardware, the manual unroll runs
c               slower
c
c      if( ngp .ne. 8 ) then
!DIR$ LOOP COUNT MAX=27
        do k = 1, ngp
           do j = 1, hist_size
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
              do i = 1, span
                 global_history(j,k,i) = local_history(i,j,k)
              end do
           end do
        end do
        return
c      end if
c
c                       number of gauss points = 8
c
c      do j = 1, hist_size
c@!DIR$ IVDEP
c        do i = 1, span
c          global_history(j,1,i) = local_history(i,j,1)
c          global_history(j,2,i) = local_history(i,j,2)
c          global_history(j,3,i) = local_history(i,j,3)
c          global_history(j,4,i) = local_history(i,j,4)
c          global_history(j,5,i) = local_history(i,j,5)
c          global_history(j,6,i) = local_history(i,j,6)
c          global_history(j,7,i) = local_history(i,j,7)
c          global_history(j,8,i) = local_history(i,j,8)
c        end do
c      end do
c
c      return
      end
c
