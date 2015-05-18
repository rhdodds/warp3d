c     ****************************************************************
c     *                                                              *
c     *                      subroutine ebeslv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 04/13/97 asg               *
c     *                                                              *
c     *     this subroutine creates the pseudo-residual by           *
c     *     multiplying the inverse of the preconditioner and the    *
c     *     residual.  The pseudo-residual allows faster convergence.*
c     *                                                              *
c     *     this is the parallel version, where the scheduling of    *
c     *     the ebe precoonditioner is determined by a balanced      *
c     *     coloring algorithm for maximum parallel execution.       *
c     *                                                              *
c     *     see routine find_ebe_order in mpi_lnpcg.f for more       *
c     *     details.                                                 *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ebeslv(r,z,local_diag,refdof,mydof, iter )
      use elem_block_data, only: pcm_blocks, edest_blocks
      use pcg_data 
      use mpi_lnpcg
      implicit integer (a-z)
c
$add common.main
c
#sgl      real
#dbl      double precision
     &     dummy(1), r(*), z(*), local_diag(*), zero, ddum 
      logical debug, ldum
      data debug /.false./
      integer, dimension(:,:), pointer :: edest
#sgl      real,
#dbl      double precision,
     &     dimension(:,:), pointer :: pcm
      logical internal_done, send_done, recv_done
      dimension node_list(mxnmbl)
      data zero /0.0/
c
c
      if ( debug ) write(*,*) '>>>>> inside ebeslv', refdof, mydof
      if ( .not.nostif .and. ebe_pre ) then
         if ( debug ) write (*,*) 'preconditioner is ebe.' 
c
c              perform initial global diagonal scaling.
c
         call vec_ops( r, local_diag, z,  refdof, 1 )
c
c              perform element forward reduction.
c
c	          the ordering of operations is pre-determined by
c	          a balanced coloring algorithm to ensure good
c                 levels of mpi-level parallelism
c
         if ( debug ) write(*,*) 'forward reduction'
c
         if ( debug ) write(*,*) '    fixed ordering'
         int_ptr = 0
         num_sends = 0
         num_list = 0
         do ptr = 1, order_ptr
            oper = order_list (1,ptr)
c
c                          if oper = 1, apply preconditioner to specified
c                          internal blocks
c
            if ( oper .eq. 1 ) then
c
               do i = 1, order_list(2,ptr)
                  blk = local_nodes%internal_blks(int_ptr + i)
                  call efrdrv(r,z,local_diag,refdof,mydof, blk)
               end do
               int_ptr = int_ptr + order_list(2,ptr)
c
c                          if oper = 2, apply preconditioner to specified
c                          group of external blocks
c
            else if ( oper .eq. 2 ) then
c
               grnode = order_list (2,ptr)
c
c                             loop over blks in group and apply the
c                             preconditioner
c
               do i = 1, local_graph(grnode)%num_blks
                  blk = local_graph(grnode)%blks(i)
                  call efrdrv(r,z,local_diag,refdof,mydof, blk)
               end do
               num_list = num_list + 1
               node_list(num_list) = grnode
c
c                          if oper = 3, we are are completing a group.
c                          Make sure all previous sends have completed,
c                          and that all receives needed for this 
c                          group are done.
c
            else if ( oper .eq. 3 ) then
               do i= 1, num_list
                 call dd_send ( z, ldum, num_sends, node_list(i), 1)
               end do
               num_list = 0
c
c                          if oper = 4, we start a new group.  Register 
c                          all required receives for this group 
c
            else if ( oper .eq. 4 ) then
c
               grp = order_list(2,ptr)
               if ( grp .gt. 1 ) then
c
                  if ( num_sends .gt. 0 ) then
                     call dd_send ( ddum, ldum, num_sends, idum, 3 )
                     call dd_send ( ddum, ldum, num_sends, idum, 4 )
                  endif
c
                  call dd_recv ( idum, ddum, ldum, 3 )
                  call dd_recv ( grp-1, z, ldum, 4 )
                  num_sends = 0 
c
               endif
c
               if ( grp .lt. num_groups + 1) then
                  call dd_recv ( grp, ddum, ldum, 1 )
               endif
c
            endif
c
         enddo
c
c              perform element diagonal scaling. 
c
 3000    continue
c
         if ( debug ) write (*,*) ' diagonal scaling'
         call vec_ops( pdiag, z, dummy, refdof, 2 )
c
c              perform element back substitution.
c
c                 to enable faster execution, ebksub is fully
c                 unrolled in seperate versions for each element
c                 type.
c
c	          the ordering of operations is pre-determined by
c	          a balanced coloring algorithm to ensure good
c                 levels of parallelism.  Do the process
c                 backwards to ensure symmetry of the preconditioner
c
c
         if ( debug ) write (*,*) ' back substitution.'
c
         num_sends = 0
         num_list = 0
         int_ptr = local_nodes%num_int_blks + 1
         do ptr = order_ptr, 1, -1
            oper = order_list (1,ptr)
c
c                          if oper = 1, apply preconditioner to specified
c                          internal blocks
c
            if ( oper .eq. 1 ) then
               do i = 1, order_list(2,ptr)
                  blk = local_nodes%internal_blks(int_ptr - i)
                  call ebkdrv(r,z,local_diag,refdof,mydof, blk)
               end do
               int_ptr = int_ptr - order_list(2,ptr)
c
c                          if oper = 2, apply preconditioner to specified
c                          group of external blocks
c
            else if ( oper .eq. 2 ) then
               grnode = order_list (2,ptr)
c
c                             loop over blks in group and apply the
c                             preconditioner
c
               do i = local_graph(grnode)%num_blks, 1, -1
                  blk = local_graph(grnode)%blks(i)
                  call ebkdrv(r,z,local_diag,refdof,mydof, blk)
               end do
               num_list = num_list + 1
               node_list(num_list) = grnode
c
c                          if oper = 3, we start a new group.  Register 
c                          all required receives for this group 
c
            else if ( oper .eq. 3 ) then
               grp = order_list(2,ptr)
c
               if ( grp .lt. num_groups + 1 ) then
c
                  if ( num_sends .ne. 0) then
                     call dd_send ( ddum, ldum, num_sends, idum, 3 )
                     call dd_send ( ddum, ldum, num_sends, idum, 4 )
                  endif
c     
                  call dd_recv ( idum, ddum, ldum, 3 )
                  call dd_recv ( grp, z, ldum, 4 )
                  num_sends = 0
               endif
c
               if ( grp .gt. 1) then
                  call dd_recv ( grp - 1, ddum, ldum, 1 )
               endif
c
c                          if oper = 4, we are are completing a group.
c                          Make sure all previous sends have completed,
c                          and that all receives needed for this 
c                          group are done. Also set up for the next
c                          set of recvs
c
            else if ( oper .eq. 4 ) then
               do  i = 1, num_list
                  call dd_send ( z, ldum, num_sends, node_list(i), 1)
               end do
               num_list = 0
c
            endif
c
         end do
c
c                       perform final global diagonal scaling.
c
         if (debug) write (*,*) ' final diagonal scaling'
         call vec_ops( local_diag, z, dummy, refdof, 2 )
c
c
      else
c
         if (debug) write (*,*) 'preconditioner is diagonal.'
c
c                       the preconditioning matrix is diagonal.
c
         call vec_ops( r, local_diag, z, refdof, 3 )
c
c
      end if
c
c
 9999 continue
      if ( debug ) write(*,*) '<<<<< leaving ebeslv'
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine efrdrv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/27/09 rhd               *
c     *                                   comments about threads     *
c     *                                                              *
c     *     this subroutine conducts the forward substituion stage   *
c     *     of the ebe preconditioner for a block of elements        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine efrdrv(r,z,local_diag,refdof,mydof, blk)
      use elem_block_data, only: pcm_blocks, edest_blocks
      use pcg_data 
      use mpi_lnpcg
      implicit integer (a-z)
c
$add common.main
c
#sgl      real
#dbl      double precision
     & dummy(1), r(*), z(*), local_diag(*), zero
      logical debug, valid_type
      data debug /.false./
      integer, dimension(:,:), pointer :: edest
#sgl      real,
#dbl      double precision,
     &     dimension(:,:), pointer :: pcm
      data zero /0.0/
c
c
      if ( debug ) write(*,*) '>>>>> inside ebeslv', refdof, mydof
c
c	 At present, this process is not executed in thread
c        parallel. The blocking within an MPI doamin is not
c        "vectorized."  Vectorized blocking is required for the
c        threads only version of WARP3D. Then the span loop in
c        the ebksub routines can be executed in thread parallel.
c
c        It may be possible to relax serial execution here for the
c        span loop in ebksub routine. Serial here or vectorized
c        blocking is "strictly" required to maintain complete
c        symmetry of the forward/backward application of the
c        ebe preconsditioner to the current residual.
c
c        Some numerical testing may show this can be relaxed.                   

      span   = elblks(0,blk)
      felem  = elblks(1,blk)
      totdof = iprops(2,felem)*iprops(4,felem)
      utsz = (totdof * ( totdof + 1 )) / 2
      type  =  iprops(1,felem)
      edest  => ledest_blocks(blk)%ptr
      pcm    => pcm_blocks(blk)%ptr
      valid_type = ( type .ge. 1 .and. type .le. 5 ) 
     &             .or. type .eq. 12
      if ( .not. valid_type ) then
        write(*,9000) type, felem
        call die_abort
      end if
c
c	 the 8-node elements use unrolled code. others use general
c        code.  Based on benchmarks in late 2009.
c
      if ( type .eq. 2 .or. type .eq. 12 ) then
        call efrwrd8( span, totdof, edest(1,1), z, pcm(1,1), utsz )
      else
        call efrwrd_gen( span, totdof, edest(1,1), z, pcm(1,1), utsz )
      endif
c
      return 
c
 9000 format('>>> FATAL SYSTEM ERROR: warp3d does not have the',
     &     /,'                        code implemented to process',
     &     /,'                        an element type. routine efrdrv',
     &     /,'                        element type, element no.:',
     &     2i8,
     &     /,'                        job terminated...')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ebkdrv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/27/09 rhd               *
c     *                                   comments about threads     *
c     *                                                              *
c     *     this subroutine conducts the back solve stage            *
c     *     of the ebe preconditioner for a block of elements        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ebkdrv(r,z,local_diag,refdof,mydof, blk)
      use elem_block_data, only: pcm_blocks, edest_blocks
      use pcg_data 
      use mpi_lnpcg
      implicit integer (a-z)
c
$add common.main
c
#sgl      real
#dbl      double precision
     & dummy(1), r(*), z(*), local_diag(*), zero
      logical debug, valid_type
      data debug /.false./
      integer, dimension(:,:), pointer :: edest
#sgl      real,
#dbl      double precision,
     &     dimension(:,:), pointer :: pcm
      data zero /0.0/
c
c
c	 At present, this process is not executed in thread
c        parallel. The blocking within an MPI doamin is not
c        "vectorized."  Vectorized blocking is required for the
c        threads only version of WARP3D. Then the span loop in
c        the ebksub routines can be executed in thread parallel.
c
c        It may be possible to relax serial execution here for the
c        span loop in ebksub routine. Serial here or vectorized
c        blocking is "strictly" required to maintain complete
c        symmetry of the forward/backward application of the
c        ebe preconsditioner to the current residual.
c
c        Some numerical testing may show this can be relaxed.                   
c
      if ( debug ) write (*,*) ' back substitution.'
c
      span   = elblks(0,blk)
      felem  = elblks(1,blk)
      totdof = iprops(2,felem)*iprops(4,felem)
      utsz = (totdof * ( totdof + 1 )) / 2
      type  =  iprops(1,felem)
      edest  => ledest_blocks(blk)%ptr
      pcm    => pcm_blocks(blk)%ptr
      valid_type = ( type .ge. 1 .and. type .le. 5 ) 
     &             .or. type .eq. 12
      if ( .not. valid_type ) then
        write(*,9000) type, felem
        call die_abort
      end if
c
c	 the 8-node elements use unrolled code. others use general
c        code.  Based on benchmarks in late 2009.
c
      if ( type .eq. 2 .or. type .eq. 12 ) then
        call ebksub8( span, totdof, edest(1,1), z, pcm(1,1), utsz )
      else 
        call ebksub_gen( span, totdof, edest(1,1), z, pcm(1,1), utsz )
      endif
c
      return
c
 9000 format('>>> FATAL SYSTEM ERROR: warp3d does not have the',
     &     /,'                        code implemented to process',
     &     /,'                        an element type. rotuine ebkdrv',
     &     /,'                        element type, element no.:',
     &     2i8,
     &     /,'                        job terminated...')
c
      end




















