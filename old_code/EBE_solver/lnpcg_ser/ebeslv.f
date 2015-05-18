c     ****************************************************************
c     *                                                              *
c     *                      subroutine ebeslv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/27/09 rhd               *
c     *                                   comments about threads     *
c     *                                                              *
c     *     this subroutine creates the pseudo-residual by           *
c     *     multiplying the inverse of the preconditioner and the    *
c     *     residual.  The pseudo-residual allows faster convergence.*
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ebeslv(r,z)
      use elem_block_data, only: pcm_blocks, edest_blocks
      use pcg_data 
      implicit integer (a-z)
c
$add common.main
c
#sgl      real
#dbl      double precision
     & dummy(1), r(*), z(*)
      logical debug, valid_type
      data debug /.false./
      integer, dimension(:,:), pointer :: edest
#sgl      real,
#dbl      double precision,
     &  dimension(:,:), pointer :: pcm

c
      if ( debug ) write(*,*) '>>>>> inside ebeslv'
      if ( .not.nostif .and. ebe_pre ) then
c
c                       perform initial global diagonal scaling.
c 
      call vec_ops( r, diag, z,  nodof, 1 )
c
c                       perform element forward reduction.
c
c                       to enable faster execution, efrwrd is fully
c                       unrolled for some element types.
c
c                       loop over elements inside the block is threaded
c
      if ( debug ) write (*,*) ' forward reduction'
      do blk = 1, nelblk
         span   = elblks(0,blk)
         felem  = elblks(1,blk)
         totdof = iprops(2,felem)*iprops(4,felem)
         utsz   = (totdof * ( totdof + 1 )) / 2
         type   = iprops(1,felem)
         edest  => edest_blocks(blk)%ptr
         pcm    => pcm_blocks(blk)%ptr 
         valid_type = ( type .ge. 1 .and. type .le. 5 ) 
     &                .or. type .eq. 12
         if ( .not. valid_type ) then
             write(*,9000) type, felem
             call die_abort
         end if
c
         if ( type .eq. 2 .or. type .eq. 12 ) then
            call efrwrd8 ( span, totdof, edest(1,1), z,
     &           pcm(1,1), utsz )
         else 
            call efrwrd_gen ( span, totdof, edest(1,1), z,
     &           pcm(1,1), utsz )
         endif 
c
      end do
c
c                       perform element diagonal scaling. 
c
         if ( debug ) write (*,*) ' diagonal scaling'
         call vec_ops( pdiag, z, dummy, nodof, 2 )
c
c                       perform element back substitution.
c                       to enable faster execution, ebksub is fully
c                       unrolled for some elements
c
         if ( debug ) write (*,*) ' back substitution.'
         do blk = nelblk, 1, -1
            span   = elblks(0,blk)
            felem  = elblks(1,blk)
            totdof = iprops(2,felem)*iprops(4,felem)
            utsz   = (totdof * ( totdof + 1 )) / 2
            type   = iprops(1,felem)
            edest  => edest_blocks(blk)%ptr
            pcm    => pcm_blocks(blk)%ptr
c
            if ( type .eq. 2 .or. type .eq. 12 ) then
               call ebksub8 ( span, totdof, edest(1,1), z,
     &              pcm(1,1), utsz )
            else 
               call ebksub_gen( span, totdof, edest(1,1), z,
     &              pcm(1,1), utsz )
            endif
         end do
c
c                       perform final global diagonal scaling.
c
         if (debug) write (*,*) ' final diagonal scaling'
         call vec_ops( diag, z, dummy, nodof, 2 )
c
c
      else
         if (debug) write (*,*) 'preconditioner is diagonal.'
c
c                       the preconditioning matrix is diagonal.
c
          call vec_ops( r, diag, z, nodof, 3 )
c
      end if
c
c
 9999 continue
      if ( debug ) write(*,*) '<<<<< leaving ebeslv' 
      return
c
 9000 format('>>> FATAL SYSTEM ERROR: warp3d does not have the',
     &     /,'                        code implemented to process',
     &     /,'                        an element type. routine ebeslv',
     &     /,'                        element type, element no.:',
     &     2i8,
     &     /,'                        job terminated...')
c
      end











