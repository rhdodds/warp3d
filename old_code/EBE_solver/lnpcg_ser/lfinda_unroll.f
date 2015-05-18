c     ****************************************************************
c     *                                                              *
c     *                      subroutine lfinda                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/27/09 rhd               *
c     *                                   added thread parallel
c     *                                                              *
c     *     this subroutine computes the step length for the         *
c     *     current iteration of the linear pcg algorithm.           *
c     *                                                              *
c     *     this version calls the unrolled versions of lmprd, the   *
c     *     matrix-vector multiplication routine.  There is an       *
c     *     unrolled version for each element in warp3d.             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lfinda( step, stpitr, iter, trcpcg, z1tr1, alpha,
     &                   sdir )
c
      use elem_block_data, only: estiff_blocks, edest_blocks 
c
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &  z1tr1, alpha, ptap, zero, sdir(*)
c
#dbl      double precision,
#sgl      real,
     & dimension(:,:), pointer :: estiff_block
       integer, dimension(:,:), pointer :: edest_block
#dbl      double precision,
#sgl      real,
     & dimension(:,:), allocatable :: ifv_threads
c
      logical trcpcg, valid_type
      data zero /0.0/
c
c                       perform linear matrix product. the lm
c                       product is stored in the internal force
c                       vector, as the ifv will not be used again
c                       before its next update.
c
c                       initialize the product vector.
c
c
      ifv(1:nodof) = zero
      allocate( ifv_threads(nodof,num_threads) )
      call zero_vector( ifv_threads, num_threads*nodof )
c
c
c                       perform the mutliplication in element blocks.
c
c                       to provide faster execution, fully unrolled 
c                       versions of lmprd for each element type are provided.
c
c
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO PRIVATE( blk, now_thread, span, felem, type,
c$OMP&                     totdof, estiff_block, edest_block,
c$OMP&                     valid_type ) 
c$OMP&         SHARED( nelblk, elblks, iprops, estiff_blocks,
c$OMP&                 edest_blocks, sdir, ifv_threads )
      do blk = 1, nelblk
         now_thread = omp_get_thread_num() + 1
         span   = elblks(0,blk)                   
         felem  = elblks(1,blk)
         type   = iprops(1,felem)
         totdof = iprops(2,felem) * iprops(4,felem)
         utsz   = (totdof * ( totdof + 1 )) / 2
         estiff_block => estiff_blocks(blk)%ptr
         edest_block  => edest_blocks(blk)%ptr
         valid_type = ( type .ge. 1 .and. type .le. 5 ) 
     &                .or. type .eq. 12
         if ( .not. valid_type ) then
             write(*,9000) type, felem
             call die_abort
         end if
c
         felem = 1
c
         if ( type . eq. 2 .or. type .eq. 12 )  then
            call lmprd_8node( sdir, ifv_threads(1,now_thread), 
     &           estiff_block(1,1), span, edest_block(1,1) )
         else 
            call lmprd_gen( sdir, ifv_threads(1,now_thread), 
     &           estiff_block(1,1), span, totdof, utsz,
     &           edest_block(1,1) )
         endif
      end do
c$OMP END PARALLEL DO
c
c                       reduction of thread ifv vectors into
c                       unique system vector.
c
      do j = 1, num_threads
       ifv(1:nodof) = ifv(1:nodof) + 
     &                       ifv_threads(1:nodof,j)
      end do
      deallocate( ifv_threads )
c
c                       compute the triple product sdir transpose * lm *
c                       sdir.
c
      ptap = zero
      do i = 1, nodof
         ptap = ptap + sdir(i)*ifv(i)
      end do
c
c                      compute the step length.
c
      alpha = z1tr1 / ptap
c
c                       output the value of alpha if the trace flag
c                       is on.
c
      if( trcpcg ) then
         call oualph( step, stpitr, iter, alpha )
      end if
c
      return
c
 9000 format('>>> FATAL SYSTEM ERROR: warp3d does not have the',
     &     /,'                        code implemented to process',
     &     /,'                        an element type. routine lfinda',
     &     /,'                        element type, element no.:',
     &     2i8,
     &     /,'                        job terminated...')
c
      end













