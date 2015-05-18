c     ****************************************************************
c     *                                                              *
c     *                      subroutine lfinda                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/27/2009 rhd             *
c     *                                   added thread parallel      *
c     *                                                              *
c     *     this subroutine computes the step length for the         *
c     *     current iteration of the linear pcg algorithm.           *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lfinda( step, stpitr, iter, trcpcg, z1tr1, alpha,
     &                   sdir, mydof, refdof )
      use elem_block_data, only: estiff_blocks
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &  z1tr1, alpha, ptap, zero, sdir(*)
      logical trcpcg, debug, valid_type
      data zero, debug /0.0, .false./
      integer, dimension(:,:), pointer :: ledest
#dbl      double precision,
#sgl      real,
     & dimension(:,:), allocatable :: ifv_threads
      external omp_get_thread_num

c
c                       perform linear matrix product. the lm
c                       product is stored in the internal force
c                       vector, as the ifv will not be used again
c                       before its next update.
c
c                       initialize the product vector.
c
c
      allocate( ifv_threads(refdof,num_threads) )
      ifv_threads(1:refdof,1:num_threads) = zero
      ifv(1:refdof) = zero
c
c                       perform the mutliplication in element blocks.
c
c     
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO PRIVATE( blk, now_thread, span, felem, type,
c$OMP&                     totdof, utsz, ledest, valid_type ) 
c$OMP&         SHARED( nelblk, elblks, myid, iprops, ledest_blocks,
c$OMP&                    sdir, ifv_threads, estiff_blocks )
      do blk = 1, nelblk
         now_thread = omp_get_thread_num() + 1
c
         if (elblks(2,blk) .ne. myid ) cycle
         span  = elblks(0,blk)                   
         felem = elblks(1,blk)
         type  =  iprops(1,felem)
         totdof = iprops(2,felem) * iprops(4,felem)
         utsz   = (totdof * ( totdof + 1 )) / 2
         ledest  => ledest_blocks(blk)%ptr
         valid_type = ( type .ge. 1 .and. type .le. 5 ) 
     &             .or. type .eq. 12
         if ( .not. valid_type ) then
            write(*,9000) type, felem
            call die_abort
         end if
c
         if ( type . eq. 2 .or. type .eq. 12 )  then
            call lmprd_8node( sdir, ifv_threads(1,now_thread), 
     &         estiff_blocks(blk)%ptr(1,1), span, ledest )
         else 
            call lmprd_gen( sdir, ifv_threads(1,now_thread), 
     &         estiff_blocks(blk)%ptr(1,1), span, totdof,
     &         utsz, ledest )
         endif
c
      end do
c$OMP END PARALLEL DO
c
c                       reduction of thread ifv vectors into
c                       unique system vector.
c
      do j = 1, num_threads
       ifv(1:refdof) = ifv(1:refdof) + 
     &                       ifv_threads(1:refdof,j)
      end do
      deallocate( ifv_threads )
c
c                       conduct the necessary communications between
c                       processors for the MPI version
c
      call lnpcg_comm ( ifv, mydof, refdof, 1 )
c
c                       compute the triple product sdir transpose * lm *
c                       sdir.
c
      ptap = zero
      do i = 1, mydof
         ptap = ptap + sdir(i)*ifv(i)
      end do
c
c                       reduce dot product among processors and rebroadcast
c
      call wmpi_dotprod (ptap)
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
c
      return
c
 9000 format('>>> FATAL SYSTEM ERROR: warp3d does not have the',
     &     /,'                        code implemented to process',
     &     /,'                        an element type. rotuine lfinda',
     &     /,'                        element type, element no.:',
     &     2i8,
     &     /,'                        job terminated...')
c
      end








