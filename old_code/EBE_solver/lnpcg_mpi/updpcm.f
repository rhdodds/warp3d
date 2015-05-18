c     ****************************************************************
c     *                                                              *
c     *                      subroutine updpcm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/27/2009 rhd             *
c     *                                   add threads                *
c     *                                                              *
c     *     this subroutine updates the preconditioning matrix       *
c     *                                                              *
c     ****************************************************************
c
c                      
c
      subroutine updpcm
c
      use elem_block_data, only : pcm_blocks, edest_blocks
      use pcg_data,        only : ebe_pre, pdiag, diag
      use main_data,       only : kdiag, local_mdiag
      use mpi_lnpcg
c
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &     one, zero, dumd, nfac, sign
      real dumr
      character *1 dums
      logical debug1
c
      integer, dimension(:,:), pointer :: edest
#dbl      double precision,
#sgl      real,
     & dimension(:,:), pointer :: pcm
c
#dbl      double precision,
#sgl      real,
     & dimension(:,:), allocatable :: pdiag_threads
c
#dbl      double precision,
#sgl      real,
     & dimension(:), allocatable :: full_diag, chk_diag
c
      data zero, one, debug1 / 0.0, 1.0, .false. /
c
c
      if ( debug1 ) write (*,*) myid,':>>>>> inside updpcm'
      call thyme( 10, 1 )
c                                   
c                       tell the slave processors to start updpcm, then
c                       allocate lnpcg vectors if not already done
c
      call wmpi_alert_slaves ( 32 ) 
c
c     
      if ( .not. allocated (diag) ) call mem_allocate( 17 )
c
c                       diagonal vector contains the current
c                       estmate of the tangent stiffness and the
c                       nodal mass scaled by newmark beta and dt.
c
      refdof = local_nodes%num_local_nodes * 3
      mydof = (local_nodes%num_private + 
     &     local_nodes%num_own_shared)*3
      nfac = one / (nbeta*dt*dt)
      do i = 1, refdof
          diag(i) =  kdiag(i) + local_mdiag(i)*nfac
      end do
c      if ( root_processor) write (*,*) '>>> Here is full kdiag:'
c      call dd_print(kdiag)
c      if ( root_processor) write (*,*) '>>> Here is full local_mdiag:'
c      call dd_print(local_mdiag)
c      if ( root_processor) write (*,*) '>>> Here is full diag:'
c      call dd_print(diag)
c      call wmpi_wait
c
c                      if the pcm is diagonal then return
c
      if ( .not. ebe_pre ) then
         go to 9999
      end if
c
c                       create the preconditioning matrix and 
c                       apply the constraints to it and the
c                       diagonal vector.
c
      call applcn  
c
c                      invert and take the square root of the
c                      diagonal vector. initialize the element
c                      diagonal scaling vector.
c
      do i = 1, refdof
	 if ( diag(i) .lt. zero ) then
	    sign = -one
	 else
	    sign = one
	 end if
         diag(i)  = sign / sqrt(abs(diag(i)))           
      end do
c
c                      set up pdiag vector for each thread, then reduce
c                      at end of block loop
c
      allocate( pdiag_threads(refdof,num_threads) )
      pdiag(1:refdof) = one
      pdiag_threads(1:refdof,1:num_threads) = one
c
c                      perform winget regularization and crout 
c                      factorization. multiply the element diagonal 
c                      scaling vector. can be threaded w/o problems
c                      since the "order" of computing the pre-conditioning
c                      matrix for each element does not matter.
c
c$OMP PARALLEL DO DEFAULT( shared )
c$OMP&            PRIVATE( blk, now_thread, span, felem, totdof, utsz, 
c$OMP&            edest, pcm ) 
      do blk = 1, nelblk
c
         now_thread = omp_get_thread_num() + 1
         if ( elblks(2,blk) .ne. myid ) cycle
         span   = elblks(0,blk)
         felem  = elblks(1,blk)
         totdof = iprops(2,felem)*iprops(4,felem)
         utsz   = ( totdof * ( totdof + 1 ) ) / 2
         edest  => ledest_blocks(blk)%ptr
         pcm    => pcm_blocks(blk)%ptr
c
         call wnregl( span, cp, dcp, totdof, edest(1,1),
     &                diag, pcm(1,1), utsz )
         call croutf( span, cp, totdof, pcm(1,1), utsz )
         call edgscl( span, dcp, totdof, edest(1,1),
     &                pdiag, pcm(1,1), utsz )
c
      end do
c$OMP END PARALLEL DO
c
c                       reduction of thread pdiag vectors into
c                       unique system vector.
c
      do j = 1, num_threads
       pdiag(1:refdof) = pdiag(1:refdof) * pdiag_threads(1:refdof,j)
      end do
      deallocate( pdiag_threads )
c
c                      sync up pdiag
c      
      call lnpcg_comm( pdiag, mydof, refdof, 2 )
c
c                      if we are having miserable lives, print out the
c                      peconditioning matrix
c
      if ( debug1 ) then
         write(*,*) '>>>>>> the pcm matrix'
         do blk = 1, nelblk
            if ( elblks(2,blk) .ne. myid ) then
               call wmpi_wait
               cycle
            endif
            span   = elblks(0,blk)
            felem  = elblks(1,blk)
            do elem = 1, span
               write (*,'("pcm for element:",i5)') elem + felem - 1
               write (*,'(2x,6e13.4)') (pcm_blocks(blk)%ptr(j,elem),
     &              j=1,300)
c               tot = 0 
c               do j= 1, 24
c                  write (*,'(i2,1x,12e13.4)') j,
c     &              (pcm_blocks(blk)%ptr(k+tot,elem),k=1,j)
c                  tot = tot + j
c               enddo
            enddo
            call wmpi_wait
         enddo
         if ( root_processor) write (*,*) '>>> Here is full pdiag:'
         call dd_print(pdiag)
         if ( root_processor) write (*,*) '>>> Here is full diag:'
         call dd_print(diag)
      endif
      call wmpi_wait
c
c                      deactivate flags necessitating a new pcm.
c    
 9999 continue
      newcns = .false.
      newmas = .false.
      newstf = .false.
c
      call thyme( 10, 2 )
      if ( debug1 ) write (*,*) '<<<<< leaving updpcm'
c
      return
      end



