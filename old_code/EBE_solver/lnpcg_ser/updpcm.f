c     ****************************************************************
c     *                                                              *
c     *                      subroutine updpcm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/27/09 (rhd)             *
c     *                                   added thread support       *
c     *                                                              *
c     *     this subroutine updates the diagonal and HW              *
c     *     preconditioning matrix                                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine updpcm
c
      use elem_block_data, only : pcm_blocks, edest_blocks
      use pcg_data,        only : ebe_pre, pdiag, diag
      use main_data,       only : kdiag, mdiag
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
      data zero, one, debug1 / 0.0, 1.0, .false. /
c
c
      if ( debug1 ) write (*,*) '>>>>> inside updpcm'
      call thyme( 10, 1 )
c                                   
c                       diagonal vector contains the current
c                       estmate of the tangent stiffness and the
c                       nodal mass scaled by newmark beta and dt.
c
      nfac = one / (nbeta*dt*dt)
      do i = 1, nodof
          diag(i) =  kdiag(i) + mdiag(i)*nfac
      end do
c
c                      if the pcm is diagonal then return
c
      if ( .not. ebe_pre ) go to 9999
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
      do i = 1, nodof
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
      pdiag = one
      allocate( pdiag_threads(nodof,num_threads) )
      pdiag_threads = one
c
c                      perform winget regularization and crout 
c                      factorization. multiply the element diagonal 
c                      scaling vector. the order of elements does
c                      not matter in computing the preconditioning
c                      matrix for each element.
c
c$OMP PARALLEL DO DEFAULT( shared )
c$OMP&            PRIVATE( blk, now_thread, span, felem, totdof, utsz, 
c$OMP&            edest, pcm ) 
      do blk = 1, nelblk
         now_thread = omp_get_thread_num() + 1
         span   = elblks(0,blk)
         felem  = elblks(1,blk)
         totdof = iprops(2,felem)*iprops(4,felem)
         utsz   = ( totdof * ( totdof + 1 ) ) / 2
         edest  => edest_blocks(blk)%ptr
         pcm    => pcm_blocks(blk)%ptr
c
         call wnregl( span, cp, dcp, totdof, edest(1,1),
     &                diag, pcm(1,1), utsz )
         call croutf( span, cp, totdof, pcm(1,1), utsz )
         call edgscl( span, dcp, totdof, edest(1,1),
     &                pdiag_threads(1,now_thread), pcm(1,1), utsz )
      end do
c$OMP END PARALLEL DO
c
c                       reduction of thread pdiag vectors into
c                       unique system vector.
c
      do j = 1, num_threads
       pdiag(1:nodof) = pdiag(1:nodof) * pdiag_threads(1:nodof,j)
      end do
      deallocate( pdiag_threads )
c
c                      optional print out the peconditioning matrix
c
      if ( debug1 ) then
         write(*,*) '>>>>>> the pcm matrix'
         do blk = 1, nelblk
            span   = elblks(0,blk)
            felem  = elblks(1,blk)
            do elem = 1, span
               write (*,'(" pcm for element:",i5)') elem + felem - 1
               tot = 0 
               do j= 1, 24
                  write (*,'(i2,1x,12e13.4)') j,
     &              (pcm_blocks(blk)%ptr(k+tot,elem),k=1,j)
                  tot = tot + j
               enddo
            enddo
         enddo
         write(*,*) '<<<<< end of the pcm matrix'
         write(*,*) '>>>>>> diag:'
         write (*,'(2x,6e13.4)') (pdiag(j),j=1, nodof)
      endif
         
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



