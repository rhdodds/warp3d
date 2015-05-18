c     ****************************************************************
c     *                                                              *
c     *                      subroutine gpifv1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/2/12 rhd                *
c     *                                                              *
c     *     computes the internal resisting force                    *
c     *     vectors for a block of similar elements in uniform       *
c           global coordinates for a integration point               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gpifv1( eleifv, nrow_ifv, gpn, weight, dj, b, 
     &                   urcs_blk_n1, local_work, element_volumes )
c
      use segmental_curves, only : max_seg_points, max_seg_curves
      implicit integer (a-z)
$add param_def 
$add include_sig_up
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &  eleifv(nrow_ifv,*), weight, dj(*), element_volumes(*),
     &  b(mxvl,mxedof,*), urcs_blk_n1(mxvl,*)
c
c                       local declarations
c
      logical geonl, bbar
c
#dbl      double precision
#sgl      real
     &   qtn1(mxvl,nstr,nstr), cs_blk_n1(mxvl,nstr),
     &   eps_bbar, w, scalar 
c
c
      span            = local_work%span
      felem           = local_work%felem
      type            = local_work%elem_type
      nnode           = local_work%num_enodes  
      geonl           = local_work%geo_non_flg
      bbar            = local_work%bbar_flg
      eps_bbar        = local_work%eps_bbar
      totdof          = local_work%totdof
      iter            = local_work%iter
c
c                       if element is triangle, wedge, tet sacle
c                       weights to correctly integrate B*sigma. 
c
#sgl      scalar = 1.0
#dbl      scalar = 1.0d00
      if ( local_work%adjust_const_elem  ) then
        call adjust_scalar_weights( type, scalar )
      end if
c
      w = weight*scalar
c
c                       compute the strain-displacement matrices
c                       for the gauss point. for large
c                       displacement analysis we are computing [B]
c                       linear using n+1 configuration and
c                       multiplying into Cauchy stresses.
c                       the cohesive materials have only 3 tractions
c                       located in 1st 3 cols of urcs..
c
      if ( local_work%is_cohes_elem ) then
        call blcmp_cohes( span, b, local_work%cohes_rot_block,
     &                    local_work%shape(1,gpn), type, nnode )
c
          do j = 1, totdof
            do i = 1, span
               eleifv(i,j) = eleifv(i,j) +
     &                       b(i,j,1) * urcs_blk_n1(i,1) * w * dj(i) +
     &                       b(i,j,2) * urcs_blk_n1(i,2) * w * dj(i) +
     &                       b(i,j,3) * urcs_blk_n1(i,3) * w * dj(i) 
            end do
          end do
          return
      end if
c
c                       solid elements
c
      call blcmp1( span, b,
     &               local_work%gama(1,1,1,gpn),
     &               local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &               local_work%nzeta(1,gpn), nnode )
c
      if ( bbar .and. type .eq. 2 ) then
            call bmod ( b, local_work%vol_block, 
     &                  span, mxvl, eps_bbar, mxedof ) 
      end if
c
c                       compute gauss point contributions to
c                       the internal force vectors for the element
c                       block. treat [B] as full even though it is
c                       not. b-bar fills upper part. eliminating
c                       zero multiplications makes code unreadable.
c                       for geonl, convert urcs -> cauchy stress at
c                       n+1 using [R,n+1] from polar decompositions
c                       during strain update.
c   
       if ( geonl ) then
          call getrm1( span, qtn1, local_work%rot_blk_n1(1,1,gpn), 2 )
          call qmply1( span, mxvl, nstr, qtn1,
     &                 local_work%urcs_blk_n1(1,1,gpn),
     &                 cs_blk_n1 )
          do j = 1, totdof
            do i = 1, span 
               eleifv(i,j) = eleifv(i,j) +
     &                       b(i,j,1)*cs_blk_n1(i,1)*w*dj(i) +
     &                       b(i,j,2)*cs_blk_n1(i,2)*w*dj(i) +
     &                       b(i,j,3)*cs_blk_n1(i,3)*w*dj(i) +
     &                       b(i,j,4)*cs_blk_n1(i,4)*w*dj(i) +
     &                       b(i,j,5)*cs_blk_n1(i,5)*w*dj(i) +
     &                       b(i,j,6)*cs_blk_n1(i,6)*w*dj(i)
            end do 
          end do 
          go to 9000
      end if
c
c                       small displacement formulation. the urcs
c                       are stresses. no transformation needed.
c                        
      do j = 1, totdof
        do i = 1, span 
           eleifv(i,j) = eleifv(i,j) +
     &                   b(i,j,1)*urcs_blk_n1(i,1)*w*dj(i) +
     &                   b(i,j,2)*urcs_blk_n1(i,2)*w*dj(i) +
     &                   b(i,j,3)*urcs_blk_n1(i,3)*w*dj(i) +
     &                   b(i,j,4)*urcs_blk_n1(i,4)*w*dj(i) +
     &                   b(i,j,5)*urcs_blk_n1(i,5)*w*dj(i) +
     &                   b(i,j,6)*urcs_blk_n1(i,6)*w*dj(i)
        end do 
      end do 
      go to 9000
c
c                       update gauss point contribution to volume
c                       of element.
c
 9000 continue
      do i = 1, span 
        element_volumes(i) = element_volumes(i)   +   w * dj(i) 
      end do
c
      return
      end


