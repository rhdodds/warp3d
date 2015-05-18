c     ****************************************************************
c     *                                                              *
c     *                      subroutine addifv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 05/25/04 rhd               *
c     *                                                              *
c     *     this subroutine assembles the internal force vectors of  *
c     *     a block of similar, non-conflicting elements into the    *
c     *     global internal force vector.                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine addifv( span, bedst, totdof, ifv, iprops, 
     &                   felem, sum_ifv, num_term_ifv,
     &                   eleifv, dam_ifv, dam_state )
      use damage_data, only : dam_ptr, growth_by_kill
      implicit integer (a-z)
$add param_def
c
#dbl      double precision
#sgl      real
     &   ifv(*), sum_ifv, eleifv(span,*), dam_ifv(mxedof,*)
c
      integer bedst(totdof,*), iprops(mxelpr,*),
     &        dam_state(*)
      logical debug

      debug = .false.
c
      if ( debug ) write (*,*) '>>>>  inside addifv'
c
      do j = 1, totdof    
         do i = 1, span
            ifv(bedst(j,i)) = ifv(bedst(j,i)) + eleifv(i,j)
            sum_ifv = sum_ifv + abs(eleifv(i,j))
         end do
      end do
      num_term_ifv = num_term_ifv + totdof * span
c
      if( debug ) then
       write(*,9200) felem
       do i = 1, span
         write(*,*) ' '
         write(*,*) 'element: ', felem+i-1
         write(*,9300) eleifv(i,1:totdof)
       end do
      end if
c
c             return if crack growth by element extinction
c             not specified. return if block does not contain
c             killable elems
c
      if ( .not. growth_by_kill ) go to 9999
      if ( iand( iprops(30,felem),2 ) .eq. 0 ) go to 9999
c
c             loop for each killable element defined in this block.
c             if we have NOT started releasing the internal forces,
c             store the element contribution to the ifv in the
c             crack growth data structures.         
c             this is not real efficient but we don't know here
c             which elements have just been killed.
c
      if( .not. allocated( dam_ptr ) ) then ! sanity check
          write(*,9100)
          call die_abort
      end if
c
      do i = 1, span
         element = felem + i - 1
         if ( dam_state(dam_ptr(element)) .eq. 0 ) then
            do j = 1, totdof
               dam_ifv(j,dam_ptr(element)) = eleifv(i,j)
            end do
         end if
      end do
c
 9999 continue
      if ( debug ) write (*,*) '<<<<  leaving  addifv'
      return
c
 9000 format(1x,'dof:',i2,' element ifv:',e14.6,' total ifv:',e14.6)
 9100 format(1x,'FATAL ERROR: in adifv. contact WARP3D developers')
 9200 format(/,2x,'... addifv for block with first element: ',i6)
 9300 format(5x,8e14.6)
c
      end
