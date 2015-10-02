c     ****************************************************************
c     *                                                              *
c     *                      subroutine addifv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 9/28/2015 rhd              *
c     *                                                              *
c     *     assembles the internal force vectors for a block of      *
c     *     similar, elements into the global internal force vector. *
c     *     the calling routine must insure that concurrent access   *
c     *     on other threads is prevented                            *
c     *     this is a classic scatter operation                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine addifv( span, bedst, totdof, ifv, iprops, 
     &                   felem, sum_ifv, num_term_ifv,
     &                   eleifv, dam_ifv, dam_state, out )
      use damage_data, only : dam_ptr, growth_by_kill
      implicit integer (a-z)
$add param_def
c
#dbl      double precision ::
#sgl      real ::
     &   ifv(*), sum_ifv, eleifv(span,*), dam_ifv(mxedof,*)
c
      integer :: bedst(totdof,*), iprops(mxelpr,*),
     &           dam_state(*)
      logical :: debug

      debug = .false.
c
      if( debug ) write (out,*) '>>>>  inside addifv'
c
      do j = 1, totdof    
@!DIR$ LOOP COUNT MAX=###  
         do i = 1, span
c$OMP ATOMIC UPDATE
            ifv(bedst(j,i)) = ifv(bedst(j,i)) + eleifv(i,j)
         end do
      end do
      do j = 1, totdof    
@!DIR$ LOOP COUNT MAX=###  
         do i = 1, span
            sum_ifv = sum_ifv + abs(eleifv(i,j))
         end do
      end do
      num_term_ifv = num_term_ifv + totdof * span
c
      if( debug ) then
       write(out,9200) felem
       do i = 1, span
         write(out,*) ' '
         write(out,*) 'element: ', felem+i-1
         write(out,9300) eleifv(i,1:totdof)
       end do
      end if
c
c             return if crack growth by element extinction
c             not specified. return if block does not contain
c             killable elems
c
      if( .not. growth_by_kill ) return
      if( iand( iprops(30,felem),2 ) .eq. 0 ) return
c
c             loop for each killable element defined in this block.
c             if we have NOT started releasing the internal forces,
c             store the element contribution to the ifv in the
c             crack growth data structures.         
c             this is not real efficient but we don't know here
c             which elements have just been killed.
c
      if( .not. allocated( dam_ptr ) ) then ! sanity check
          write(out,9100)
          call die_abort
      end if
c
@!DIR$ LOOP COUNT MAX=###  
      do i = 1, span
         element = felem + i - 1
         if( dam_state(dam_ptr(element)) .eq. 0 ) then
            do j = 1, totdof
               dam_ifv(j,dam_ptr(element)) = eleifv(i,j)
            end do
         end if
      end do
c
      if( debug ) write (out,*) '<<<<  leaving  addifv'
      return
c
 9000 format(1x,'dof:',i2,' element ifv:',e14.6,' total ifv:',e14.6)
 9100 format(1x,'FATAL ERROR: in adifv. contact WARP3D developers')
 9200 format(/,2x,'... addifv for block with first element: ',i6)
 9300 format(5x,8e14.6)
c
      end
