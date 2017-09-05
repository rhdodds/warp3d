c     ****************************************************************
c     *                                                              *
c     *                      subroutine prcsel                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 8/11/2017 rhd              *
c     *                                                              *
c     *     sets up all the elements                                 *
c     *     to have their properties placed in global                *
c     *     storage, and calls the appropriate processing subroutine *
c     *     based on element type.                                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine prcsel
c
      use global_data ! old common.main
      use main_data, only : elstor
c
      implicit none
c
      integer :: elem, type, now_thread
      logical, parameter :: local_debug = .false.
      double precision, external ::  omp_get_wtime
      double precision :: start_time
c
c                       for each element in the structure, store
c                       element properties in permanent storage.
c
      call omp_set_dynamic( .false. )
      if( local_debug ) start_time = omp_get_wtime()
c
c$OMP PARALLEL DO  PRIVATE( elem, type, now_thread )
       do elem = 1, noelem
c
       type = elstor(1,elem)
       select case( type )
c
c                       type 1: q3disop
c
         case( 1 )
         call elprp1( elem, type )
c
c                       type 2: l3disop
c
         case( 2 )
         call elprp2( elem, type )
c
c                       type 3: ts12disop
c
         case( 3 )
         call elprp3( elem, type )
c
c                       type 4: ts15disop
c
         case( 4 )
         call elprp4( elem, type )
c
c                       type 5: ts9disop
c
         case( 5 )
         call elprp5( elem, type )
c
c                       type 6: tet10
c
         case( 6 )
         call elprp6( elem, type )
c
c                       type 7: wedge15
c
         case( 7 )
         call elprp7( elem, type )
c
c                       type 8: tri6
c
         case( 8 )
         call elprp8( elem, type )
c
c                       type 9: quad8
c
         case( 9 )
         call elprp9( elem, type )
c
c                       type 10: axiquad8
c
         case( 10 )
         call elprp10( elem, type )
c
c                       type 11: axitri6
c
         case( 11 )
         call elprp11( elem, type )
c
c                        type 12: inter_8
c
         case( 12 )
         call elprp12( elem, type  )
c
c                        type 13: tet4
c
         case( 13 )
         call elprp13( elem, type  )
c
c                        type 14: trint6
c
         case( 14 )
         call elprp14( elem, type  )
c
c                        type 15: trint12
c
         case( 15 )
         call elprp15( elem, type  )
c
c                        type 18: bar2
c
         case( 18 )
         call elprp18( elem, type  )
c
c                        type 19: link2
c
         case( 19 )
         call elprp19( elem, type  )
c         
         case default
             write(out,9000) type, elem
             call die_abort
c
       end select
c
      end do
c$OMP END PARALLEL DO
c
      return
c
 9000 format(/,'>>>> FATAL ERROR: prcsel. elem, type: ', 2i8,
     &       /,'                  job terminated' )
      end
