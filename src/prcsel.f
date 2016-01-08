c     ****************************************************************
c     *                                                              *
c     *                      subroutine prcsel                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 08/22/97                   *
c     *                                 : 12/23/00 sushovan          *
c     *                                                              *
c     *     this subroutine sets up all the elements                 *
c     *     to have their properties placed in global                * 
c     *     storage, and calls the appropriate processing subroutine *
c     *     based on element type.                                   * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine prcsel
      use main_data, only : elstor
      implicit integer (a-z)     
$add common.main
c
c                       for each element in the structure, store
c                       element properties in permanent storage.
c                                   
      call omp_set_dynamic( .false. )
      if( local_debug ) start_time = omp_get_wtime()
c
c$OMP PARALLEL DO  PRIVATE( elem, type, now_thread ) 
       do 10 elem = 1, noelem
c
       type = elstor(1,elem) 
c
c                       branch on element type.
c
         go to ( 100, 200, 300, 400, 500, 600, 700, 800, 900,
     &          1000, 1100, 1200, 1300, 1400, 1500 ), type
c
c                       type 1: q3disop
c
 100     continue
         call elprp1( elem, type )
         go to 10
c
c                       type 2: l3disop
c
 200     continue
         call elprp2( elem, type )
         go to 10
c
c                       type 3: ts12disop
c
 300     continue
         call elprp3( elem, type )
         go to 10
c
c                       type 4: ts15disop
c
 400     continue
         call elprp4( elem, type )
         go to 10
c
c                       type 5: ts9disop
c
 500     continue
         call elprp5( elem, type )
         go to 10
c
c                       type 6: tet10
c
 600     continue
         call elprp6( elem, type )
         go to 10
c
c                       type 7: wedge15
c
 700     continue
         call elprp7( elem, type )
         go to 10
c
c                       type 8: tri6
c
 800     continue
         call elprp8( elem, type )
         go to 10
c
c                       type 9: quad8
c
 900     continue
         call elprp9( elem, type )
         go to 10
c
c                       type 10: axiquad8
c
 1000    continue
         call elprp10( elem, type )
         go to 10
c
c                       type 11: axitri6
c
 1100    continue
         call elprp11( elem, type )
         go to 10
c
c                        type 12: inter_8
c
 1200    continue
         call elprp12( elem, type  )
         go to 10
c
c                        type 13: tet4
c
 1300    continue
         call elprp13( elem, type  )
         go to 10
c
c                        type 14: trint6
c
 1400    continue
         call elprp14( elem, type  )
         go to 10
c
c                        type 15: trint12
c
 1500    continue
         call elprp15( elem, type  )
         go to 10
c
 10   continue
c$OMP END PARALLEL DO
c
c
 9999 return
      end
