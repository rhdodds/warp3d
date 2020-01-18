c     ****************************************************************
c     *                                                              *
c     *                      subroutine patch_data                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/19/2019 rhd             *
c     *                                                              *
c     *     patch command to support development mods of data during *
c     *     execution                                                *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine patch_data( sbflg1, sbflg2 )
      use global_data ! old common.main
c
      use main_data
c
      implicit none
c
      logical ::sbflg1, sbflg2
      logical, external :: matchs, numd, integr, endcrd, true, label,
     &                     string, realn
c
c
      if( matchs('props',4) ) then
         call patch_data_props
         return
      else
         write(out,9000)
         call die_abort
      end if
c
 9000 format(2x,"Fatal error: invalid patch command...")
c
      contains
c     ========
c
      subroutine patch_data_props
      implicit none
c
      real :: realval
      integer :: intval, row, col, itype
c      

      if( .not. matchs('row',3) ) then
          write(out,9000) 
          call die_abort
      end if
      if( .not. integr( row ) ) then
          write(out,9000) 
          call die_abort
      end if
      if( .not. matchs('column',3) ) then
          write(out,9000) 
          call die_abort
      end if
      if( .not. integr( col ) ) then
          write(out,9000) 
          call die_abort
      end if
c
      itype = 1
      if( .not. integr( intval ) ) then
          if( .not. realn( realval ) ) then
            write(out,9000) 
            call die_abort
          end if
          itype = 2
      end if
c
      if( itype == 1 ) iprops(row,col) = intval
      if( itype == 2 ) props(row,col)  = realval
c
      return
c
 9000 format(2x,"Fatal error: invalid patch command...")
c
      end subroutine patch_data_props
      end subroutine patch_data
