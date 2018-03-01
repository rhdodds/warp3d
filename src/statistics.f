c     ****************************************************************
c     *                                                              *
c     *                                                              *
c     *        subroutine to calculate necessary statistics for      *
c     *        direct solver                                         *
c     *                                                              *
c     *                  (1)   total number of equations             *
c     *                  (2)   storage for profile (megawords)       *
c     *                  (3)   max and avg col heights               *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 03/8/94                    *
c     *                                   06/20/94 ASG--put in zero  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine statistics ( stats, col_max, col_heigth, diag_pos, out)
      implicit integer (a-z)
      double precision
     &    words, mega_fact, avg_col_hgth, d1024, zero
      data d1024, zero /1024.0, 0.0/
      dimension  col_heigth(*), diag_pos(*)
c
      logical stats
       word_fact = 8
       stats        = .false.
       mega_fact    = d1024 * d1024
       words        = real(diag_pos(col_max) * word_fact) / mega_fact
       sum          = 0
       max_col_hgth = 0
c
       do i = 1, col_max
          sum          = sum + col_heigth(i)
          max_col_hgth = max( max_col_hgth, col_heigth(i) )
       end do
c
       if ( col_max .ne. 0 ) then
           avg_col_hgth = real( sum ) / real( col_max )
       else
           avg_col_hgth = zero
       end if
c
       write (out,9000) col_max, words, max_col_hgth, avg_col_hgth
c
       return
c
 9000  format (
     & /10x, '>> direct solver statistics:'
     & /15x, 'total number of equations:       ',i12,
     & /15x, 'total megabytes needed for [k]:  ',f12.1,
     & /15x, 'maximum column heigth: ',i10,' average: ',f12.1,/  )
c
       end
