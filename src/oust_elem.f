c     ****************************************************************
c     *                                                              *
c     *                      subroutine oust_elem                    *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 2/10/2017 rhd              *
c     *                                                              *
c     *     output element (center) stress or strain results         *
c     *     to (1) patran file in either binary or formatted forms   *
c     *     or (2) flat text or stream file                          *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine oust_elem( stress, oubin, ouasc, num_vals,
     &                      elem_values, nrowd, nrow, first_block,
     &                      felem, last_block, flat_file,
     &                      stream_file, text_file, compressed )
      use global_data ! old common.main
      implicit none
c
      integer :: num_vals, nrowd, nrow, felem
      logical :: oubin, ouasc, stress, first_block, last_block,
     &           flat_file, stream_file, text_file, compressed
c
      double precision :: elem_values(nrowd,num_vals)
c
c                       local declarations
c
      integer :: data_type, relelem, elem, i, quantity
      integer, save :: fileno, flat_file_number
      integer :: titl(80), titl1(80)
      integer, external ::  warp3d_get_device_number
      logical ::  connected
      logical ::  patran_file
      double precision, parameter :: small_tol = 1.d-80, zero = 0.d0

      character(len=4) :: title(80), title1(80)
      character(len=1) dummy_char
      character(len=80) :: string, strng1, stepstring*7
      equivalence (title,titl), (title1,titl1)
c
      patran_file = oubin .or. ouasc
      data_type = 1
      if( stress ) data_type = 2
c
c                       close patran file(s) or flat file
c
      if( last_block ) then
        call ouocst_elem( data_type, ltmstp, oubin, ouasc, fileno, 2,
     &                    use_mpi, myid, flat_file,
     &                    stream_file, text_file, compressed,
     &                    flat_file_number, dummy_char )
        return
      end if
c
c                       Patran: open either binary or formatted files
c                               or both for output. assign names to
c                               the files based on type of output,
c                               step number (ltmstp), and possibly
c                               MPI rank
c
c                       Flat: open either text or stream file
c                             assign name to the file based on type
c                             of output, step number (ltmstp), and
c                             possibly MPI rank
c
c
c                       Logic here mingles patran and flat file
c                       writing of element results since not much
c                       difference.
c
      if( first_block ) then
c
        fileno = warp3d_get_device_number()
        if( fileno .eq. -1 ) then
          write(out,9000)
          call die_abort
        end if
        call ouocst_elem( data_type, ltmstp, oubin, ouasc, fileno, 1,
     &                    use_mpi, myid, flat_file,
     &                    stream_file, text_file, compressed,
     &                    flat_file_number, dummy_char )
c
c                        header info files
c
        if( patran_file ) call oust_elem_patran_header
        if( flat_file .and. text_file ) then
          quantity = 1
          if( stress ) quantity = 2
          call ouddpa_flat_header( 3, quantity, fileno )
        end if
c
        if( use_mpi .and. flat_file ) then
          if( text_file )   write(fileno,*) "  ", num_vals
          if( stream_file ) write(fileno) num_vals
        end if
c
      end if ! on first_block test
c
c                        write element data records into Patran or
c                        flat file. For Patran file(s), the element
c                        type is hard coded into this line as a
c                        hex element
c
c                        zero small values to prevent 3-digit exponents
c                        in formatted files
c
      where( abs(elem_values) .lt. small_tol ) elem_values = zero
c
      if( patran_file .and. oubin ) then
         do relelem = 1, nrow
             elem = felem + relelem - 1
             write(fileno) elem, 8,
     &           (sngl(elem_values(relelem,i)),i=1,num_vals)
         end do
      end if
c
      if( patran_file .and. ouasc ) then
         do relelem = 1, nrow
            elem = felem + relelem - 1
            write(fileno,920) elem, 8,
     &              (elem_values(relelem,i),i=1,num_vals)
         end do
      end if
c
      if( flat_file .and. stream_file  ) then
         do relelem = 1, nrow
           elem = felem + relelem - 1
           if( use_mpi ) then
             write(fileno) elem, elem_values(relelem,1:num_vals)
           else
             write(fileno) elem_values(relelem,1:num_vals)
           end if
         end do
      end if
c
      if( flat_file .and. text_file ) then
         do relelem = 1, nrow
           elem = felem + relelem - 1
           if( use_mpi ) then
             write(fileno,9100) elem,
     &                          elem_values(relelem,1:num_vals)
           else
             write(fileno,9200) elem_values(relelem,1:num_vals)
           end if
         end do
      end if
c
      return
c
 900  format(80a1)
 916  format(i10)
 915  format(i5)
 920  format(2i8,/,(6e13.6))
 9000 format('>> WARNING: could not find a unit number to write',
     & /     '            patran or flat file. action skipped...',/)
 9100 format(i8,30d15.6)
 9200 format(30e15.6)
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oust_elem_patran_header               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/16/2014 (rhd)            *
c     *                                                              *
c     ****************************************************************
c
      subroutine oust_elem_patran_header
      implicit none
c
      string = ' '
      strng1 = ' '
      if( stress ) then
          string = 'element stress results for structure '
      else
          string = 'element strain results for structure '
      end if
      strng1 = 'loading '
c
      string(38:45) = stname
      strng1(9:16)  = lsldnm
      write(stepstring,'(i7)') ltmstp
      string(46:) = ', step '//stepstring
c
      do i = 1, 80
         title(i)  = string(i:i)
         title1(i) = strng1(i:i)
      end do
c
c                     write header records for each Patran
c                     file type.
c
c                     flat files do not have header records
c                     at this time.
c
      if( patran_file .and. oubin ) then
       write(fileno) titl,num_vals
       write(fileno) titl1
       write(fileno) titl1
       if ( use_mpi ) write(fileno) noelem
      end if
c
      if( patran_file .and. ouasc ) then
       write(fileno,900) titl
       write(fileno,915) num_vals
       write(fileno,900) titl1
       write(fileno,900) titl1
       if ( use_mpi ) write(fileno,916) noelem
      end if
c
      return
c
 900  format(80a1)
 915  format(i5)
 916  format(i10)
c
      end subroutine oust_elem_patran_header

      end subroutine oust_elem
