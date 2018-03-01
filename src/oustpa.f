c     ****************************************************************
c     *                                                              *
c     *                      subroutine oustpa                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/30/2017 rhd              *
c     *                                                              *
c     *     output stress or strain nodal results to (1) patran file *
c     *     in either binary or formatted forms, or (2) flat file    *
c     *     in stream or text formats                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oustpa( stress, oubin, ouasc, num_vals,
     &                   nodal_values, num_struct_nodes, flat_file,
     &                   stream_file, text_file, compressed )
      use global_data ! old common.main
      implicit integer (a-z)
c
      logical oubin, ouasc, stress, flat_file, stream_file, text_file,
     &       compressed
c
      type :: node_entry
         integer :: count
         double precision, dimension(:), pointer :: node_values
      end type node_entry
c
      type (node_entry), dimension (num_struct_nodes) :: nodal_values
         double precision, dimension(:), pointer :: snode_values
c
c                       local declarations
c
      logical :: patran_file
      real :: dum2
      character(len=4) :: title(80), title1(80)
      character(len=80) :: string, strng1, stepstring*7
      dimension titl(80), titl1(80)
      equivalence (title,titl), (title1,titl1)
      double precision :: small_tol, zero
      data small_tol, zero / 1.0d-80, 0.0d00 /
c
      patran_file = oubin .or. ouasc
      stepno = ltmstp ! more descriptive name for use here
c
c                       open Patran or flat files
c
      call ouocst( stress, stepno, oubin, ouasc, bnfile, fmfile, 1,
     &             out, use_mpi, myid, flat_file, stream_file,
     &             text_file, compressed, flat_file_number )
c
c                        set Patran title records.
c
      if( patran_file ) then
        string = ' '
        strng1 = ' '
        string = 'nodal strain results for structure '
        if( stress ) string = 'nodal stress results for structure '
        strng1 = 'loading '
c
        string(36:43) = stname
        strng1(9:16)  = lsldnm
        write(stepstring,'(i7)') stepno
        string(44:) = ', step '// stepstring
c
        do i = 1, 80
           title(i)  = string(i:i)
           title1(i) = strng1(i:i)
        end do
      end if   ! patran_file
c
c                       set dummy values not needed since the
c                       results are not displacement results.
c
      dum1 = 0
      dum2 = 0.0
      dum3 = 0
c
c                       write to a Patran binary file
c
      if( patran_file .and. oubin ) then
c
c                       records 1-3
c
         write(bnfile) titl, nonode, dum1, dum2, dum3, num_vals
         write(bnfile) titl1
         write(bnfile) titl1
c
c                       nodal results records. only those with non-zero
c                       counts are in this domain. we zero very small
c                       values to prevent 3 digit exponents in
c                       formatted output.
c
         do nod = 1, nonode
           if ( nodal_values(nod)%count .eq. 0 ) cycle
           snode_values => nodal_values(nod)%node_values
           if ( use_mpi .and. numprocs .gt.1 ) then
               write(bnfile) nod, nodal_values(nod)%count,
     &                       (snode_values(i),i=1,num_vals)
           else
         write(bnfile) nod, (sngl(snode_values(i)),i=1,num_vals)
           end if
         end do
      end if  ! patran & oubin
c
c                       write to a Patran formatted file
c
      if( patran_file .and. ouasc ) then
c
c                       records 1-4
c
         write(fmfile,900) titl
         write(fmfile,910) nonode, dum1, dum2, dum3, num_vals
         write(fmfile,900) titl1
         write(fmfile,900) titl1
c
c                       nodal results records. only those with non-zero
c                       counts are in this domain
c
         do nod = 1, nonode
          if ( nodal_values(nod)%count .eq. 0 ) cycle
          snode_values => nodal_values(nod)%node_values
          where( abs(snode_values) .lt. small_tol ) snode_values = zero
          if ( use_mpi .and. numprocs .gt.1 ) then
            write(fmfile,930) nod, nodal_values(nod)%count
            write(fmfile,940) snode_values(1:num_vals)
          else
            write(fmfile,920) nod, snode_values(1:num_vals)
          end if
         end do
      end if  ! patran  & ascii
c
      if( patran_file ) then  ! close and leave
         call ouocst( stress, stepno, oubin, ouasc, bnfile, fmfile,
     &                2, out, use_mpi, myid, flat_file, stream_file,
     &                text_file, compressed, flat_file_number )
         return
      end if
c
c                       write to a flat text or stream file
c
      if( .not. flat_file ) then ! big error somewhere
         write(out,9200)
         call die_abort
      end if
c
      if( text_file ) then
        quantity = 1
        if( stress ) quantity = 2
        call ouddpa_flat_header( 2, quantity, flat_file_number )
      end if
c
      if( use_mpi .and. numprocs .gt.1 ) then
         if( stream_file ) write(flat_file_number) num_vals
         if( text_file ) write(flat_file_number,*) "  ", num_vals
         do nod = 1, nonode
           if ( nodal_values(nod)%count .eq. 0 ) cycle
           snode_values => nodal_values(nod)%node_values
           where( abs(snode_values) .lt. small_tol )
     &            snode_values = zero
           if( stream_file) write(flat_file_number) nod,
     &               nodal_values(nod)%count,
     &               snode_values(1:num_vals)
           if( text_file) write(flat_file_number,9100) nod,
     &               nodal_values(nod)%count,
     &               snode_values(1:num_vals)
         end do
       else    !  just threaded or mpi w/ num_procs = 1
         do nod = 1, nonode
            if ( nodal_values(nod)%count .eq. 0 ) cycle
            snode_values => nodal_values(nod)%node_values
            where( abs(snode_values) .lt. small_tol )
     &             snode_values = zero
             if( stream_file) write(flat_file_number)
     &                    snode_values(1:num_vals)
            if( text_file) write(flat_file_number,9300)
     &               snode_values(1:num_vals)
         end do
      end if
c
c                       close flat file
c
      call ouocst( stress, stepno, oubin, ouasc, bnfile, fmfile, 2,
     &             out, use_mpi, myid, flat_file,
     &             stream_file, text_file, compressed,
     &             flat_file_number )
c
      return
c
 900  format(80a1)
 910  format(2i9,e15.6,2i9)
 920  format(i8,5(e13.6))
 930  format(2i8)
 940  format(e23.15)
 9100 format(2i9,30d15.6)
 9200 format(/1x,'>>>>> Fatal Error: routine oustpa.',
     &   /16x,   'should have not reach this point @ 1.',
     &   /16x,   'job terminated' )
 9300 format(30e15.6)
c
      end



