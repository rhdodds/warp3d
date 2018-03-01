c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouddpa                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 4/11/2017 (rhd)            *
c     *                                                              *
c     *     write nodal displacements, velocities, accelerations,    *
c     *     reactions, temperatures to (1) patran file in            *
c     *     either binary or formatted forms, or (2) a flat text     *
c     *     on stream (binary) file                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine ouddpa( dva, oubin, ouasc, nodmax, defmax, nwidth,
     &                 flat_file, stream_file, text_file, compressed )
      use global_data ! old common.main
c
      use main_data, only : trn, trnmat, mdiag, pbar, rload,
     &                      inverse_incidences, temper_nodes,
     &                      temper_nodes_ref
      implicit none
c
      integer :: dva, nodmax, nwidth
      logical :: oubin, ouasc, flat_file, stream_file, text_file,
     &           compressed
      double precision :: defmax
c
c                       local declarations
c
      integer :: step_num, bnfile, fmfile, flat_file_number, nod, elem,
     &           dof, ndof, sdof, i
      integer :: titl(80), titl1(80)
      double precision, allocatable :: edva(:,:), trnmte(:,:,:)
      double precision :: nfac
      double precision, parameter :: one = 1.0d0
      real :: sgl_defmax, sgl_vals(10)
      logical :: patran_file
      logical, allocatable :: trne(:,:)

      character(len=4) :: title(80), title1(80)
      character(len=80) :: string, strng1, stepstring*7
      equivalence (title,titl), (title1,titl1)
c
c                       make file name and open file
c                       patran:  binary or formatted files or
c                       both for output.
c
c                       flat file: text or stream. no header lines
c                       in file and no node numbers written
c
c                       for mpi, all nodal results have been collected
c                       to rank 0 -- only 1 patran or flat file
c                       is needed.
c
      allocate( trne(mxvl,mxndel), edva(mxvl,mxndof), ! keep off stack
     &          trnmte(mxvl,mxedof,mxndof) )
c
      step_num = ltmstp
      call ouocdd( dva, step_num, oubin, ouasc, bnfile, fmfile, 1,
     &             .false., myid, flat_file, stream_file, text_file,
     &             compressed, flat_file_number )
c
c                        writing of patran file(s) and flat files
c                        mingled together since same logic is followed
c                        for both.
c
      patran_file = oubin .or. ouasc
c
c                        write info at head of file
c
      if( patran_file ) call ouddpa_patran_header
      if( flat_file .and. text_file )
     &   call ouddpa_flat_header( 1, dva, flat_file_number )
c
c                       loop over the all nodes & write records
c
      do nod = 1, nonode
c
        elem = inverse_incidences(nod)%element_list(1)
        ndof = iprops(4,elem)
c
c                       set the vector of dis/vel/acc/reactions to be
c                       printed. for reactions, it is a bit complicated
c                       since we have to remove internal forces and
c                       inertia effects. nfac is from newmark beta
c                       integration scheme. transform to global from
c                       node local if needed
c
        if( dva .le. 4 ) then
          nfac = one / (nbeta*dt*dt)
          do dof = 1, ndof
            sdof = dstmap(nod)+dof-1
            if( dva .eq. 1 ) then
               edva(1,dof) = u(sdof)
            else if( dva .eq. 2 ) then
               edva(1,dof) = v(sdof)
            else if( dva .eq. 3 ) then
               edva(1,dof) = a(sdof)
            else if( dva .eq. 4 ) then
               edva(1,dof) = pbar(sdof) - ifv(sdof) -
     &                       mdiag(sdof)*du(sdof)*nfac
               if( cstmap(sdof) .ne. 0 ) edva(1,dof) =
     &                       -rload(sdof)+ifv(sdof)+
     &                        mdiag(sdof)*a(sdof)
            end if
          end do
          trne(1,1) = trn(nod)
          if ( trne(1,1) ) then
            trnmte(1,1:3,1:3) = trnmat(nod)%mat(1:3,1:3)
            call trnvec( edva, trnmte, trne, ndof, 1, 1, 2 )
          end if
        end if     !  on dva <= 4
c
        if( dva .eq. 5 ) then   ! temperature only
         edva(1,1) = temper_nodes(nod)
         edva(1,2) = temper_nodes(nod) - temper_nodes_ref(nod)
        end if
c
c                       write the results records for this node.
c                       patran binary results must be real*4. use
c                       real*4 for ASCII to keep compiler from
c                       using a "d" format.
c
c                       flat files: just use double precision
c                       values
c
        sgl_vals(1:nwidth) = edva(1,1:nwidth)
        if( patran_file .and. oubin ) write(bnfile) nod,
     &                                sgl_vals(1:nwidth)
        if( patran_file .and. ouasc ) write(fmfile,920) nod,
     &                                sgl_vals(1:nwidth)
c
        if( flat_file .and. stream_file )
     &    write(flat_file_number) edva(1,1:nwidth)
        if( flat_file .and. text_file )
     &    write(flat_file_number,930) edva(1,1:nwidth)
c
      end do  ! on all nodes
c
c                       close patran or flat file. done.
c
      call ouocdd( dva, ltmstp, oubin, ouasc, bnfile, fmfile, 2,
     &             .false., myid, flat_file, stream_file, text_file,
     &             compressed, flat_file_number )
c
      return
c
 920  format(i8,5(e13.6))
 930  format(3e15.6)
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine ouddpa_patran_header                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/16/2014 (rhd)            *
c     *                                                              *
c     *     Write Patran file header info                            *
c     *                                                              *
c     ****************************************************************

      subroutine ouddpa_patran_header
      implicit none
c
      string = ' '
      strng1 = ' '
c
      if( dva.eq.1 ) then
         string = 'nodal displacement results for structure '
      else if( dva.eq.2) then
         string = 'nodal velocity results for structure     '
      else if(dva.eq.3) then
         string = 'nodal acceleration results for structure '
      else if( dva.eq.4 ) then
         string = 'nodal reactions for structure            '
      else if( dva.eq.5 ) then
         string = 'nodal temperatures for structure         '
      end if
      strng1 = 'loading '
c
c                     step number is ltmstp
c
      string(42:49) = stname
      strng1(9:16)  = lsldnm
      write(stepstring,'(i7)') step_num
      string(50:) = ', step '//stepstring
c
      do i = 1, 80
       title(i)  = string(i:i)
       title1(i) = strng1(i:i)
      end do
c
c                       write initial records to either binary or
c                       formatted file
c
      sgl_defmax = sngl( defmax )
c
      if( patran_file .and. oubin ) then
        write(bnfile) titl, nonode, nonode, sgl_defmax, nodmax, nwidth
        write(bnfile) titl1
        write(bnfile) titl1
      end if
c
c                       write to a formatted file, if necessary.
c
      if( patran_file .and. ouasc ) then
        write(fmfile,900) titl
        write(fmfile,910) nonode, nonode, sgl_defmax, nodmax, nwidth
        write(fmfile,900) titl1
        write(fmfile,900) titl1
      end if
c
      return
c
 900  format(80a1)
 910  format(2i9,e15.6,2i9)
c
c
      end subroutine ouddpa_patran_header
      end subroutine ouddpa

c     ****************************************************************
c     *                                                              *
c     *             subroutine ouddpa_flat_header                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/16/2014 (rhd)            *
c     *                                                              *
c     *     Write header lines for flat text files                   *
c     *                                                              *
c     ****************************************************************

      subroutine ouddpa_flat_header( type, quantity,
     &                               flat_file_number )
      use global_data ! old common.main
      implicit none

      integer :: type, quantity, flat_file_number
c
c                       local declarations
c
      integer :: step_num
      character(len=15) :: vector_types(5), tensor_types(2),
     &                     scalar_types(1)
      character(len=24) :: sdate_time_tmp
c
      data vector_types / 'displacements', 'velocities',
     &                    'accelerations', 'reactions',
     &                    'temperatures' /
      data tensor_types / 'strains', 'stresses' /
      data scalar_types / 'states' /
c
c                       type:
c                         = 1 a nodal quantity of vector type
c                         = 2 a nodal tensor type
c                         = 3 an element tensor type
c                         = 4 scalar type
c                       quantity for type = 1, see entries in
c                         vector_types above
c                       quantity for type = 2,3 see entries in
c                         tensor_types
c                       quantity for type = see scalar_types
c
      step_num = ltmstp  ! from common.main
c
      write(flat_file_number,9000)
      if( type .eq. 1 ) then
        write(flat_file_number,9010) vector_types(quantity)
      end if
      if( type .eq. 2 ) then
        write(flat_file_number,9010) tensor_types(quantity)
      end if
      if( type .eq. 3 ) then
        write(flat_file_number,9012) tensor_types(quantity)
      end if
      if( type .eq. 4 ) then
        write(flat_file_number,9014) scalar_types(quantity)
      end if
      write(flat_file_number,9020) stname
      write(flat_file_number,9030) nonode, noelem
      call fdate( sdate_time_tmp )
      write(flat_file_number,9040) sdate_time_tmp
      write(flat_file_number,9050) step_num
      write(flat_file_number,9000)
c
      return
c
 9000 format('#')
 9010 format('#  WARP3D nodal results: ',a15)
 9012 format('#  WARP3D element results: ',a15)
 9014 format('#  WARP3D element results: ',a15)
 9020 format('#  Structure name: ',a8 )
 9030 format('#  Model nodes, elements: ',2i8)
 9040 format('#  ',a24)
 9050 format('#  Load(time) step: ',i8 )
c
      end
