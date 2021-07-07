      implicit none
c
      double precision, allocatable :: x_coord(:), j_values(:)
      double precision :: zero, crack_extension, jv, angle
      integer, allocatable :: step_list(:), master_node_list(:)
      integer :: i, step, snode, nvalues_growth, nsteps,
     &           num_comments, ifile, term
      character(len=30) :: scratch_string
      logical :: debug

      term = 6 
      write(term,*) '... startup...'
      zero = 0.0d0
      debug = .false.
      ifile = 10
c
      open(unit=ifile,file='./j_values.inp')
      write(term,*) '... J_values file opened...' 
      read(ifile,*) num_comments
      write(term,*) '    ... num_comments: ', num_comments
      do i = 1, num_comments
          read(ifile,fmt="(a10)") scratch_string
      end do
c
      read(ifile,*) nsteps
      allocate( j_values(nsteps) )
      do i = 1, nsteps
        read(ifile,*) j_values(i)
      end do
      close(unit=ifile)
      write(term,*) '... J-values read...'
c
      open(unit=ifile,file='./crack_lengths.inp')
c
      read(ifile,*) num_comments
      write(term,*) '    ... num_comments: ', num_comments
      do i = 1, num_comments
          read(ifile,fmt="(a10)") scratch_string
      end do
c
      read(ifile,*) nvalues_growth
      allocate( step_list(nvalues_growth), 
     &          master_node_list(nvalues_growth),
     &          x_coord(nvalues_growth) )
      do i = 1, nvalues_growth
        read(ifile,*)  master_node_list(i), angle, step_list(i),
     &                 x_coord(i)
      end do
      close(unit=ifile)
      write(term,*) '... front data read...'
c

      do i = 1, nvalues_growth
        step = step_list(i)
        snode = master_node_list(i)
        crack_extension = x_coord(i) 
        jv = j_values(step)
        write(term,9000) i, step, snode, crack_extension, jv
      end do



      
 


      stop
 9000 format(i6, i6, i6, f10.4, f10.4)
      end
