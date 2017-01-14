c     Checklist

c
c     ***********************************************************************
c     *                                                                     *
c     *     This file contains all subroutines and modules needed to        *
c     *     assemble the linearized equations in parallel on all ranks.     *
c     *                                                                     *
c     *     Created by MCM                                                  *
c     *                                                                     *
c     ***********************************************************************
c
c     ***********************************************************************
c     *                                                                     *
c     *     Module distributed_stiffness_data                               *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *                                                                     *
c     ***********************************************************************
c
      module distributed_stiffness_data
            use local_stiffness_mod
            implicit none
c                       Parameters (user input or otherwise)
            logical :: parallel_assembly_allowed
            logical :: parallel_assembly_used
            logical :: distributed_stiffness_used
            integer :: initial_map_type
            integer :: final_map_type

            integer :: mxepn
            parameter (mxepn = 2)

c                       Actual data
            type(local_stiffness_dist) :: local_k 
c                       Whether we need to recreate the mapping (1) or just
c                       redistribute the new coefs
            integer :: assembly_type
c                       Distributed (full) dof_eqn_map and eqn_node_map
            integer, allocatable :: dist_dof_eqn_map(:),
     &                        dist_eqn_node_map(:)
c                       Local sparsity structure
            integer              :: g_n, g_nnz
            integer              :: my_n, my_nnz
            integer, allocatable :: my_eqns(:), my_eqns_inv(:)
            integer, allocatable :: my_row_ptrs(:),my_col_indices(:)
            double precision, allocatable :: my_coefs(:)
            integer, allocatable :: global_nnz_vec(:,:)
c                       Initial map : just for sparse structure so we
c                       can assemble the global sparse structure and 
c                       load balance.  Get locality only.
            integer, allocatable :: initial_map(:)
c                       Your portion of the global sparsity structure
c                       we need this because we're still using the WARP
c                       communicator
            integer :: g_my_n, g_my_nnz
            integer, allocatable :: g_row_ptrs(:), g_col_indices(:)
            integer, allocatable :: int_my_eqns(:),int_my_eqns_inv(:)
c                       Map to actual distributed stiffness
            integer, allocatable :: final_map(:)
            integer              :: final_nprocs
            integer, allocatable :: final_my_eqns(:), final_eqn_map(:)
            integer, allocatable :: final_my_eqns_inv(:)
            integer, allocatable :: final_eqn_map_inv(:)
c
            double precision, allocatable :: load_vec(:)
c                 This is rather annoying, but we don't have much choice
            integer, allocatable :: old_col_indexes(:)
      end module
c
c     ***********************************************************************
c     *                                                                     *
c     *     Module qsort                                                    *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *     Code outline courtesy of rosettacode.org                        *
c     *                                                                     *
c     ***********************************************************************
c
      module quicksort
      implicit none
      contains
            subroutine qsort_interface(a,lpo)
                  integer :: a(:)
                  integer :: lpo
c
                  integer :: s,i,j
c
                  s = lpo - 1
c                 qsort the partial array
                  call qsort(a(1:s))
c                 remove duplicates
                  j = 1
                  do i=2,s
                        if (a(i) .ne. a(j)) then
                              j=j+1
                              a(j) = a(i)
                        end if
                  end do
                  lpo = j+1

            end subroutine
c
            recursive subroutine qsort(a)
                  integer :: a(:)
                  integer :: split
c
                  if(size(a) > 1) then
                        call partition(a, split)
                        call qsort(a(:split-1))
                        call qsort(a(split:))
                  end if
c 
            end subroutine
c 
            subroutine partition(a, splitter)
                  integer :: a(:)
                  integer :: splitter
                  integer :: left, right, pivot, temp
c 
                  pivot = (a(1) + a(size(a))) / 2 
                  left = 0
                  right = size(a) + 1
c 
                  do while (left < right)
                        right = right - 1
                        do while (a(right) > pivot)
                              right = right-1
                        end do
                        left = left + 1
                        do while (a(left) < pivot)
                              left = left + 1
                        end do
                        if (left < right) then
                              temp = a(left)
                              a(left) = a(right)
                              a(right) = temp
                        end if
                  end do
c 
                  if (left == right) then
                        splitter = left + 1
                  else
                        splitter = left
                  end if
            end subroutine

 
            recursive subroutine qsort_col_coef(a,b)
                  integer :: a(:)
                  double precision :: b(:)
                  integer :: split
c
                  if(size(a) > 1) then
                        call partition_two(a,b, split)
                        call qsort_col_coef(a(:split-1),b(:split-1))
                        call qsort_col_coef(a(split:),b(split:))
                  end if
c 
            end subroutine
c 
            subroutine partition_two(a, b, splitter)
                  integer :: a(:)
                  double precision :: b(:)
                  integer :: splitter
                  integer :: left, right, pivot, temp
                  double precision :: tempd
c 
                  pivot = (a(1) + a(size(a))) / 2 
                  left = 0
                  right = size(a) + 1
c 
                  do while (left < right)
                        right = right - 1
                        do while (a(right) > pivot)
                              right = right-1
                        end do
                        left = left + 1
                        do while (a(left) < pivot)
                              left = left + 1
                        end do
                        if (left < right) then
                              temp = a(left)
                              tempd = b(left)
                              a(left) = a(right)
                              b(left) = b(right)
                              a(right) = temp
                              b(right) = tempd
                        end if
                  end do
c 
                  if (left == right) then
                        splitter = left + 1
                  else
                        splitter = left
                  end if
            end subroutine
      end module quicksort
c
c     ***********************************************************************
c     *                                                                     *
c     *     Module vararray                                                 *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *     Basically a variable length array of variable length arrays     *
c     *     The data structure comes up a few times so I put it here to be  *
c     *     commonly accessible.                                            *
c     *                                                                     *
c     ***********************************************************************
c
      module vararray_mod
            implicit none
c           
            type subarray
                  integer :: length
                  integer, allocatable :: sarray(:)
            end type subarray

            type vararray
                  integer :: length
                  type (subarray), pointer :: array(:)
            end type vararray

            contains

                  subroutine alloc_vararray(this,s)
                        implicit none
                        type (vararray) :: this
                        integer :: s

                        this%length = s
                        allocate(this%array(s))

                  end subroutine

                  subroutine alloc_subarray(this,i,s)
                        implicit none
                        type (vararray) :: this
                        integer :: i,s
                        
                        this%array(i)%length = s
                        allocate(this%array(i)%sarray(s))
                  end subroutine

                  subroutine free_vararray(this)
                        implicit none
                        type(vararray) :: this
                        
                        integer :: i

                        do i=1,this%length
                              if (allocated(this%array(i)%sarray)) then
                                    deallocate(this%array(i)%sarray)
                              end if
                        end do

                        deallocate(this%array)

                  end subroutine

      end module
c
c     ***********************************************************************
c     *                                                                     *
c     *     distribute_from_assembled                                       *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *     This does assembly "the old way."  We have a fully assembled    *
c     *     matrix on the root process and we distribute it to all the      *
c     *     worker ranks.                                                   *
c     *                                                                     *
c     *     At the present, options are only full = true (assemble full,    *
c     *     not LT matrix) and blocking = 'rows' but I leave other options  *
c     *                                                                     *
c     *     Other inputs are the arrays containing the VSR assembled matrix *
c     *                                                                     *
c     ***********************************************************************
c
      subroutine distribute_from_assembled( neq, ncoeff, k_diag, rhs,
     &                  sol_vec, eqn_coeffs, k_pointers, k_indices,
     &                  full,mapping, itype, out)
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
            implicit none
c                 Input
            logical :: full
            character (len=*) :: mapping
c
            integer :: neq, ncoeff, k_pointers, k_indices
            integer :: out, itype
            double precision :: k_diag, rhs, sol_vec, eqn_coeffs
            dimension :: k_diag(*), rhs(*), eqn_coeffs(*),
     &            k_pointers(*), k_indices(*), sol_vec(*)
c                 Some local variables
            real :: wcputime
            external :: wcputime
            integer :: ncoeff_copy,ierr
c
c           Set basic data
            ncoeff_copy = ncoeff+neq
            call setup_comm(local_k,MPI_COMM_WORLD)

            assembly_type = itype
            call MPI_Bcast(assembly_type,1,MPI_INTEGER,0,
     &            local_k%comm,ierr)
c     
c           Check the input (if root)
            if (local_k%my_rank .eq. 0) then
              if (.not. full) then
                write(out,*) "Error: Invalid option for flag full - "
                write(out,*) "must be true, requiring assembly of full"
                write(out,*) "(not LT) stiffness)."
                call MPI_Abort(local_k%comm,-1,ierr)
              end if
c
              if ( mapping .ne. "blockrow" ) then
                write(out,*) "Error: Invalid mapping type -"
                write(out,*) "must be blockrow, indicating"
                write(out,*) "contiguous row blocking"
                call MPI_Abort(local_k%comm,-1,ierr)
              end if
            end if
c           Convert format (if root)
            if (local_k%my_rank .eq. 0) then
              call t_start_assembly(start_assembly_step)
              if (assembly_type .le. 2) then
                call pardiso_symmetric_map( neq, ncoeff, k_diag,
     &                            rhs,eqn_coeffs,
     &                            k_pointers, k_indices )
                call convert_to_full( neq, ncoeff_copy, 
     &                   eqn_coeffs, k_pointers,k_indices)
              end if
              call t_end_assembly(assembly_total,start_assembly_step)
              write(out,
     &      '(15x,"convert to (full) CSR done", "    @ ", f10.2)')
     &            wcputime(1)
            end if

c           Actually distribute
            if (local_k%my_rank .eq. 0) then
                  call t_start_assembly(start_assembly_step)
            end if
            if (assembly_type .eq. 1) then
                  call clear_dist(local_k)
                  call setup_dist(local_k,neq,ncoeff_copy,k_pointers)
                  call distribute_full(local_k,k_pointers,k_indices,
     &                  eqn_coeffs,rhs,sol_vec)
            elseif (assembly_type .eq. 2) then
                  call distribute_coefs(local_k,k_pointers,eqn_coeffs,
     &                  rhs)
            else
                  return
            end if
            if (local_k%my_rank .eq. 0) then
                  call t_end_assembly(assembly_total,
     &                  start_assembly_step)
            end if
            if (local_k%my_rank .eq. 0) then
                  write(out,
     &      '(15x,"distribution finished         @ ",  f10.2)')
     &            wcputime(1)
            end if
            return

      end subroutine
c
c
c     ***********************************************************************
c     *                                                                     *
c     *     determine_local_sparsity                                        *
c     *                                                                     *
c     *     Determine the sparsity structure of your local blocks of        *
c     *     elements                                                        *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *                                                                     *
c     ***********************************************************************
c
      subroutine determine_local_sparsity()
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
            use main_data, only : inverse_incidences, elems_to_blocks
            use quicksort
            use vararray_mod
            implicit integer(a-z)
      include 'common.main'
            integer :: num_enode_dof
            integer :: num_struct_dof
            integer,allocatable :: scratch(:), edest_vec(:)
            integer :: blk, scratch_ptr, felem, span, i,j,eq,ierr
            real :: wcputime
            external :: wcputime
            type(vararray) :: storage_array
            integer :: mineqns
            logical :: ok, gok
c
c           Setup some basic maps (on each proc)
            num_enode_dof = iprops(4,1)
            num_struct_dof = nonode * num_enode_dof
c           
           
            if ( .not. allocated ( dist_dof_eqn_map ) ) then
              allocate( dist_dof_eqn_map(num_struct_dof) )
              allocate( dist_eqn_node_map(num_struct_dof) ) 
            end if
            call dof_map( dist_dof_eqn_map, cstmap, nonode,
     &                  num_enode_dof,dist_eqn_node_map, neqns )
c
c           This is because cstmap doesn't get updated on workers and I
c           don't think I should broadcast the whole thing every time.
            call MPI_Allreduce(neqns, mineqns, 1, MPI_INTEGER,MPI_MAX,
     &            MPI_COMM_WORLD,ierr)
            if (neqns .ne. mineqns) then
                  ok = .false.
            else
                  ok = .true.
            end if
            call MPI_Allreduce(ok,gok,1,MPI_LOGICAL,MPI_LAND,
     &            MPI_COMM_WORLD,ierr)
            if ( .not. gok) then
                  call MPI_Bcast(cstmap,nonode*num_enode_dof,
     &                  MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
                  call dof_map( dist_dof_eqn_map, cstmap, nonode,
     &                  num_enode_dof,dist_eqn_node_map, neqns )
            end if
c           Set basic data
            g_n = neqns
c
c           Allocate the scratch space
c           Note: worst case we could have each column in a row with an entry
c           (g_n) times a number of duplicated entries equal to that times the
c           maximum number of elements connected to a single node.  In practice
c           overallocate by a factor and warn.
            allocate(scratch(mxepn*g_n))
            do i=1,g_n
                  scratch(i) = 0
            end do
            scratch_ptr = 1
c           Loop on equations to determine those we have a part of
c$OMP PARALLEL DO PRIVATE(i,node,numele,j,elem,blk)
            do i=1,g_n
                  node = dist_eqn_node_map(i)
                  numele = inverse_incidences(node)%element_count
                  do j=1,numele
                        elem = inverse_incidences(node)%element_list(j)
                        blk = elems_to_blocks(elem,1)
                        if (elblks(2,blk) .ne. myid) cycle
c                         We established that we own elem elem, therefore
c                         we own a member of eqn i
                        scratch(i)=1
                  end do
            end do
c$OMP END PARALLEL DO
c           Pick up our local n
            my_n = sum(scratch(1:g_n))
c           So the issue here is that I assumed g_n is always the same
c           this is not true.
            if (size(my_eqns) .ne. my_n) then
                  if (allocated(my_eqns)) deallocate(my_eqns)
                  allocate(my_eqns(my_n))
            end if
            if (size(my_row_ptrs) .ne. (my_n+1)) then
                  if (allocated(my_row_ptrs)) deallocate(my_row_ptrs)
                  allocate(my_row_ptrs(my_n+1))
            end if
c           Loop on scratch, inserting into my_eqns
            do i=1,g_n
                  if (scratch(i) .eq. 1) then
                        my_eqns(scratch_ptr) = i
                        scratch_ptr=scratch_ptr+1
                  end if
            end do
c           Allocate our vararray
            call alloc_vararray(storage_array,my_n)

c           Loop on my equations to get row_ptrs and nnz
            my_nnz=0
            my_row_ptrs(1) = 1
            do i=1,my_n
              scratch_ptr = 1
              eqn = my_eqns(i)
              node = dist_eqn_node_map(eqn)
              numele = inverse_incidences(node)%element_count
              do j=1,numele
                elem = inverse_incidences(node)%element_list(j)
                if (elem .le. 0) cycle
                blk = elems_to_blocks(elem,1)
                if (elblks(2,blk) .ne. myid) cycle
                totdof = iprops(2,elem)*iprops(4,elem)
                allocate(edest_vec(totdof))
                call get_single_edest_terms(edest_vec,elem)
                do erow = 1, totdof
                  srow = dist_dof_eqn_map(edest_vec(erow)) 
                  if (srow .ne. eqn) cycle
                  do ecol = 1,totdof
                    scol = dist_dof_eqn_map(edest_vec(ecol))
                    if (scol .eq. 0) cycle
                    scratch(scratch_ptr)=scol
                    scratch_ptr = scratch_ptr+1
                    if (scratch_ptr .gt. mxepn*g_n) then
      write(out,*) ""
      write(out,*) "Error: Increase parameter mxepn found in the module"
      write(out,*) "       local_stiffness_mod in distributed_assembly_"
      write(out,*) "       real in the linux_packages/source directory."
      write(out,*) ""
                        call die_gracefully
                    end if
                  end do
                end do
                deallocate(edest_vec)
              end do
c             Now we need to sort and remove the duplicates
              call qsort_interface(scratch,scratch_ptr)
c             scratch_ptr is row nnz+1
c             scratch contains the sorted col indices for the row
              my_nnz = my_nnz+scratch_ptr-1
              my_row_ptrs(i+1) = my_row_ptrs(i) + scratch_ptr-1

c             Allocate the local subarray
              call alloc_subarray(storage_array,i,scratch_ptr-1)
c             Insert into that subarray
              storage_array%array(i)%sarray(:) = 
     &                  scratch(1:(scratch_ptr-1))
            end do
            deallocate(scratch)
c           Allocate col_indices
            if (size(my_col_indices) .ne. my_nnz) then
                  if (allocated(my_col_indices))
     &                  deallocate(my_col_indices)
                  allocate(my_col_indices(my_nnz))
            end if
c$OMP PARALLEL DO PRIVATE(i,startv,endv)
            do i=1,my_n
                  startv = my_row_ptrs(i)
                  endv   = my_row_ptrs(i+1)-1
                  my_col_indices(startv:endv)=
     &                  storage_array%array(i)%sarray(:)
            end do
c$OMP END PARALLEL DO
c           We now have the local non-zero structure done
            call free_vararray(storage_array)
c           Create the inverse my_eqns map (which will tell us, given a 
c           equation number, whether or not we own the equation and where it is)
            if (size(my_eqns_inv) .ne. g_n) then
                  if (allocated(my_eqns_inv)) deallocate(my_eqns_inv)
                  allocate(my_eqns_inv(g_n))
            end if
            do i=1,g_n
                  my_eqns_inv(i) = 0
            end do
            do i=1,my_n
                  my_eqns_inv(my_eqns(i)) = i
            end do

c           Create a large nnz vec, which we'll use many times hereafter
            if (size(global_nnz_vec,2) .ne. g_n) then
                  if (allocated(global_nnz_vec))
     &                deallocate(global_nnz_vec)
                  allocate(global_nnz_vec(numprocs,g_n))
            end if
c$OMP PARALLEL DO PRIVATE(i,eq)
            do i=1,g_n
                  eq = my_eqns_inv(i)
                  if (eq .eq. 0) then
                        global_nnz_vec(myid+1,i) = 0
                  else
                        global_nnz_vec(myid+1,i) = 
     &                   my_row_ptrs(eq+1)-my_row_ptrs(eq)
                  end if
            end do
c$OMP END PARALLEL DO

            do i=1,numprocs
                  call MPI_Bcast(global_nnz_vec(i,:),g_n,MPI_INTEGER,
     &                  i-1,MPI_COMM_WORLD,ierr)
            end do

            if (myid .eq. 0) then
           write (out,'(15x,"Local sparsity assembled      @ ",f10.2)')
     &            wcputime(1)       
            end if

      end subroutine
c
c     ***********************************************************************
c     *                                                                     *
c     *     determine_initial_map                                           *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *     Determine an initial map of the assembly of the sparse          *
c     *     structure.  Right now, we only care about locality in terms of  *
c     *     nnz that remain local.  Later on (after we have the actual      *
c     *     structure) we will load balance                                 *
c     *                                                                     *
c     *     We allow two options (switch on initial_map_type):              *
c     *                                                                     *
c     *     1 = simple block row - allocate rows to procs based on simple   *
c     *         blocking.  Poor load balance, but fast mapping.             *
c     *                                                                     *
c     *     2 = load balance rows - Assign each row to the processors which *
c     *         "naturally" owns the most entries in that row               *
c     *                                                                     *
c     ***********************************************************************
c
      subroutine determine_initial_map
            use performance_data
            use distributed_stiffness_data
            use local_stiffness_mod
            implicit integer (a-z)
      include 'common.main'
            integer :: i, ierr
            real :: wcputime
            external :: wcputime
c
c                 Match up parameters and set defaults
            call MPI_Bcast(initial_map_type,1,MPI_INTEGER,0,
     &            MPI_COMM_WORLD,ierr)
            if (size(initial_map) .ne. g_n) then
                  if (allocated(initial_map)) deallocate(initial_map)
                  allocate(initial_map(g_n))
            end if
            if (initial_map_type .eq. 2) then
c                 Loop on eqns, making map
c                 As a first guess, assign the row to the proc with the
c                 most nnz in it.
            g_my_n = 0
c$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:G_MY_N)
            do i=1,g_n
                  initial_map(i) = maxloc(global_nnz_vec(:,i),DIM=1)-1
                  if (initial_map(i) .eq. myid) g_my_n=g_my_n+1
            end do
c$OMP END PARALLEL DO

            else if (initial_map_type .eq. 1) then
            g_my_n = 0
            do i=1,g_n
                  initial_map(i) = block_owner(i,numprocs,g_n)
                  if (initial_map(i) .eq. myid) g_my_n=g_my_n+1
            end do

            else
                  if (myid .eq. 0) then
                        write (*,*) "Error: unknown mapping type."
                        call die_gracefully
                  end if

            end if

            if (myid .eq. 0) then
                  write (out,
     &             '(15x,"Initial map created           @ ",f10.2)')
     &             wcputime(1)  
            end if

      end subroutine
c
c     ***********************************************************************
c     *                                                                     *
c     *     assemble_sparsity                                               *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *     Assemble the global sparsity pattern across processors based    *
c     *     on our initial map.                                             *
c     *                                                                     *
c     ***********************************************************************
c
      subroutine assemble_sparsity
            use performance_data
            use distributed_stiffness_data
            use local_stiffness_mod
            use vararray_mod
            use quicksort
            implicit integer (a-z)
      include 'common.main'
c            
            integer, allocatable :: buffer(:)
            type(vararray) :: storage_array
            integer :: sendcount,sendn,s,e,sz,geq,ierr,bptr,p,rptr
            integer :: leq,eq,eq_ptr
            integer, allocatable :: sends(:)
            integer, allocatable :: rcount(:),receives(:)
            real :: wcputime
            external :: wcputime
c
            call alloc_vararray(storage_array,g_my_n)
            if (size(int_my_eqns) .ne. g_my_n) then
                  if (allocated(int_my_eqns)) deallocate(int_my_eqns)
                  allocate(int_my_eqns(g_my_n))
            end if
c           Determine your equations
            eq_ptr=1
            do eq=1,g_n
                  if (initial_map(eq) .eq. myid) then
                        int_my_eqns(eq_ptr) = eq
                        eq_ptr = eq_ptr+1
                  end if
            end do
c           Count our sends
            sendcount = 0
            do eq=1,my_n
                  g_e = my_eqns(eq)
                  if (initial_map(g_e) .ne. myid) then
                        sendcount = sendcount + 1
                  end if
            end do
c           Send nonblocking
            allocate(sends(sendcount))
            sendn = 1
            do eq=1,my_n
                  geq = my_eqns(eq)
                  if (initial_map(geq) .ne. myid) then
                        s=my_row_ptrs(eq)
                        e=my_row_ptrs(eq+1)-1
                        sz=e-s+1
                        call MPI_Issend(my_col_indices(s:e),sz,
     &                   MPI_INTEGER, initial_map(geq), geq,
     &                   MPI_COMM_WORLD,sends(sendn),ierr)
                        sendn = sendn+1
                  end if
            end do
c           Count our receives
            allocate(rcount(g_my_n))
            do i=1,g_my_n
                  geq = int_my_eqns(i)
                  rcount(i) = count(global_nnz_vec(:,geq) .ne. 0)
                  if (global_nnz_vec(initial_map(geq)+1,geq) .ne. 0) 
     &                        then
                        rcount(i) = rcount(i) - 1
                  end if
            end do
c           Receive, sort, and insert
            do eq=1,g_my_n
                  geq = int_my_eqns(eq)
                  leq = my_eqns_inv(geq)
                  szb=sum(global_nnz_vec(:,geq))
                  allocate(buffer(szb))
                  bptr = 1
                  rptr = 1
                  allocate(receives(rcount(eq)))
                  do p=1,numprocs
                        sz = global_nnz_vec(p,geq)
                        if (( sz .ne. 0) .and. 
     &                        ((p-1) .ne. myid)) then
                              call MPI_Irecv(buffer(bptr:(bptr+sz-1)),
     &                              sz, MPI_INTEGER,p-1,geq,
     &                              MPI_COMM_WORLD,receives(rptr),ierr)
                              rptr = rptr+1
                              bptr = bptr+sz
                        end if
                  end do
                  call MPI_Waitall(rcount(eq),receives,
     &                  MPI_STATUSES_IGNORE,ierr)
                  deallocate(receives)
                  if (leq .ne. 0) then
                        s = my_row_ptrs(leq)
                        e = my_row_ptrs(leq+1)-1
                        sz = e-s+1
                        buffer(bptr:(bptr+sz-1))=my_col_indices(s:e)
                        bptr = bptr+sz
                  end if
                  call qsort_interface(buffer,bptr)
                  call alloc_subarray(storage_array,eq,bptr-1)
                  storage_array%array(eq)%sarray(:) =
     &                  buffer(1:(bptr-1))                  

                  deallocate(buffer)
            end do
            deallocate(rcount)
c                 Now we have everything in the vararray, we can extract
c                 our local nnz, and get g_row_ptrs and g_col_indices
            g_my_nnz = 0
            do i=1,g_my_n
                  g_my_nnz = g_my_nnz + storage_array%array(i)%length
            end do
            if (size(g_row_ptrs) .ne. (g_my_n+1)) then
                  if (allocated(g_row_ptrs)) deallocate(g_row_ptrs)
                  allocate(g_row_ptrs(g_my_n+1))
            end if
            if (size(g_col_indices) .ne. g_my_nnz) then
                  if (allocated(g_col_indices))
     &                  deallocate(g_col_indices)
                  allocate(g_col_indices(g_my_nnz))
            end if

           g_row_ptrs(1) = 1
            do i=2,g_my_n+1
                  g_row_ptrs(i) = g_row_ptrs(i-1)+
     &                        storage_array%array(i-1)%length
                  g_col_indices(g_row_ptrs(i-1):(g_row_ptrs(i)-1))=
     &                 storage_array%array(i-1)%sarray(:)
            end do

            call free_vararray(storage_array)

            call MPI_Waitall(sendcount,sends,MPI_STATUSES_IGNORE,ierr)
c                 We can get the global nnz! (That took long enough)
            call MPI_Allreduce(g_my_nnz,g_nnz,1,MPI_INTEGER,MPI_SUM,
     &                  MPI_COMM_WORLD,ierr)

            if (myid .eq. 0) then
                  write (out,
     &             '(15x,"Sparsity structure assembled  @ ",f10.2)')
     &             wcputime(1)  
                  write (out,
     &             '(20x,"n                        =     ", i8)') g_n
                  write (out,
     &             '(20x,"nnz                      =   ", i12)') g_nnz
            end if
c
            return
c           
      end subroutine

c
c     ***********************************************************************
c     *                                                                     *
c     *     determine_ordering                                              *
c     *                                                                     *
c     *     created by: mcm 10/11                                           *
c     *     last modified:                                                  *
c     *                                                                     *
c     *     Determine the final ordering of the sparse matrix.  Relevant    *
c     *     parameter is final_map_type:                                    *
c     *                                                                     *
c     *     1     Same ordering as initial_map_type                         *
c     *     2     Simple block row                                          *
c     *     3     Efficient order for mat-vec product                       *
c     *                                                                     *
c     ***********************************************************************
c
      subroutine determine_ordering
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
            implicit integer (a-z)
      include 'common.main'
c
            integer :: ierr
            real :: wcputime
            external :: wcputime
c
c           Sync parameters
            call MPI_Bcast(final_map_type,1,MPI_INTEGER,0,
     &            MPI_COMM_WORLD,ierr)

c           This will be an issue in the future...
            final_nprocs = numprocs
            if (size(final_map) .ne. g_n) then
                  if (allocated(final_map)) deallocate(final_map)
                  allocate(final_map(g_n))
            end if
c           Switch on the type
            if (final_map_type .eq. 1) then
             final_map(:) = initial_map(:)
            else if (final_map_type .eq. 2) then
             do i=1,g_n
                   final_map(i) = block_owner(i,final_nprocs,g_n)
             end do
            else if (final_map_type .eq. 3) then
             if (myid .eq. 0) then
                  write (*,*) "Error: Feature not supported yet."
                  call die_gracefully
             end if
            else
             if (myid .eq. 0) then
                  write (*,*) "Error: Unknown mapping type."
                  call die_gracefully
             end if
            end if

            if (myid .eq. 0) then
                  write (out,
     &             '(15x,"Final map created             @ ",f10.2)')
     &             wcputime(1)  
            end if
c
      end subroutine
c
c     *******************************************************************
c     *                                                                 *
c     *     move_sparsity                                               *
c     *                                                                 *
c     *     Move sparsity pattern from initial_map to final_map         *
c     *     Note: this will be the first complicated routine when I     *
c     *     resetup rank spawning.  The ranks will need to be spawned   *
c     *     by now.                                                     *
c     *                                                                 *
c     *     Written by: mcm 10/11                                       *
c     *                                                                 *
c     *******************************************************************
c
      subroutine move_sparsity
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
c
            implicit integer (a-z)
      include 'common.main'
c
            integer :: ierr,s,e,sz
            integer,allocatable :: temp_offset(:),temp_recv_cnts(:)

            integer :: sendcount, geq, sendn, receivecount
            integer :: receiven,l_s,l_e
            integer, allocatable :: sends(:),receives(:),sendb(:)
            real :: wcputime
            external :: wcputime
c
c           This will be a pain when we have rank spawning
            local_k%comm=MPI_COMM_WORLD
            local_k%global_n=g_n
            local_k%global_nnz=g_nnz
c
            local_k%num_procs = final_nprocs
c           Obviously this is cheating
            local_k%my_rank = myid
c           We need yet another inverse map...
            if (size(int_my_eqns_inv) .ne. g_n) then
                  if (allocated(int_my_eqns_inv))
     &                  deallocate(int_my_eqns_inv)
                  allocate(int_my_eqns_inv(g_n))
            end if
            do eq=1,g_n
                  int_my_eqns_inv(eq) = 0
            end do
            do eq=1,g_my_n
                  int_my_eqns_inv(int_my_eqns(eq)) = eq
            end do
c
c                 Do the next bit in 2 rounds.  In round one you figure
c                 out your n and nnz, and allocate the local structures
c                 appropriately.  Also note which rows you own (in absolute
c                 reference).  This becomes the map we need to get back to
c                 the correct ordering.
            local_k%local_n=0
            do eq=1,g_n
                  if (final_map(eq) .eq. local_k%my_rank) then
                        local_k%local_n=local_k%local_n+1
                  end if
            end do
c
            local_k%local_nnz=0
            if (size(final_my_eqns) .ne. local_k%local_n) then
                  if (allocated(final_my_eqns))
     &              deallocate(final_my_eqns)
                  allocate(final_my_eqns(local_k%local_n))
            end if
            if (size(local_k%row_ptrs) .ne. (local_k%local_n+1)) then
                  if (allocated(local_k%row_ptrs))
     &                  deallocate(local_k%row_ptrs)
                  allocate(local_k%row_ptrs(local_k%local_n+1))
            end if
c
c           Count the number of sends, and non-blocking send the nnzs
            sendcount = 0
            do eq=1,g_my_n
                  geq = int_my_eqns(eq)
                  if (final_map(geq) .ne. local_k%my_rank) then
                        sendcount = sendcount + 1
                  end if
            end do
            allocate(sends(sendcount))
            allocate(sendb(sendcount))
            sendn=1
            do eq=1,g_my_n
                  geq = int_my_eqns(eq)
                  if (final_map(geq) .ne. local_k%my_rank) then
                        sendb(sendn) = g_row_ptrs(eq+1)-g_row_ptrs(eq)
                        call MPI_Issend(sendb(sendn),1,MPI_INTEGER,
     &                       final_map(geq),geq,local_k%comm,
     &                       sends(sendn),ierr)
                        sendn = sendn+1
                  end if
            end do
c            
c           Count the number of receives and make your equation list,
c           receive or copy all your data
            eqn_ptr = 1
            receivecount = 0
            do eq=1,g_n
                  if (final_map(eq) .eq. local_k%my_rank) then
                        final_my_eqns(eqn_ptr) = eq
                        if (initial_map(eq) .ne. myid) then
                              receivecount = receivecount+1
                        end if
                        eqn_ptr = eqn_ptr + 1
                  end if
            end do
c           
            allocate(receives(receivecount))
            receiven=1
            do eq=1,local_k%local_n
                  geq = final_my_eqns(eq)
                  if (initial_map(geq) .ne. myid) then
                        call MPI_Irecv(local_k%row_ptrs(eq+1),1,
     &                   MPI_INTEGER,initial_map(geq),geq,
     &                   local_k%comm,receives(receiven),ierr)
                        receiven = receiven + 1
                  else
                        ieq = int_my_eqns_inv(geq)
                        local_k%row_ptrs(eq+1) = g_row_ptrs(ieq+1)-
     &                        g_row_ptrs(ieq)
                  end if
            end do
c           Wait for the data
            call MPI_Waitall(receivecount,receives,MPI_STATUSES_IGNORE,
     &            ierr)
c            
c           Now put 1 in the first spot in your local_k%row_ptrs and add the
c           nnzs to it.  Also sum our nnz
            local_k%local_nnz=0
            local_k%row_ptrs(1) = 1
            do eq=2,local_k%local_n+1
                  local_k%local_nnz=local_k%local_nnz+
     &                  local_k%row_ptrs(eq)
                  local_k%row_ptrs(eq) = local_k%row_ptrs(eq)+
     &                  local_k%row_ptrs(eq-1)
            end do
c           Wait for the sends
            call MPI_Waitall(sendcount,sends,MPI_STATUSES_IGNORE,ierr)
                  deallocate(sendb)
c
c                 I can now allocate the col_indexes
            if (size(local_k%col_indexes) .ne. local_k%local_nnz) then
                  if (allocated(local_k%col_indexes))
     &                  deallocate(local_k%col_indexes)
                  allocate(local_k%col_indexes(local_k%local_nnz))
            end if
c            
c           Send the col indexes
            sendn=1
            do eq=1,g_my_n
                  geq = int_my_eqns(eq)
                  if (final_map(geq) .ne. local_k%my_rank) then
                        s=g_row_ptrs(eq)
                        e=g_row_ptrs(eq+1)-1
                        sz=e-s+1
                        call MPI_Issend(g_col_indices(s:e),sz,
     &                       MPI_INTEGER,
     &                       final_map(geq),geq,local_k%comm,
     &                       sends(sendn),ierr)
                        sendn = sendn+1
                  end if
            end do
c            
c           Recv or copy your col indexes
            receiven=1
            do eq=1,local_k%local_n
                  geq = final_my_eqns(eq)
                  s = local_k%row_ptrs(eq)
                  e = local_k%row_ptrs(eq+1)-1
                  sz = e-s+1
                  if (initial_map(geq) .ne. myid) then
                        call MPI_Irecv(local_k%col_indexes(s:e),sz,
     &                   MPI_INTEGER,initial_map(geq),geq,
     &                   local_k%comm,receives(receiven),ierr)
                        receiven = receiven + 1
                  else
                        ieq = int_my_eqns_inv(geq)
                        l_s = g_row_ptrs(ieq)
                        l_e = g_row_ptrs(ieq+1)-1
                        local_k%col_indexes(s:e)=g_col_indices(l_s:l_e)
                  end if
            end do
c           Wait for the data
            call MPI_Waitall(receivecount,receives,MPI_STATUSES_IGNORE,
     &            ierr)
            deallocate(receives)
c            
c           Wait for the sends
            call MPI_Waitall(sendcount,sends,MPI_STATUSES_IGNORE,ierr)
            deallocate(sends)
c           
c           All done
c
c                 Assemble the global map which takes new eqn numbers to
c                 old eqn numbers
            allocate(temp_recv_cnts(local_k%num_procs))
            allocate(temp_offset(local_k%num_procs))
            if (size(final_eqn_map) .ne. local_k%global_n) then
                  if (allocated(final_eqn_map))
     &                  deallocate(final_eqn_map)
                  allocate(final_eqn_map(local_k%global_n))
            end if
            call MPI_Allgather(local_k%local_n,1,MPI_INTEGER,
     &            temp_recv_cnts,1,MPI_INTEGER,local_k%comm,ierr) 
            temp_offset(1) = 0
            do i=2,local_k%num_procs
                  temp_offset(i)=temp_offset(i-1)+temp_recv_cnts(i-1)
            end do
            call MPI_Allgatherv(final_my_eqns,local_k%local_n,
     &            MPI_INTEGER,final_eqn_map,temp_recv_cnts,
     &            temp_offset,MPI_INTEGER,local_k%comm,ierr)
c
            local_k%local_start = temp_offset(local_k%my_rank+1)+1
            deallocate(temp_offset)
            deallocate(temp_recv_cnts)
c
            if (myid .eq. 0) then
                  write (out,
     &             '(15x,"Sparsity moved                @ ",f10.2)')
     &             wcputime(1)  
            end if
c
      end subroutine
c
c     *******************************************************************
c     *                                                                 *
c     *     assemble_coefs                                              *
c     *                                                                 *
c     *     Assemble the coefficients into our sparse matrix (can go    *
c     *     directly from local part to final map now)                  *
c     *                                                                 *
c     *     Written by: mcm 10/11                                       *
c     *                                                                 *
c     *******************************************************************
c
      subroutine assemble_coefs(new_size_in)
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
            use elem_block_data, only: estiff_blocks
            use main_data, only : inverse_incidences,elems_to_blocks
            implicit integer (a-z)
      include 'common.main'
            integer :: ierr,eqn,node,numele,relcol
            integer :: elem,blk,totdof,erow,srow,scol,k
            integer, allocatable :: edest_vec(:)
            integer :: linind
            double precision, dimension(:,:), pointer :: emat
            double precision :: coef
            integer :: l_eq,f_eq
            integer :: t
            real :: wcputime
            external :: wcputime
            logical :: new_size_in, new_size

            integer :: sendcount,eq,geq,s,e,sz,leq,bufsz,recvc,bptr
            integer :: recvn,p,i,j
            integer, allocatable :: colsends(:),coefsends(:),colrecvs(:)
            integer, allocatable :: coefrecvs(:), colbuffer(:)
            double precision, allocatable :: coefbuffer(:)

c                 VERY IMPORTANT
c                 If we are re-using an old sparsity pattern, we need to
c                 revert to the old column indexes
            if (local_k%my_rank .eq. 0) new_size = new_size_in
            call MPI_Bcast(new_size,1,MPI_LOGICAL,0,local_k%comm,ierr)
            if ( .not. new_size) then
                  local_k%col_indexes = old_col_indexes
            end if

c           First assemble the local part...
            if (.not. allocated(my_coefs)) then
                  allocate(my_coefs(my_nnz))
            end if
            do i=1,my_nnz
                  my_coefs(i) = 0
            end do
c           Loop on my equations to get row_ptrs and nnz
            do i=1,my_n
              eqn = my_eqns(i)
              node = dist_eqn_node_map(eqn)
              numele = inverse_incidences(node)%element_count
              do j=1,numele
                elem = inverse_incidences(node)%element_list(j)
                if (elem .le. 0) cycle
                blk = elems_to_blocks(elem,1)
                if (elblks(2,blk) .ne. myid) cycle
                relcol = elems_to_blocks(elem,2)
                totdof = iprops(2,elem)*iprops(4,elem)
                emat => estiff_blocks(blk)%ptr
                allocate(edest_vec(totdof))
                call get_single_edest_terms(edest_vec,elem)
                do erow = 1, totdof
                  srow = dist_dof_eqn_map(edest_vec(erow))    
                  if (srow .ne. eqn) cycle
                  do ecol = 1,totdof
                    scol = dist_dof_eqn_map(edest_vec(ecol))
                    if (scol .eq. 0) cycle
c                       We are on row i (locally), column scol
c                       This corresponds to erow,ecol
c                       Need to convert erow, ecol to the linear
c                       indexing, find our position in coefs
c                       and insert  
c
                        linind = dcp(max(ecol,erow))-abs(ecol-erow)
                        coef = emat(linind,relcol)
c                       Now find the correct place and insert
                        do k=my_row_ptrs(i),(my_row_ptrs(i+1)-1)
                              if (scol .eq. my_col_indices(k)) exit
                        end do
                        my_coefs(k) = my_coefs(k)+coef
                  end do
                end do
                deallocate(edest_vec)
              end do
            end do
c
c           Allocate and set up maps
c
c            if (size(final_my_eqns_inv) .ne. g_n) then
c                  if (allocated(final_my_eqns_inv))
c     &                  deallocate(final_my_eqns_inv)
c                  allocate(final_my_eqns_inv(g_n))
c            end if
            allocate(final_my_eqns_inv(g_n))
            do eq=1,g_n
                  final_my_eqns_inv(eq) = 0
            end do
            do eq=1,local_k%local_n
                  final_my_eqns_inv(final_my_eqns(eq)) = eq
            end do

            if ( size(local_k%coefs) .ne. local_k%local_nnz) then
                  if (allocated(local_k%coefs)) 
     &                  deallocate(local_k%coefs)
                  allocate(local_k%coefs(local_k%local_nnz))
            end if
            do i=1,local_k%local_nnz
                  local_k%coefs(i)=0.0
            end do

c           Now we need to get the coefs to their proper owner
 
c           Count the number of sends
            sendcount = 0
            do eq=1,my_n
                  geq = my_eqns(eq)
                  if (final_map(geq) .ne. local_k%my_rank) then
                        sendcount = sendcount + 1
                  end if
            end do
            allocate(colsends(sendcount))
            allocate(coefsends(sendcount))
c           Send the appropriate col_indexes and coefs
            sendn = 1
            do eq=1,my_n
                  geq = my_eqns(eq)
                  if (final_map(geq) .ne. local_k%my_rank) then
                        s = my_row_ptrs(eq)
                        e = my_row_ptrs(eq+1)-1
                        sz = e-s+1
                        call MPI_Issend(my_col_indices(s:e),sz,
     &                        MPI_INTEGER,final_map(geq),geq,
     &                        local_k%comm,colsends(sendn),ierr)
                        call MPI_Issend(my_coefs(s:e),sz,
     &                        MPI_DOUBLE_PRECISION,final_map(geq),
     &                        geq+g_n,local_k%comm,coefsends(sendn),
     &                        ierr)
                        sendn = sendn+1
                  end if
            end do
c           Loop on your final equations, receiving and processing as appropriate
            do eq=1,local_k%local_n
                  geq = final_my_eqns(eq)
                  leq = my_eqns_inv(geq)
                  bufsz = sum(global_nnz_vec(:,geq))
                  recvc = count((global_nnz_vec(:,geq) .ne. 0), DIM=1)
                  if (leq .ne. 0) recvc=recvc-1
                  allocate(colbuffer(bufsz))
                  allocate(coefbuffer(bufsz))
                  allocate(colrecvs(recvc))
                  allocate(coefrecvs(recvc))
                  bufptr = 1
                  if (leq .ne. 0) then
                        s=my_row_ptrs(leq)
                        e=my_row_ptrs(leq+1)-1
                        sz = e-s+1
                        colbuffer(bufptr:(bufptr+sz-1))=
     &                   my_col_indices(s:e)
                        coefbuffer(bufptr:(bufptr+sz-1))=
     &                   my_coefs(s:e)
                        bufptr = bufptr+sz
                  end if
                  recvn=1
                  do p=0,local_k%num_procs-1
                   if ((p .ne. myid) .and. 
     &                        (global_nnz_vec(p+1,geq) .ne. 0)) then
                     sz=global_nnz_vec(p+1,geq)
                     call MPI_Irecv(colbuffer(bufptr:(bufptr+sz-1)),
     &                  sz, MPI_INTEGER, p, geq, local_k%comm,
     &                  colrecvs(recvn),ierr)
                     call MPI_Irecv(coefbuffer(bufptr:(bufptr+sz-1)),
     &                  sz,MPI_DOUBLE_PRECISION,p,geq+g_n,local_k%comm,
     &                  coefrecvs(recvn),ierr)
                     recvn=recvn+1
                     bufptr = bufptr+sz
                   end if
                  end do
                  call MPI_Waitall(recvc,colrecvs,MPI_STATUSES_IGNORE,
     &                  ierr)
                  call MPI_Waitall(recvc,coefrecvs,MPI_STATUSES_IGNORE,
     &                  ierr)
                  deallocate(colrecvs)
                  deallocate(coefrecvs)
                  do i=1,bufsz
                    do j=local_k%row_ptrs(eq),
     &                 (local_k%row_ptrs(eq+1)-1)
                      if (local_k%col_indexes(j) .eq.
     &                    colbuffer(i)) then
                         local_k%coefs(j)=local_k%coefs(j)+
     &                       coefbuffer(i)
                         exit
                      end if
                     end do
                  end do
                  deallocate(colbuffer)
                  deallocate(coefbuffer)
            end do
            
c           Wait for your sends to be done
            call MPI_Waitall(sendcount,colsends,MPI_STATUSES_IGNORE,
     &            ierr)
            call MPI_Waitall(sendcount,coefsends,MPI_STATUSES_IGNORE,
     &            ierr)
            deallocate(colsends)
            deallocate(coefsends)

c           Well, going to have to realloc anyway, may as well free the memory
            deallocate(my_coefs)

            if (myid .eq. 0) then
                  write (out,
     &             '(15x,"Coefficients distributed      @ ",f10.2)')
     &             wcputime(1)  
            end if

            return

      end subroutine
c
c     *******************************************************************
c     *                                                                 *
c     *     assem_load_vec                                              *
c     *                                                                 *
c     *     Assemble the distributed load vector                        *
c     *                                                                 *
c     *     written by: mcm 10/11                                       *
c     *                                                                 *
c     *******************************************************************
c
      subroutine assem_load_vec
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
            use elem_block_data, only: edest_blocks
            implicit integer (a-z)
      include 'common.main'
c
            integer :: blk,felem,totdof,span
            integer :: ierr,eq,meqn,i
            real :: wcputime
            external :: wcputime

c           NOTE THIS PART IS NOT DISTRIBUTED
            if (myid .eq. 0) then
c           We'll just do this in full on each processor
c            if (size(load_vec) .ne. g_n) then
c                  if (allocated(load_vec)) deallocate(load_vec)
c                  allocate(load_vec(g_n))
c            end if
            if (.not. allocated(load_vec)) then
                  allocate(load_vec(g_n))
            end if
            do i=1,g_n
                  load_vec(i) = 0.0
            end do
c           Assemble your part
            do blk=1,nelblk
                  felem = elblks(1,blk)
                  totdof = iprops(2,felem) * iprops(4,felem)
                  span = elblks(0,blk)
                  call load_vector(load_vec,res,
     &                  edest_blocks(blk)%ptr(1,1),dist_dof_eqn_map,
     &                  totdof,span)
            end do
            end if
c           Get it to the correct processors
            if (size(local_k%rhs) .ne. local_k%local_n) then
                  if (allocated(local_k%rhs))
     &                  deallocate(local_k%rhs)
                  allocate(local_k%rhs(local_k%local_n))
            end if
c            if (.not. allocated(local_k%rhs)) then
c                  allocate(local_k%rhs(local_k%local_n))
c            end if
            do i=1,local_k%local_n
                  local_k%rhs(i) = 0.0
            end do

            do eq=1,local_k%global_n
                  meqn=final_my_eqns_inv(eq)
                  if (myid .eq. 0) then
                        if (meqn .ne. 0) then
                              local_k%rhs(meqn) = load_vec(eq)
                        else
                              call MPI_Send(load_vec(eq),1,
     &                         MPI_DOUBLE_PRECISION,final_map(eq),
     &                         1,local_k%comm,ierr)
                        end if
                  else
                        if (meqn .ne. 0) then
                              call MPI_Recv(local_k%rhs(meqn),1,
     &                         MPI_DOUBLE_PRECISION,0,1,
     &                         local_k%comm,mpi_status,ierr)
                        end if
                  end if
            end do
            if (myid .eq. 0) then
                  deallocate(load_vec)
                  write (out,
     &             '(15x,"Load vector distributed       @ ",f10.2)')
     &             wcputime(1) 
            end if

            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            
            deallocate(final_my_eqns_inv)

                        

            return

      end subroutine
c
c     *******************************************************************
c     *                                                                 *
c     *     dist_final_setup                                            *
c     *                                                                 *
c     *     Perform some final setup on struts.  Most importantly, we   *
c     *     may no longer have a symmetric matrix, so permute our col   *
c     *     and coef vectors to get a symmetric matrix                  *
c     *                                                                 *
c     *     written by: mcm 11/11                                       *
c     *                                                                 *
c     *******************************************************************
c
      subroutine dist_final_setup
            use local_stiffness_mod
            use distributed_stiffness_data
            use performance_data
            use quicksort
            implicit integer (a-z)
      include 'common.main'
            integer :: i,j,sz,ptr,s,e,k
            integer, allocatable :: temp_cols(:)
            double precision, allocatable :: temp_coefs(:)
            integer :: pos
            real :: wcputime
            external :: wcputime
c
c                 First back up old col_indexes so we can restore later
            if (size(old_col_indexes) .ne. local_k%local_nnz) then
                  if (allocated(old_col_indexes))
     &                  deallocate(old_col_indexes)
                  allocate(old_col_indexes(local_k%local_nnz))
            end if
            old_col_indexes = local_k%col_indexes
c
            if (size(final_eqn_map_inv) .ne. g_n) then
                  if (allocated(final_eqn_map_inv)) 
     &                  deallocate(final_eqn_map_inv)
                  allocate(final_eqn_map_inv(g_n))
            end if
            do i=1,g_n
                  final_eqn_map_inv(final_eqn_map(i)) = i
            end do

c           Loop on your local equations, fixing up the pointers
c$OMP PARALLEL DO PRIVATE(i,s,e,sz,temp_cols,temp_coefs,ptr,j)
            do i=1,local_k%local_n
                  s = local_k%row_ptrs(i)
                  e = local_k%row_ptrs(i+1)-1
                  sz = local_k%row_ptrs(i+1)-local_k%row_ptrs(i)
                  allocate(temp_cols(sz))
                  allocate(temp_coefs(sz))
                  ptr = 1
                  do j=s,e
                        temp_cols(ptr) = final_eqn_map_inv(
     &                        local_k%col_indexes(j))
                        temp_coefs(ptr) = local_k%coefs(j)
                        ptr = ptr+1
                  end do
c                 Sort temp_cols into order and match local_k%coefs(s:e)
c                 to the new order.  Copy over.
                  call qsort_col_coef(temp_cols,temp_coefs)
                  local_k%col_indexes(s:e) = temp_cols
                  local_k%coefs(s:e) = temp_coefs
                  deallocate(temp_cols)
                  deallocate(temp_coefs)
            end do
c$OMP END PARALLEL DO

c           Probably should pass in and check solver number
            if (size(local_k%cols_per_row) .ne. local_k%local_n) then
                  if (allocated(local_k%cols_per_row))
     &                  deallocate(local_k%cols_per_row)
                  allocate(local_k%cols_per_row(local_k%local_n))
            end if
            call get_cols_per_row(local_k)
            if (size(local_k%vec_indexes) .ne. local_k%local_n) then
                  if (allocated(local_k%vec_indexes)) 
     &                  deallocate(local_k%vec_indexes)
                  allocate(local_k%vec_indexes(local_k%local_n))
            end if
            call get_vec_indexes(local_k)
            if (size(local_k%soln) .ne. local_k%local_n) then
                  if (allocated(local_k%soln))
     &                  deallocate(local_k%soln)
                  allocate(local_k%soln(local_k%local_n))
            end if
            call get_guess(local_k)
c            call stiffness_check_sums(local_k)
c
            if (myid .eq. 0) then
                  write (out,
     &             '(15x,"Final setup completed         @ ",f10.2)')
     &             wcputime(1) 
            end if

      end subroutine
c
c     *******************************************************************
c     *                                                                 *
c     *     curr_neqns                                                  *
c     *                                                                 *
c     *     Function to return neqns, needed for root to alloc structs  *
c     *                                                                 *
c     *     written by: mcm 11/11                                       *
c     *                                                                 *
c     *******************************************************************
c
      integer function curr_neqns
            use local_stiffness_mod
            use distributed_stiffness_data
            implicit none

            curr_neqns = local_k%global_n

            return

      end function
c
c     *******************************************************************
c     *                                                                 *
c     *     reorder_soln_vec                                            *
c     *                                                                 *
c     *     Puts our solution vector back in the order WARP expect it   *
c     *     to be in.                                                   *
c     *                                                                 *
c     *     written by: mcm 11/11                                       *
c     *                                                                 *
c     *******************************************************************
c
      subroutine reorder_soln_vec(u_vec)
            use distributed_stiffness_data
            implicit none

            double precision :: u_vec(*)

            integer :: i
            double precision, allocatable :: temp_vec(:)

            allocate(temp_vec(g_n))
            do i=1,g_n
                  temp_vec(final_eqn_map(i)) = u_vec(i)
            end do

            u_vec(1:g_n) = temp_vec(1:g_n)

            deallocate(temp_vec)

            return

      end subroutine
c
c     *******************************************************************
c     *                                                                 *
c     *     convert_to_full                                             *
c     *                                                                 *
c     *     Converts our symmetric (LT) CSR matrices to full CSR format *
c     *                                                                 *
c     *     written by: mcm 2/11                                        *
c     *                                                                 *
c     *******************************************************************
c
      subroutine convert_to_full( neq, ncoeff, 
     &       eqn_coeffs, k_pointers,k_indices)
      implicit none    
c           Input
      INTEGER :: neq, ncoeff, k_pointers, k_indices
      double precision ::  eqn_coeffs
      dimension eqn_coeffs(*), k_pointers(*),
     &       k_indices(*)
c           Local variables
      integer :: new_nnz, ierr, i, ltsum, j, offsetv, temp
      integer, allocatable :: work(:)
c
c    The new number of non zeros, this space should already be allocated
      new_nnz = ncoeff * 2 - neq
c    Allocate for temp storage
      allocate(work(neq+1))
c    Zero it out 
      do i=1,neq+1
            work(i) = 0
      end do
c
c    Isn't algorithm design fun?
c    Step 1 - transverse col_indices figuring out how many LT element there will be
      do i=1,ncoeff
            work(k_indices(i))=work(k_indices(i))+1
      end do
c
c    Please remember that the real nnz in each row LT is actually -1 this (diagonal)
      do i=1,neq
            work(i) = work(i)-1
      end do
c
c    step 2 - sum up number of LT terms
      ltsum = sum(work)
      offsetv = ltsum
c
c    Step 3 - shift elements in indices and coefs to correct position
c    Do this backwards as to avoid overwrite
      do i=neq,1,-1
c
c          Do the actual shift (backwards again...)
            do j=k_pointers(i+1)-1,k_pointers(i),-1
                  k_indices(j+offsetv) = k_indices(j)
                  eqn_coeffs(j+offsetv) = eqn_coeffs(j)
            end do
c
c          Adjust the shift amount
            offsetv = offsetv - work(i)
c
      end do
c
c    Step 4 - adjust row_ptrs to the new correct amounts (each row gets shifted by nnz last row)
      offsetv = ltsum
      do i=neq+1,2,-1
            k_pointers(i) = k_pointers(i) + offsetv
            offsetv = offsetv-work(i-1)
      end do
c
c    Step 5 - now transverse our entire UT minus diagonal and insert the transpose of each element
c    Rows are rows (i.e. don't change)
      do i=neq,1,-1
c          Our column offsets are now wrong.  we want "new" offset+LT offset
            temp = k_pointers(i)+work(i)
            do j=temp,k_pointers(i+1)-1
c                Now we can index into column and coefs
c                Get the column from i and the row from column_indexes(j)
c                Note that we are looping on rows, so we know that if we insert
c                into a row it will be in column order
c
c                Rem: row is k_indices(j)
                  k_indices(k_pointers(k_indices(j))+work(
     &                  k_indices(j)))=i
                  eqn_coeffs(k_pointers(k_indices(j))+work(k
     &                  _indices(j)))= eqn_coeffs(j)
                  work(k_indices(j)) = work(k_indices(j))-1
            end do
      end do
c    Clear scratch
      deallocate(work)
c    Change the value of nnz
      ncoeff = new_nnz
      return
      end
      
c     ****************************************************************
c     *                                                              *
c     *  build effective load vector for solution of global          *
c     *  constrained equations                                       *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 03/8/94                    *
c     *                                                              *
c     ****************************************************************
c
c
c     include constraints into loading vector by global loading vector
c     (dofs) being passed to constrained data structure (eqns)
c
      subroutine load_vector( p, res, edst_block, dof_eqn_map, nedof,
     &                        span )
      implicit integer (a-z)
      include 'param_def'
      dimension edst_block(nedof,*), dof_eqn_map(*)
      double precision
     &   res(*), p(*)
c
c                 locally allocated
c
      logical local_debug
      data local_debug / .false. /
c
      if ( local_debug ) write (*,*) 'calc p ',nedof,span
      do i = 1, nedof
        do ispan = 1, span
           if (dof_eqn_map(edst_block(i,ispan)) .ne. 0 )
     &       p(dof_eqn_map(edst_block(i,ispan)))=
     &                        res(edst_block(i,ispan))
        end do
      end do
c
      return
      end
      
