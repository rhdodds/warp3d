c
c     ***********************************************************************
c     *                                                                     *
c     *     Module local_stiffness_mod                                      *
c     *                                                                     *
c     *     created by: mcm 2/11                                            *
c     *     last modified: mcm 5/11                                         *
c     *                                                                     *
c     *     Contains a data structure for holding the local part of         *
c     *     a row-block distributed matrix-vector system and functions      *
c     *     to take a fully  assembled stiffness + vector system            *
c     *     on 1 processor and distribute it (by contiguous row blocks)     *
c     *     to all the processors on a communicator                         *
c     *                                                                     *
c     *     Notes: Will fail if a diagonal entry is not represented in the  *
c     *     sparsity structure.  Of course this should never happen         *
c     *     in WARP (b/c we store the full diagonal) but it is something    *
c     *     to keep in mind if you want to reuse the module.                *
c     *                                                                     *
c     ***********************************************************************
c     
      module local_stiffness_mod
            implicit none
            include 'mpif.h'
c ***********************************************************************
c           Data structure for the local rows
c ***********************************************************************
            type local_stiffness_dist
c                Communication/global info
                  integer :: global_n
                  integer :: global_nnz
                  integer :: my_rank
                  integer :: num_procs
                  integer :: comm
c                 Local part of stiffness
                  integer :: local_n
                  integer :: local_nnz
                  integer :: local_start
                  integer, allocatable :: row_ptrs(:) 
                  integer, allocatable :: col_indexes(:)
                  double precision, allocatable :: coefs(:)
                  double precision, allocatable :: rhs(:)
                  double precision, allocatable :: soln(:)
c                 Data for hypre
                  integer, allocatable :: cols_per_row(:)
                  integer, allocatable :: vec_indexes(:) 
            end type local_stiffness_dist
c ***********************************************************************
c
      contains
c                 Various helper functions
c                 For whatever wacky compiler reason, functions must come
c                 before subroutines.
c
c                 A helper function to determine the block_low of any proc
            integer function block_low(p,procs,n)
                  implicit none
                  integer :: p,procs,n
                  block_low = (p*n/procs)+1
                  return
            end function block_low
c
c                A helper to determine the last row a proc owns
            integer function block_high(p,procs,n)
                  implicit none
                  integer :: p, procs,n
                  block_high = block_low(p+1,procs,n)-1
                  return
            end function block_high
c
c                A helper to determine the size of a block
            integer function block_size(p,procs,n)
                  implicit none
                  integer :: p,procs,n
                  block_size = block_high(p,procs,n)-
     &                  block_low(p,procs,n)+1
                  return
            end function block_size
c
c                A helper to determine which proc owns a row
            integer function block_owner(row,procs,n)
                  implicit none
                  integer :: row,procs,n
                  block_owner = (procs*row-1)/n
                  return
            end function block_owner
c
c                 Now the actual subroutines.
c
c                 Subroutine to reassemble complete solution vector
c                 from partials.
            subroutine merge_soln(this, soln)
                  implicit none
c                 Input
                  TYPE(local_stiffness_dist) :: this
                  double precision :: soln(*)
c                 local
                  integer, allocatable :: recv_counts(:), recv_ptrs(:)
                  integer :: ierr,i,mpi_status,fr
                  double precision, allocatable :: test(:)
c                       Allocate and setup our recv arrays for the MPI
c                       gather (see documentation)
                  allocate(recv_counts(this%num_procs))
                  allocate(recv_ptrs(this%num_procs))
                  call MPI_Gather(this%local_n,1,MPI_INTEGER,
     &             recv_counts,1,MPI_INTEGER,0,this%comm,ierr)
                  recv_ptrs(1) = 0
                  do i=2,this%num_procs
                        recv_ptrs(i) = recv_ptrs(i-1)+recv_counts(i-1)
                  end do
c                       Gather and free the recv arrays
                  call MPI_Gatherv(this%soln,this%local_n,
     &                  MPI_DOUBLE_PRECISION,soln, recv_counts, 
     &                  recv_ptrs, MPI_DOUBLE_PRECISION,0, 
     &                  this%comm, ierr)
                  deallocate(recv_counts)
                  deallocate(recv_ptrs)
            end subroutine merge_soln
c
c                 Subroutine to count the number of columns per row (need
c                 for hypre).
            subroutine get_cols_per_row(this)
                  implicit none
c                       Input
                  Type(local_stiffness_dist) :: this
c                       Local
                  integer :: row
c                       Each row has row_ptrs(i+1)-row_ptrs(i) cols
                  do row=1,this%local_n
                        this%cols_per_row(row) = this%row_ptrs(row+1)
     &                        - this%row_ptrs(row)
                  end do
            end subroutine get_cols_per_row
c
c                 Get a vector indicating the global vector indices
            subroutine get_vec_indexes(this)
                  implicit none
c                       Input
                  TYPE(local_stiffness_dist) :: this
c                       Local
                  integer :: row
                  do row=1,this%local_n
                        this%vec_indexes(row) = row+this%local_start-1
                  end do
            end subroutine get_vec_indexes
c
c                 Subroutine to get a guess at the soln of the system
c                 For now zero, but we could also use the last converged soln
c                 if it is available.
            subroutine get_guess(this)
            implicit none
c                 Input
            type(local_stiffness_dist) :: this
c                 Local
            integer :: row
            do row=1,this%local_n
                  this%soln(row) = 0.0
            end do
            end subroutine get_guess
c
c                       Subroutine to setup distribution to ranks
c                       and allocate data structures
            subroutine setup_dist(this,n,nnz,k_ptrs)
                  implicit none
c                       Input
                  TYPE(local_stiffness_dist) :: this
                  integer :: n, nnz, k_ptrs
                  dimension :: k_ptrs(*)
c                       Local
                  integer :: ierr, temp, iter, mpi_status
c                       Root sets n and broadcasts.  Everyone then uses
c                       the global n to figure their local number of rows
c                       and the start of their row block.
                  if (this%my_rank .eq. 0) then
                        this%global_n = n
                  end if
                  call MPI_Bcast(this%global_n,1,MPI_INTEGER,
     &                  0,this%comm,ierr)
                  if (this%my_rank .eq. 0) then
                        this%global_nnz = nnz
                  end if
                  call MPI_Bcast(this%global_nnz,1,MPI_INTEGER,
     &                  0,this%comm,ierr)
                  this%local_start = block_low(this%my_rank, 
     &                  this%num_procs, this%global_n)
                  this%local_n = block_size(this%my_rank,
     &                  this%num_procs,this%global_n)
c                       Root figures out everyone's local_nnz and sends it.
c                       Do this via the usual MPI if root then deal with your
c                       local block, then send everyone elses, else recieve.
                  if (this%my_rank .eq. 0) then
                        this%local_nnz=k_ptrs(block_high(this%my_rank, 
     &                 this%num_procs, this%global_n)+1)-k_ptrs(
     &                 block_low(this%my_rank,this%num_procs,
     &                 this%global_n))
                        do iter=1,this%num_procs-1
                              temp=k_ptrs(block_high(
     &                       iter, this%num_procs, this%
     &                       global_n)+1)-k_ptrs(block_low(
     &                       iter,this%num_procs,this%global_n))
                              call MPI_Send(temp,1,MPI_INTEGER,
     &                       iter,1,this%comm,ierr)
                        end do
                  else
                        call MPI_Recv(this%local_nnz,1,MPI_INTEGER,
     &                        0, 1, this%comm, mpi_status, ierr)
                  end if
c                       Allocate the actual data structs
                  allocate(this%row_ptrs(this%local_n+1))
                  allocate(this%col_indexes(this%local_nnz))
                  allocate(this%coefs(this%local_nnz))
                  allocate(this%rhs(this%local_n))
                  allocate(this%soln(this%local_n))
                  allocate(this%cols_per_row(this%local_n))
                  allocate(this%vec_indexes(this%local_n))
c
                  return
            end subroutine setup_dist
c
c           Subroutine to deallocate the local data structure.
            subroutine clear_dist(this)
                  implicit none
c                 Input
                  type(local_stiffness_dist) ::this
c                 
                  if (allocated(this%row_ptrs)) then
                        deallocate(this%row_ptrs)
                        deallocate(this%col_indexes)
                        deallocate(this%coefs)
                        deallocate(this%rhs)
                        deallocate(this%soln)
                        deallocate(this%cols_per_row)
                        deallocate(this%vec_indexes)
                  end if
                  return
            end subroutine clear_dist
c
c                 Distribute everything (i.e. call this on code 1)
c                 We use a subroutine to distribute the coefs and vectors
c                 so that if we want to do a partial distribution (same
c                 sparsity structure) we can
            subroutine distribute_full(this,k_ptrs,k_indexes,k_coefs,
     &            rhs,soln)
                  implicit none
c                       Input
                  type(local_stiffness_dist) :: this
                  integer :: k_ptrs, k_indexes
                  double precision :: k_coefs, rhs, soln
                  dimension :: k_ptrs(*), k_indexes(*), k_coefs(*),
     &                  rhs(*), soln(*)
c                       Local
                  integer :: ierr, mpi_status,i, temp
                  integer, allocatable :: send_counts(:), send_ptrs(:)
c                       Generally, for each item we need to send (k_ptrs,
c                       col_indices) create the MPI send arrays on root,                       root deal with its own local block, then send out the
c                       and then all MPI_Gatherv which takes care of the rest.
c                  if (this%my_rank .eq. 0) then
                        allocate(send_counts(this%num_procs))
                        allocate(send_ptrs(this%num_procs))
c                  end if
c                       k_ptrs
                  if (this%my_rank .eq. 0) then
                  do i=1,this%num_procs
                        send_counts(i) = block_size(i-1,this%num_procs, 
     &                        this%global_n)+1
                        send_ptrs(i) = block_low(i-1,this%num_procs,
     &                        this%global_n)-1
                  end do
                  end if
                  call MPI_Scatterv(k_ptrs,send_counts,send_ptrs,
     &                  MPI_INTEGER, this%row_ptrs, this%local_n+1,
     &                  MPI_INTEGER, 0, this%comm, ierr)
c                       Adjust to 1 (fortran...)
                  temp = this%row_ptrs(1)
                  do i=1,this%local_n+1
                        this%row_ptrs(i) = this%row_ptrs(i)-
     &                  temp+1
                  end do
c                       k_indexes
                  if (this%my_rank .eq. 0) then
                  do i=1,this%num_procs
                        send_counts(i)=k_ptrs(block_high(
     &                       i-1, this%num_procs, this%
     &                       global_n)+1)-k_ptrs(block_low(
     &                       i-1,this%num_procs,this%global_n))
                        send_ptrs(i) = k_ptrs(block_low(i-1,
     &                 this%num_procs,this%global_n))-1
                  end do
                  end if
                  call MPI_Scatterv(k_indexes,send_counts,send_ptrs,
     &                  MPI_INTEGER, this%col_indexes, this%local_nnz,
     &                  MPI_INTEGER,0,this%comm,ierr)
c                       Free our send arrays (all done).
                  if (this%my_rank .eq. 0) then
                        deallocate(send_counts)
                        deallocate(send_ptrs)
                  end if
c                       get hypre allocation data
                  call get_cols_per_row(this)
                  call get_vec_indexes(this)
                  call get_guess(this)
c                       Call the coef distribute function to do this rest
                  call distribute_coefs(this,k_ptrs,k_coefs,rhs)
                  return
            end subroutine distribute_full
c
c                 Subroutine to distribute just the stiffness coefs
c                 and RHS.  This (obviously) assumes that the structs
c                 were already allocated and we are keeping the same
c                 sparsity structure.
            subroutine distribute_coefs(this,k_ptrs,k_coefs,rhs)
                  implicit none
c                       Input
                  type(local_stiffness_dist) :: this
                  integer :: k_ptrs
                  double precision :: k_coefs, rhs
                  dimension :: k_coefs(*), rhs(*), k_ptrs(*)
c                       Local
                  integer :: ierr, mpi_status,i, temp
                  integer, allocatable :: send_counts(:), send_ptrs(:)
c                       Basically this is going to follow the same pattern
c                       as above.  Root will create the MPI send arrays and
c                       then everyone calls MPI_Gatherv which nicely deals
c                       with all the rest (including copying root's portion
c                       over to the allocated memory).
c                  if (this%my_rank .eq. 0) then
                        allocate(send_counts(this%num_procs))
                        allocate(send_ptrs(this%num_procs))
c                  end if
c                       rhs
                  if (this%my_rank .eq. 0) then
                  do i=1,this%num_procs
                        send_counts(i) = block_size(i-1,this%num_procs, 
     &                        this%global_n)
                        send_ptrs(i) = block_low(i-1,this%num_procs,
     &                        this%global_n)-1
                  end do
                  end if
                  call MPI_Scatterv(rhs,send_counts,send_ptrs,
     &                  MPI_DOUBLE_PRECISION, this%rhs, this%local_n,
     &                  MPI_DOUBLE_PRECISION, 0, this%comm, ierr)
c                       coefs
                  if (this%my_rank .eq. 0) then
                  do i=1,this%num_procs
                        send_counts(i)=k_ptrs(block_high(
     &                       i-1, this%num_procs, this%
     &                       global_n)+1)-k_ptrs(block_low(
     &                       i-1,this%num_procs,this%global_n))
                        send_ptrs(i) = k_ptrs(block_low(i-1,
     &                 this%num_procs,this%global_n))-1
                  end do
                  end if
                  call MPI_Scatterv(k_coefs,send_counts,send_ptrs,
     &                  MPI_DOUBLE_PRECISION, this%coefs, 
     &                  this%local_nnz,MPI_DOUBLE_PRECISION,0,
     &                  this%comm,ierr)
c                       Free our send arrays
                  if (this%my_rank .eq. 0) then
                        deallocate(send_counts)
                        deallocate(send_ptrs)
                  end if
                  return
            end subroutine distribute_coefs
c                 Determine which communicator to use and set your rank, etc.
c                 This is probably where we should move the icky business about
c                 spawning ranks.
            subroutine setup_comm(this, incomm)
                  implicit none
                  integer :: incomm
                  type(local_stiffness_dist) :: this
                  integer :: ierr
c
                  this%comm = incomm
                  call MPI_Comm_size(this%comm,this%num_procs,ierr)
                  call MPI_Comm_rank(this%comm,this%my_rank,ierr)
                  return
            end subroutine
c
c           Get some checksums
            subroutine stiffness_check_sums(this)
                  implicit none
                  type (local_stiffness_dist) :: this
c
                  integer :: l_col_sum, g_col_sum, ierr
                  integer :: l_ptr_sum, g_ptr_sum
                  double precision :: l_coef_sum,g_coef_sum
                  double precision :: l_rhs_sum,g_rhs_sum
c
                  l_ptr_sum = sum(this%row_ptrs)
                  l_col_sum = sum(this%col_indexes)
                  l_coef_sum = sum(this%coefs)
                  l_rhs_sum=sum(this%rhs)
c
                  call MPI_Reduce(l_ptr_sum,g_ptr_sum,1,MPI_INTEGER,
     &                  MPI_SUM,0,this%comm,ierr)
                  call MPI_Reduce(l_col_sum,g_col_sum,1,MPI_INTEGER,
     &                  MPI_SUM,0,this%comm,ierr)
                  call MPI_Reduce(l_coef_sum,g_coef_sum,1,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,this%comm,ierr)
                  call MPI_Reduce(l_rhs_sum,g_rhs_sum,1,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,this%comm,ierr)
c                  
                  if (this%my_rank .eq. 0) then
                        write (*,*) "Checksums: rows ", g_ptr_sum,
     &                        " columns ",g_col_sum,
     &                        " coefs ", g_coef_sum, " rhs ", g_rhs_sum
                  end if
            end subroutine

      end module local_stiffness_mod
