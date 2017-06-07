c     Checklist                                                                 
c           1) Fix element killing issues                                       
c           2) Get mat-vec map working?                                         
c           4) Fix whatever issue we're having with the adaptive stepping       
                                                                                
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
c            use local_stiffness_mod                                            
            implicit none                                                       
c                       Parameters (user input or otherwise)                    
            logical :: parallel_assembly_allowed                                
            logical :: parallel_assembly_used                                   
            logical :: distributed_stiffness_used                               
            integer :: initial_map_type                                         
            integer :: final_map_type                                           
                                                                                
c                       Actual data                                             
c            type(local_stiffness_dist) :: local_k                              
c                       Whether we need to recreate the mapping (1) or just     
c                       redistribute the new coefs                              
c            integer :: assembly_type                                           
c                       Distributed (full) dof_eqn_map and eqn_node_map         
c            integer, allocatable :: dist_dof_eqn_map(:),                       
c     &                        dist_eqn_node_map(:)                             
c                       Local sparsity structure                                
c            integer              :: g_n, g_nnz                                 
c            integer              :: my_n, my_nnz                               
c            integer, allocatable :: my_eqns(:), my_eqns_inv(:)                 
c            integer, allocatable :: my_row_ptrs(:),my_col_indices(:)           
c            double precision, allocatable :: my_coefs(:)                       
c            integer, allocatable :: global_nnz_vec(:,:)                        
c                       Initial map : just for sparse structure so we           
c                       can assemble the global sparse structure and            
c                       load balance.  Get locality only.                       
c            integer, allocatable :: initial_map(:)                             
c                       Your portion of the global sparsity structure           
c                       we need this because we're still using the WARP         
c                       communicator                                            
c            integer :: g_my_n, g_my_nnz                                        
c            integer, allocatable :: g_row_ptrs(:), g_col_indices(:)            
c            integer, allocatable :: int_my_eqns(:),int_my_eqns_inv(:)          
c                       Map to actual distributed stiffness                     
c            integer, allocatable :: final_map(:)                               
c            integer              :: final_nprocs                               
c            integer, allocatable :: final_my_eqns(:), final_eqn_map(:)         
c            integer, allocatable :: final_my_eqns_inv(:)                       
c            integer, allocatable :: final_eqn_map_inv(:)                       
c                                                                               
c            double precision, allocatable :: load_vec(:)                       
                                                                                
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
c     *     2 = load balance rows - balance the rows based on a simple      *   
c     *         switching algorithm.  Subsidiary parameters are             *   
c     *         load_bal_tol and max_load_bal_iters                         *   
c     *                                                                     *   
c     ***********************************************************************   
c                                                                               
      subroutine determine_initial_map                                          
                                                                                
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
      subroutine assemble_coefs                                                 
                                                                                
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
            implicit none                                                       
                                                                                
            curr_neqns = -1                                                     
                                                                                
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
                                                                                
      end                                                                       
