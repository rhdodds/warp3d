c     ****************************************************************
c     *                                                              *
c     *  consider the upper triangle of the assembled equilibrium    *
c     *  equations. this routine computes the number of non-zero     *
c     *  terms that appear above the diagonal. the process employs   *
c     *  a psuedo-assembly algorithm sequentially by structural      *
c     *  equation number. this is opposite of the conventional       *
c     *  element-by-element assembly. for a given row of the         *
c     *  equilibrium equations, we simulate the actual assembly      *
c     *  of the row to count the final number of columns with        *
c     *  non-zero entries.                                           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/25/97                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine count_profile_terms( neqns, eqn_node_map, ncoeff,
     &                                iprops, mxelpr,
     &                                edest, mxedof, dof_eqn_map,
     &                                scol_flags, scol_list,
     &                                repeat_incid)
      use main_data, only : inverse_incidences
      implicit integer (a-z)
c
c                    parameter declarations
c
      dimension   eqn_node_map(*), dof_eqn_map(*), scol_flags(*),
     &            scol_list(*), iprops(mxelpr,*),
     &            edest(mxedof,*)
c
c
      logical repeat_incid(*)
c
c                    local declarations
c
c                 simulate assembly of the equilibrium equations directly
c                 in sparse matrix format, i.e., only the non-zero
c                 coefficients are stored on each row of the upper
c                 triangle. the assembly process
c                 proceeds equation-by-equation, rather than the
c                 traditional element-by-element manner. the outline
c                 of the algorithm is:
c
c                 loop over all equations:
c
c                  (a) find the corresponding structure node number
c
c                  (b) get the number of elements incident on the node
c
c                  (b.1) get the displacement indexes for elements
c                        attached to the node. (skip if this node
c                        same as last node).
c
c                  (c) loop over elements connected to the node
c                       (1) skip null elements
c                       (2) we need to know what row of the element
c                           stiffness corresponds to the
c                           current equilibrium equation. the
c                           best current algorithm uses the simple
c                           search over the element destination vector.
c                           this resolves to erow as the current element
c                           row.
c                       (3) loop over all columns of the element [ke]
c                           for row erow. mark the corresponding
c                           strutural column numbers in which the
c                           element terms appear (scol_flags)
c                       (4) we update two integer (row) vectors in this
c                           process.
c                           scol_flags() has neqns entries. initially
c                              zeroed. we put a 1 in each structure column
c                              that has a [ke] entry.
c                           scol_list() has entries listing columns
c                              used in scol_flags. there will be many
c                              duplicates here since multiple elements
c                              can contribute stiffness terms at the
c                              same scol. allocated space of scol_list
c                              must be a multiple of the number of element
c                              dofs (mxconn * # ele dofs * a safety
c                              number of say 10)
c                           the first time we mark the structure column
c                           as having an element contribution, we increment
c                           the count of non-zero terms above the diagonal.
c
c                  (d) zero out the entries of scol_flags used for
c                      this row to set up for the next equation.
c
      do scol = 1, neqns
         scol_flags(scol) = 0
      end do
      num_non_zero_terms = 0
c
c                 loop over all structural equations (skip last one
c                 since it has no terms right of diagonal).
c
      count         = 0
      previous_node = 0
      do srow = 1, neqns-1
       node           = eqn_node_map(srow)
       numele         = inverse_incidences(node)%element_count
       next_scol_list = 0
       if ( node .ne. previous_node ) then
         call get_edest_terms(
     &        edest, inverse_incidences(node)%element_list(1), numele )
         previous_node = node
       end if
c
c                 process all elements which contribute stiffness terms
c                 to this structural equation.
c
       do 300 j = 1, numele
          felem = inverse_incidences(node)%element_list(j)
          if ( felem .le. 0 ) go to 300
          totdof = iprops(2,felem) * iprops(4,felem)
          do 200 erow = 1, totdof
            if ( dof_eqn_map(edest(erow,j)) .ne. srow ) go to 200
            do 100 ecol = 1, totdof
              scol = dof_eqn_map(edest(ecol,j))
              if ( scol .le. srow ) go to 100
              if ( scol_flags(scol) .eq. 0 ) then
                scol_flags(scol)   = 1
                num_non_zero_terms = num_non_zero_terms + 1
                next_scol_list     = next_scol_list + 1
                scol_list(next_scol_list) = scol
              end if
 100        continue
            if ( .not. repeat_incid(felem)) go to 300
 200      continue
 300    continue
c
c                 zero the flags used by this equation to set up next one.
c
        do j = 1, next_scol_list
          scol_flags(scol_list(j)) = 0
        end do
c
      end do
c
c                 set the total number of non-zero terms in the upper
c                 triangle for the calling routine.
c
      ncoeff = num_non_zero_terms
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *  create the key sparse data structures that describe the     *
c     *  location of non-zero terms in the upper triangle of the     *
c     *  assembled equilibrium equations. these are k_ptrs and       *
c     *  k_indexes (see comments at end of this file)                *
c     *  the req'd length of k_indexes has been previously computed  *
c     *  with the space allocated by the calling routine.            *
c     *  the algortihm here uses an equation-by-equation process to  *
c     *  build the terms in k_indexes and the corresponding k_ptrs   *
c     *  entry.                                                      *
c     *                                                              *
c     *                   written by : rhd                           *
c     *                                                              *
c     *                   last modified : 04/25/97                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine build_col_sparse( neqns, eqn_node_map, last_k_index,
     &                             dof_eqn_map, scol_list, k_indexes,
     &                             k_ptrs, iprops, mxelpr, edest,
     &                             mxedof, repeat_incid )
      use main_data, only : inverse_incidences
      implicit integer (a-z)
c
c                    parameter declarations
c
      dimension   eqn_node_map(*), dof_eqn_map(*), k_ptrs(*),
     &            k_indexes(*), scol_list(*),
     &            iprops(mxelpr,*), edest(mxedof,*)
c
      logical repeat_incid(*)
c
c                    local declarations
c
c
      data one / 1 /
c
c                 simulate assembly of the equilibrium equations directly
c                 in sparse matrix format, i.e., only the non-zero
c                 coefficients are stored on each row of the upper
c                 triangle. the assembly process
c                 proceeds equation-by-equation, rather than the
c                 traditional element-by-element manner. the outline
c                 of the algorithm is:
c
c                 loop over all equations:
c
c                  (a) find the corresponding structure node number
c
c                  (b) get the number of elements on the node
c
c                  (b.1) get the displacent indexes for elements
c                        attached to the node. (skip if this node
c                        same as last node).
c
c                  (c) zero a vector that can hold a full row
c                      of the equilibrum equations.
c
c                  (d) loop over elements connected to the node
c                       (1) skip null elements
c                       (2) we need to know what row of the element
c                           stiffness corresponds to the
c                           current equilibrium equation. the
c                           best current algorithm uses the simple
c                           search over the element destination vector.
c                           this resolves to erow as the current element
c                           row.
c                       (3) loop over all columns of the element [ke]
c                           for row erow. add the structure column
c                           number into which the element tern would be
c                           added to the list of used column numbers
c                           on this equation. creates scol_list.
c
c                  (e) sort the list of used column numbers on this
c                      equation (to right of diagonal) not including
c                      the diagonal. there are duplicate column
c                      numbers in the list
c
c                  (f) copy the list of used column numbers in
c                      ascending order (skip duplicates) into k_indexes.
c                      number used columns for the equation is the
c                      k_ptrs. diagonal terms are not included.
c
      next_k_index = 0
c
c                 loop over all structural equations (skip last one
c                 since it has no terms right of diagonal).
c
      previous_node = 0
      do srow = 1, neqns-one
       node           = eqn_node_map(srow)
       numele         = inverse_incidences(node)%element_count
       num_scols_used = 0
       if ( node .ne. previous_node ) then
         call get_edest_terms(
     &       edest, inverse_incidences(node)%element_list(1), numele )
         previous_node = node
       end if
c
c                 process all elements which contribute stiffness terms
c                 to this structural equation.
c
       do 300 j = 1, numele
          felem = inverse_incidences(node)%element_list(j)
          if ( felem .le. 0 ) go to 300
          totdof = iprops(2,felem) * iprops(4,felem)
          do 200 erow = 1, totdof
            if ( dof_eqn_map(edest(erow,j)) .ne. srow ) go to 200
            do 100 ecol = 1, totdof
              scol = dof_eqn_map(edest(ecol,j))
              if ( scol .le. srow ) go to 100
              num_scols_used = num_scols_used + one
              scol_list(num_scols_used) = scol
 100        continue
            if ( .not. repeat_incid(felem)) go to 300
 200      continue
 300    continue
c
c               sort the list of non-zero columns on this row. use an
c               in-place heap sort algorithm. duplicates in the list
c               of columns are still present.
c
        if ( num_scols_used .lt. 2 ) go to 600
        k  = num_scols_used/2 + one
        ir = num_scols_used
400     continue
          if( k .gt. one ) then
            k   = k - one
            rra = scol_list(k)
          else
            rra           = scol_list(ir)
            scol_list(ir) = scol_list(one)
            ir            = ir - one
            if( ir .eq. one ) then
              scol_list(one) = rra
              go to 600
            end if
          end if
          i = k
          j = k + k
500       if( j .le. ir ) then
            if( j .lt. ir ) then
              if( scol_list(j) .lt. scol_list(j+1) ) j = j + one
            end if
            if( rra .lt. scol_list(j) ) then
              scol_list(i) = scol_list(j)
              i = j
              j = j + j
            else
              j = ir + one
            end if
          go to 500
          end if
          scol_list(i) = rra
        go to 400
c
c               sort for this equation completed. copy the used
c               (unique) column numbers from the sorted list into
c               k_indexes (skip duplicates). count of moved
c               column numbers is k_ptrs for this equations. we
c               can have no terms to right of diagonal or only
c               1 term. treat as special cases.
c
 600    continue
        if ( num_scols_used .eq. 0 ) then
          k_ptrs(srow) = 0
          go to 700
        end if
        num_unique_cols         = one
        next_k_index            = next_k_index + one
        k_indexes(next_k_index) = scol_list(one)
        if ( num_scols_used .eq. one ) then
           k_ptrs(srow) = one
           go to 700
        end if
        do j = 2, num_scols_used
           if ( scol_list(j) .ne. scol_list(j-1) ) then
               num_unique_cols         = num_unique_cols + one
               next_k_index            = next_k_index + one
               k_indexes(next_k_index) = scol_list(j)
           end if
        end do
        k_ptrs(srow) = num_unique_cols
 700    continue
c
c               do next equation...
c
      end do
c
      k_ptrs(neqns) = 0
      last_k_index  = next_k_index
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *  build mapping between dof at all structure nodes to the     *
c     *  equilibrium equations generated for unconstrained dof       *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 04/15/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine dof_map( dof_eqn_map, cstmap, nonode, ndof,
     &                    eqn_node_map, neqns  )
      implicit none
c
c         parameters
c
      integer :: cstmap(*), dof_eqn_map(*), eqn_node_map(*)
      integer :: nonode, ndof, neqns
c
c         local
c
      integer :: eqn_counter, dof 
c
c         assign the equation number for each structure dof.
c         dof with absolute constraints do not appear in the equations.
c
c         dof_eqn_map(sdof) -> equation number for structure sdof
c
c         eqn_node_map(equation) -> structure node number for
c                                   given equation number
c
      eqn_counter = 0
      do dof = 1, nonode*ndof ! ndof is always = 3
        eqn_node_map(dof) = 0
        dof_eqn_map(dof)  = 0
      end do
c
      do dof =  1, ndof*nonode
       if( cstmap(dof) .eq. 0 ) then ! dof has no abs constraint
           eqn_counter      = eqn_counter + 1
           dof_eqn_map(dof) = eqn_counter
           eqn_node_map(eqn_counter) = (dof-1)/ndof + 1
        end if
      end do
c      
c         set the number of equilibrium equations to be solved in each 
c         iteration of the load(time) step.
c         
      neqns = eqn_counter
c      
      return
      end
c ------------------------------------------------------------------------
c                        NASA - VSS
c
c               sparse matrix storage format
c
c   example:
c
c                1    2    3    4     5     6
c
c          1  | 100   1    2                5  |  | x1 |     | 201 |
c          2  |     200    6    7           9  |  | x2 |     | 202 |
c          3  |          300   10    11    12  |  | x3 |     | 203 |
c      a = 4  |                400   13    14  |  | x4 |  =  | 204 |
c          5  |                     500    15  |  | x5 |     | 205 |
c          6  |                           600  |  | x6 |     | 206 |
c
c     number of equations    = 6
c
c     number of coefficients = 12
c
c
c     k_ptrs    = { 3, 3, 3, 2, 1, 0}
c
c     k_indesxs = { 2, 3, 6,  3, 4, 6,  4, 5, 6,  5, 6,  6}
c
c     k_coefs   = { 1, 2, 5,  6, 7, 9, 10,11,12, 13,14, 15}
c
c     k_diag    = { 100, 200, 300, 400, 500, 600}
c
c     k_rhs     = { 201, 202, 203, 204, 205, 206}
c
c ------------------------------------------------------------------------
c     ****************************************************************
c     *                                                              *
c     *  assembly of the equilibrium equations in sparse format.     *
c     *  we use a structure row-by-row assembly procedure rather     *
c     *  than conventional element-by-element. this structure        *
c     *  supports threaded parallel assembly: a thread builds        *
c     *  a full row of the quations. multiple threads can build      *
c     *  multiple rows simultaneously.                               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *            last modified : 02/22/10 rhd fix threads          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine assem_by_row( neqns, num_threads, eqn_node_map,
     &                         dof_eqn_map, k_diag, k_coeffs,
     &                         k_indexes, k_ptrs, iprops, dcp,
     &                         noelem )
      implicit integer (a-z)
$add param_def
c
c                    parameter declarations
c
      dimension   eqn_node_map(*), dof_eqn_map(*), k_ptrs(*),
     &            k_indexes(*), iprops(mxelpr,*), dcp(*)
#dbl      double precision
#sgl      real
     & k_diag(*), k_coeffs(*)
c
c                    local declarations
c
#dbl      double precision
#sgl      real
     & zero, coeff_row
      allocatable coeff_row(:,:)
      allocatable row_start_index(:)
      allocatable edest(:,:,:)
      integer thread_previous_node(max_threads)

c
      data zero / 0.0 /
c
c                 assemble the equilibrium equations directly in
c                 sparse matrix format, i.e., only the non-zero
c                 coefficients are stored. the assembly process
c                 proceeds equation-by-equation, rather than the
c                 traditional element-by-element manner. the outline
c                 of the algorithm is:
c
c                 loop over all equations (in threaded parallel)
c
c                  (a) find the corresponding structure node number
c
c                  (b) get the number of elements on the node
c
c                  (b.1) get the displacement indexes for elements
c                        attached to the node. (skip if this node
c                        same as last node). We keep track of previous
c                        node processed by each thread. Each thread has its
c                        own block of indexes.
c
c                  (c) zero a vector that can hold a full row
c                      of the equilibrum equations. each thread has
c                      its own allocated row for this.
c
c                  (d) loop over elements connected to the node
c                       (1) skip null elements (killed elems)
c                       (2) we need to know what row of the element
c                           stiffness corresponds to the
c                           current equilibrium equation. the
c                           best current algorithm uses the simple
c                           search over the element destination vector.
c                           this resolves to erow as the current element
c                           row.
c                       (3) loop over all columns of the element [ke]
c                           for row erow. add the [ke] terms into
c                           correct position on the structure row,
c                           using coeff_row to support random access.
c                           (the diagonal term is copied as well but
c                            never used.)
c
c                  (e) update the diagonal term for the current eqn.
c                      copy the non-zero terms in coeff_row into the
c                      packed k_coeffs vector. use k_indexes to know
c                      which terms on row are non-zero. re-zero the
c                      terms in coeff_row to set up next equation.
c
c                 row_start_index(i) gives the starting index in the
c                 packed vector of assembled equations for equation i.
c
      allocate( coeff_row(neqns,num_threads) )
      allocate( row_start_index(neqns) )
c
      coeff_row(1:neqns,1:num_threads) = zero
      row_start_index(1) = 1
c
      do i = 2, neqns
       row_start_index(i) = row_start_index(i-1) + k_ptrs(i-1)
      end do
c
      allocate( edest(mxedof,mxconn,num_threads) )
      thread_previous_node(1:num_threads) = 0
c
c                 loop over all structural equations. the row by row
c                 assembly makes threaded parallelism reasonably simple.
c                 Code above handles making private copies
c                 of arrays for threads above.
c
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO PRIVATE( srow, now_thread )
c
      do srow = 1, neqns
c
         now_thread = omp_get_thread_num() + 1
         call assem_a_row(
     &     srow, neqns, eqn_node_map, coeff_row(1,now_thread),
     &     dof_eqn_map, k_diag, k_coeffs, k_indexes, k_ptrs,
     &     iprops, dcp, noelem, row_start_index,
     &     edest(1,1,now_thread),
     &     thread_previous_node(now_thread) )
c
      end do
c$OMP END PARALLEL DO
c
      deallocate( coeff_row, row_start_index, edest )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *  assembly a single row of the equilibrium equations          *
c     *  in sparse format. this code is designed to be run inside    *
c     *  a threaded loop over all equations.                         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/29/09 rdh add threads   *
c     *                                                              *
c     ****************************************************************
c
      subroutine assem_a_row( srow, neqns, eqn_node_map, coeff_row,
     &                        dof_eqn_map, k_diag, k_coeffs,
     &                        k_indexes, k_ptrs, iprops,  dcp,
     &                        noelem, row_start_index, edest,
     &                        previous_node )
      use elem_block_data, only : estiff_blocks
      use main_data,       only : elems_to_blocks, repeat_incid,
     &                            inverse_incidences
      implicit integer (a-z)
$add param_def
c
c                    parameter declarations
c
      dimension   eqn_node_map(*), dof_eqn_map(*), k_ptrs(*),
     &             k_indexes(*), iprops(mxelpr,*),
     &             dcp(*), row_start_index(*), edest(mxedof,*)
#dbl      double precision
#sgl      real
     & k_diag(*), k_coeffs(*), coeff_row(*)
c
c                    local declarations
c
#dbl      double precision
#sgl      real
     & zero, ekterm
#dbl      double precision, dimension(:,:), pointer :: emat
#sgl      real, dimension(:,:), pointer :: emat
c
      logical repeated
      data zero / 0.0 /
c
c                 the structure node number corresponding to this equation
c                 pull the number of elements connected to the node.
c                 get the destination numbers for elements connected to
c                 this node. note we keep track of previous
c                 structure node processed for this thread. if the
c                 same we don't need to get new destination indexes.
c                 (see edest and previous_node made thread
c                 private by our code in calling routine).
c
      node   = eqn_node_map(srow)
      numele = inverse_incidences(node)%element_count
      if ( node .ne. previous_node ) then
        call get_edest_terms(
     &      edest, inverse_incidences(node)%element_list(1), numele )
        previous_node = node
      end if
c
c                 process all elements which contribute stiffness terms
c                 to this structural equation.
c
       do 300 j = 1, numele
          felem = inverse_incidences(node)%element_list(j)
          repeated = repeat_incid(felem)
          if ( felem .le. 0 ) go to 300
          totdof  = iprops(2,felem) * iprops(4,felem)
          blk     = elems_to_blocks(felem,1)
          rel_col = elems_to_blocks(felem,2)
          emat    => estiff_blocks(blk)%ptr
          do 200 erow = 1, totdof
             if ( dof_eqn_map(edest(erow,j)) .ne. srow ) go to 200
             do 100 ecol = 1, totdof
                scol = dof_eqn_map(edest(ecol,j))
                if ( scol .lt. srow ) go to 100
                kk = dcp(max(ecol,erow))-abs(ecol - erow)
                ekterm = emat(kk,rel_col)
                coeff_row(scol) = coeff_row(scol) + ekterm
 100         continue
             if ( .not. repeated) go to 300
 200      continue
 300  continue
c
c                 copy the diagonal term from coeff_row into k_diag.
c
      k_diag(srow) = k_diag(srow) + coeff_row(srow)
      coeff_row(srow) = zero
c
c                 pack the non-zero stiffness coefficients on this
c                 equation (right of diagonal) into the vector
c                 of coefficients. zero the used entries in the
c                 work row (coeff_row) which supports the random access
c                 above. the k_indexes tell us which terms in coeff_row
c                 have been filled. k_ptrs tells us number of non-zero
c                 terms right of diagonal. these key sparse data
c                 structures have been built earlier in assembly process.
c
      if ( srow .eq. neqns ) return
      start_loc = row_start_index(srow)-1
      do j = 1, k_ptrs(srow)
         k_coeffs(start_loc+j) = coeff_row(k_indexes(start_loc+j))
         coeff_row(k_indexes(start_loc+j)) = zero
      end do
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *  assembly of the equilibrium equations in sparse format.     *
c     *  we use a structure row-by-row assembly procedure rather     *
c     *  than conventional element-by-element. this structure        *
c     *  supports threaded parallel assembly: a thread builds        *
c     *  a full row of the quations. multiple threads can build      *
c     *  multiple rows simultaneously.                               *
c     *                                                              *
c     *  This version is for asymmetric assembly with the full CSR   *
c     *  sparsity.
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *            last modified : 12/10/13                          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine assem_by_row_a( neqns, nnz, num_threads, eqn_node_map,
     &                         dof_eqn_map, k_ptrs, k_indexes, iprops,
     &                         dcp, noelem, k_coeffs)
      implicit integer (a-z)
$add param_def
c
c                    parameter declarations
c
      dimension   eqn_node_map(*), dof_eqn_map(*),
     &            iprops(mxelpr,*), dcp(*)
      integer :: neqns, nnz, num_threads, k_ptrs(neqns+1), 
     &            k_indexes(nnz)
#dbl      double precision
#sgl      real
     & k_coeffs(nnz)
c     Internal use
      integer, allocatable :: edest(:,:,:)
      integer :: srow, now_thread, nnzr
c 
c       Algorithm is a bit easier than the original one
c             1) Loop on each equation in parallel
c             2) Pass each thread (via a subroutine) the row number,
c                the row column indexes, and the section of the coeffs
c                array where we want out data and all the other 
c                mapping info.
c             3) For each row then:
c                   a) Find the node number
c                   b) Look up all elements connected to this node
c                   c) Loop over those elements
c                   d) Skip killed elements
c                   e) Find the row of the element corresponding to
c                      the node/dof we're looking at
c                   f) Loop over all columns in this local row and
c                      insert/add their values where they belong
c
      allocate(edest(mxedof,mxconn,num_threads))
c
c$OMP PARALLEL DO PRIVATE(srow, now_thread, nnzr)
      do srow = 1, neqns
        now_thread = omp_get_thread_num() + 1
        nnzr = k_ptrs(srow+1) - k_ptrs(srow)
        call assem_a_row_a(srow, neqns, nnz, eqn_node_map, dof_eqn_map,
     &    iprops, dcp, noelem, nnzr,
     &    k_indexes(k_ptrs(srow):k_ptrs(srow+1)-1),
     &    edest(:,:,now_thread), 
     &    k_coeffs(k_ptrs(srow):k_ptrs(srow+1)-1))
      end do
c$OMP END PARALLEL DO

      deallocate(edest)

      return

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *  assembly a single row of the equilibrium equations          *
c     *  in sparse format. this code is designed to be run inside    *
c     *  a threaded loop over all equations.                         *
c     *                                                              *
c     *  This version is for full, asymmetric assembly               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 12/10/13                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine assem_a_row_a(srow, neqns, nnz, eqn_node_map, 
     &    dof_eqn_map, iprops, dcp, noelem, nnzr, k_indexes,
     &    edest, k_coeffs)
      use elem_block_data, only : estiff_blocks
      use main_data, only: elems_to_blocks, repeat_incid,
     &                  inverse_incidences
      implicit integer (a-z)
$add param_def
      integer :: srow, neqns, nnz, noelem, nnzr
      integer :: eqn_node_map(*), dof_eqn_map(*), dcp(*),
     &           iprops(mxelpr,*), edest(mxedof,*)
      integer :: k_indexes(nnzr)
#dbl      double precision
#sgl      real
     & k_coeffs(nnzr)
c
#dbl      double precision, dimension(:,:), pointer :: emat
#sgl      real, dimension(:,:), pointer :: emat
#dbl      double precision :: val
#sgl      real :: val
      integer :: node, numele, e, felem, totdof, blk, rel_col,
     &            erow, ecol, scol, k
      logical :: repeated
c
c     Get the node number, element connections, and edest maps for
c     this row
c
      node = eqn_node_map(srow)
      numele = inverse_incidences(node)%element_count
      call get_edest_terms(edest, 
     &      inverse_incidences(node)%element_list(1), numele)
c
c     Loop on elements
c
      do e = 1, numele
        felem = inverse_incidences(node)%element_list(e)
        repeated = repeat_incid(felem)
        if ( felem .le. 0) cycle
        totdof = iprops(2,felem) * iprops(4,felem)
        blk = elems_to_blocks(felem,1)
        rel_col = elems_to_blocks(felem,2)
        emat => estiff_blocks(blk)%ptr
c           Loop through the element to find the row.  This is the
c           inefficiency mentioned above
        do erow = 1, totdof
          if (dof_eqn_map(edest(erow,e)) .ne. srow) cycle
c           Now loop through the columns
          do ecol = 1, totdof
            scol = dof_eqn_map(edest(ecol,e))
            ! Okay now we're trying to insert the reshaped 
            ! emat(erow,ecol,rel_col)
            ! into global position srow, scol
            val = emat((ecol-1)*totdof+erow,rel_col)
            do k=1,nnzr
              if (k_indexes(k) .eq. scol) then
                k_coeffs(k) = k_coeffs(k) + val
              end if
            end do
          end do
          if (.not. repeated) exit
        end do
      end do
c
      return

      end subroutine
