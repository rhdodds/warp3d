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
c     *                   last modified : 6/6/2017 rhd               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine count_profile_symmetric( neqns, ncoeff, num_threads,           
     &                                    eqn_node_map, iprops,                 
     &                                    dof_eqn_map, 
     &                                    start_kindex_locs,
     &                                    edest, scol_flags, 
     &                                    scol_lists, nrow_lists )                         
      use main_data, only : inverse_incidences                                 
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer :: neqns, ncoeff, num_threads, nrow_lists                                 
      integer :: eqn_node_map(*), iprops(mxelpr,*), dof_eqn_map(*),
     &           start_kindex_locs(*),
     &           edest(mxedof,mxconn,num_threads),
     &           scol_flags(neqns,num_threads),      
     &           scol_lists(nrow_lists,num_threads)       
c                                                                               
c                    locals                                                     
c                                                                               
      integer :: srow, now_thread                                
      integer, external :: omp_get_thread_num                                   
      logical, parameter :: local_debug = .false.                               
      integer :: thread_previous_node(max_threads),                             
     &           num_non_zero_terms(max_threads)  
c                                                                               
c            simulate assembly of the equilibrium equations directly            
c            in sparse matrix format, i.e., only the non-zero                   
c            coefficients are stored on each row of the upper                   
c            triangle. the assembly process proceeds equation-by-equation,      
c            rather than the traditional element-by-element manner.             
c            algorithm outline:                                                 
c                                                                               
c            loop over all equations: (thread parallel)                         
c                                                                               
c             (a) find the corresponding structure node number                  
c                                                                               
c             (b) get the number of elements incident on the node               
c                                                                               
c             (b.1) get the displacement indexes for elements                   
c                   attached to the node. (skip if this node                    
c                   same as last node).                                         
c                                                                               
c             (c) loop over elements connected to the node                      
c                  (1) skip null elements                                       
c                  (2) we need to know what row of the element                  
c                      stiffness corresponds to the                             
c                      current equilibrium equation. the                        
c                      best current algorithm uses the simple                   
c                      search over the element destination vector.              
c                      this resolves to erow as the current element             
c                      row.                                                     
c                  (3) loop over all columns of the element [ke]                
c                      for row erow. mark the corresponding                     
c                      strutural column numbers in which the                    
c                      element terms appear (scol_flags)                        
c                  (4) we update two integer (row) vectors in this              
c                      process.                                                 
c                      scol_flags() has neqns entries. initially                
c                         zeroed. we put a 1 in each structure column           
c                         that has a [ke] entry.                                
c                      scol_list() has entries listing columns                  
c                         used in scol_flags. there will be many                
c                         duplicates here since multiple elements               
c                         can contribute stiffness terms at the                 
c                         same scol. allocated space of scol_list               
c                         must be a multiple of the number of element           
c                         dofs (mxconn * # ele dofs * a safety                  
c                         number of say 10)                                     
c                      the first time we mark the structure column              
c                      as having an element contribution, we increment          
c                      the count of non-zero terms above the diagonal.          
c                                                                               
c             (d) zero out the entries of scol_flags used for                   
c                 this row to set up for the next equation.                     
c                                                                               
c             repeat_incid = .true. if the same structure node appears          
c             more than once in the incid list for element, i.e.,               
c             crack front collapsed elements. Then have to process              
c             all element [Ke] rows since more than 1 erow will                 
c             correspond to srow being processed.                               
c 
      if( local_debug ) write(*,9010) 
      scol_flags           = 0 ! 2d array                                       
      scol_lists           = 0 ! 2d array                                       
      num_non_zero_terms   = 0 ! vector                                         
      thread_previous_node = 0 ! vector      
c                                                                               
c             loop over structural equations. omit last eqn since it            
c             has no terms right of diagonal. with a few temps based on         
c             number of threads, the srow loop runs in parallel                 
c                                                                               
      call omp_set_dynamic( .false. )                                           
c$OMP PARALLEL DO PRIVATE( srow, now_thread ) ! all else shared                 
      do srow = 1, neqns-1                                                      
       now_thread = omp_get_thread_num() + 1                                    
       call count_profile_symmetric_srow( srow, neqns,                          
     &      thread_previous_node(now_thread),                                   
     &      num_non_zero_terms(now_thread), eqn_node_map,                       
     &      edest(1,1,now_thread), iprops, dof_eqn_map,                         
     &      scol_flags(1,now_thread), scol_lists(1,now_thread),
     &      start_kindex_locs(srow), nrow_lists )                
       end do                                                                   
c$OMP END PARALLEL DO                                                           
c                                                                               
c                 set the total number of non-zero terms in the upper           
c                 triangle                                                      
c                                                                               
      ncoeff = sum( num_non_zero_terms )                                        
      if( local_debug ) write(*,9000) ncoeff 
c                                                                               
      return                                                                    
c                                                                               
 9000 format(10x,
     & '... number of non-zero K-coeff above diagonal (ncoeff): ',i10,                             
     & /,1x,'..... leaving count_profile_symmetric .....')
 9010 format(1x,'..... entering count_profile_symmetric .....')                                                 
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                    count_profile_symmetric_srow              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 8/10/2017 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine count_profile_symmetric_srow(                                  
     &  srow, neqns, previous_snode, num_non_zero_terms, eqn_node_map,          
     &  edest, iprops, dof_eqn_map, scol_flags, scol_list, temp_kptr,
     &  nrow_lists )                     
c                                                                               
      use main_data, only : repeat_incid, inverse_incidences                    
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
      integer :: srow, neqns, previous_snode, num_non_zero_terms,
     &           nrow_lists                
      integer :: eqn_node_map(*), edest(mxedof,mxconn), 
     &           iprops(mxelpr,*),            
     &           dof_eqn_map(*), scol_flags(*), scol_list(nrow_lists)                    
c                                                                               
      integer :: snode, num_ele_on_snode, next_scol_list, ele_on_snode,         
     &           totdof, erow, ecol, scol, j  
     
      integer :: temp_kptr    
      logical, parameter :: ldebug = .false.                              
c                                                                               
      snode             = eqn_node_map(srow) ! structure node                   
      num_ele_on_snode  = inverse_incidences(snode)%element_count               
      next_scol_list    = 0                                                     
      if( snode .ne. previous_snode ) then                                      
         call get_edest_terms_assemble( edest,                                  
     &         inverse_incidences(snode)%element_list(1),                       
     &         num_ele_on_snode, iprops )                                       
         previous_snode = snode                                                 
      end if 
c      
      if( ldebug ) then 
         write(*,*) '    ... inside  count_profile_symmetric_srow ...'
         write(*,*) '          srow, snode, num_ele_on_snode: ',
     &      srow, snode, num_ele_on_snode
         write(*,*) '          edest:'
         do j = 1, 24
            write(*,9000) j, edest(j,1:num_ele_on_snode)
         end do
      end if   
c      
      temp_kptr = 0                                                                  
c                                                                               
c                 process all elements which contribute stiffness terms         
c                 to this structural equation.                                  
c                                                                               
      do j = 1, num_ele_on_snode                                                
          ele_on_snode = inverse_incidences(snode)%element_list(j)      
c          if( ldebug ) write(*,*) '... doing element: ',ele_on_snode
          if( ele_on_snode .le. 0 ) cycle                                       
          totdof = iprops(2,ele_on_snode) * iprops(4,ele_on_snode)  
c          if( ldebug ) write(*,*) '     totdof: ', totdof
          do erow = 1, totdof ! which erow matches srow  
c            if( ldebug ) write(*,*)
c     &        '    erow,  dof_eqn_map(edest(erow,j))', erow,
c     &            dof_eqn_map(edest(erow,j))                        
            if( dof_eqn_map(edest(erow,j)) .ne. srow ) cycle                    
            do ecol = 1, totdof                                                 
              scol = dof_eqn_map(edest(ecol,j))  
c              if( ldebug) write(*,*) '      ecol, scol: ', ecol, scol
c              if( scol .gt. neqns ) then
c                  write(*,*) '.... bad scol: ', scol
c                  call die_abort
c              end if                               
              if( scol .le. srow ) cycle ! would be on lower-triangle           
              if( scol_flags(scol) .eq. 0 ) then                                
                scol_flags(scol)   = 1                                          
                num_non_zero_terms = num_non_zero_terms + 1 
                temp_kptr = temp_kptr + 1                    
                next_scol_list     = next_scol_list + 1                         
                scol_list(next_scol_list) = scol                                
              end if                                                            
            end do ! on ecol                                                    
            if( .not. repeat_incid(ele_on_snode) ) exit ! erow loop             
          end do  ! on erow                                                     
      end do  ! on j                                                            
c                                                                               
c                 zero the flags used by this equation to set                   
c                 up next one. faster than zeroing entire vector                
c                                                                               
      do j = 1, next_scol_list                                                  
          scol_flags(scol_list(j)) = 0                                          
      end do                                                                    
c                                                                               
      return
c                                                                          
 9000 format(15x,i2,3i8)  
c      
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
c     *                   last modified : 8/10/2017 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine build_col_sparse_symm(                                         
     &  neqns, num_threads, eqn_node_map, dof_eqn_map,               
     &  k_indexes, k_ptrs, start_kindex_locs, edest, scol_lists,
     &  nrow_lists, ncoeff_above_diagonal )  
c                                                             
      use main_data, only : inverse_incidences
      use global_data, only : out, iprops, mxedof, mxconn, max_threads                        
      implicit none   
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer :: neqns, num_threads, nrow_lists, ncoeff_above_diagonal                          
      integer :: eqn_node_map(*), dof_eqn_map(*), k_ptrs(*),                    
     &           k_indexes(*), start_kindex_locs(*), 
     &           edest(mxedof,mxconn,num_threads),
     &           scol_lists(nrow_lists,num_threads)       
c                                                                               
c                    local declarations                                         
c                                                                               
      integer :: i, max_k_index, previous_node, srow, scol, snode,             
     &           num_scols_used, j, totdof, ele_on_snode, erow, ecol,           
     &           num_unique_cols, num_ele_on_snode, now_thread
      integer :: thread_previous_node(max_threads),
     &           last_k_index(max_threads)  ! for checking
      integer, external :: omp_get_thread_num                                   
c                                   
      logical, parameter :: local_debug = .false.    
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
c                 loop over all structural equations (skip last one             
c                 since it has no terms right of diagonal).                     
c 
      scol_lists           = 0 ! 2d array                                       
      thread_previous_node = 0 ! vector   
      last_k_index         = 0 ! vector for checking consistency  
c                                                         
      call omp_set_dynamic( .false. )                                           
c$OMP PARALLEL DO PRIVATE( srow, now_thread ) ! all else shared                 
      do srow = 1, neqns-1                                                      
       now_thread = omp_get_thread_num() + 1    
       call build_col_sparse_symm_srow( 
     &     srow, neqns, thread_previous_node(now_thread), 
     &     start_kindex_locs(srow), eqn_node_map(srow), dof_eqn_map, 
     &     k_indexes, k_ptrs, edest(1,1,now_thread), 
     &     scol_lists(1,now_thread), last_k_index(now_thread)  )  
      end do !  srow                                                            
c$OMP END PARALLEL DO                                                           
c                                                                               
      k_ptrs(neqns) = 0   
      max_k_index = maxval( last_k_index )
      if( max_k_index .ne. ncoeff_above_diagonal ) then
         write(out,9020)
         call die_abort
      end if                                                         
c                                                                               
      if( .not. local_debug ) return
      write(out,9000) max_k_index, neqns                                       
      write(out,*) '   k_ptrs ...'   
      do i = 1, neqns
         write(out,9010) i, k_ptrs(i)
      end do                                 
      write(out,*) '   k_indexes ...'                                           
      write(out,9010) (i,k_indexes(i),i=1,max_k_index)                         
c                                                                               
      return                                                                    
c                                                                               
 9000 format(2x,'... in build_col_sparse_symm: ',                               
     & /,10x,'max_k_index, neqns: ',2i8)                                       
 9010 format(10x,2i8)   
 9020 format(1x,                                                                
     &'>> FATAL ERROR: Job Aborted.',                                           
     &5x,' inconsistent data structure in build_col_sparse_symm' ) 
c                                                                               
      end     
      
c     ****************************************************************          
c     *                                                              *          
c     *                 build_col_sparse_symm_srow                   *   
c     *               (routine runs thread parallel)                 *       
c     *                                                              *          
c     *                   written by : rhd                           *          
c     *                                                              *          
c     *                   last modified : 8/10/2017 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine build_col_sparse_symm_srow(                                       
     &  srow, neqns, previous_snode,
     &  start_kindex_loc, snode, dof_eqn_map,               
     &  k_indexes, k_ptrs,
     &  edest, scol_list, last_k_index )    
c                                                           
      use main_data, only : inverse_incidences, repeat_incid     
      use global_data, only : mxelpr, mxedof, iprops, out          
      implicit none                                                             
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer :: srow, neqns, previous_snode, start_kindex_loc, snode,
     &           last_k_index      ! for checking                   
      integer :: dof_eqn_map(*), k_ptrs(*), k_indexes(*), scol_list(*),  
     &           edest(mxedof,*)                              
c                                                                               
c                    local declarations                                         
c                                                                               
      integer :: i, j, num_ele_on_snode, num_scols_used,
     &           ele_on_snode, totdof, erow, ecol, scol, 
     &           num_unique_cols, next_k_index
      logical, parameter :: local_debug = .false.    
c                                                                               
c                 for a single equation processed here (srow):                                      
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
c                                                                               
      num_ele_on_snode = inverse_incidences(snode)%element_count               
      num_scols_used   = 0                                                     
      if( snode .ne. previous_snode ) then                                      
         call get_edest_terms_assemble( edest,                                  
     &        inverse_incidences(snode)%element_list(1),                        
     &        num_ele_on_snode, iprops )                                        
         previous_snode = snode                                                  
      end if                                                                   
c                                                                               
c                 process all elements which contribute stiffness terms         
c                 to this structural equation.                                  
c                                                                               
      do j = 1, num_ele_on_snode                                               
          ele_on_snode = inverse_incidences(snode)%element_list(j)              
          if( ele_on_snode .le. 0 ) cycle                                       
          totdof = iprops(2,ele_on_snode) * iprops(4,ele_on_snode)              
          do erow = 1, totdof                                                   
            if( dof_eqn_map(edest(erow,j)) .ne. srow ) cycle                    
            do ecol = 1, totdof                                                 
              scol = dof_eqn_map(edest(ecol,j))                                 
              if( scol .le. srow ) cycle                                        
              num_scols_used = num_scols_used + 1                               
              scol_list(num_scols_used) = scol                                  
            end do                                                              
            if( .not. repeat_incid(ele_on_snode)) exit ! erow loop              
          end do ! erow                                                         
      end do  !                                                                
c                                                                               
c               sort the list of non-zero columns on this row.                  
c               duplicates in the list of columns are still present.            
c                                                                               
      if( num_scols_used > 1 )                                                
     &       call build_col_sparse_sort( num_scols_used, scol_list )            
c                                                                               
c               copy the used (unique) column numbers from the sorted           
c               list into k_indexes (skip duplicates). count of moved           
c               column numbers is k_ptrs for this equations. we                 
c               can have zero terms to right of diagonal or only                
c               1 term. treat as special cases.                                 
c                                                                               
      if( num_scols_used .eq. 0 ) then                                        
          k_ptrs(srow) = 0                                                      
          return                                                                
      end if                                                                  
c
      num_unique_cols         = 1                                             
      next_k_index            = start_kindex_loc   
      if( next_k_index == 0 ) then
          write(*,9300) '... next_k_index =0, srow: ', srow
          call die_abort
      end if                           
      k_indexes(next_k_index) = scol_list(1)      
      last_k_index = next_k_index                         
      if( num_scols_used .eq. 1 ) then                                        
           k_ptrs(srow) = 1                                                     
           return                                                              
      end if     
c                                                                   
      do j = 2, num_scols_used                                                
           if( scol_list(j) .eq. scol_list(j-1) ) cycle                         
           num_unique_cols         = num_unique_cols + 1                        
           next_k_index            = next_k_index + 1                           
           k_indexes(next_k_index) = scol_list(j)        
           last_k_index            = next_k_index  ! for checking                   
      end do                                                                  
      k_ptrs(srow) = num_unique_cols                                          
c                                                                               
      return                                                                    
c
 9300 format(1x,                                                                
     & '>> FATAL ERROR: Job Aborted.',                                           
     & 5x,'inconsistent data structure. build_col_sparse_symm_srow')
c
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *  heap sort the list of used column numbers for the equation  *          
c     *    -- this routine is inlined by the compilers --            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  7/20/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
                                                                                
      subroutine build_col_sparse_sort( n, vec )                                
      implicit none                                                             
      integer :: n, vec(n)                                                      
c                                                                               
      integer :: l, ir,rra, i, j                                                
c                                                                               
      l = n/2 + 1                                                               
      ir = n                                                                    
c                                                                               
c               the index l will be decremented from its initial                
c               value during the "hiring" (heap creation) phase. once           
c               it reaches 1, the index ir will be decremented from             
c               its initial value down to 1 during the                          
c               "retirement-and-promotion" (heap selection) phase.              
c                                                                               
10    continue                                                                  
      if( l > 1 )then                                                           
        l = l - 1                                                               
        rra = vec(l)                                                            
      else                                                                      
        rra = vec(ir)                                                           
        vec(ir) = vec(1)                                                        
        ir = ir-1                                                               
        if( ir .eq. 1 )then                                                     
          vec(1) = rra                                                          
          return                                                                
        end if                                                                  
      end if                                                                    
      i = l                                                                     
      j = l + l                                                                 
20    if( j .le. ir )then                                                       
        if( j < ir )then                                                        
           if( vec(j) < vec(j+1) )  j = j + 1                                   
        end if                                                                  
        if( rra < vec(j) )then                                                  
          vec(i) = vec(j)                                                       
          i = j; j = j + j                                                      
        else                                                                    
          j = ir + 1                                                            
        end if                                                                  
        go to 20                                                                
      end if                                                                    
      vec(i) = rra                                                              
      go to 10                                                                  
c                                                                               
      end       
c                                                                               
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
      logical, parameter :: local_debug = .false.                                              
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
      do dof = 1, ndof*nonode                                                   
       if( cstmap(dof) .ne. 0 ) cycle ! no abs constraint                       
       eqn_counter      = eqn_counter + 1                                       
       dof_eqn_map(dof) = eqn_counter                                           
       eqn_node_map(eqn_counter) = (dof-1)/ndof + 1                             
      end do                                                                    
c                                                                               
c         set the number of equilibrium equations to be solved in each          
c         iteration of the load(time) step.                                     
c                                                                               
      neqns = eqn_counter   
      if( local_debug ) then
         write(*,*)'... leaving dof_map ....'
         write(*,*)'      dof,  dof_eqn_map    eqn_node_map'
         do dof = 1,  ndof*nonode
          write(*,9000) dof, dof_eqn_map(dof), eqn_node_map(dof)
         end do
      end if
c
      return         
 9000 format(10x,i3,2i8 )                                                             
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
c     *  assembly of the symmetric equilibrium equations in sparse   *          
c     *  format. we use a structure row-by-row assembly procedure    *          
c     *  rather than conventional element-by-element. this enables   *          
c     *  threaded parallel assembly: a thread builds a complete row  *          
c     *  of the equations. multiple threads build rows concurrently  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *            last modified : 6/20/2016 rhd                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine assem_by_row( neqns, num_threads, eqn_node_map,                
     &                         dof_eqn_map, k_diag, k_coeffs,                   
     &                         k_indexes, k_ptrs, iprops, dcp,                  
     &                         noelem )                                         
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer ::  eqn_node_map(*), dof_eqn_map(*), k_ptrs(*),                   
     &            k_indexes(*), iprops(mxelpr,*), dcp(*)                        
      integer :: neqns, num_threads, noelem                                     
      double precision :: k_diag(*), k_coeffs(*)                                
c                                                                               
c                    local declarations                                         
c                                                                               
      double precision, parameter :: zero = 0.d0                                                 
      double precision, allocatable ::  coeff_row(:,:)                          
      integer :: i, srow, now_thread                                            
      integer, external :: omp_get_thread_num                                   
      integer, allocatable :: row_start_index(:), edest(:,:,:)                  
      integer :: thread_previous_node(max_threads)   
      logical, parameter :: check_values = .false.                           
c                                                                               
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
      coeff_row = zero                                                          
      row_start_index(1) = 1                                                    
c                                                                               
      do i = 2, neqns                                                           
       row_start_index(i) = row_start_index(i-1) + k_ptrs(i-1)                  
      end do                                                                    
c                                                                               
      allocate( edest(mxedof,mxconn,num_threads) )                              
      thread_previous_node = 0 ! vector                                         
c                                                                               
c                 loop over all structural equations. the row by row            
c                 assembly makes threaded parallelism reasonably simple.        
c                 Code above handles making private copies                      
c                 of arrays for threads above.                                  
c     
      call omp_set_dynamic( .false. )    
      if( check_values ) call estiff_allocate( 4 )                                                                           
c$OMP PARALLEL DO PRIVATE( srow, now_thread ) ! all else shared                 
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
c     *  assembly a single row (symmetric) of the equilibrium eqns   *          
c     *  in sparse format. this code is designed to be run inside    *          
c     *  a threaded loop over all equations.                         *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *           last modified : 05/03/2016 rhd                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine assem_a_row( srow, neqns, eqn_node_map, coeff_row,             
     &                        dof_eqn_map, k_diag, k_coeffs,                    
     &                        k_indexes, k_ptrs, iprops,  dcp,                  
     &                        noelem, row_start_index, edest,                   
     &                        previous_snode )         
      use global_data, only : out                         
      use elem_block_data, only : estiff_blocks                                 
      use main_data,       only : elems_to_blocks, repeat_incid,                
     &                            inverse_incidences                            
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer :: srow, neqns, noelem, previous_snode                            
      integer :: eqn_node_map(*), dof_eqn_map(*),                               
     &           k_ptrs(*), k_indexes(*),                                       
     &           iprops(mxelpr,*), dcp(*),                                      
     &           row_start_index(*), edest(mxedof,*)                            
      double precision :: k_diag(*), k_coeffs(*), coeff_row(*)                                     
c                                                                               
c                    local declarations                                         
c                                                                               
      integer :: local_scol(mxedof)                                             
      integer :: snode, num_ele_on_snode, j, ele_on_snode, totdof, blk,         
     &           rel_col, start_loc                                             
      double precision, parameter :: zero = 0.d0, k_tol = 1.d-30                                                  
      double precision, dimension(:,:), pointer :: emat                         
c                                                                               
      logical :: repeated                                                       
c                                                                               
c                 get the structure node number corresponding to this           
c                 equation.                                                     
c                 pull the number of elements connected to snode.               
c                 get the destination numbers for elements connected to         
c                 snode. note we keep track of previous                         
c                 structure snode processed for this thread. if the             
c                 same we don't need to get new destination indexes.            
c                 (see edest and previous_snode made thread                     
c                 private by our code in calling routine).                      
c                                                                               
      snode   = eqn_node_map(srow)                                              
      num_ele_on_snode = inverse_incidences(snode)%element_count                
      if( snode .ne. previous_snode ) then                                      
        call get_edest_terms_assemble( edest,                                   
     &     inverse_incidences(snode)%element_list(1), num_ele_on_snode,         
     &      iprops )                                                            
        previous_snode = snode                                                  
      end if                                                                    
c                                                                               
c                 process all elements which contribute stiffness terms         
c                 to this structural equation.                                  
c                                                                               
      do j = 1, num_ele_on_snode                                                
          ele_on_snode = inverse_incidences(snode)%element_list(j)              
          if( ele_on_snode .le. 0 ) cycle                                       
          repeated = repeat_incid(ele_on_snode)                                 
          totdof  = iprops(2,ele_on_snode) * iprops(4,ele_on_snode)             
          blk     = elems_to_blocks(ele_on_snode,1)                             
          rel_col = elems_to_blocks(ele_on_snode,2)    
          if( .not. associated( estiff_blocks(blk)%ptr ) ) then
            write(out,9100) srow, ele_on_snode
            call die_abort
          end if                        
          emat    => estiff_blocks(blk)%ptr                                     
          if( totdof .eq. 24 ) then                                             
             call assem_a_row_24                                                
             cycle                                                              
          end if                                                                
          if( totdof .eq. 30 ) then                                             
             call assem_a_row_30                                                
             cycle                                                              
          end if                                                                
          if( totdof .eq. 36 ) then                                             
             call assem_a_row_36                                                
             cycle                                                              
          end if                                                                
          if( totdof .eq. 60 ) then                                             
             call assem_a_row_60                                                
             cycle                                                              
          end if                                                                
          call assem_a_row_gen                                                  
      end do                                                                    
c                                                                               
c                 copy the diagonal term from coeff_row into k_diag.            
c                                                                               
      k_diag(srow)    = k_diag(srow) + coeff_row(srow)  
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
      if( srow .eq. neqns ) return                                              
c                                                                               
      start_loc = row_start_index(srow)-1                                       
!DIR$ IVDEP                                                                     
      do j = 1, k_ptrs(srow)                                                    
         k_coeffs(start_loc+j) = coeff_row(k_indexes(start_loc+j))              
         coeff_row(k_indexes(start_loc+j)) = zero                               
      end do                                                                    
c                                                                               
      return     
c
 9000 format('>> FATAL ERROR: assem_a_row. < = 0 diagonal for node #,',        
     &  /,   '                equation #: ',2i10,' value: ',d20.8,                                     
     &  /,   '                job terminated' )  
 9100 format('>> FATAL ERROR: assem_a_row. bad block ptr. srow,',        
     &  /,   '                element: ',2i10,                                  
     &  /,   '                job terminated' )  

                                                                     
                                                                                
      contains                                                                  
c     ****************************************************************          
c     *                                                              *          
c     *  support routines for symmetric assembly of a single         *          
c     *  row. specific for most common elements to get a bit more    *          
c     *  efficiency in vectorization.                                *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
      subroutine assem_a_row_24                                                 
      implicit none                                                             
c                                                                               
      integer :: erow, ecol, scol, kk                                           
      double precision :: ekterm                                                
c                                                                               
      local_scol(1:24) =  dof_eqn_map(edest(1:24,j))                            
c                                                                               
      do erow = 1, 24                                                           
       if( local_scol(erow) .ne. srow ) cycle                                   
!DIR$ IVDEP                                                                     
       do ecol = 1, 24                                                          
         scol = local_scol(ecol) ! dof_eqn_map(edest(ecol,j))                   
         if ( scol .lt. srow ) cycle ! lower triange                            
         kk = dcp(max0(ecol,erow))-iabs(ecol - erow)                            
         ekterm = estiff_blocks(blk)%ptr(kk,rel_col)    
         coeff_row(scol) = coeff_row(scol) + ekterm                             
       end do                                                                   
       if( .not. repeated ) return                                              
      end do                                                                    
      return                                                                    
      end subroutine assem_a_row_24                                             
                                                                                
      subroutine assem_a_row_30                                                 
      implicit none                                                             
c                                                                               
      integer :: erow, ecol, scol, kk                                           
      double precision :: ekterm                                                
c                                                                               
      local_scol(1:30) =  dof_eqn_map(edest(1:30,j))                            
c                                                                               
      do erow = 1, 30                                                           
       if( local_scol(erow) .ne. srow ) cycle                                   
!DIR$ IVDEP                                                                     
       do ecol = 1, 30                                                          
         scol = local_scol(ecol) ! dof_eqn_map(edest(ecol,j))                   
         if( scol .lt. srow ) cycle                                             
         kk = dcp(max0(ecol,erow))-iabs(ecol - erow)                            
         ekterm = emat(kk,rel_col)                                              
         coeff_row(scol) = coeff_row(scol) + ekterm                             
       end do                                                                   
       if( .not. repeated ) return                                              
      end do                                                                    
      return                                                                    
      end subroutine assem_a_row_30                                             
                                                                                
      subroutine assem_a_row_36                                                 
      implicit none                                                             
c                                                                               
      integer :: erow, ecol, scol, kk                                           
      double precision :: ekterm                                                
c                                                                               
      local_scol(1:36) =  dof_eqn_map(edest(1:36,j))                            
c                                                                               
      do erow = 1, 36                                                           
       if( local_scol(erow) .ne. srow ) cycle                                   
!DIR$ IVDEP                                                                     
       do ecol = 1, 36                                                          
         scol = local_scol(ecol) ! dof_eqn_map(edest(ecol,j))                   
         if( scol .lt. srow ) cycle                                             
         kk = dcp(max0(ecol,erow))-iabs(ecol - erow)                            
         ekterm = emat(kk,rel_col)                                              
         coeff_row(scol) = coeff_row(scol) + ekterm                             
       end do                                                                   
       if( .not. repeated ) return                                              
      end do                                                                    
      return                                                                    
      end subroutine assem_a_row_36                                             
                                                                                
      subroutine assem_a_row_60                                                 
      implicit none                                                             
c                                                                               
      integer :: erow, ecol, scol, kk                                           
      double precision :: ekterm                                                
c
      local_scol(1:60) =  dof_eqn_map(edest(1:60,j))                     
c                                                                                
      do erow = 1, 60                                                           
       if( local_scol(erow) .ne. srow ) cycle                                   
!DIR$ IVDEP                                                                     
       do ecol = 1, 60                                                          
         scol = local_scol(ecol) ! dof_eqn_map(edest(ecol,j))                   
         if( scol .lt. srow ) cycle                                             
         kk = dcp(max0(ecol,erow))-iabs(ecol - erow)                            
         ekterm = emat(kk,rel_col)                                              
         coeff_row(scol) = coeff_row(scol) + ekterm                             
       end do                                                                   
       if( .not. repeated ) return                                              
      end do                                                                    
      return                                                                    
      end subroutine assem_a_row_60                                             
                                                                                
      subroutine assem_a_row_gen                                                
      implicit none                                                             
c                                                                               
      integer :: erow, ecol, scol, kk                                           
      double precision :: ekterm                                                
c                                                                               
      local_scol(1:totdof) =  dof_eqn_map(edest(1:totdof,j))                    
c                                                                               
      do erow = 1, totdof                                                       
       if( local_scol(erow) .ne. srow ) cycle                                   
!DIR$ IVDEP                                                                     
        do ecol = 1, totdof                                                     
          scol = local_scol(ecol) ! dof_eqn_map(edest(ecol,j))                  
          if( scol .lt. srow ) cycle                                            
          kk = dcp(max0(ecol,erow))-iabs(ecol - erow)                           
          ekterm = emat(kk,rel_col)                                             
          coeff_row(scol) = coeff_row(scol) + ekterm                            
        end do                                                                  
        if( .not. repeated ) return                                             
      end do                                                                    
      return                                                                    
      end subroutine assem_a_row_gen                                            
                                                                                
      end subroutine assem_a_row                                                
                                                                                
                                                                                
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
c     *  sparsity.                                                   *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *            last modified : 6/23/2016 rhd                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine assem_by_row_asymmetric(                                       
     &  neqns, num_threads, eqn_node_map, dof_eqn_map, k_ptrs,                  
     &  k_indexes, iprops, k_coeffs )                                           
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                    parameter declarations                                     
c                                                                               
      integer :: neqns, num_threads, eqn_node_map(*), dof_eqn_map(*),           
     &           k_ptrs(*), k_indexes(*), iprops(mxelpr,*)                      
      double precision :: k_coeffs(*)                                           
c                                                                               
c                    local declarations                                         
c                                                                               
      integer, allocatable :: edest(:,:,:)                                      
      integer :: srow, now_thread, num_indexes_for_srow                         
      integer :: thread_previous_node(max_threads)                              
      integer, external :: omp_get_thread_num                                   
      logical, parameter :: print_cpu_time = .false.                            
      real :: start_time, stop_time                                             
c                                                                               
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
      allocate( edest(mxedof,mxconn,num_threads) )                              
      thread_previous_node = 0 ! vector                                         
c                                                                               
      if( print_cpu_time ) call cpu_time( start_time )                          
c                                                                               
c$OMP PARALLEL DO  PRIVATE(srow, now_thread, num_indexes_for_srow)              
      do srow = 1, neqns                                                        
        now_thread = omp_get_thread_num() + 1                                   
        num_indexes_for_srow = k_ptrs(srow+1) - k_ptrs(srow)                    
        call assem_a_row_asym( srow, eqn_node_map,                              
     &    dof_eqn_map, iprops, num_indexes_for_srow,                            
     &    k_indexes(k_ptrs(srow)), edest(1,1,now_thread),                       
     &    k_coeffs(k_ptrs(srow)), thread_previous_node(now_thread) )            
      end do                                                                    
c$OMP END PARALLEL DO                                                           
      if( print_cpu_time ) then                                                 
          call cpu_time( stop_time )                                            
          write(*,* ) '.... assembly time: ', stop_time - start_time            
      end if                                                                    
      deallocate( edest )                                                       
                                                                                
      return                                                                    
                                                                                
      end subroutine                                                            
c     ****************************************************************          
c     *                                                              *          
c     *  assemble a single row of the equilibrium equations          *          
c     *  in sparse format. this code is designed to be run inside    *          
c     *  a threaded loop over all equations.                         *          
c     *                                                              *          
c     *  This version is for full, asymmetric assembly               *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 6/23/2016 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine assem_a_row_asym( srow, eqn_node_map, dof_eqn_map,             
     &                             iprops, num_k_indexes, k_indexes,            
     &                             edest, k_coeffs, previous_snode)             
      use elem_block_data, only : estiff_blocks                                 
      use main_data, only: elems_to_blocks, repeat_incid,                       
     &                     inverse_incidences                                   
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                 parameter declarations                                        
c                                                                               
      integer :: srow, num_k_indexes, previous_snode                            
      integer :: dof_eqn_map(*), eqn_node_map(*), iprops(mxelpr,*),             
     &           edest(mxedof,*), k_indexes(num_k_indexes)                      
      double precision :: k_coeffs(num_k_indexes)                               
c                                                                               
c                 local declarations                                            
c                                                                               
      integer :: snode, num_ele_on_snode, e, ele_on_snode, totdof, blk,         
     &           rel_col, erow, ecol, scol, k, ekrow, bs_start,                 
     &           bs_finish, bs_range                                            
      double precision, dimension(:,:), pointer :: emat                         
c                                                                               
c                 - get structure node number for this equation                 
c                 - number of element connected to snode                        
c                 - get structure dof numbers for each element                  
c                   connected to snode (sdof, not eqn #s)                       
c                                                                               
      snode            = eqn_node_map(srow)                                     
      num_ele_on_snode = inverse_incidences(snode)%element_count                
      if( snode .ne. previous_snode ) then                                      
        call get_edest_terms( edest,                                            
     &      inverse_incidences(snode)%element_list(1),                          
     &      num_ele_on_snode, iprops)                                           
        previous_snode = snode                                                  
      end if                                                                    
c                                                                               
c                 - loop over element connected to this snode                   
c                   - skip killed elements                                      
c                   - get corresponding element blok and column                 
c                     within that block                                         
c                   - use local pointer to block of element                     
c                     stffness matrices for more efficent access                
c                   - find which row of element K corresponds to the            
c                     equation number being processed. uses linear              
c                     search since no ordered index available                   
c                   - process each column of Ke on row erow.                    
c                   - find which col of equations to add the Ke term            
c                                                                               
c                 we use two algortihms here. for an element with               
c                 repeated structure nodes, a linear search                     
c                 proves to be fastest.                                         
c                                                                               
c                 otherwise we binary search the sorted k_indexes               
c                 for the equation column that matches element                  
c                 column.                                                       
c                                                                               
      do e = 1, num_ele_on_snode                                                
        ele_on_snode = inverse_incidences(snode)%element_list(e)                
        if( ele_on_snode .le. 0) cycle                                          
        totdof = iprops(2,ele_on_snode) * iprops(4,ele_on_snode)                
        blk = elems_to_blocks(ele_on_snode,1)                                   
        rel_col = elems_to_blocks(ele_on_snode,2)                               
        emat => estiff_blocks(blk)%ptr                                          
c                                                                               
        if( repeat_incid(ele_on_snode) ) then                                   
         call assem_a_row_asym_repeated                                         
         cycle                                                                  
        end if                                                                  
c                                                                               
        do erow = 1, totdof ! which row of Ke for this eqn #                    
          if( dof_eqn_map(edest(erow,e)) .ne. srow) cycle                       
          ekrow = ( erow - 1 ) * totdof                                         
          do ecol = 1, totdof                                                   
            scol  = dof_eqn_map(edest(ecol,e))                                  
            ekrow = ekrow + 1                                                   
            if( scol .eq. 0 ) cycle                                             
            bs_start  =  1   ! start in-place binary search                     
            bs_finish = num_k_indexes                                           
            bs_range  = bs_finish - bs_start                                    
            k = (bs_start + bs_finish) / 2                                      
            do while( k_indexes(k) /= scol .and. bs_range >  0)                 
               if( scol > k_indexes(k) ) then                                   
                  bs_start = k + 1                                              
               else                                                             
                  bs_finish = k - 1                                             
               end if                                                           
               bs_range = bs_finish - bs_start                                  
               k = (bs_start + bs_finish) / 2                                   
            end do                                                              
            if( k_indexes(k) /= scol ) call assem_a_row_asym_error              
            k_coeffs(k) = k_coeffs(k) + emat(ekrow,rel_col)                     
          end do ! on ecol                                                      
        end do ! on erow                                                        
      end do ! on e over elements connected to snode                            
c                                                                               
      return                                                                    
c                                                                               
      contains                                                                  
c     ========                                                                  
c                                                                               
      subroutine assem_a_row_asym_error                                         
      use global_data ! old common.main
      implicit none                                                             
c                                                                               
      write(out,9000)                                                           
      write(out,9010) erow, ecol, num_k_indexes, scol                           
      write(out,9020)                                                           
      call die_abort                                                            
      return                                                                    
c                                                                               
 9000 format(/1x,'>>>>> Fatal Error: invalid state in assem_a_row_asym')        
 9010 format(10x,'erow, ecol, num_k_indexes. scol: ',2i3, i7, i7 )              
 9020 format(/1x,'>>>>> Job terminated',//)                                     
c                                                                               
      end subroutine  assem_a_row_asym_error                                    
c                                                                               
      subroutine assem_a_row_asym_repeated  ! and small num_k_indexes           
      implicit none                                                             
c                                                                               
      do erow = 1, totdof ! which row of Ke for this eqn #                      
          if( dof_eqn_map(edest(erow,e)) .ne. srow) cycle                       
          ekrow = ( erow - 1 ) * totdof                                         
          do ecol = 1, totdof                                                   
            scol  = dof_eqn_map(edest(ecol,e))                                  
            ekrow = ekrow + 1                                                   
            if( scol .eq. 0 ) cycle                                             
            do k = 1, num_k_indexes ! which ecol if any. may repeat             
              if( k_indexes(k) .eq. scol ) then                                 
                k_coeffs(k) = k_coeffs(k) + emat(ekrow,rel_col)                 
              end if                                                            
            end do ! on k                                                       
          end do ! on ecol                                                      
      end do ! on erow                                                          
c                                                                               
      return                                                                    
      end subroutine  assem_a_row_asym_repeated                                 
      end subroutine  assem_a_row_asym                                          
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *           extract indexes into the global                    *          
c     *           displacmenent, velocity, acceleration vectors      *          
c     *           for equations for a list of elements. this         *          
c     *           essentially de-blocks the indexes for a set of     *          
c     *           elements                                           *          
c     *                                                              *          
c     *                       written by  : rhd                      *          
c     *                   last modified : 8/10/2017 rhd              *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
                                                                                
      subroutine get_edest_terms_assemble( table, elem_list,                    
     &                                     list_length, iprops )                
      use elem_block_data, only :  edest_blocks                                 
      use main_data,       only :  elems_to_blocks                              
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
      integer :: table(mxedof,mxconn), elem_list(*), list_length,                    
     &           iprops(mxelpr,*)                                               
      integer, dimension (:,:), pointer :: edest                                
c                                                                               
      integer :: i, elem, totdof, blk, rel_elem, dof                            
c        
c              could zero table for clarity since new values may
c              not overwrite all old values - but accesses into 
c              table know their sizes.
c
      do i = 1, list_length                                                     
       elem = elem_list(i)   ! absolute element number                          
       if( elem .le. 0 ) cycle                                                  
       totdof   = iprops(2,elem) * iprops(4,elem)                               
       blk      = elems_to_blocks(elem,1)                                       
       rel_elem = elems_to_blocks(elem,2)                                       
       edest    => edest_blocks(blk)%ptr                                                                                                         
       table(1:totdof,i) = edest(1:totdof,rel_elem)    
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
