c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *   this routine modifies the stiffness matrix to enforce the  *          
c     *  multi-point constraint equations either entered by the user *          
c     *          or resulting from the tied-contact processor        *          
c     *                                                              *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                      last modifed : 2/19/2016 rhd            *          
c     *                                                              *          
c     * edit: add optional checks for uninitialized variables in     *          
c     * k_coeffs (through its various resizings)                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_insert_terms(neqns, k_ptrs, k_diag, dstmap,               
     &                             dof_eqn_map)                                 
c                                                                               
      use mod_mpc, only : nmpc, num_tied_con_mpc, num_user_mpc                  
      implicit integer (a-z)                                                    
      integer, allocatable, dimension(:) :: abs_ptr, abs_trm                    
      real  dumr                                                                
      double precision  dumd                                                    
      double precision                                                          
     &          k_diag                                                          
      character(len=1) :: dums                                                  
      dimension  k_ptrs(*), k_diag(*), dstmap(*), dof_eqn_map(*)                
      logical ldebug                                                            
      data ldebug / .false. /                                                   
c                                                                               
c        allocate local variables                                               
c                                                                               
      if( ldebug ) write(*,*) '...@ 1 mpc_insert_terms ...'                     
      nmpc = num_tied_con_mpc + num_user_mpc                                    
      allocate ( abs_ptr(neqns),                                                
     &           abs_trm(nmpc),                                                 
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
      if( ldebug ) write(*,*) '...@ 2 mpc_insert_terms ...'                     
c                                                                               
c        1. Create new data structures for mpcs and expanded eqns               
c                                                                               
      call mpc_chk_nan( 10 )                                                    
      call mpc_new_structures(neqns, max_dep, abs_trm, dstmap,                  
     &                        dof_eqn_map)                                      
      call mpc_chk_nan( 11 )                                                    
      if( ldebug ) write(*,*) '...@ 3 mpc_insert_terms ...'                     
c                                                                               
c        2. Find and re-organize the full dependent equations and               
c              all the new non-zero terms to be added to [K]                    
c                                                                               
      call mpc_find_dep_terms(neqns, k_ptrs, abs_trm, max_dep, max_len)         
      call mpc_chk_nan( 12 )                                                    
      if( ldebug ) write(*,*) '...@ 4 mpc_insert_terms ...'                     
c                                                                               
c        3. Sort and copy new terms back into old data structures               
c                                                                               
      call mpc_copy_new_terms(neqns, k_ptrs, abs_ptr)                           
      call mpc_chk_nan( 13 )                                                    
      if( ldebug ) write(*,*) '...@ 5 mpc_insert_terms ...'                     
c                                                                               
c        4. Find the locations in k_coeffs for modification routine             
c                                                                               
      call mpc_find_locations(neqns, k_ptrs, k_diag, abs_ptr, max_len)          
      call mpc_chk_nan( 14 )                                                    
      if( ldebug ) write(*,*) '...@ 6 mpc_insert_terms ...'                     
c                                                                               
c        deallocate temporary memory                                            
c                                                                               
      deallocate (abs_ptr,abs_trm)                                              
      if( ldebug ) write(*,*) '... leaving mpc_insert_terms ...'                
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *   run a check on contents of k_coeffs for NaN entry.         *          
c     *   useful for debugging                                       *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                   last modified : 01/19/04                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mpc_chk_nan( iwhere )                                          
      use stiffness_data, only : k_coeffs, dep_locations, ncoeff,               
     &                           ind_locations, diag_locations                  
      implicit none                                                             
      integer :: i, iwhere                                                      
c                                                                               
      return                                                                    
      write(*,*) '... checking NaN @ ', iwhere                                  
      do i = 1, size(k_coeffs)                                                  
        if( .not. isnan( k_coeffs(i) ) ) cycle                                  
          write(*,*) '  NaN found. i: ', i                                      
          call die_abort                                                        
      end do                                                                    
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *   this routine modifies the stiffness matrix to enforce the  *          
c     *  multi-point constraint equations either entered by the user *          
c     *          or resulting from the tied-contact processor        *          
c     *                                                              *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                   last modified : 11/12/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_new_structures( neqns, max_dep, abs_trm, dstmap,          
     &                                dof_eqn_map )                             
c                                                                               
      use mod_mpc, only : num_tied_con_mpc, tied_con_mpc_table, nmpc,           
     &                    num_user_mpc, user_mpc_table, dep_check,              
     &                    dep_dof, ind_dof, num_terms, multi_list               
      implicit integer (a-z)                                                    
      real  dumr                                                                
      double precision  dumd                                                    
      character(len=1) :: dums                                                  
      logical  last_good                                                        
      logical, parameter :: local_debug = .false.                               
      dimension  abs_trm(*), dstmap(*), dof_eqn_map(*)                          
c                                                                               
c        allocate module variables                                              
c                                                                               
      if (allocated(dep_check))   deallocate(dep_check)                         
      if (allocated(dep_dof)) deallocate(dep_dof)                               
      if (allocated(ind_dof)) deallocate(ind_dof)                               
      if (allocated(num_terms)) deallocate(num_terms)                           
      allocate ( dep_check(neqns),                                              
     &           dep_dof(nmpc),                                                 
     &           num_terms(nmpc),                                               
     &           ind_dof(nmpc*8),                                               
     &           multi_list(nmpc*8),                                            
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        reorganize tied_mesh and user mpc equations into simpler               
c        structures                                                             
c                                                                               
      dep_check = 0  ! vector                                                   
      max_dep = 0                                                               
      skip = 0                                                                  
      ptr = 0                                                                   
      abs_trm(1) = 1                                                            
      last_good = .true.                                                        
c                                                                               
      nxt_mpc: do mpc = 1, nmpc                                                 
         if (mpc .le. num_tied_con_mpc) then                                    
            ntrms = tied_con_mpc_table(mpc)%num_terms                           
            do trm = 1, ntrms                                                   
               node = tied_con_mpc_table(mpc)%node_list(trm)                    
               dof  = tied_con_mpc_table(mpc)%dof_list(trm)                     
               sdof = dstmap(node) + dof - 1                                    
               eqn  = dof_eqn_map(sdof)                                         
               if (trm .eq. 1) then                                             
                  dep = eqn                                                     
                  if (dep .gt. max_dep)  max_dep = dep                          
                  dep_dof(mpc) = dep                                            
                  dep_check(dep) = mpc                                          
                  num_terms(mpc) = ntrms - 1                                    
                  if (mpc .gt. 1) then                                          
                     abs_trm(mpc) = abs_trm(mpc-1) +                            
     &                              num_terms(mpc-1) + skip                     
                  end if                                                        
                  cycle                                                         
               end if                                                           
               if (eqn .eq. dep) then                                           
                  dep_dof(mpc)   = 0                                            
                  dep_check(dep) = 0                                            
                  num_terms(mpc) = 0                                            
                  if (mpc .gt. 1) then                                          
                     abs_trm(mpc) = abs_trm(mpc-1)                              
                  end if                                                        
                  ptr = ptr - (trm - 2)                                         
                  if (last_good)  skip = num_terms(mpc-1)                       
                  last_good = .false.                                           
                  call errmsg2(65,mpc,dums,dumr,dumd)                           
                  cycle nxt_mpc                                                 
               end if                                                           
               if (trm .eq. ntrms) then                                         
                  skip = 0                                                      
                  last_good = .true.                                            
               end if                                                           
               ptr = ptr + 1                                                    
               ind_dof(ptr) = eqn                                               
               multi_list(ptr) =                                                
     &            tied_con_mpc_table(mpc)%multiplier_list(trm)                  
            end do                                                              
         else                                                                   
            pnt = mpc - num_tied_con_mpc                                        
            ntrms = user_mpc_table(pnt)%num_terms                               
            do trm = 1, ntrms                                                   
               node = user_mpc_table(pnt)%node_list(trm)                        
               dof  = user_mpc_table(pnt)%dof_list(trm)                         
               sdof = dstmap(node) + dof - 1                                    
               eqn  = dof_eqn_map(sdof)                                         
               if (trm .eq. 1) then                                             
                  dep = eqn                                                     
                  if (dep .gt. max_dep)  max_dep = dep                          
                  dep_dof(mpc) = dep                                            
                  dep_check(dep) = mpc                                          
                  num_terms(mpc) = ntrms - 1                                    
                  if (mpc .gt. 1) then                                          
                     abs_trm(mpc) = abs_trm(mpc-1) +                            
     &                              num_terms(mpc-1) + skip                     
                  end if                                                        
                  cycle                                                         
               end if                                                           
               if (eqn .eq. dep) then                                           
                  dep_dof(mpc)   = 0                                            
                  dep_check(dep) = 0                                            
                  num_terms(mpc) = 0                                            
                  if (mpc .gt. 1) then                                          
                     abs_trm(mpc) = abs_trm(mpc-1)                              
                  end if                                                        
                  ptr = ptr - (trm - 2)                                         
                  if (last_good)  skip = num_terms(mpc-1)                       
                  last_good = .false.                                           
                  call errmsg2(65,mpc,dums,dumr,dumd)                           
                  cycle nxt_mpc                                                 
               end if                                                           
               if (trm .eq. ntrms) then                                         
                  skip = 0                                                      
                  last_good = .true.                                            
               end if                                                           
               ptr = ptr + 1                                                    
               ind_dof(ptr) = eqn                                               
               multi_list(ptr) =                                                
     &            user_mpc_table(pnt)%multiplier_list(trm)                      
            end do                                                              
         end if                                                                 
      end do  nxt_mpc                                                           
c                                                                               
      if( local_debug ) then                                                    
        write(*,*) " .... mpc_new_structures ...."                              
        write(*,*) " "                                                          
        write(*,*) "       nmpc, max_dep, ptr: ",nmpc,max_dep,ptr               
        write(*,*)                                                              
        write(*,*) "       ... dep_check ..."                                   
        write(*,9000) ( i, dep_check(i), i = 1, neqns)                          
        write(*,*) " "                                                          
        write(*,*) "       ... dep_dof ..."                                     
        write(*,9000) ( i, dep_dof(i), i = 1, nmpc)                             
        write(*,*) " "                                                          
        write(*,*) "       ... num_terms ..."                                   
        write(*,9000) ( i, num_terms(i), i = 1, nmpc)                           
        write(*,*) " "                                                          
        write(*,*) "       ... ind_dof ..."                                     
        write(*,9000) ( i, ind_dof(i), i = 1, 8*nmpc)                           
        write(*,*) " "                                                          
        write(*,*) "       ... multi_list ..."                                  
        write(*,9010) ( i, multi_list(i), i = 1, 8*nmpc)                        
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(2x,8i10)                                                           
 9010 format(2x,8(i5,d14.6))                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *    find all terms in the dependent equations that must be    *          
c     *               moved to the independent equations             *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                      last modifed : rhd                      *          
c     *                                                              *          
c     *                   last modified : 11/12/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_find_dep_terms( neqns, k_ptrs, abs_trm, max_dep,          
     &                                max_len )                                 
c                                                                               
      use mod_mpc, only : tied_con_mpc_table, user_mpc_table, nmpc,             
     &                    eqn_row, dep_check, dep_dof, ind_dof,                 
     &                    num_terms, dep_ptr, num_dep_trms, dep_trms,           
     &                    abs_dep_ptr, dep_trms_len                             
      use stiffness_data, only : new_ptrs, k_indexes                            
      implicit integer (a-z)                                                    
      integer, allocatable, dimension(:) :: tmp_trms, eqn_tmp                   
      real  dumr                                                                
      double precision  dumd                                                    
      character(len=1) :: dums                                                  
      dimension  k_ptrs(*), abs_trm(*)                                          
      intrinsic size                                                            
      logical, parameter :: local_debug = .false.                               
c                                                                               
c                                                                               
      if( local_debug ) write(*,9000)                                           
c                                                                               
c        allocate local variables                                               
c                                                                               
      if(allocated(new_ptrs) )      deallocate(new_ptrs)                        
      if(allocated(num_dep_trms) )  deallocate(num_dep_trms)                    
      if(allocated(dep_trms) )      deallocate(dep_trms)                        
c                                                                               
      if( allocated(eqn_tmp) ) deallocate(eqn_tmp)                              
      if( allocated(dep_ptr) ) deallocate(dep_ptr)                              
      if( allocated(abs_dep_ptr) ) deallocate(abs_dep_ptr)                      
      if( allocated(tmp_trms) ) deallocate(tmp_trms)                            
c                                                                               
c        number of equations could be changing in                               
c        model. discard data vectors from ptrs then                             
c        allocated vector itself. (neqns will be number of                      
c        new equations).                                                        
c                                                                               
      if ( allocated(eqn_row) ) then                                            
           jnrows = size(eqn_row,1)                                             
           do i = 1, jnrows                                                     
             if( associated( eqn_row(i)%loc_list ) ) then                       
               deallocate( eqn_row(i)%loc_list, stat=err)                       
               if (err .ne. 0) then                                             
                 call errmsg2(48,dumi,dums,dumr,dumd)                           
                 call die_abort                                                 
               end if                                                           
               nullify(eqn_row(i)%loc_list)                                     
             end if                                                             
           end do                                                               
           deallocate(eqn_row,stat=err)                                         
           if (err .ne. 0) then                                                 
             call errmsg2(48,dumi,dums,dumr,dumd)                               
             call die_abort                                                     
           end if                                                               
      end if                                                                    
c                                                                               
      if( local_debug ) write(*,*) '   @ 2. neqns: ', neqns                     
c                                                                               
      dep_trms_len = neqns                                                      
      allocate( eqn_row(neqns),                                                 
     &          eqn_tmp(neqns),                                                 
     &          dep_ptr(neqns),                                                 
     &          new_ptrs(neqns),                                                
     &          num_dep_trms(neqns),                                            
     &          abs_dep_ptr(neqns),                                             
     &          dep_trms(dep_trms_len),                                         
     &          tmp_trms(8),                                                    
     &          stat=err)                                                       
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        size new data structures for eqns to include new non-zero terms        
c        an initial guess is made for the size of these new structures,         
c        but a resizer routine is called if needed to prevent overflow          
c        note that the 'eqn_row' structure stores the entire equation           
c        (row, diag, col) for a dep dof but stores only the row section         
c        of every other eqn                                                     
c                                                                               
      new_ptrs = 0  ! vector                                                    
      dep_ptr  = 0  ! vector                                                    
      do eqn = 1, neqns                                                         
         jsize = max(10,k_ptrs(eqn))                                            
         eqn_row(eqn)%length = jsize                                            
         allocate (eqn_row(eqn)%loc_list(jsize), stat=err)                      
         if (err .ne. 0) then                                                   
            call errmsg2(48,dumi,dums,dumr,dumd)                                
            call die_abort                                                      
         end if                                                                 
      end do                                                                    
c                                                                               
c        find and store the entire dep eqn and new terms                        
c                                                                               
      max_len = 0                                                               
      tmp_idx = 1                                                               
      dep_idx = 1                                                               
c                                                                               
      nxt_row: do row = 1, max_dep                                              
            len = k_ptrs(row)                                                   
            eqn_tmp(1:len) = k_indexes(tmp_idx:tmp_idx+len-1)                   
            tmp_idx = tmp_idx + len                                             
            mpc = dep_check(row)                                                
         if (mpc .gt. 0) then                                                   
            call mpc_store_term(row, row)                                       
            ntrms = num_terms(mpc)                                              
            tmp_ptr = 0                                                         
            trm_ptr = abs_trm(mpc)                                              
            do trm = trm_ptr, trm_ptr+ntrms-1                                   
               ind = ind_dof(trm)                                               
               tmp_ptr = tmp_ptr + 1                                            
               tmp_trms(tmp_ptr) = ind                                          
               call mpc_store_term(row, ind)                                    
            end do                                                              
            do ptr = 1, len                                                     
               col = eqn_tmp(ptr)                                               
               call mpc_store_term(row, col)                                    
               mpc = dep_check(col)                                             
               if (mpc .gt. 0) then                                             
                  ntrms = num_terms(mpc)                                        
                  trm_ptr = abs_trm(mpc)                                        
                  do trm = trm_ptr, trm_ptr+ntrms-1                             
                     ind = ind_dof(trm)                                         
                     call mpc_store_term(row, ind)                              
                  end do                                                        
                  call mpc_store_term(col, row)                                 
                  do trm = 1, tmp_ptr                                           
                     ind = tmp_trms(trm)                                        
                     call mpc_store_term(col, ind)                              
                  end do                                                        
               end if                                                           
            end do                                                              
            dlen = dep_ptr(row)                                                 
            call mpc_heapsort(dlen, eqn_row(row)%loc_list(1) )                  
            if (dlen .gt. max_len)  max_len = dlen                              
            do dptr = 1, dlen                                                   
               eqn = eqn_row(row)%loc_list(dptr)                                
               if (eqn .eq. row)  dep_ptr(row) = dptr                           
               do trm = 1, tmp_ptr                                              
                  ind = tmp_trms(trm)                                           
                  if (eqn .gt. ind) then                                        
                     if (dep_check(ind) .gt. 0) then                            
                        cycle                                                   
                     end if                                                     
                     call mpc_store_term(ind, eqn)                              
                  else if (eqn .lt. ind) then                                   
                     if (dep_check(eqn) .gt. 0) then                            
                        cycle                                                   
                     end if                                                     
                     call mpc_store_term(eqn, ind)                              
                  end if                                                        
               end do                                                           
            end do                                                              
            fin = dep_idx + dlen - 1                                            
            if (fin .gt. dep_trms_len)  call mpc_resize_vector(4)               
c                                                                               
c              store dep eqn in simpler form, deallocate structure              
c                                                                               
            dep_trms(dep_idx:fin) = eqn_row(row)%loc_list(1:dlen)               
            deallocate(eqn_row(row)%loc_list)                                   
            num_dep_trms(row) = dlen                                            
            abs_dep_ptr(row) = dep_idx                                          
            dep_idx = dep_idx + dlen                                            
         else                                                                   
            do ptr = 1, len                                                     
               col = eqn_tmp(ptr)                                               
               if (col .gt. max_dep)  cycle nxt_row                             
               if (dep_check(col) .gt. 0) then                                  
                  call mpc_store_term(col, row)                                 
               end if                                                           
            end do                                                              
         end if                                                                 
      end do  nxt_row                                                           
c                                                                               
c        deallocate temp space                                                  
c                                                                               
      deallocate(tmp_trms,eqn_tmp)                                              
                                                                                
      if( local_debug ) then                                                    
        write(*,*) " .... updated data ...."                                    
        write(*,*) " "                                                          
        write(*,*) "  max_len, tmp_idx, dep_idx: ", max_len, tmp_idx,           
     &                   dep_idx                                                
        write(*,*) "       ... eqn_row data structure..."                       
        do eqn = 1, neqns                                                       
         jsize = eqn_row(eqn)%length                                            
         write(*,9020) eqn, jsize                                               
         if( associated( eqn_row(eqn)%loc_list) ) then                          
            write(*,9025) eqn_row(eqn)%loc_list(1:jsize)                        
         else                                                                   
            write(*,*) '       already deleted'                                 
         end if                                                                 
        end do                                                                  
                                                                                
        write(*,*) "       ... dep_ptr ..."                                     
        write(*,9005) ( i, dep_ptr(i), i = 1, neqns)                            
        write(*,*) " "                                                          
        write(*,*) "       ... new_ptrs ..."                                    
        write(*,9005) ( i, new_ptrs(i), i = 1, neqns)                           
        write(*,*) " "                                                          
        write(*,*) "       ... num_dep_trms ..."                                
        write(*,9005) ( i, num_dep_trms(i), i = 1,neqns)                        
        write(*,*) " "                                                          
        write(*,*) "       ... abs_dep_ptr..."                                  
        write(*,9005) ( i, abs_dep_ptr(i), i = 1, neqns)                        
        write(*,*) " "                                                          
        write(*,*) "       ... dep_trms..."                                     
        write(*,9005 ) ( i, dep_trms(i), i = 1, dep_trms_len)                   
      end if                                                                    
                                                                                
      return                                                                    
c                                                                               
 9000 format(/,2x,'   .... entered  mpc_find_dep_terms ....',//)                
 9005 format(2x,8i10)                                                           
 9010 format(2x,8(i5,d14.6))                                                    
 9020 format(10x,'eqn #, size: ',2i8)                                           
 9025 format(15x,10i6)                                                          
      end                                                                       
c                                                                               
c                                                                               
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *     stores in new data structures the new nonzero terms      *          
c     *                   found in 'find_dep_eqns'                   *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 11/11/2016 rhd            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_store_term( eqn, index )                                  
c                                                                               
      use mod_mpc, only :  eqn_row, dep_check, dep_ptr                          
      use stiffness_data, only : new_ptrs                                       
      implicit none                                                             
c                                                                               
      integer :: eqn, index, mpc, ptr, len                                      
c                                                                               
c        increment counter, check current size of structure                     
c                                                                               
      mpc = dep_check(eqn)                                                      
      if( mpc .gt. 0 ) then                                                     
         ptr = dep_ptr(eqn) + 1                                                 
      else                                                                      
         ptr = new_ptrs(eqn) + 1                                                
      end if                                                                    
      len = eqn_row(eqn)%length                                                 
      if( mpc .gt. 0 ) then                                                     
         dep_ptr(eqn) = ptr                                                     
      else                                                                      
         new_ptrs(eqn) = ptr                                                    
      end if                                                                    
c                                                                               
c        if structure too small, call resizer routine. also                     
c        updates stored length                                                  
c                                                                               
      if( ptr .gt. len ) call mpc_resize_eqn_row( eqn, len )                    
c                                                                               
c        add new term to eqn                                                    
c                                                                               
      eqn_row(eqn)%loc_list(ptr) = index                                        
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *           increases the size of eqn_row(eqn)%loc_list        *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 11/11/2016 rhd            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_resize_eqn_row( eqn, old_len )                            
c                                                                               
      use mod_mpc, only : eqn_row                                               
      implicit none                                                             
                                                                                
      integer :: eqn, old_len, new_len                                          
                                                                                
      integer :: err, dumi                                                      
      integer, allocatable, dimension (:) :: temp_loc                           
      real  dumr                                                                
      double precision  dumd                                                    
      character(len=1) :: dums                                                  
c                                                                               
c        allocate temp space                                                    
c                                                                               
      allocate( temp_loc(old_len), stat=err )                                   
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data to temp space, deallocate, double size, reallocate           
c                                                                               
      temp_loc(1:old_len) = eqn_row(eqn)%loc_list(1:old_len)                    
      deallocate( eqn_row(eqn)%loc_list )                                       
      nullify( eqn_row(eqn)%loc_list )                                          
      new_len = old_len*2                                                       
      allocate( eqn_row(eqn)%loc_list(new_len), stat=err )                      
      if( err .ne. 0 ) then                                                     
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data back into structure, deallocate temp space. no need          
c        to zero end of new space.                                              
c                                                                               
      eqn_row(eqn)%loc_list(1:old_len) = temp_loc(1:old_len)                    
      eqn_row(eqn)%length = new_len                                             
c                                                                               
      deallocate(temp_loc)                                                      
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *      uses a heapsort algorithm to sort an integer array      *          
c     *              then all the duplicates are removed             *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 7/22/03                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_heapsort(n, ra)                                           
c                                                                               
      implicit integer (a-z)                                                    
      dimension  ra(*)                                                          
c                                                                               
      if (n .lt. 2) then                                                        
         return                                                                 
      end if                                                                    
      l  = n/2 + 1                                                              
      ir = n                                                                    
 10   continue                                                                  
      if (l .gt. 1) then                                                        
         l   = l - 1                                                            
         rra = ra(l)                                                            
      else                                                                      
         rra = ra(ir)                                                           
         ra(ir) = ra(1)                                                         
         ir = ir - 1                                                            
         if (ir .eq. 1) then                                                    
            ra(1) = rra                                                         
            goto 30                                                             
         end if                                                                 
      end if                                                                    
      i = l                                                                     
      j = l + l                                                                 
 20   if (j .le. ir) then                                                       
         if (j .lt. ir) then                                                    
            if (ra(j) .lt. ra(j+1))  j = j + 1                                  
         end if                                                                 
         if (rra .lt. ra(j)) then                                               
            ra(i) = ra(j)                                                       
            i = j                                                               
            j = j + j                                                           
         else                                                                   
            j = ir + 1                                                          
         end if                                                                 
         goto 20                                                                
      end if                                                                    
      ra(i) = rra                                                               
      goto 10                                                                   
c                                                                               
c     sort completed, now get rid of duplicates and <= 0                        
c                                                                               
 30   continue                                                                  
      count = 1                                                                 
      do i = 2, n                                                               
         if (ra(i) .ne. ra(i-1)) then                                           
            count = count + 1                                                   
            ra(count) = ra(i)                                                   
         end if                                                                 
      end do                                                                    
      n = count                                                                 
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     * inserts the new non-zero terms into the old data structures  *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 02/19/2016 rhd            *          
c     *                                                              *          
c     * edit: zero end of new k_coeffs so no uninitializeds          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_copy_new_terms( neqns, k_ptrs, abs_ptr )                  
c                                                                               
      use mod_mpc, only : eqn_row, dep_check, dep_ptr, abs_dep_ptr,             
     &                    num_dep_trms, dep_trms                                
      use stiffness_data, only : k_coeffs, k_indexes, ncoeff,                   
     &                           new_locations, new_indexes, temp_len,          
     &                           newcount, new_ptrs, ind_temp,                  
     &                           cof_temp, big_ncoeff,                          
     &                           new_len, new_loc, new_ind                      
      implicit integer (a-z)                                                    
c                                                                               
      integer :: k_ptrs(*), abs_ptr(*)                                          
c                                                                               
      logical, parameter :: local_debug = .false.                               
      integer, intrinsic :: size                                                
      integer, allocatable, dimension(:) :: eqn_tmp, old_ind                    
      real  dumr                                                                
      double precision  dumd                                                    
      double precision,                                                         
     &          allocatable, dimension(:) :: old_cof                            
      character(len=1) :: dums                                                  
c                                                                               
c        allocate temp storage space, copy old data to temp space               
c                                                                               
      if( local_debug ) write(*,*) ' ... entered mpc_copy_new_terms '           
c                                                                               
      eqn_len = max(100,neqns/10)                                               
      new_len = max(300,neqns*2)                                                
      allocate ( eqn_tmp(eqn_len),                                              
     &           new_loc(new_len),                                              
     &           new_ind(new_len),                                              
     &           old_ind(ncoeff+neqns),                                         
     &           old_cof(ncoeff+neqns),                                         
     &           ind_temp(ncoeff+neqns),                                        
     &           cof_temp(ncoeff+neqns),                                        
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
      temp_len = ncoeff + neqns                                                 
c                                                                               
      if( local_debug ) then  ! some sanity cehcks                              
        write(*,*) '     neqns, ncoeff: ',neqns, ncoeff                         
        write(*,*) '     eqn_len, new_len, temp_len: ',eqn_len, new_len,        
     &                     temp_len                                             
        write(*,*) '     size k_indexes: ', size(k_indexes)                     
        write(*,*) '     size k_coeffs:  ', size(k_coeffs)                      
      end if                                                                    
c                                                                               
      old_ind = k_indexes ! vector                                              
      old_cof = k_coeffs  ! vector                                              
      call mpc_chk_nan( 15 )                                                    
c                                                                               
c        initialize counters and structures, move data now in 'eqn_row'         
c        structures into stiffness matrix structures being sure to note         
c        which terms are new terms--the new terms are stored to make            
c        subsequent runs considerably faster                                    
c        'eqn_row' no longer needed, so deallocate                              
c                                                                               
      ncoeff   = 0                                                              
      old_idx  = 1                                                              
      ptr_idx  = 1                                                              
      cof_idx  = 1                                                              
      new_idx  = 0                                                              
      ind_temp = 0       ! vector                                               
      cof_temp = 0.0d00  ! vector                                               
      abs_ptr(1) = 1                                                            
c                                                                               
      do eqn = 1, neqns                                                         
         if (dep_check(eqn) .gt. 0) then                                        
            dlen = num_dep_trms(eqn)                                            
            dptr = dep_ptr(eqn)                                                 
            len = dlen-dptr                                                     
            if (len .gt. eqn_len) then                                          
               deallocate(eqn_tmp)                                              
               allocate(eqn_tmp(len), stat=err)                                 
               if (err .ne. 0) then                                             
                  call errmsg2(48,dumi,dums,dumr,dumd)                          
                  call die_abort                                                
               end if                                                           
               eqn_len = len                                                    
            end if                                                              
            new_ptrs(eqn) = len                                                 
            beg = abs_dep_ptr(eqn) + dptr                                       
            end = beg + len                                                     
            eqn_tmp(1:len) = dep_trms(beg:end)                                  
            old_idx = old_idx + k_ptrs(eqn)                                     
         else                                                                   
            len = new_ptrs(eqn)                                                 
            beg = len + 1                                                       
            num = k_ptrs(eqn)                                                   
            len = len + num                                                     
            if (len .gt. eqn_len) then                                          
               deallocate(eqn_tmp)                                              
               allocate(eqn_tmp(len), stat=err)                                 
               if (err .ne. 0) then                                             
                  call errmsg2(48,dumi,dums,dumr,dumd)                          
                  call die_abort                                                
               end if                                                           
               eqn_len = len                                                    
            end if                                                              
            eqn_tmp(1:beg-1) = eqn_row(eqn)%loc_list(1:beg-1)                   
            deallocate( eqn_row(eqn)%loc_list )                                 
            nullify(eqn_row(eqn)%loc_list)                                      
            eqn_tmp(beg:len) = old_ind(old_idx:old_idx+num-1)                   
            call mpc_heapsort(len, eqn_tmp)                                     
            new_ptrs(eqn) = len                                                 
            old_idx = old_idx + num                                             
         end if                                                                 
         ncoeff = ncoeff + new_ptrs(eqn)                                        
         tmp_idx = 0                                                            
         max_cof = cof_idx + k_ptrs(eqn) - 1                                    
         fin = ptr_idx+len-1                                                    
         if(fin .gt. temp_len) call mpc_resize_tempk(neqns, eqn, len)           
         do ind_idx = ptr_idx, fin                                              
            tmp_idx = tmp_idx + 1                                               
            ind_temp(ind_idx) = eqn_tmp(tmp_idx)                                
            if ( (ind_temp(ind_idx).eq.old_ind(cof_idx)) .and.                  
     &           (cof_idx .le. max_cof) ) then                                  
               cof_temp(ind_idx) = old_cof(cof_idx)                             
               cof_idx = cof_idx + 1                                            
            else                                                                
               new_idx = new_idx + 1                                            
               if (new_idx .gt. new_len) then                                   
                  call mpc_resize_vector(5)                                     
                  call mpc_resize_vector(6)                                     
               end if                                                           
               new_loc(new_idx) = ind_idx                                       
               new_ind(new_idx) = ind_temp(ind_idx)                             
            end if                                                              
         end do                                                                 
         k_ptrs(eqn) = len                                                      
         if (eqn .gt. 1)  abs_ptr(eqn) = abs_ptr(eqn-1) + k_ptrs(eqn-1)         
         ptr_idx = ptr_idx + len                                                
      end do                                                                    
      big_ncoeff = ncoeff                                                       
c                                                                               
c        allocate structures to keep new terms, copy data from temp space,      
c        deallocate temp space                                                  
c        'new_locations' contains the exact locations in the k_indexes &        
c        k_coeffs vectors of all the new terms added                            
c        'new_indexes' contains the index to be added to k_indexes for          
c        each new location                                                      
c        there is no need to keep the coefficient to be added to k_coeffs       
c        because all new coefficients are 0.0                                   
c                                                                               
      deallocate(k_indexes,k_coeffs)                                            
      allocate(k_indexes(ncoeff+neqns),                                         
     &         k_coeffs(ncoeff+neqns),                                          
     &         stat=err)                                                        
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
      if( local_debug ) then  ! some sanity checks                              
        write(*,*) '     neqns, ncoeff: ',neqns, ncoeff                         
        write(*,*) '     eqn_len, new_len, temp_len: ',eqn_len, new_len,        
     &                     temp_len                                             
        write(*,*) '     size k_indexes: ', size(k_indexes)                     
        write(*,*) '     size k_coeffs:  ', size(k_coeffs)                      
      end if                                                                    
c                                                                               
      k_indexes(1:ncoeff) = ind_temp(1:ncoeff)                                  
      k_indexes(ncoeff) = neqns      ! critical. was long standing bug          
c                                                                               
      if( local_debug ) then                                                    
         write(*,*) '    ..... new k_indexes from ind_temp .....'               
         write(*,9000) (i,k_indexes(i), i = 1, ncoeff)                          
      end if                                                                    
c                                                                               
      do i = 1, ncoeff ! check for critical error                               
         if( k_indexes(i) .eq. 0 ) then                                         
             write(*,9010) i                                                    
             call die_abort                                                     
         end if                                                                 
      end do                                                                    
c                                                                               
      k_coeffs(1:ncoeff)  = cof_temp(1:ncoeff)                                  
      k_coeffs(ncoeff+1:) = 0.0d00                                              
c                                                                               
      newcount = new_idx                                                        
      if (allocated(new_locations))  deallocate(new_locations)                  
      if (allocated(new_indexes))    deallocate(new_indexes)                    
      allocate ( new_locations(newcount),                                       
     &           new_indexes(newcount),                                         
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
      new_locations(1:newcount) = new_loc(1:newcount)                           
      new_indexes(1:newcount)   = new_ind(1:newcount)                           
c                                                                               
      deallocate( eqn_tmp, old_ind, old_cof, new_loc, new_ind,                  
     &            ind_temp, cof_temp )                                          
c                                                                               
      return                                                                    
c                                                                               
 9000 format(2x,8i8)                                                            
 9010 format('>> FATAL ERROR: routine mpc_copy_new_terms',                      
     & /,     '                zero entry in k_indexes @ ',i8,                  
     & /,     '                job terminated....')                             
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              increases the size of temp [K] vectors          *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 11/12/2016 rhd            *          
c     *                                                              *          
c     * edit: zero to end of new cof_temp to prevent uninitializeds  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_resize_tempk( neqns, eqn, len )                           
c                                                                               
      use stiffness_data, only : ind_temp, cof_temp, temp_len                   
      implicit none                                                             
                                                                                
      integer :: neqns, eqn, len                                                
                                                                                
      integer :: add, err, dumi, old_len                                        
      integer, allocatable, dimension (:) :: itmp                               
      real  :: dumr                                                             
      double precision ::  dumd                                                 
      double precision,                                                         
     &          allocatable, dimension (:) :: dtmp                              
      character(len=1) :: dums                                                  
      logical, parameter :: local_debug = .false.                               
c                                                                               
c        this routine is used to resize                                         
c                                                                               
c        get length of indicated vector, allocate temp space                    
c                                                                               
      allocate( itmp(temp_len),                                                 
     &          dtmp(temp_len),                                                 
     &          stat=err )                                                      
      if( err .ne. 0 ) then                                                     
         call errmsg2( 48,dumi,dums,dumr,dumd )                                 
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data from indicated vector, deallocate vector                     
c                                                                               
      if( local_debug ) then   ! sanity checks on sizes                         
        write(*,*) '   .... mpc_resize_tempk ....'                              
        write(*,*) ' '                                                          
        write(*,*) '     neqns, eqn, len, temp_len: ',                          
     &                   neqns, eqn, len, temp_len                              
        write(*,*) '     size ind_temp: ', size(ind_temp)                       
        write(*,*) '     size cof_temp: ', size(cof_temp)                       
      end if                                                                    
c                                                                               
      itmp = ind_temp  ! vector                                                 
      dtmp = cof_temp  ! vector                                                 
      deallocate( ind_temp, cof_temp )                                          
c                                                                               
c        increase size of vector, reallocate vector                             
c                                                                               
      add = (neqns-eqn)*len/2                                                   
      old_len = temp_len                                                        
      temp_len = temp_len + add                                                 
      allocate( ind_temp(temp_len+1), ! not sure why Barron has +1              
     &          cof_temp(temp_len+1),                                           
     &          stat=err )                                                      
      if( err .ne. 0 ) then                                                     
         call errmsg2( 48,dumi,dums,dumr,dumd )                                 
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data from temp space, deallocate temp space                       
c                                                                               
      ind_temp(1:old_len) = itmp(1:old_len)                                     
      cof_temp(1:old_len) = dtmp(1:old_len)                                     
      cof_temp(old_len+1:) = 0.0d00                                             
c                                                                               
      deallocate( itmp, dtmp )                                                  
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *            creates data structures to be passed to           *          
c     *                    modify_stiffness routine                  *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 11/04/03                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_find_locations( neqns, k_ptrs, k_diag, abs_ptr,           
     &                                max_len)                                  
c                                                                               
      use mod_mpc, only : dep_ptr, dep_dof, ind_dof, num_terms,                 
     &                    nmpc, dep_trms, dep_coef, num_dep_trms,               
     &                    abs_dep_ptr, dep_trms_len, dep_check                  
      use stiffness_data, only : ncoeff, k_indexes, k_coeffs,                   
     &                           dep_locations, ind_locations,                  
     &                           diag_locations, dep_loc, dep_len,              
     &                           ind_loc, ind_len, dia_loc, dia_len             
      implicit none                                                             
c                                                                               
c                      parameter declarations                                   
c                                                                               
      integer :: neqns, max_len                                                 
      integer :: k_ptrs(*), abs_ptr(*)                                          
      double precision :: k_diag(*)                                             
                                                                                
c                                                                               
c                      local declarations                                       
c                                                                               
      integer, allocatable, dimension(:) :: dep_eqn_tmp                         
      real :: dumr                                                              
      double precision ::  dumd                                                 
      character(len=1) :: dums                                                  
c                                                                               
      integer :: err, mpc, dep, len, dptr, beg, eqn, idx, dep_idx,              
     &           ptr, ntrms,trm, ind, dia_idx, ind_idx, ind_ptr,                
     &           dumi, cnt                                                      
                                                                                
c                                                                               
c                                                                               
c        allocate temp storage space                                            
c        an initial guess is made for the sizes of the temp spaces              
c        for the dep term locations, ind term locations, and the                
c        locations of the terms that update the ind diagonals, but              
c        a resizer routine is used as needed to prevent overflow                
c                                                                               
      if( allocated(dep_coef) )  deallocate( dep_coef )                         
      dep_len = neqns*10                                                        
      ind_len = neqns*10                                                        
      dia_len = neqns                                                           
      allocate( dep_eqn_tmp(max_len),                                           
     &          dep_loc(dep_len),                                               
     &          ind_loc(ind_len),                                               
     &          dia_loc(dia_len),                                               
     &          dep_coef(dep_trms_len),                                         
     &          stat=err )                                                      
      if( err .ne. 0 ) then                                                     
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        intialize counters, loop through mpc equations                         
c                                                                               
      dep_idx = 0                                                               
      ind_idx = 0                                                               
      dia_idx = 0                                                               
      ind_ptr = 1                                                               
c                                                                               
      do mpc = 1, nmpc                                                          
c                                                                               
c           find all the dep term locations for current mpc                     
c                                                                               
         dep  = dep_dof(mpc)                                                    
         if( dep .eq. 0 )  cycle                                                
         len  = num_dep_trms(dep)                                               
         dptr = dep_ptr(dep)                                                    
         beg  = abs_dep_ptr(dep)                                                
         dep_eqn_tmp(1:len) = dep_trms(beg:beg+len-1)                           
         nxt_eqn: do ptr = 1, dptr-1                                            
            eqn = dep_eqn_tmp(ptr)                                              
            do idx = abs_ptr(eqn), abs_ptr(eqn)+k_ptrs(eqn)-1                   
               if( k_indexes(idx) .eq. dep ) then                               
                  dep_idx = dep_idx + 1                                         
                  if( dep_idx .gt. dep_len ) call mpc_resize_vector(1)          
                  dep_loc(dep_idx) = idx                                        
                  dep_coef(dep_idx) = k_coeffs(idx)                             
                  cycle nxt_eqn                                                 
               end if                                                           
            end do                                                              
         end do  nxt_eqn                                                        
         dep_idx = dep_idx + 1                                                  
         if( dep_idx .gt. dep_len )  call mpc_resize_vector(1)                  
         dep_loc(dep_idx) = 0                                                   
         dep_coef(dep_idx) = k_diag(dep)                                        
         do ptr = abs_ptr(dep), abs_ptr(dep)+k_ptrs(dep)-1                      
            dep_idx = dep_idx + 1                                               
            if( dep_idx .gt. dep_len )  call mpc_resize_vector(1)               
            dep_loc(dep_idx) = ptr                                              
            dep_coef(dep_idx) = k_coeffs(ptr)                                   
         end do                                                                 
c                                                                               
c           find all the ind term locations for each ind dof for                
c           current mpc                                                         
c           the diagonal terms are also found in this loop                      
c                                                                               
         ntrms = num_terms(mpc)                                                 
         nxt_trm: do trm = ind_ptr, ind_ptr+ntrms-1                             
            ind = ind_dof(trm)                                                  
            do cnt = 1, len                                                     
               eqn = dep_eqn_tmp(cnt)                                           
               if( ind .gt. eqn ) then                                          
                  do ptr = abs_ptr(eqn), abs_ptr(eqn)+k_ptrs(eqn)-1             
                     if( k_indexes(ptr) .eq. ind ) then                         
                        if (eqn .eq. dep) then                                  
                         dia_idx = dia_idx + 1                                  
                         if( dia_idx.gt.dia_len )                               
     &                     call mpc_resize_vector(3)                            
                         dia_loc(dia_idx) = ptr                                 
                        end if                                                  
                        ind_idx = ind_idx + 1                                   
                        if( ind_idx .gt. ind_len )                              
     &                    call mpc_resize_vector(2)                             
                        ind_loc(ind_idx) = ptr                                  
                        exit                                                    
                     end if                                                     
                  end do                                                        
               else if( ind .eq. eqn ) then                                     
                  ind_idx = ind_idx + 1                                         
                  if( ind_idx .gt. ind_len )  call mpc_resize_vector(2)         
                  ind_loc(ind_idx) = 0                                          
               else if( ind .lt. eqn ) then                                     
                  do ptr = abs_ptr(ind), abs_ptr(ind)+k_ptrs(ind)-1             
                     if( k_indexes(ptr) .eq. eqn ) then                         
                        if( eqn .eq. dep ) then                                 
                         dia_idx = dia_idx + 1                                  
                         if( dia_idx.gt.dia_len )                               
     &                         call mpc_resize_vector(3)                        
                         dia_loc(dia_idx) = ptr                                 
                        end if                                                  
                        ind_idx = ind_idx + 1                                   
                        if (ind_idx .gt. ind_len)                               
     &                  call mpc_resize_vector(2)                               
                        ind_loc(ind_idx) = ptr                                  
                        exit                                                    
                     end if                                                     
                  end do                                                        
               end if                                                           
            end do                                                              
         end do  nxt_trm                                                        
         ind_ptr = ind_ptr + ntrms                                              
      end do                                                                    
c                                                                               
c        allocate storage space of exact size, copy data from temp space,       
c        deallocate temp space                                                  
c                                                                               
      if( allocated(dep_locations) )   deallocate(dep_locations)                
      if( allocated(ind_locations) )   deallocate(ind_locations)                
      if( allocated(diag_locations) )  deallocate(diag_locations)               
      allocate(  dep_locations(dep_idx),                                        
     &           ind_locations(ind_idx),                                        
     &           diag_locations(dia_idx),                                       
     &           stat=err )                                                     
      if( err .ne. 0 ) then                                                     
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
      dep_locations(1:dep_idx)  = dep_loc(1:dep_idx)                            
      ind_locations(1:ind_idx)  = ind_loc(1:ind_idx)                            
      diag_locations(1:dia_idx) = dia_loc(1:dia_idx)                            
c                                                                               
      deallocate(dep_eqn_tmp,dep_ptr,dep_loc,dia_loc,ind_loc,                   
     &           abs_dep_ptr)                                                   
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              increases the size of indicated vector          *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 07/30/2016 rhd            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mpc_resize_vector(vec_ptr)                                     
c                                                                               
      use mod_mpc, only : dep_trms, dep_trms_len                                
      use stiffness_data, only : dep_loc, dep_len, ind_loc, ind_len,            
     &                           dia_loc, dia_len, new_loc, new_ind,            
     &                           new_len                                        
      implicit integer (a-z)                                                    
      integer, allocatable, dimension (:) :: temp_loc                           
      real  dumr                                                                
      double precision  dumd                                                    
      character(len=1) :: dums                                                  
c                                                                               
c        this routine is used to resize several different vectors               
c        the vectors are referenced through the 'vec_ptr' flag                  
c        a '1' indicates the 'dep_loc'  vector                                  
c        a '2' indicates the 'ind_loc'  vector                                  
c        a '3' indicates the 'dia_loc'  vector                                  
c        a '4' indicates the 'dep_trms' vector                                  
c        a '5' indicates the 'new_loc'  vector                                  
c        a '6' indicates the 'new_ind'  vector                                  
c                                                                               
c        get length of indicated vector, allocate temp space                    
c                                                                               
      if (vec_ptr .eq. 1) then                                                  
         len = dep_len                                                          
      else if (vec_ptr .eq. 2) then                                             
         len = ind_len                                                          
      else if (vec_ptr .eq. 3) then                                             
         len = dia_len                                                          
      else if (vec_ptr .eq. 4) then                                             
         len = dep_trms_len                                                     
      else if ((vec_ptr.eq.5) .or. (vec_ptr.eq.6)) then                         
         len = new_len                                                          
      end if                                                                    
      allocate (temp_loc(len), stat=err)                                        
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data from indicated vector, deallocate vector                     
c                                                                               
      if (vec_ptr .eq. 1) then                                                  
         temp_loc = dep_loc                                                     
         deallocate (dep_loc)                                                   
      else if (vec_ptr .eq. 2) then                                             
         temp_loc = ind_loc                                                     
         deallocate (ind_loc)                                                   
      else if (vec_ptr .eq. 3) then                                             
         temp_loc = dia_loc                                                     
         deallocate (dia_loc)                                                   
      else if (vec_ptr .eq. 4) then                                             
         temp_loc = dep_trms                                                    
         deallocate (dep_trms)                                                  
      else if (vec_ptr .eq. 5) then                                             
         temp_loc = new_loc                                                     
         deallocate (new_loc)                                                   
      else if (vec_ptr .eq. 6) then                                             
         temp_loc = new_ind                                                     
         deallocate (new_ind)                                                   
      end if                                                                    
c                                                                               
c        double size of vector, reallocate vector                               
c                                                                               
      len = len*2                                                               
      if (vec_ptr .eq. 1) then                                                  
         allocate (dep_loc(len), stat=err)                                      
      else if (vec_ptr .eq. 2) then                                             
         allocate (ind_loc(len), stat=err)                                      
      else if (vec_ptr .eq. 3) then                                             
         allocate (dia_loc(len), stat=err)                                      
      else if (vec_ptr .eq. 4) then                                             
         allocate (dep_trms(len), stat=err)                                     
      else if (vec_ptr .eq. 5) then                                             
         allocate (new_loc(len), stat=err)                                      
      else if (vec_ptr .eq. 6) then                                             
         allocate (new_ind(len), stat=err)                                      
      end if                                                                    
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data from temp space, deallocate temp space                       
c                                                                               
      if (vec_ptr .eq. 1) then                                                  
         dep_loc(1:len/2) = temp_loc(1:len/2)                                   
         dep_len = len                                                          
      else if (vec_ptr .eq. 2) then                                             
         ind_loc(1:len/2) = temp_loc(1:len/2)                                   
         ind_len = len                                                          
      else if (vec_ptr .eq. 3) then                                             
         dia_loc(1:len/2) = temp_loc(1:len/2)                                   
         dia_len = len                                                          
      else if (vec_ptr .eq. 4) then                                             
         dep_trms(1:len/2) = temp_loc(1:len/2)                                  
         dep_trms_len = len                                                     
      else if (vec_ptr .eq. 5) then                                             
         new_loc(1:len/2) = temp_loc(1:len/2)                                   
      else if (vec_ptr .eq. 6) then                                             
         new_ind(1:len/2) = temp_loc(1:len/2)                                   
         new_len = len                                                          
      end if                                                                    
      deallocate (temp_loc)                                                     
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                modifies the stiffness matrix for             *          
c     *                   tied contact implementation                *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 11/06/03                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_modify_stiffness( neqns, k_diag, p_vec )                  
c                                                                               
      use mod_mpc, only : nmpc, dep_dof, ind_dof, num_terms, multi_list,        
     &                    dep_rhs, num_dep_trms                                 
      use stiffness_data, only : k_coeffs, dep_locations, ncoeff,               
     &                           ind_locations, diag_locations                  
      implicit integer (a-z)                                                    
      real  dumr, mlt                                                           
      double precision  dumd                                                    
      double precision                                                          
     &         dep_trm, k_diag, p_vec                                           
      character(len=1) :: dums                                                  
      dimension  k_diag(*), p_vec(*)                                            
                                                                                
c                                                                               
c        intialize counters, use the terms in the dep_locations to              
c        modify those in the ind_locations                                      
c        use the terms in the diag_locations to modify the ind diag terms       
c                                                                               
      if (allocated(dep_rhs))  deallocate(dep_rhs)                              
      allocate(dep_rhs(neqns), stat=err)                                        
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
      trm_ptr = 1                                                               
      dia_ptr = 1                                                               
      dep_loc_idx = 1                                                           
      ind_loc_idx = 1                                                           
c                                                                               
      do mpc = 1, nmpc                                                          
         dep = dep_dof(mpc)                                                     
         if (dep .eq. 0)  cycle                                                 
         len = num_dep_trms(dep)                                                
         ntrms = num_terms(mpc)                                                 
         do trm = trm_ptr, trm_ptr+ntrms-1                                      
            ind = ind_dof(trm)                                                  
            mlt = multi_list(trm)                                               
            p_vec(ind) = p_vec(ind) + mlt*p_vec(dep)                            
            if (ind .gt. dep) then                                              
               dia_trm = diag_locations(dia_ptr)                                
               k_diag(ind) = k_diag(ind) + mlt*k_coeffs(dia_trm)                
               dia_ptr = dia_ptr + 1                                            
            end if                                                              
            do dep_loc_ptr = dep_loc_idx, dep_loc_idx+len-1                     
               ind_ptr = ind_locations(ind_loc_idx)                             
               dep_ptr = dep_locations(dep_loc_ptr)                             
               if (dep_ptr .eq. 0) then                                         
                  dep_trm = k_diag(dep)                                         
               else                                                             
                  dep_trm = k_coeffs(dep_ptr)                                   
               end if                                                           
               if (ind_ptr .eq. 0) then                                         
                  k_diag(ind) = k_diag(ind) + mlt*dep_trm                       
               else                                                             
c                  write(*,*) ' at critical line'                               
c                  write(*,*) ind_ptr, mlt, dep_trm, size(k_coeffs)             
c                  write(*,*) k_coeffs(ind_ptr)                                 
                  k_coeffs(ind_ptr) = k_coeffs(ind_ptr) + mlt*dep_trm           
               end if                                                           
               ind_loc_idx = ind_loc_idx + 1                                    
            end do                                                              
            if (ind .lt. dep) then                                              
               dia_trm = diag_locations(dia_ptr)                                
               k_diag(ind) = k_diag(ind) + mlt*k_coeffs(dia_trm)                
               dia_ptr = dia_ptr + 1                                            
            end if                                                              
         end do                                                                 
         dep_rhs(dep) = p_vec(dep)                                              
         p_vec(dep) = 0.0d00                                                    
         k_diag(dep) = 1.0d00                                                   
         dep_loc_idx = dep_loc_idx + len                                        
         trm_ptr = trm_ptr + ntrms                                              
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *      inserts the new terms created by the mpc equations      *          
c     *                    for subsequent solves                     *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 01/22/04                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_reinsert_terms(neqns, k_ptrs, dstmap, dof_eqn_map)        
c                                                                               
      use mod_mpc, only : nmpc, num_tied_con_mpc, num_user_mpc,                 
     &                    tied_con_mpc_table, user_mpc_table, dep_dof,          
     &                    ind_dof, num_terms, multi_list                        
      use stiffness_data, only : ncoeff, k_indexes, k_coeffs, newcount,         
     &                           new_locations, new_indexes, new_ptrs           
      implicit integer (a-z)                                                    
      integer, allocatable, dimension (:) :: ind_tmp                            
      real  dumr                                                                
      double precision  dumd                                                    
      double precision,                                                         
     &          allocatable, dimension (:) :: cof_tmp                           
      character(len=1) :: dums                                                  
      logical  last_good                                                        
      dimension  k_ptrs(*), dstmap(*), dof_eqn_map(*)                           
c                                                                               
c        allocate temp storage space, copy data into temp space                 
c                                                                               
      allocate ( ind_tmp(ncoeff),                                               
     &           cof_tmp(ncoeff),                                               
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
      k_ptrs(1:neqns)   = new_ptrs(1:neqns)                                     
      ind_tmp(1:ncoeff) = k_indexes(1:ncoeff)                                   
      cof_tmp(1:ncoeff) = k_coeffs(1:ncoeff)                                    
c                                                                               
c        insert each new_index into its corresponding new_location              
c        deallocate temp space                                                  
c                                                                               
      ptr = 0                                                                   
      do cnt = 1, newcount                                                      
         new_loc = new_locations(cnt)                                           
         new_ind = new_indexes(cnt)                                             
         if (new_loc .eq. ptr+1) then                                           
            k_indexes(new_loc) = new_ind                                        
            k_coeffs(new_loc)  = 0.0d00                                         
            ptr = new_loc                                                       
            cycle                                                               
         end if                                                                 
         do loc = ptr+1, new_loc-1                                              
            k_indexes(loc) = ind_tmp(loc-cnt+1)                                 
            k_coeffs(loc)  = cof_tmp(loc-cnt+1)                                 
         end do                                                                 
         k_indexes(new_loc) = new_ind                                           
         k_coeffs(new_loc)  = 0.0d00                                            
         ptr = new_loc                                                          
      end do                                                                    
      do loc = ptr+1, ncoeff                                                    
         k_indexes(loc) = ind_tmp(loc-cnt+1)                                    
         k_coeffs(loc)  = cof_tmp(loc-cnt+1)                                    
      end do                                                                    
      deallocate(ind_tmp,cof_tmp)                                               
c                                                                               
c        allocate module variables                                              
c                                                                               
      allocate ( dep_dof(nmpc),                                                 
     &           num_terms(nmpc),                                               
     &           ind_dof(nmpc*8),                                               
     &           multi_list(nmpc*8),                                            
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        reorganize tied_mesh and user mpc equations into f77 structures        
c                                                                               
      ptr = 0                                                                   
      nxt_mpc: do mpc = 1, nmpc                                                 
         if (mpc .le. num_tied_con_mpc) then                                    
            ntrms = tied_con_mpc_table(mpc)%num_terms                           
            do trm = 1, ntrms                                                   
               node = tied_con_mpc_table(mpc)%node_list(trm)                    
               dof  = tied_con_mpc_table(mpc)%dof_list(trm)                     
               sdof = dstmap(node) + dof - 1                                    
               eqn  = dof_eqn_map(sdof)                                         
               if (trm .eq. 1) then                                             
                  dep = eqn                                                     
                  dep_dof(mpc) = dep                                            
                  num_terms(mpc) = ntrms - 1                                    
                  cycle                                                         
               end if                                                           
               if (eqn .eq. dep) then                                           
                  dep_dof(mpc)   = 0                                            
                  num_terms(mpc) = 0                                            
                  ptr = ptr - (trm - 2)                                         
                  call errmsg2(65,mpc,dums,dumr,dumd)                           
                  cycle nxt_mpc                                                 
               end if                                                           
               ptr = ptr + 1                                                    
               ind_dof(ptr) = eqn                                               
               multi_list(ptr) =                                                
     &            tied_con_mpc_table(mpc)%multiplier_list(trm)                  
            end do                                                              
         else                                                                   
            pnt = mpc - num_tied_con_mpc                                        
            ntrms = user_mpc_table(pnt)%num_terms                               
            do trm = 1, ntrms                                                   
               node = user_mpc_table(pnt)%node_list(trm)                        
               dof  = user_mpc_table(pnt)%dof_list(trm)                         
               sdof = dstmap(node) + dof - 1                                    
               eqn  = dof_eqn_map(sdof)                                         
               if (trm .eq. 1) then                                             
                  dep = eqn                                                     
                  dep_dof(mpc) = dep                                            
                  num_terms(mpc) = ntrms - 1                                    
                  cycle                                                         
               end if                                                           
               if (eqn .eq. dep) then                                           
                  dep_dof(mpc)   = 0                                            
                  num_terms(mpc) = 0                                            
                  ptr = ptr - (trm - 2)                                         
                  call errmsg2(65,mpc,dums,dumr,dumd)                           
                  cycle nxt_mpc                                                 
               end if                                                           
               ptr = ptr + 1                                                    
               ind_dof(ptr) = eqn                                               
               multi_list(ptr) =                                                
     &            user_mpc_table(pnt)%multiplier_list(trm)                      
            end do                                                              
         end if                                                                 
      end do  nxt_mpc                                                           
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *      inserts the new terms created by the mpc equations      *          
c     *                    for subsequent solves                     *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                                                              *          
c     *                    last modified : 02/19/2016 rhd            *          
c     *                                                              *          
c     *  edit: zero to end of new k_coeffs to prevent uninitialized  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine  mpc_remove_dep_eqns(neqns, k_ptrs)                            
c                                                                               
      use mod_mpc, only : dep_check                                             
      use stiffness_data, only : ncoeff, k_indexes, k_coeffs                    
      implicit integer (a-z)                                                    
      integer, allocatable, dimension(:) :: ind_tmp                             
      real  dumr                                                                
      double precision  dumd                                                    
      double precision,                                                         
     &          allocatable, dimension(:) :: cof_tmp                            
      character(len=1) :: dums                                                  
      logical new_size                                                          
      dimension  k_ptrs(*)                                                      
c                                                                               
c        allocate temp storage space                                            
c                                                                               
      allocate ( ind_tmp(ncoeff),                                               
     &           cof_tmp(ncoeff),                                               
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        remove all references to dep equations--entire row and column          
c        store terms to be kept in temp spaces, calculate new ncoeff            
c                                                                               
      cnt = 0                                                                   
      ptr = 1                                                                   
      do eqn = 1, neqns                                                         
         if (dep_check(eqn) .gt. 0) then                                        
            ptr = ptr + k_ptrs(eqn)                                             
            k_ptrs(eqn) = 0                                                     
         else                                                                   
            num = k_ptrs(eqn)                                                   
            do loc = ptr, ptr+num-1                                             
               idx = k_indexes(loc)                                             
               if (dep_check(idx) .gt. 0) then                                  
                  k_ptrs(eqn) = k_ptrs(eqn) - 1                                 
               else                                                             
                  cnt = cnt + 1                                                 
                  ind_tmp(cnt) = idx                                            
                  cof_tmp(cnt) = k_coeffs(loc)                                  
               end if                                                           
            end do                                                              
            ptr = ptr + num                                                     
         end if                                                                 
      end do                                                                    
      ncoeff = cnt                                                              
c                                                                               
c        deallocate structures, reallocate in new size, initialize              
c                                                                               
      deallocate(k_indexes,k_coeffs)                                            
      allocate ( k_indexes(ncoeff+neqns),                                       
     &           k_coeffs(ncoeff+neqns),                                        
     &           stat=err)                                                      
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        copy data from temp space, deallocate temp space                       
c                                                                               
      k_indexes(1:ncoeff) = ind_tmp(1:ncoeff)                                   
      k_coeffs(1:ncoeff)  = cof_tmp(1:ncoeff)                                   
      k_coeffs(ncoeff+1:) = 0.0d00                                              
      deallocate(ind_tmp,cof_tmp)                                               
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *   modify displacement vector to account for mpc equations    *          
c     *             and calculate lagrange multipliers               *          
c     *                                                              *          
c     *                       written by  : bjb                      *          
c     *                       modified by : rhd                      *          
c     *                                                              *          
c     *                   last modified : 04/19/2015                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mpcs_apply(x, neqns, nodof, cstmap)                            
c                                                                               
      use mod_mpc, only : num_tied_con_mpc, tied_con_mpc_table,                 
     &                    num_user_mpc, user_mpc_table, nmpc,                   
     &                    num_dep_trms, dep_coef, dep_trms,                     
     &                    dep_dof, ind_dof, num_terms, multi_list,              
     &                    dep_rhs                                               
      use stiffness_data, only : i_lagrange_forces                              
      implicit integer (a-z)                                                    
      real :: dumr, mlt                                                         
      double precision :: dumd, zero                                            
      double precision                                                          
     &          x, cof                                                          
      double precision,                                                         
     &          allocatable, dimension (:) :: lagmlt                            
      character(len=1) :: dums                                                  
      dimension  x(*), cstmap(*)                                                
      data zero / 0.0d00 /                                                      
                                                                                
c                                                                               
      allocate( lagmlt(neqns), stat=err)                                        
      if (err .ne. 0) then                                                      
         call errmsg2(48,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
c        calculate displacements for dep dofs                                   
c                                                                               
      ptr = 1                                                                   
      do mpc = 1, nmpc                                                          
         dep = dep_dof(mpc)                                                     
         if (dep .eq. 0)  cycle                                                 
         ntrms = num_terms(mpc)                                                 
         do trm = ptr, ptr+ntrms-1                                              
            ind = ind_dof(trm)                                                  
            mlt = multi_list(trm)                                               
            x(dep) = x(dep) + mlt*x(ind)                                        
         end do                                                                 
         ptr = ptr + ntrms                                                      
      end do                                                                    
c                                                                               
c        calculate lagrange multipliers for dep dofs                            
c                                                                               
      lagmlt = zero                                                             
      ptr = 1                                                                   
      do mpc = 1, nmpc                                                          
         dep = dep_dof(mpc)                                                     
         if (dep .eq. 0)  cycle                                                 
         ntrms = num_dep_trms(dep)                                              
         do trm = ptr, ptr+ntrms-1                                              
            eqn = dep_trms(trm)                                                 
            cof = dep_coef(trm)                                                 
            lagmlt(dep) = lagmlt(dep) + cof*x(eqn)                              
         end do                                                                 
         lagmlt(dep) = lagmlt(dep) - dep_rhs(dep)                               
         ptr = ptr + ntrms                                                      
      end do                                                                    
c                                                                               
c        calculate lagrange multipliers for ind dofs                            
c                                                                               
      ptr = 1                                                                   
      do mpc = 1, nmpc                                                          
         dep = dep_dof(mpc)                                                     
         if (dep .eq. 0)  cycle                                                 
         ntrms = num_terms(mpc)                                                 
         do trm = ptr, ptr+ntrms-1                                              
            ind = ind_dof(trm)                                                  
            mlt = multi_list(trm)                                               
            lagmlt(ind) = lagmlt(ind) - mlt*lagmlt(dep)                         
         end do                                                                 
         ptr = ptr + ntrms                                                      
      end do                                                                    
c                                                                               
c        map lagrange multipliers to global dofs                                
c                                                                               
      eqn = 1                                                                   
      do dof = 1, nodof                                                         
         i_lagrange_forces(dof) = zero                                          
         if( cstmap(dof) .eq. 0 ) then ! no abs constraint                      
            i_lagrange_forces(dof) = lagmlt(eqn)                                
            eqn = eqn + 1                                                       
         end if                                                                 
      end do                                                                    
c                                                                               
c        deallocate space no longer needed                                      
c                                                                               
      deallocate (dep_dof,ind_dof,num_terms,multi_list,dep_rhs,                 
     &            lagmlt)                                                       
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
