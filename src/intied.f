c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine intied                       *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 04/09/03                   *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *               the tied contact pair information              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine intied()                                                       
      use global_data ! old common.main
c                                                                               
      use mod_mpc, only : tied_contact_table                                    
      parameter (max_prs=10)                                                    
      integer  dumi, npairs, err                                                
      integer, allocatable, dimension (:) :: mstr_lst, slv_lst                  
      real     dumr, tied_tol                                                   
      double precision  dumd                                                    
      character(len=16) :: setid, mstrid, slvid                                 
      character(len=1) :: dums                                                  
      logical  label, matchs, numr, true, adj_flag                              
c                                                                               
      allocate (tied_contact_table(max_tied_sets),                              
     &          mstr_lst(max_prs), slv_lst(max_prs), stat=err)                  
      if (err .ne. 0) then                                                      
         call errmsg2(46,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
      new_set: do                                                               
c                                                                               
c        check for name of set, if none found, skip set data                    
c                                                                               
         if (matchs('mesh',4)) call splunj                                      
         if (.not. label(dumi)) then                                            
            call errmsg2(41,dumi,dums,dumr,dumd)                                
            do                                                                  
               call readsc                                                      
               if (matchs('tolerance',5))  cycle                                
               if (matchs('master',6))     cycle                                
               if (matchs('tie',3)) then                                        
                  cycle  new_set                                                
               else                                                             
                  call reset                                                    
                  if (true(dumi))  call splunj                                  
                  deallocate (mstr_lst, slv_lst)                                
                  return                                                        
               end if                                                           
            end do                                                              
         end if                                                                 
c                                                                               
         setid = ' '                                                            
         call entits(setid,len)                                                 
         call check_setid(setid)                                                
c                                                                               
c        set has a name, initialize flags and counters, read set data           
c                                                                               
         adj_flag = .true.                                                      
         npairs   = 0                                                           
         tied_tol = 0.05                                                        
c                                                                               
         set_data: do                                                           
            call readsc                                                         
            if (matchs('tie',3)) then                                           
c                                                                               
c              a new set found, store current set data (if any)                 
c              cycle to look for new set                                        
c                                                                               
               if (npairs .gt. 0) then                                          
                  call intied_store(setid,tied_tol,adj_flag,npairs,             
     &                              mstr_lst,slv_lst)                           
                  cycle  new_set                                                
               else                                                             
                  call errmsg2(53,dumi,setid,dumr,dumd)                         
                  cycle  new_set                                                
               end if                                                           
            end if                                                              
c                                                                               
c           check for tolerance and adjust flag                                 
c           both have default values                                            
c                                                                               
            if (matchs('tolerance',5)) then                                     
               if (.not.numr(tied_tol)) then                                    
                  call errmsg2(49,dumi,dums,tied_tol,dumd)                      
               end if                                                           
               if (matchs('adjust',3)) then                                     
                  if (matchs('off',3)) then                                     
                     adj_flag = .false.                                         
                  end if                                                        
               end if                                                           
               cycle  set_data                                                  
            end if                                                              
c                                                                               
c           read tied contact pairs                                             
c           either slave or master can come first, but must be one of each      
c                                                                               
            if (matchs('master',6)) then                                        
               if (label(dumi)) then                                            
                  mstrid = ' '                                                  
                  call entits(mstrid,len)                                       
               else                                                             
                  call errmsg2(50,dumi,dums,dumr,dumd)                          
                  cycle set_data                                                
               end if                                                           
               if (matchs('slave',5)) then                                      
                  if (label(dumi)) then                                         
                     slvid = ' '                                                
                     call entits(slvid,len)                                     
                  else                                                          
                     call errmsg2(50,dumi,dums,dumr,dumd)                       
                     cycle set_data                                             
                  end if                                                        
                  npairs = npairs + 1                                           
               else                                                             
                  call errmsg2(51,dumi,dums,dumr,dumd)                          
                  cycle set_data                                                
               end if                                                           
               if (npairs .gt. max_prs) then                                    
                  call errmsg2(62,max_prs,dums,dumr,dumd)                       
                  call die_abort                                                
               end if                                                           
               call find_surfs(mstrid,slvid,npairs,mstr_lst,slv_lst)            
               cycle  set_data                                                  
            end if                                                              
            if (matchs('slave',5)) then                                         
               if (label(dumi)) then                                            
                  slvid = ' '                                                   
                  call entits(slvid,len)                                        
               else                                                             
                  call errmsg2(50,dumi,dums,dumr,dumd)                          
                  cycle set_data                                                
               end if                                                           
               if (matchs('master',6)) then                                     
                  if (label(dumi)) then                                         
                     mstrid = ' '                                               
                     call entits(mstrid,len)                                    
                  else                                                          
                     call errmsg2(50,dumi,dums,dumr,dumd)                       
                     cycle set_data                                             
                  end if                                                        
                  npairs = npairs + 1                                           
               else                                                             
                  call errmsg2(52,dumi,dums,dumr,dumd)                          
                  cycle set_data                                                
               end if                                                           
               if (npairs .gt. max_prs) then                                    
                  call errmsg2(62,max_prs,dums,dumr,dumd)                       
                  call die_abort                                                
               end if                                                           
               call find_surfs(mstrid,slvid,npairs,mstr_lst,slv_lst)            
               cycle  set_data                                                  
            end if                                                              
c                                                                               
c           if none of these commands have been found, something other          
c           than a tied contact set has been encountered                        
c           store what has been read (if anything) and return                   
c                                                                               
            if (npairs .gt. 0) then                                             
               call intied_store(setid,tied_tol,adj_flag,npairs,                
     &                           mstr_lst,slv_lst)                              
               return                                                           
            else                                                                
               call errmsg2(53,dumi,setid,dumr,dumd)                            
               return                                                           
            end if                                                              
c                                                                               
         end do  set_data                                                       
      end do  new_set                                                           
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine find_surfs                    *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 05/19/03                   *          
c     *                                                              *          
c     *       this subroutine relates the master and slave id's      *          
c     *          to their index number in the surface table          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine find_surfs(mstrid,slvid,npairs,mstr_lst,slv_lst)               
c                                                                               
      use mod_mpc, only : surface_table, num_surfaces                           
      integer  dumi, npairs, mstr_lst(*), slv_lst(*)                            
      real  dumr                                                                
      double precision  dumd                                                    
      character(len=16) :: mstrid, slvid                                        
      character(len=1) :: dums                                                  
      logical  badid                                                            
c                                                                               
      mstr_lst(npairs) = 0                                                      
      slv_lst(npairs)  = 0                                                      
      badid = .false.                                                           
c                                                                               
      do i = 1, num_surfaces                                                    
         if (surface_table(i)%id .eq. mstrid)  mstr_lst(npairs) = i             
         if (surface_table(i)%id .eq. slvid)   slv_lst(npairs)  = i             
         if (mstr_lst(npairs).ne.0 .and. slv_lst(npairs).ne.0) exit             
      end do                                                                    
c                                                                               
c     an error occurs if:                                                       
c           master = slave                                                      
c           master = 0       no match                                           
c           slave  = 0       no match                                           
c     write an appropriate error message for each, skip the pair                
c     set npairs = npairs - 1 so it is overwritten or doesn't get saved         
c                                                                               
      if ((mstr_lst(npairs) .eq. slv_lst(npairs)) .and.                         
     &    (mstr_lst(npairs) .ne. 0))               then                         
         call errmsg2(54,dumi,dums,dumr,dumd)                                   
         badid = .true.                                                         
      end if                                                                    
      if (mstr_lst(npairs) .eq. 0) then                                         
         call errmsg2(55,dumi,dums,dumr,dumd)                                   
         badid = .true.                                                         
      end if                                                                    
      if (slv_lst(npairs) .eq. 0) then                                          
         call errmsg2(56,dumi,dums,dumr,dumd)                                   
         badid = .true.                                                         
      end if                                                                    
c                                                                               
      if (badid)  npairs = npairs - 1                                           
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine intied_store                   *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 04/09/03                   *          
c     *                                                              *          
c     *       this subroutine stores the tied contact set data       *          
c     *                into the proper data structures               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine intied_store(setid,tied_tol,adj_flag,npairs,                   
     &                        mstr_lst,slv_lst)                                 
      use global_data ! old common.main
c                                                                               
      use mod_mpc, only : num_tied_sets, tied_contact_table,                    
     &                    tied_sets_exist                                       
      integer  dumi, err, npairs, mstr_lst(*), slv_lst(*)                       
      real  dumr, tied_tol                                                      
      double precision  dumd                                                    
      character(len=16) :: setid                                                
      character(len=1) :: dums                                                  
      logical       adj_flag                                                    
c                                                                               
      tied_sets_exist  = .true.                                                 
      num_tied_sets = num_tied_sets + 1                                         
      if (num_tied_sets .gt. max_tied_sets) then                                
         call errmsg2(42,max_tied_sets,dums,dumr,dumd)                          
         call die_abort                                                         
      end if                                                                    
c                                                                               
      tied_contact_table(num_tied_sets)%id         = setid                      
      tied_contact_table(num_tied_sets)%tolerance  = tied_tol                   
      tied_contact_table(num_tied_sets)%adjust_gap = adj_flag                   
      tied_contact_table(num_tied_sets)%num_pairs  = npairs                     
c                                                                               
      allocate (tied_contact_table(num_tied_sets)%master_list(npairs),          
     &          tied_contact_table(num_tied_sets)%slave_list(npairs),           
     &          stat=err)                                                       
      if (err .ne. 0) then                                                      
         call errmsg2(46,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
      do i = 1, npairs                                                          
         tied_contact_table(num_tied_sets)%master_list(i) = mstr_lst(i)         
         tied_contact_table(num_tied_sets)%slave_list(i)  = slv_lst(i)          
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine check_setid                   *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 01/06/04                   *          
c     *                                                              *          
c     *   this subroutine checks the tied set name input to see if   *          
c     *   it has already been used, then clears the previous entry   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine check_setid(setid)                                             
c                                                                               
      use mod_mpc, only : tied_contact_table, num_tied_sets                     
      integer  dumi, set                                                        
      real  dumr                                                                
      double precision  dumd                                                    
      character(len=16) :: setid                                                
      character(len=1) :: dums                                                  
c                                                                               
      nxt_set: do set = 1, num_tied_sets                                        
         if (setid .eq. tied_contact_table(set)%id) then                        
            call errmsg2(60,dumi,setid,dumr,dumd)                               
            tied_contact_table(set)%id         = ' '                            
            tied_contact_table(set)%tolerance  = 0.0                            
            tied_contact_table(set)%adjust_gap = .true.                         
            tied_contact_table(set)%num_pairs  = 0                              
            num_tied_sets = num_tied_sets - 1                                   
            exit nxt_set                                                        
         end if                                                                 
      end do  nxt_set                                                           
c                                                                               
      return                                                                    
      end                                                                       
