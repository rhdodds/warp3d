c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine allocate_temp_load           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 11/15/11 jcs               *          
c     *                                                              *          
c     *     this subroutine handles the allocation and deallocation  *          
c     *     of the temporary element load input structures.          *          
c     *                                                              *          
c     *     These temporary structures hold the element loadings     *          
c     *     as specified for a specific loading case. Following the  *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine allocate_temp_load( status )                                   
      use global_data ! old common.main
      use elem_load_data                                                        
      implicit integer (a-z)                                                    
c                                                                               
      logical :: debug                                                          
      real :: zero     !  note:  single precision                               
      data zero, debug /0.0, .false./                                           
c                                                                               
      if ( debug ) write (*,*) ' >>> inside allocate_temp_load'                 
      go to (10,20), status                                                     
c                                                                               
      go to 9999                                                                
c                                                                               
c                allocate and zero temp_body                                    
c                                                                               
 10   continue                                                                  
      temp_allocated = .true.                                                   
      if (allocated(temp_body)) goto 9998                                       
      allocate(temp_body(noelem,mxndof))                                        
      call zero_vector1( temp_body, noelem*mxndof )                             
c                                                                               
c                allocate and zero temp_face                                    
c                                                                               
      if (allocated(temp_face)) goto 9998                                       
      allocate(temp_face(noelem,numfaces,mxndof))                               
      call zero_vector1( temp_face, noelem*mxndof*numfaces )                    
c                                                                               
c                allocate and zero temp_press                                   
c                                                                               
      if (allocated(temp_press)) goto 9998                                      
      allocate(temp_press(noelem,numfaces))                                     
      call zero_vector1( temp_press, noelem*numfaces )                          
c                                                                               
c                allocate and zero temp_temper                                  
c                                                                               
      if (allocated(temp_temper)) go to 9998                                    
      allocate(temp_temper(noelem))                                             
      do k = 1, noelem                                                          
       temp_temper(k) = zero                                                    
      end do                                                                    
c                                                                               
c                allocate and zero temp_piston                                  
c                                                                               
      if (allocated(temp_piston)) go to 9998                                    
      allocate(temp_piston(noelem,numfaces))                                    
      call zero_vector2( temp_piston, noelem*numfaces )                         
      go to 9999                                                                
c                                                                               
c                deallocate all                                                 
c                                                                               
 20   continue                                                                  
      if (allocated(temp_body)) deallocate(temp_body)                           
      if (allocated(temp_face)) deallocate(temp_face)                           
      if (allocated(temp_press)) deallocate(temp_press)                         
      if (allocated(temp_temper)) deallocate(temp_temper)                       
      if (allocated(temp_piston)) deallocate(temp_piston)                       
      temp_allocated = .false.                                                  
      goto 9999                                                                 
c                                                                               
 9998 if ( debug ) write(*,*) 'skipping'                                        
 9999 continue                                                                  
      if ( debug ) write (*,*) ' <<< leaving allocate_temp_load'                
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine allocate_perm_load           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 11/15/11 jcs               *          
c     *                                                              *          
c     *     this subroutine handles the allocation and deallocation  *          
c     *     of the permanent body and face load input structures.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine allocate_perm_load( status, loadnum )                          
      use global_data ! old common.main
      use elem_load_data                                                        
      implicit integer (a-z)                                                    
c                                                                               
      logical debug                                                             
c                                                                               
      real zero     !  note:  single precision                                  
      data zero, debug /0.0, .false./                                           
c                                                                               
      if ( debug ) write (*,*) ' >>> allocate_perm_load'                        
      go to (10,20), status                                                     
c                                                                               
      go to 9999                                                                
c                                                                               
c                  allocate the permanent space -- skip if temp space           
c                  has not been allocated. ( Passing 0 in as the number         
c                  of loads into do_perm_allo will skip the allocation )        
c                                                                               
 10   continue                                                                  
      num_loads = 0                                                             
      if ( .not. temp_allocated ) go to 15                                      
c                                                                               
c                       first find number of entries needed                     
c                                                                               
      if ( debug ) write (*,*) ' counting loadings...'                          
      do elem = 1, noelem                                                       
c         if ( debug ) write (*,*) '   ** element:',elem                        
         do dof = 1, mxndof                                                     
            if ( temp_body(elem,dof).ne.zero ) then                             
               num_loads = num_loads + 1                                        
c               if (debug) write (*,*) '    body'                               
            endif                                                               
         end do                                                                 
         do face = 1, numfaces                                                  
            if ( temp_press(elem,face).ne.zero ) then                           
               num_loads = num_loads + 1                                        
c               if (debug) write (*,*) '     pressure'                          
            elseif ( temp_piston(elem,face).ne.0 ) then                         
               num_loads = num_loads + 1                                        
c               if (debug) write(*,*) '     piston'                             
            else                                                                
               do dof = 1, mxndof                                               
                  if ( temp_face(elem,face,dof).ne.zero ) then                  
                     num_loads = num_loads + 1                                  
c                     if ( debug ) write (*,*) '     face'                      
                  end if                                                        
               end do                                                           
            end if                                                              
         end do                                                                 
         if ( temp_temper(elem) .ne. zero ) num_loads = num_loads + 1           
      end do                                                                    
c                                                                               
 15   continue                                                                  
      if ( debug ) write (*,'(" for loading ",i5," numloads=",i3)')             
     &     loadnum, num_loads                                                   
c                                                                               
c                       now allocate and fill permanent variables               
c                                                                               
      call do_perm_allo( 1, loadnum, num_loads, .true. )                        
      go to 9999                                                                
c                                                                               
c                      deallocate all the space                                 
c                                                                               
 20   continue                                                                  
      call do_perm_allo( 2, loadnum, 0, .false. )                               
      go to 9999                                                                
c                                                                               
c                                                                               
 9999 continue                                                                  
      if (debug) write (*,*) ' <<< leaving allocate_perm_load'                  
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine do_perm_allo                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 11/15/11 jcs               *          
c     *                                                              *          
c     *     this subroutine does the allocation                      *          
c     *     of the permanent body,face load and element temperature  *          
c     *     input data structures.                                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine do_perm_allo (status, loadnum, num_loads, fill)                
      use global_data ! old common.main
      use elem_load_data                                                        
      implicit integer (a-z)                                                    
c                                                                               
      logical debug, fill                                                       
      data debug / .false./                                                     
c                                                                               
      if ( debug ) write (*,*) ' >>> inside do_perm_allo'                       
      go to (1, 1000), status                                                   
c                                                                               
c                  Allocation of variables.  First check if the general         
c                  elem_loads structure has been allocated; allocate it         
c                  if needed.  Then set the size of the variables               
c                  from num_loads and allocate the loading structures           
c                  inside elem_loads for the current loading condition          
c                  number.  If num_loads is 0, skip allocation.                 
c                                                                               
 1    continue                                                                  
c                                                                               
      if (.not. allocated(elem_loads)) then                                     
         write (*,*) ' elem_loads vector not allocated. stop'                   
         call die_gracefully                                                    
         stop                                                                   
      endif                                                                     
c                                                                               
      if (num_loads.eq.0) then                                                  
         elem_loads(loadnum)%size = 0                                           
         if (debug) write (*,*) '  skipping allocation:'                        
         goto 9999                                                              
      endif                                                                     
c                                                                               
      if (elem_loads(loadnum)%size .ne. 0) then                                 
         deallocate (elem_loads(loadnum)%data)                                  
         deallocate (elem_loads(loadnum)%vals)                                  
         deallocate (elem_loads(loadnum)%piston_tabnum)                         
         deallocate (elem_loads(loadnum)%thread_number)                         
      endif                                                                     
c                                                                               
      allocate (elem_loads(loadnum)%data(num_loads,3))                          
      allocate (elem_loads(loadnum)%vals(num_loads))                            
      allocate (elem_loads(loadnum)%piston_tabnum(num_loads))                   
      allocate (elem_loads(loadnum)%thread_number(num_loads))                   
      elem_loads(loadnum)%size = num_loads                                      
c                                                                               
      if (fill) call store_perm_allo (                                          
     &     elem_loads(loadnum)%data(1,1),                                       
     &     elem_loads(loadnum)%vals(1),                                         
     &     elem_loads(loadnum)%piston_tabnum(1),                                
     &     elem_loads(loadnum)%size,                                            
     &     elem_loads(loadnum)%thread_number(1) )                               
c                                                                               
      goto 9999                                                                 
c                                                                               
c                  Deallocation of variables                                    
c                                                                               
 1000 continue                                                                  
c                                                                               
      if (elem_loads(loadnum)%size .ne. 0) then                                 
         deallocate (elem_loads(loadnum)%data)                                  
         deallocate (elem_loads(loadnum)%vals)                                  
         deallocate (elem_loads(loadnum)%piston_tabnum)                         
         deallocate (elem_loads(loadnum)%thread_number)                         
      endif                                                                     
c                                                                               
      goto 9999                                                                 
c                                                                               
c                                                                               
 9999 continue                                                                  
      if (debug) write (*,*) ' <<< leaving do_perm_allo'                        
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine store_perm_allo              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 11/15/11 jcs               *          
c     *                                                              *          
c     *     this subroutine stores the loading information into      *          
c     *     the permanent body and face load input structures.       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c           Here is how the face/body/pressure/temperature loading              
c           storage strategy works:                                             
c                                                                               
c                 There is one 3 by 'size' integer array and one                
c                 'size' long real array that store the loading                 
c                 information:                                                  
c                                                                               
c                   eload_data(i,1) -- the element with the loading             
c                                       description for entry i                 
c                   eload_data(i,2) -- the type of the loading for entry i:     
c                               1 to  6 : face number for traction loading      
c                                  0    : body force loading                    
c                              -1 to -6 : pressure loading                      
c                              -7 to -12: -face # for piston loading -6         
c                                -100   : element temperature                   
c                   eload_data(i,3)  -- the dof for the loading for entry i     
c                               (zero for pressure loading and                  
c                                temperature)                                   
c                   eload_val(i)     -- the value of the force for entry i      
c                                       value of element temperature            
c                   eload_pist(i)   -- piston table number for entry i          
c                   thread_number(i) -- thread number processing entry i        
c                                                                               
      subroutine store_perm_allo(                                               
     &     eload_data, eload_val, eload_pist, size, thread_number )             
      use global_data ! old common.main
      use elem_load_data                                                        
      implicit integer (a-z)                                                    
c                                                                               
      dimension eload_data(size,3), thdnum(noelem),                             
     &          thread_number(size), eload_pist(size)                           
      real eload_val(size)                                                      
      real zero    !  note: single precision                                    
      double precision                                                          
     &     mx, my, mz, norm                                                     
      logical debug                                                             
      data zero, debug / 0.0, .false. /                                         
c                                                                               
      if ( debug ) write (*,*) ' >>> inside store_perm_allo'                    
      force_entry = 0                                                           
c                                                                               
      thdnum(1:noelem) = -1                                                     
      count = 1                                                                 
c                                                                               
      do elem = 1, noelem                                                       
c                                                                               
c                      body loadings for the element                            
c                                                                               
         if ( debug ) write (*,*) '--------- body loadings:'                    
         do dof = 1, mxndof                                                     
            if ( temp_body(elem,dof).ne.zero ) then                             
               call eloads_thread_count(count,thdnum(elem),num_threads)         
               if(debug)write(*,'("elem,dof:",2i6,"temp_body:",e14.6)')         
     &              elem,dof,temp_body(elem,dof)                                
               force_entry               = force_entry + 1                      
               eload_data(force_entry,1) = elem                                 
               eload_data(force_entry,2) = 0                                    
               eload_data(force_entry,3) = dof                                  
               eload_val(force_entry)    = temp_body(elem,dof)                  
               eload_pist(force_entry)   = 0                                    
            end if                                                              
         end do                                                                 
c                                                                               
c                      face loadings for the element                            
c                                                                               
         if ( debug )write (*,*) '-------- face loadings:'                      
         do face = 1, numfaces                                                  
            if ( temp_press(elem,face).ne.zero ) then                           
               call eloads_thread_count(count,thdnum(elem),num_threads)         
               if(debug)write (*,'("press;el,fa:",2i6,"temp:",e14.6)')          
     &              elem,face,temp_press(elem,face)                             
               force_entry               = force_entry + 1                      
               eload_data(force_entry,1) = elem                                 
               eload_data(force_entry,2) = -face                                
               eload_data(force_entry,3) = 0                                    
               eload_val(force_entry)    = temp_press(elem,face)                
               eload_pist(force_entry)   = 0                                    
               if(debug)write (*,'("press;el,fa:",2i6,"temp:",e14.6)')          
     &              elem,eload_data(force_entry,2),                             
     &              eload_val(force_entry)                                      
            elseif (temp_piston(elem,face).ne.0 ) then                          
               if (debug) write (*,'("pistn;el,fa:",2i6)') elem,face            
               call eloads_thread_count(count,thdnum(elem),num_threads)         
               force_entry = force_entry + 1                                    
               eload_data(force_entry,1) = elem                                 
               eload_data(force_entry,2) = -face-6                              
               eload_data(force_entry,3) = 0                                    
               eload_val(force_entry)    = zero                                 
               eload_pist(force_entry)   = temp_piston(elem,face)               
            else                                                                
               do dof = 1, mxndof                                               
                if ( temp_face(elem,face,dof).ne.zero ) then                    
                   call eloads_thread_count(count,thdnum(elem),                 
     &                                      num_threads)                        
                   if(debug)write (*,'("face;el,fa,do:",3i6,"val:",             
     &                e14.6)')elem,face,dof,temp_face(elem,face,dof)            
                   force_entry               = force_entry + 1                  
                   eload_data(force_entry,1) = elem                             
                   eload_data(force_entry,2) = face                             
                   eload_data(force_entry,3) = dof                              
                   eload_val(force_entry)    = temp_face(elem,face,dof)         
                   eload_pist(force_entry)   = 0                                
                endif                                                           
               end do                                                           
            end if                                                              
         end do                                                                 
c                                                                               
c                      element temperature change                               
c                                                                               
         if ( debug ) write (*,*) '-------- temperature loadings:'              
         if ( temp_temper(elem).ne.zero ) then                                  
            call eloads_thread_count(count,thdnum(elem),num_threads)            
            force_entry               = force_entry + 1                         
            eload_data(force_entry,1) = elem                                    
            eload_data(force_entry,2) = -100                                    
            eload_data(force_entry,3) = 0                                       
            eload_val(force_entry)    = temp_temper(elem)                       
            eload_pist(force_entry)   = 0                                       
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
c                                                                               
c                                                                               
      if ( force_entry .ne. size ) then                                         
         write (*,*) ' >>>> Error: a memory fault in the element'               
         write (*,*) '     has been detected.  Stopping execution.'             
         call die_gracefully                                                    
         stop                                                                   
      end if                                                                    
c                                                                               
c                     set thread number                                         
c                                                                               
      do erow = 1, size                                                         
         elem = eload_data(erow,1)                                              
         thread_number(erow) = thdnum(elem)                                     
      end do                                                                    
c                                                                               
      if (debug) then                                                           
        call dump_load (eload_data, eload_val,                                  
     &     thread_number, size)                                                 
        write (*,*) ' <<< leaving store_perm_allo'                              
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dump_load                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/03/11 jcs               *          
c     *                                                              *          
c     *     this subroutine prints out the face and body loadings    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dump_load(eload_data, eload_val, thread_number, size)          
      use global_data ! old common.main
      implicit integer (a-z)                                                    
c                                                                               
      real eload_val(size)                                                      
      dimension eload_data(size,3), thread_number(size)                         
c                                                                               
c                                                                               
      write (*,*) ' ==================='                                        
      write (*,*) '   DUMPING LOADING  '                                        
      write (*,*) ' ==================='                                        
c                                                                               
c                                                                               
      write (*,1000)                                                            
      if (size .gt. 0) then                                                     
         do i = 1, size                                                         
            write (*,1010) i,(eload_data(i,k),k=1,3),                           
     &           thread_number(i),eload_val(i)                                  
         end do                                                                 
      else                                                                      
         write (*,*) 'size is zero.'                                            
      endif                                                                     
c                                                                               
      return                                                                    
 1000 format ( 'entry    elem    type  dof    thd#  val')                       
 1010 format (1x,i4,2x,i6,4x,i4,4x,i1,i5,'   ',e14.6)                           
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine thread_count                 *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *                   last modified : 10/26/11                   *          
c     *                                                              *          
c     *     this subroutine updates the thread count                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine eloads_thread_count(count,thdnum,num_threads)                  
      implicit integer (a-z)                                                    
c                                                                               
c                     thdnum represents the thread assigned to this             
c                     element, by default (-1), no thread is assigned           
c                     check if a thread has been assigned, if so, skip          
c                                                                               
      if (thdnum .ne. -1) return                                                
c                                                                               
c                     thdnum has not been assigned, assign thread               
c                     based on count. count provides the current                
c                     thread to assign to an element.                           
c                     each time a new element is assigned to count,             
c                     count is updated to the next thread.                      
c                     this process balances the loads processed on              
c                     each thread                                               
c                                                                               
      thdnum = count                                                            
      count  = count + 1                                                        
      if ( count .gt. num_threads ) count = 1                                   
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine eloads_rebuild_thread_list          *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *                   last modified : 10/26/11                   *          
c     *                                                              *          
c     *     this subroutine rebuilds the thread list after           *          
c     *     a restart                                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine eloads_rebuild_thread_list( eload_data, size,                  
     &                                       thread_number )                    
      use global_data ! old common.main
      implicit integer (a-z)                                                    
c                                                                               
      dimension eload_data(size,3), thdnum(noelem), thread_number(size)         
c                                                                               
      logical debug                                                             
      data debug / .false. /                                                    
c                                                                               
      if (debug) write (*,*) ' >>> inside eloads_rebuild_thread_list'           
c                                                                               
c                     set thdnum and count to initial values                    
c                                                                               
      thdnum(:) = -1                                                            
      count = 1                                                                 
c                                                                               
c                     cycle through each element loading and determine          
c                     the element                                               
c                                                                               
      do erow = 1, size                                                         
c                                                                               
c                     check if the element for this row has a thread            
c                     assigned to it already. if not, set thread                
c                                                                               
         element = eload_data( erow, 1 )                                        
c         if (debug) write(*,*) 'erow, element, thread number',                 
c     &        erow, element, thdnum(element)                                   
c                                                                               
         call eloads_thread_count(count,thdnum(element),num_threads)            
c            if (debug) write(*,*), 'new thread, count, #threads',              
c     &           thdnum(element), count, num_threads                           
c                                                                               
c                     set thread number of the element row                      
c                                                                               
         thread_number(erow) = thdnum(element)                                  
c                                                                               
      end do                                                                    
c                                                                               
      if (debug) write (*,*) ' >>> leaving rebuild_thread_list'                 
      return                                                                    
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine zero_vector1                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 11/1/11 rhd                *          
c     *                                                              *          
c     *     zero a "real" vector of specified length w/ floating     *          
c     *     zero. this local version for this .f to enable inlining  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine zero_vector1( vec, n )                                         
      real vec(*), zero    !   real not dp                                      
      data zero / 0.0 /                                                         
c                                                                               
      vec(1:n) = zero                                                           
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine zero_vector2                 *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *     zero a "integer" vector of specified length w/ integer   *          
c     *     zero. this local version for this .f to enable inlining  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine zero_vector2( vec, n )                                         
      integer vec(*), zero   ! integer, not dp                                  
      data zero / 0 /                                                           
c                                                                               
      vec(1:n) = zero                                                           
c                                                                               
      return                                                                    
      end                                                                       
