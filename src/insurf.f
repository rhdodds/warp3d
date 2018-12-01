c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine insurf                       *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 11/26/2018                 *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *         surfaces that define regions of mesh tieing          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine insurf()                                                       
      use global_data ! old common.main
c                                                                               
      use mod_mpc, only : surface_table                                         
c      parameter (max_ele=mxel/10)                                               
      integer  ele, nelem, face, dumi, errnum, icn, iplist, count,              
     &         len, err                                                         
      integer, allocatable, dimension (:) :: faces, elems                       
      real     dumr                                                             
      double precision  dumd                                                    
      character(len=16) :: surfid                                               
      character(len=1) ::   dums                                                
      logical  label, matchs, integr, true, bad_surf,                           
     &         abaqus_face_flag                                                 
      dimension  intlst(mxlsz)                                                  
c  
      max_ele = max( 5000, mxel/10 )                                                                             
      allocate (surface_table(max_surfaces),                                    
     &          faces(max_ele), elems(max_ele), stat=err)                       
      if (err .ne. 0) then                                                      
         call errmsg2(45,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
      new_surface: do                                                           
c                                                                               
c        check for name of surface, if none found, skip surface data            
c                                                                               
         if (.not. label(dumi)) then                                            
            call errmsg2(36,dumi,dums,dumr,dumd)                                
            do                                                                  
               call readsc                                                      
               if (matchs('element',7))  cycle                                  
               if (integr(ele))          cycle                                  
               if (matchs('surface',3)) then                                    
                  cycle  new_surface                                            
               else                                                             
                  call reset                                                    
                  if (true(dumi))  call splunj                                  
                  deallocate (faces, elems)                                     
                  return                                                        
               end if                                                           
            end do                                                              
         end if                                                                 
c                                                                               
         surfid = ' '                                                           
         call entits(surfid,len)                                                
         call check_surfid(surfid)                                              
c                                                                               
c        surface has name, reset counter, read surface data                     
c                                                                               
         nelem = 0                                                              
         surface_data: do                                                       
            call readsc                                                         
            if (matchs('surface',3)) then                                       
c                                                                               
c              a new surface found, store current surface data (if any)         
c                                                                               
               if (nelem .gt. 0) then                                           
                  call insurf_store(nelem,elems,faces,surfid)                   
                  cycle  new_surface                                            
               else                                                             
                  call errmsg2(57,dumi,surfid,dumr,dumd)                        
                  cycle  new_surface                                            
               end if                                                           
            end if                                                              
            if (matchs('element',7))  cycle surface_data                        
            call trlist(intlst,mxlsz,noelem,lenlst,errnum)                      
c                                                                               
c           branch on the errnum returned from trlist.                          
c           a value of 1 indicates no error.                                    
c           a value of 2 indicates that the parse rules failed in               
c                         the list.                                             
c           a value of 3 indicates that the list overflowed its                 
c                         maximum length of mxlsz.                              
c               in these last two cases, the illegal list will be               
c                ignored and a new node list will be sought.                    
c           a value of 4 indicates that no list was found.                      
c               in this case, surface input has ceased.                         
c                                                                               
            if (errnum .eq. 1) then                                             
                  iplist   = 1                                                  
                  icn      = 0                                                  
                  count    = 0                                                  
                  bad_surf = .false.                                            
                  call backsp(1)                                                
                  if (true(dumi))  call splunj                                  
               else if (errnum .eq. 2) then                                     
                  param = 1                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  cycle  surface_data                                           
               else if (errnum .eq. 3) then                                     
                  param = 2                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  cycle  surface_data                                           
               else if (errnum .eq. 4) then                                     
c                                                                               
c                 something other than a list found                             
c                 store current surface data (if any)                           
c                                                                               
                  if (nelem .gt. 0) then                                        
                     call insurf_store(nelem,elems,faces,surfid)                
                     call reset                                                 
                     if (true(dumi))  call splunj                               
                     deallocate (faces, elems)                                  
                     return                                                     
                  else                                                          
                     call errmsg2(57,dumi,surfid,dumr,dumd)                     
                     call reset                                                 
                     if (true(dumi))  call splunj                               
                     deallocate (faces, elems)                                  
                     return                                                     
                  end if                                                        
               else                                                             
                  param = 3                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  cycle  surface_data                                           
            end if                                                              
c                                                                               
c           integer list found and read, check for face data                    
c                                                                               
            if (.not. matchs("face",4)) then                                    
               call errmsg2(37,dumi,dums,dumr,dumd)                             
               cycle  surface_data                                              
            end if                                                              
c                                                                               
            abaqus_face_flag = .false.                                          
            if( matchs('abaqus',4) ) abaqus_face_flag = .true.                  
            if (.not. integr(face)) then                                        
               call errmsg2(37,dumi,dums,dumr,dumd)                             
               cycle  surface_data                                              
            end if                                                              
            if( matchs('abaqus',4) ) abaqus_face_flag = .true.                  
c                                                                               
c                      the abaqus flag applies only to hex                      
c                      elements, faces (4,5) in WARP3D are (5,4)                
c                      in abaqus. this option is convenient to allow            
c                      input of abaqus face numbers by user. should             
c                      not be used for tet elements                             
c                                                                               
            if( abaqus_face_flag ) then                                         
               if( face .eq. 4 ) then                                           
                   face = 5                                                     
               elseif( face .eq. 5 ) then                                       
                   face = 4                                                     
               endif                                                            
            endif                                                               
            abaqus_face_flag = .false.                                          
c                                                                               
c                                                                               
c           check surface data before storing                                   
c           skip surface list if it contains bad data by                        
c           cycling without incrementing the counter so                         
c           bad data is overwritten or not saved                                
c                                                                               
            do                                                                  
               count = count + 1                                                
                                                                                
               call trxlst(intlst,lenlst,iplist,icn,ele)                        
                                                                                
               if ((ele .gt. noelem).or.(ele .le. 0)) then                      
                  call errmsg2(38,ele,dums,dumr,dumd)                           
                  bad_surf = .true.                                             
               end if                                                           
                                                                                
               if (nelem+count .gt. max_ele) then                               
                  call errmsg2(44,max_ele,dums,dumr,dumd)                       
                  call die_abort                                                
               end if                                                           
                                                                                
               elems(nelem+count) = ele                                         
               faces(nelem+count) = face                                        
                                                                                
               if (iplist .eq. 0)  exit                                         
            end do                                                              
                                                                                
            if (face .le. 0) then                                               
               call errmsg2(39,ele,dums,dumr,dumd)                              
               bad_surf = .true.                                                
            end if                                                              
                                                                                
            if (bad_surf)  cycle surface_data                                   
                                                                                
            nelem = nelem + count                                               
         end do  surface_data                                                   
      end do  new_surface                                                       
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine insurf_store                   *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 04/06/03                   *          
c     *                                                              *          
c     *   this subroutine stores the surface data for mesh tieing    *          
c     *                into the proper data structures               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine insurf_store(nelem, elems, faces, surfid)                      
      use global_data ! old common.main
      use mod_mpc, only : num_surfaces, surface_table                           
      integer  err, nelem, elems(*), faces(*)                                   
      character(len=16) :: surfid                                               
c                                                                               
      num_surfaces = num_surfaces + 1                                           
      if (num_surfaces .gt. max_surfaces) then                                  
         call errmsg2(40,max_surfaces,dums,dumr,dumd)                           
         call die_abort                                                         
      end if                                                                    
      surface_table(num_surfaces)%id        = surfid                            
      surface_table(num_surfaces)%num_elems = nelem                             
      allocate (surface_table(num_surfaces)%elem_list(nelem),                   
     &          surface_table(num_surfaces)%face_list(nelem),stat=err)          
c                                                                               
      if (err .ne. 0) then                                                      
         call errmsg2(45,dumi,dums,dumr,dumd)                                   
         call die_abort                                                         
      end if                                                                    
c                                                                               
      do i = 1, nelem                                                           
         surface_table(num_surfaces)%elem_list(i) = elems(i)                    
         surface_table(num_surfaces)%face_list(i) = faces(i)                    
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine check_surfid                  *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 01/06/04                   *          
c     *                                                              *          
c     *   this subroutine checks the surface name input to see if    *          
c     *   it has already been used, then clears the previous entry   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine check_surfid(surfid)                                           
      use mod_mpc, only : surface_table, num_surfaces                           
      integer  surf                                                             
      character(len=16) :: surfid                                               
      nxt_surf: do surf = 1, num_surfaces                                       
         if (surfid .eq. surface_table(surf)%id) then                           
            call errmsg2(61,dumi,surfid,dumr,dumd)                              
            surface_table(surf)%id        = ' '                                 
            surface_table(surf)%num_elems = 0                                   
            num_surfaces = num_surfaces - 1                                     
            exit nxt_surf                                                       
         end if                                                                 
      end do  nxt_surf                                                          
      return                                                                    
      end                                                                       
