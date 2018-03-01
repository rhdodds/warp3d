c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine con_dump                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 11/10/94                   *          
c     *                                                              *          
c     *     This routine dumps the constraints table for a model     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine con_dump ( olddof )                                            
      use global_data ! old common.main
c                                                                               
      use main_data, only : cnstrn, cnstrn_in, inverse_incidences               
      use mod_mpc, only : num_tied_con_mpc, tied_con_mpc_table,                 
     &                    num_user_mpc, user_mpc_table, mpcs_exist,             
     &                    tied_con_mpcs_constructed                             
c                                                                               
      implicit integer (a-z)                                                    
      real  mult                                                                
c                                                                               
c                       locally allocated                                       
c                                                                               
c                                                                               
      logical no_con, con_set(3)                                                
c                                                                               
      character(len=1) :: con_label(3), con_string*14                           
      dimension con_string(3)                                                   
      data con_label(1),con_label(2),con_label(3)  /'u','v','w'/                
c                                                                               
      write (out,*) '>>>  Dumping constraints:'                                 
      write (out,*)                                                             
      write (out,*) '>>>    Absolute constraints:'                              
c                                                                               
      do node = 1, nonode                                                       
c                                                                               
c               find the constrains on a node, if any                           
c                                                                               
         ndof= iprops(4,inverse_incidences(node)%element_list(1))               
         no_con = .true.                                                        
         do i = 1,3                                                             
            con_string(i) = "             "                                     
            con_set(i) = .false.                                                
            dof = dstmap(node) + i - 1                                          
            if (cstmap(dof).ne.0) then                                          
               con_set(i) = .true.                                              
               no_con = .false.                                                 
            endif                                                               
         enddo                                                                  
c                                                                               
c               now print them                                                  
c                                                                               
         if (.not. no_con) then                                                 
            do i=1, 3                                                           
               dof = dstmap(node) + i - 1                                       
               if (con_set(i) .or.(dof.eq.olddof)) then                         
                  write (con_string(i),9010) cnstrn_in(dof)                     
               endif                                                            
            enddo                                                               
            write (out,9020)node,(con_label(k),con_string(k),k=1,3)             
         end if                                                                 
      enddo                                                                     
c                                                                               
      if (mpcs_exist) then                                                      
         write (out,*)                                                          
         write (out,*) '>>>    User-defined multi-point constraints:'           
         write (out,*) '>>>      ',num_user_mpc,' equations'                    
         do mpc = 1, num_user_mpc                                               
            ntrms = user_mpc_table(mpc)%num_terms                               
            write(out,*)                                                        
            do trm = 1, ntrms                                                   
               node = user_mpc_table(mpc)%node_list(trm)                        
               dof  = user_mpc_table(mpc)%dof_list(trm)                         
               mult = user_mpc_table(mpc)%multiplier_list(trm)                  
               if (trm .gt. 1) then                                             
                  write(out,9030) ' + ', node, mult, con_label(dof),','         
               else                                                             
                  write(out,9030) '   ', node, mult, con_label(dof),','         
               end if                                                           
            end do                                                              
            write(out,9040) ' = ', user_mpc_table(mpc)%constant                 
         end do                                                                 
      end if                                                                    
                                                                                
      if (tied_con_mpcs_constructed) then                                       
         write (out,*)                                                          
         write (out,*) '>>>    Tied-mesh multi-point constraints:'              
         write (out,*) '>>>      ',num_tied_con_mpc,' equations'                
         do mpc = 1, num_tied_con_mpc                                           
            ntrms = tied_con_mpc_table(mpc)%num_terms                           
            write(out,*)                                                        
            do trm = 1, ntrms                                                   
               node = tied_con_mpc_table(mpc)%node_list(trm)                    
               dof  = tied_con_mpc_table(mpc)%dof_list(trm)                     
               mult = tied_con_mpc_table(mpc)%multiplier_list(trm)              
               if (trm .gt. 1) then                                             
                  write(out,9030) ' + ', node, mult, con_label(dof),','         
               else                                                             
                  write(out,9030) '   ', node, mult, con_label(dof),','         
               end if                                                           
            end do                                                              
            write(out,9040) ' = ', tied_con_mpc_table(mpc)%constant             
         end do                                                                 
      end if                                                                    
                                                                                
      write (out,*)                                                             
      write (out,*) '<<<  Finished dumping constraints'                         
      return                                                                    
c                                                                               
 9010 format (e13.6)                                                            
 9020 format (4x,i7,3(3x,a1,1x,a13))                                            
 9030 format (8x,a,i7,es16.7,a2,a)                                              
 9040 format (8x,a,f15.4,a2,a)                                                  
c                                                                               
      end                                                                       
