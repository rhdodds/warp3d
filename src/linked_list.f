c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine add_to_list                  *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 08/29/95                   *          
c     *                                                              *          
c     *     this subroutine adds a new entry to the crack front      *          
c     *     nodes linked list.                                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine add_to_list ( data )                                           
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : crack_front_nodes                           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      logical debug                                                             
      data debug /.false./                                                      
c                                                                               
      if (debug) then                                                           
         write (*,*) '>>>> entering add_to_list'                                
         call write_list                                                        
         write (*,*) 'adding node ',data                                        
      endif                                                                     
c                                                                               
c            If this is the first entry in the list (the start pointer is not   
c            set), then take the first entry from the garbage list and make     
c            it the first entry in the list.                                    
c                                                                               
      if ( crack_front_start .eq. -1 ) then                                     
         crack_front_start = crkfrnt_garbage_start                              
         crack_front_end   = crkfrnt_garbage_start                              
c                                                                               
c                    move garbage head pointer to next garbage entry            
c                                                                               
         crkfrnt_garbage_start =                                                
     &        crack_front_nodes ( crkfrnt_garbage_start , 2 )                   
c                                                                               
      else                                                                      
c                                                                               
c            This is not the first entry.  Take the top entry from the          
c            garbage list and make it part of the real list.                    
c                                                                               
c                     set temp to head of garbage list, move                    
c                     garbage head pointer to next entry down.                  
c                                                                               
         temp = crkfrnt_garbage_start                                           
         crkfrnt_garbage_start =                                                
     &        crack_front_nodes ( crkfrnt_garbage_start , 2 )                   
c                                                                               
c                     set pointer of the last entry in data list to             
c                     point to temp, then set list end to point to temp.        
c                                                                               
         crack_front_nodes(crack_front_end,2) = temp                            
         crack_front_end = temp                                                 
      endif                                                                     
c                                                                               
c                     set values of new entry -- data to 'data' and             
c                     pointer to -1                                             
c                                                                               
      crack_front_nodes(crack_front_end,1) = data                               
      crack_front_nodes(crack_front_end,2) = -1                                 
c                                                                               
c                     check to see if garbage list is empty.  if this           
c                     happens, something is seriously wrong -- the              
c                     list is probably hosed.                                   
c                                                                               
      if ( crkfrnt_garbage_start .eq. -1 ) then                                 
         write (*,*) '>>>> Fatal Error: all the crack plane nodes are'          
         write (*,*) '        now crack front nodes. Cannot continue.'          
         call die_gracefully                                                    
         stop                                                                   
      endif                                                                     
c                                                                               
      if ( debug ) then                                                         
         call write_list                                                        
         write (*,*) '<<<< leaving add_to_list'                                 
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine rm_from_list                 *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 08/29/95                   *          
c     *                                                              *          
c     *     this subroutine, given the pointer to a node that needs  *          
c     *     to be deleted from the list and the pointer to the node  *          
c     *     in the list above it, removes the entry that             *          
c     *     contains the node from the data list and adds it to      *          
c     *     the garbage list.                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine rm_from_list( node_ptr, prev_node_ptr )                        
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : crack_front_nodes                           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      logical debug                                                             
      data debug /.false./                                                      
c                                                                               
      if ( debug ) then                                                         
         write (*,*) '>>>> entering rm_from_list'                               
         call write_list                                                        
      endif                                                                     
c                                                                               
c                  find where entry is in the list -- top, bottom, or middle.   
c                  deal with list pointers appropriately                        
c                                                                               
      if ( node_ptr .eq. crack_front_start ) then                               
c                                                                               
c                        entry is at head of list. check if end also points     
c                        to this.  if so , set both pointers to -1.  If not,    
c                        move start pointer down one entry.                     
c                                                                               
         if ( crack_front_start .eq. crack_front_end ) then                     
            crack_front_start = -1                                              
            crack_front_end   = -1                                              
         else                                                                   
            crack_front_start = crack_front_nodes(crack_front_start,2)          
         endif                                                                  
c                                                                               
      else if ( node_ptr .eq. crack_front_end ) then                            
c                                                                               
c                        entry is at bottom of list. move end pointer up        
c                        one entry and set pointer of previous entry to         
c                   -1                                                          
c                                                                               
         crack_front_end = prev_node_ptr                                        
       if (crack_front_end .ne. -1)                                             
     &          crack_front_nodes(crack_front_end,2) = -1                       
c                                                                               
      else                                                                      
c                                                                               
c                        entry is in middle of list.  move pointer of           
c                        entry before the node to the entry past                
c                        the node.                                              
c                                                                               
         crack_front_nodes(prev_node_ptr,2) =                                   
     &        crack_front_nodes(node_ptr,2)                                     
      endif                                                                     
c                                                                               
c                  set this deleted entry data to -1, and pointer to -1         
c                                                                               
      crack_front_nodes(node_ptr,1) = -1                                        
      crack_front_nodes(node_ptr,2) = -1                                        
c                                                                               
c                   now get the last garbage entry and the garbage tail         
c                   pointer to point to the newly deleted entry                 
c                                                                               
      crack_front_nodes(crkfrnt_garbage_end,2) = node_ptr                       
      crkfrnt_garbage_end = node_ptr                                            
c                                                                               
      if ( debug ) then                                                         
         call write_list                                                        
         write (*,*) '<<<<< leaving rm_from_crk_list'                           
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      function find_in_list                   *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 08/31/95                   *          
c     *                                                              *          
c     *     this function traverses the list looking for 'node'      *          
c     *     in the data portion of the list.  If it finds it,        *          
c     *     the function returns the entry number where the data     *          
c     *     was found.  If not, it returns -1.  routine also provides*          
c     *     the pointer above as an argument (used for deleting an   *          
c     *     entry).                                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      integer function find_in_list (node, entry_above)                         
c                                                                               
      use node_release_data, only : crack_front_nodes                           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      logical debug                                                             
      data debug / .false. /                                                    
c                                                                               
      if (debug) then                                                           
         write (*,*) '>>>> entering find_in_list'                               
         write (*,*) ' looking for node:',node                                  
      endif                                                                     
c                                                                               
      function = -1                                                             
      entry_above = -1                                                          
c                                                                               
c            start at beginning of list. we may not have a list yet !           
c                                                                               
      pointer = crack_front_start                                               
      if ( pointer .eq. -1 ) then                                               
        find_in_list = -1                                                       
        if ( debug ) write (*,*) ' Node not in list.'                           
        goto 9999                                                               
      end if                                                                    
                                                                                
c                                                                               
 10   continue                                                                  
      if ( crack_front_nodes(pointer,1) .eq. node ) then                        
c                                                                               
c            a match has been found. set function to entry# and leave.          
c                                                                               
         if ( debug ) write (*,*) ' Node is found!'                             
         find_in_list = pointer                                                 
         goto 9999                                                              
c                                                                               
      else                                                                      
c                                                                               
c            no match found.  go to next entry.                                 
c                                                                               
         entry_above = pointer                                                  
         pointer     = crack_front_nodes(pointer,2)                             
         if ( pointer .eq.-1 ) then                                             
c                                                                               
c            we have reached end of list and not found node.  set               
c            function to -1 and leave.                                          
c                                                                               
            find_in_list = -1                                                   
            if ( debug ) write (*,*) ' Node not in list.'                       
            goto 9999                                                           
         endif                                                                  
      endif                                                                     
      goto 10                                                                   
c                                                                               
 9999 continue                                                                  
      if ( debug ) write (*,*) '<<<< leaving find_in_list'                      
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine write_list                   *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 08/29/95                   *          
c     *                                                              *          
c     *   this subroutine writes out the contents of the crack front *          
c     *   linked list.                                               *          
c     *                                                              *          
c     *****************************************`***********************         
c                                                                               
c                                                                               
c                                                                               
      subroutine write_list                                                     
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : crack_front_nodes                           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      write (*,*) '   >>>> printing crack_front_nodes'                          
c                                                                               
      write (*,*) '      entry#, entry, pointer):'                              
      do i = 1, num_crack_plane_nodes                                           
         write (*,'(5x,3i7)')i,crack_front_nodes(i,1),                          
     &        crack_front_nodes(i,2)                                            
      enddo                                                                     
c                                                                               
      write (*,*) '      Values for begin and end pointers:'                    
      write (*,'(7x,2i10)') crack_front_start, crack_front_end                  
      write (*,*) '      Values for garbage begin and end pointers:'            
      write (*,'(7x,2i10)') crkfrnt_garbage_start, crkfrnt_garbage_end          
c                                                                               
      write (*,*) '   <<<< finshed printing list'                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
