c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouhdva                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   modified :    10/29/2012 rhd               *          
c     *                   packet mods : 11/17/00 jp                  *          
c     *                                                              *          
c     *     drive hardcopy and binary packet output for the          *          
c     *     displacements, velocities, accelerations, reactions,     *          
c     *     temperatures for a given list of nodes or nodes          *          
c     *     connected to a specified list of elements                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine ouhdva( lsttyp, dva, wide, eform, prec, intlst,                
     &                   lenlst, noheader, react_totals_only,                   
     &                   out_packet_now  )                                      
      use global_data ! old common.main
c                                                                               
      use main_data, only : incmap, incid, output_packets,                      
     &                      packet_file_no, inverse_incidences,                 
     &                      inverse_dof_map                                     
c                                                                               
      implicit integer (a-z)                                                    
      logical wide, eform, prec, noheader, react_totals_only,                   
     &        out_packet_now                                                    
      dimension intlst(*)                                                       
c                                                                               
c                        local declarations                                     
c                                                                               
      character :: doflbl(mxndof)*8, dums*1                                     
      character(len=20) :: hedtyp                                               
      logical newel, write_to_packet                                            
      real dumr                                                                 
      double precision                                                          
     &     dumd, zero, react_sums(3)                                            
      data zero / 0.0d00 /                                                      
c                                                                               
c                       initialize parameters controlling output.               
c                                                                               
      lnum        = 56                                                          
      pgnum       = 0                                                           
      lbltyp      = 0                                                           
      local_count = 0                                                           
      write_to_packet = .false.                                                 
c                                                                               
      doflbl(1:mxndof) = ' '                                                    
      react_sums(1:3) = zero                                                    
      icn    = 0                                                                
      iplist = 1                                                                
c                                                                               
c                       lsttyp = 1 => user gave a list of nodes                 
c                              = 2 => user gave a list of elements              
c                                                                               
      if( lsttyp .eq. 2 ) go to 25                                              
c                                                                               
c                                                                               
c                                                                               
c                       node list specified                                     
c                       ===================                                     
c                                                                               
c                       extract and process all entries                         
c                                                                               
c                                                                               
c                                                                               
 15   continue  ! next node entry in list                                       
      call trxlst( intlst, lenlst, iplist, icn, node )                          
c                                                                               
c                       check validity of node id                               
c                                                                               
      if( node .gt. nonode ) then                                               
        call errmsg( 16, node, dums, dumr, dumd )                               
        go to 20 ! next entry                                                   
      end if                                                                    
      if( node .lt. 0 ) then                                                    
        call errmsg( 58, node, dums, dumr, dumd )                               
        go to 20 ! next entry                                                   
      end if                                                                    
c                                                                               
c                       set the dof labels for the nodes in the list,           
c                       if they have not already been established. if           
c                       they have been established, make sure that              
c                       the current node is compatable. if not, print           
c                       a warning and change the dof labels.                    
c                                                                               
      elem   = inverse_incidences(node)%element_list(1)                         
      elnod  = inverse_dof_map(node)%edof_table(1,1)                            
      type   = iprops(1,elem)                                                   
      call oulbdd( dva, hedtyp, lbltyp, type, elem, doflbl )                    
c                                                                               
c                       print the desired results for the current node.         
c                       keep a count on the number of lines that                
c                       will be output to the packet file when needed.          
c                                                                               
      ndof  = iprops(4,elem)                                                    
      newel = .false.                                                           
      if( .not. write_to_packet ) local_count = local_count + 1                 
c                                                                               
c                        ouhnod will be called here when there is               
c                        output to the main output file.                        
c                                                                               
      if( .not. out_packet_now )                                                
     &  call ouhnod( dva, node, lnum, pgnum, doflbl, wide, eform,               
     &               prec, ndof, hedtyp, elem, newel, react_sums,               
     &               noheader, react_totals_only,                               
     &               write_to_packet, lsttyp )                                  
c                                                                               
c                        ouhnod will be called here when output is              
c                        directed to the packet file, out_packet_now            
c                        controls the output to the packet, but                 
c                        write_to_packet will only become true when             
c                        the correct local_count has been calculated.           
c                                                                               
      if( out_packet_now .and. write_to_packet )                                
     &    call ouhnod( dva, node, lnum, pgnum, doflbl, wide, eform,             
     &                 prec, ndof, hedtyp, elem, newel, react_sums,             
     &                 noheader, react_totals_only,                             
     &                 write_to_packet, lsttyp )                                
c                                                                               
c                       process the next node in user list                      
c                                                                               
 20   continue                                                                  
      if( iplist .ne. 0 ) go to 15 ! get next entry                             
c                                                                               
c                       if the packet output is finished, continue.             
c                       if not, check to see if it is needed.                   
c                                                                               
      if( write_to_packet ) go to 500                                           
c                                                                               
c                       packet output will re-execute the previous              
c                       section, directing output to the packet                 
c                       file.  skip this section when packet output             
c                       is not valid for the current warp analysis or           
c                       this specific output command.                           
c                                                                               
      if ( out_packet_now ) then ! packet header record                         
         if ( dva .eq. 4) then                                                  
           if ( .not. react_totals_only )                                       
     &        write(packet_file_no) 4, local_count+1, ltmstp, 0                 
         else                                                                   
           pkttype = dva                                                        
           if( dva .eq. 5 ) pkttype = 29 ! temperatures                         
           write(packet_file_no) pkttype, local_count, ltmstp, 0                
         end if                                                                 
c                                                                               
c                       zero reaction sums, and set write_to_packet             
c                       to true when output is directed to the packet           
c                       file.                                                   
c                                                                               
         react_sums(1:3) = zero                                                 
         write_to_packet = .true.                                               
         noheader = .true.                                                      
         icn    = 0                                                             
         iplist = 1                                                             
         go to 15  ! process node list again. write packet records              
      end if                                                                    
c                                                                               
      go to 500   ! finish updone with processing node list                     
c                                                                               
c                                                                               
c                                                                               
c                       element list specified                                  
c                       ======================                                  
c                                                                               
c                       extract all entries and process                         
c                                                                               
c                                                                               
c                                                                               
 25   continue  ! next element in user list                                     
      call trxlst( intlst, lenlst, iplist, icn, elem )                          
c                                                                               
      if( elem .gt. noelem ) then                                               
        call errmsg( 35, elem, dums, dumr, dumd )                               
        go to 35   ! next entry                                                 
      end if                                                                    
      if( elem .lt. 0 ) then                                                    
        call errmsg( 86, elem, dums, dumr, dumd )                               
        go to 35 ! next entry                                                   
      end if                                                                    
c                                                                               
c                       loop over the nodes on element                          
c                                                                               
      incptr = incmap(elem)-1                                                   
      type   = iprops(1,elem)                                                   
      nnode  = iprops(2,elem)                                                   
      newel  = .true.                                                           
c                                                                               
      do elnod = 1, nnode                                                       
        node = incid(incptr+elnod)                                              
        call oulbdd( dva, hedtyp, lbltyp, type, elem, doflbl )                  
        ndof = iprops(4,elem)                                                   
c                                                                               
c                        keep a count on the number of lines that               
c                        will be output to the packet file when                 
c                        needed.                                                
c                                                                               
        if( .not. write_to_packet ) local_count = local_count + 1               
c                                                                               
c                        ouhnod will be called here when there is               
c                        output to the main output file.                        
c                                                                               
        if( .not. out_packet_now )                                              
     &     call ouhnod( dva, node, lnum, pgnum, doflbl, wide, eform,            
     &                  prec, ndof, hedtyp, elem, newel, react_sums,            
     &                  noheader, react_totals_only,                            
     &                  write_to_packet, lsttyp )                               
c                                                                               
c                        ouhnod will be called here when output is              
c                        directed to the packet file, out_packet_now            
c                        controls the output to the packet, but                 
c                        write_to_packet will only become true when             
c                        the correct local_count has been calculated.           
c                                                                               
        if( out_packet_now .and. write_to_packet )                              
     &    call ouhnod( dva, node, lnum, pgnum, doflbl, wide, eform,             
     &                 prec, ndof, hedtyp, elem, newel, react_sums,             
     &                 noheader, react_totals_only,                             
     &                 write_to_packet, lsttyp )                                
c                                                                               
      end do                                                                    
c                                                                               
 35   continue  ! next element in list                                          
      if( iplist .ne. 0 ) go to 25                                              
c                                                                               
c                       completed element list. if the packet                   
c                       output is finished, continue.                           
c                       if not, check to see if it is needed.                   
c                                                                               
      if( write_to_packet ) go to 500                                           
c                                                                               
c                       for packet element output,                              
c                       packet output will re-execute the previous              
c                       section, directing output to the packet                 
c                       file.  skip this section when packet output             
c                       is not valid for the current output command.            
c                                                                               
      if ( out_packet_now ) then                                                
c                                                                               
         if ( dva .eq. 4) then                                                  
           if ( .not. react_totals_only )                                       
     &        write(packet_file_no) 14, local_count+1, ltmstp, 0                
         else                                                                   
           pkttype = dva + 10                                                   
           if( dva .eq. 5 ) pkttype = 30 ! temperatures                         
           write(packet_file_no) pkttype, local_count, ltmstp, 0                
         end if                                                                 
c                                                                               
         react_sums(1:3) = zero                                                 
         write_to_packet = .true.                                               
         noheader = .true.                                                      
         icn    = 0                                                             
         iplist = 1                                                             
         go to 25  ! process element list again, write packet records           
      end if                                                                    
c                                                                               
c                       for reaction output, write totals                       
c                       for the set of nodes at which results were              
c                       just printed.                                           
c                                                                               
 500  continue                                                                  
      if ( dva .eq. 4 ) then                                                    
         if ( .not. out_packet_now ) then                                       
           write(out,9000) react_sums(1), react_sums(2), react_sums(3)          
         end if                                                                 
      end if                                                                    
c                                                                               
c                       output totals to binary packet file when                
c                       required                                                
c                                                                               
      if ( .not. output_packets ) return                                        
      if ( dva .eq. 4 .and. out_packet_now ) then                               
         if ( react_totals_only .and. lsttyp .eq. 1 ) then                      
            write(packet_file_no)4,1,ltmstp,0                                   
         else if (react_totals_only .and. lsttyp .ne. 1 ) then                  
            write(packet_file_no)14,1,ltmstp,0                                  
         end if                                                                 
         write ( packet_file_no ) react_sums(1),react_sums(2),                  
     &                            react_sums(3)                                 
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/,5x,'Totals:',3(8x,e12.5))                                        
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
