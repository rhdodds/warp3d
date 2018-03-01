c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine indom                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/7/10 RHD                 *          
c     *                                   add tangent vector         *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *     the parameters defining a domain for computation of      *          
c     *     j-integrals                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine indom( sbflg1, sbflg2 )                                        
      use global_data ! old common.main
      use j_data                                                                
      implicit integer (a-z)                                                    
      real dumr, qnow, rword, rzero                                             
      character dums*1, ltitle*80                                               
      double precision                                                          
     &   dumd, forval, zero                                                     
      dimension intlst(mxlsz)                                                   
      logical sbflg1, sbflg2, stress_flag, strain_flag                          
      logical matchs, endcrd, numd, realn, true,                                
     &        label, string, numi, user_def_ct                                  
      equivalence (iword,rword)                                                 
      data rzero, zero / 0.0, 0.0d0 /                                           
c                                                                               
c                       locally allocated                                       
c                                                                               
c                       if sub flag 1 is on, indom has been re-                 
c                       entered and the last command was in error.              
c                                                                               
      if ( sbflg1 ) then                                                        
         call errmsg(205,dum,dums,dumr,dumd)                                    
      end if                                                                    
c                                                                               
c                                                                               
      call initdm                                                               
c                                                                               
c                       get the domain id if there is one.                      
c                                                                               
      ltitle(1:80) = ' '                                                        
      if ( label(dummy)) then                                                   
        call entits( ltitle, nchar )                                            
        domain_id(1:) = ltitle(1:max(24,nchar))                                 
      else if ( string(dummy) ) then                                            
        call entits( ltitle, nchar )                                            
        domain_id(1:) = ltitle(1:max(24,nchar))                                 
      end if                                                                    
c                                                                               
 10   call readsc                                                               
c                                                                               
c                       branch on sub-command to define aspects of              
c                       a domain                                                
c                                                                               
 20   if ( matchs('normal',5)    ) go to 100                                    
      if ( matchs('symmetric',5) ) go to 200                                    
      if ( matchs('front',5)     ) go to 300                                    
      if ( matchs('q',1)         ) go to 400                                    
      if ( matchs('print',5)     ) go to 500                                    
      if ( matchs('function',4)  ) go to 600                                    
      if ( matchs('elements',4)  ) go to 700                                    
      if ( matchs('debug',4)     ) go to 800                                    
      if ( matchs('use',3)       ) go to 900                                    
      if ( matchs('ignore',3)    ) go to 1000                                   
      if ( matchs('neglect',3)   ) go to 1000                                   
      if ( matchs('omit',4)      ) go to 1100                                   
      if ( matchs('output',6)    ) go to 1200                                   
      if ( matchs('crack',5)     ) go to 1300                                   
      if ( matchs('plane',5)     ) go to 1400                                   
      if ( matchs('j',1)         ) go to 1500                                   
      if ( matchs('curvature',4) ) go to 1600                                   
      if ( matchs('dump',4)      ) go to 2000                                   
      if ( matchs('node',4)      ) go to 2500                                   
      if ( matchs('tangent',4)   ) go to 2600                                   
c                                                                               
c                       no match with domain definition command.                
c                       return to driver subroutine to look for high            
c                       level command.                                          
c                                                                               
      call reset                                                                
      if ( true( dummy ) ) call splunj                                          
      go to 9999                                                                
c                                                                               
c                                                                               
c               direction cosines for normal to crack plane.                    
c               ============================================                    
c                                                                               
 100  continue                                                                  
      if ( matchs('plane',4) ) call splunj                                      
 110  if ( endcrd(dummy) ) go to 10                                             
      if ( matchs('nx',2) ) then                                                
          if ( numd(crack_plane_normal(1)) ) go to 110                          
          call errmsg(199,dum,dums,dumr,dumd)                                   
          go to 10                                                              
      end if                                                                    
      if ( matchs('ny',2) ) then                                                
          if ( numd(crack_plane_normal(2)) ) go to 110                          
          call errmsg(199,dum,dums,dumr,dumd)                                   
          go to 10                                                              
      end if                                                                    
      if ( matchs('nz',2) ) then                                                
          if ( numd(crack_plane_normal(3)) ) go to 110                          
          call errmsg(199,dum,dums,dumr,dumd)                                   
          go to 10                                                              
      end if                                                                    
      call errmsg(200,dum,dums,dumr,dumd)                                       
      go to 10                                                                  
c                                                                               
c                                                                               
c               symmetric flag (multiply j-values  by 2)                        
c               ========================================                        
c                                                                               
 200  continue                                                                  
      symmetric_domain = .true.                                                 
      go to 10                                                                  
c                                                                               
c                                                                               
c               front nodes                                                     
c               ===========                                                     
c               list of front nodes, front order and verify                     
c               coincident node verify flag.                                    
c                                                                               
 300  continue                                                                  
      if ( matchs('nodes',4) ) call splunj                                      
      if ( matchs('sets',3) ) then                                              
         user_def_ct = .true.                                                   
      else                                                                      
         user_def_ct = .false.                                                  
         if ( true( dummy ) ) call splunj                                       
         call backsp( 1 )                                                       
      end if                                                                    
      call scan                                                                 
      call trlist(intlst,mxlsz,0,lenlst,errnum)                                 
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in the last three cases, command will be ig-            
c                       nored and a new domain parameter command                
c                       will be sought.                                         
c                                                                               
            if ( errnum .eq. 2 ) then                                           
               param = 1                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 3 ) then                                      
               param = 2                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 4 ) then                                      
               param = 4                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else                                                                
               if ( errnum .ne. 1 ) then                                        
                  param = 3                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  go to 10                                                      
               end if                                                           
            end if                                                              
c                                                                               
c                      we got a valid integer list. expand the                  
c                      list of nodes into the vector of front                   
c                      nodes and get a count.                                   
c                                                                               
        icn    = 0                                                              
        iplist = 1                                                              
 310    call trxlst(intlst,lenlst,iplist,icn,node)                              
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
        if ( node .gt. nonode) then                                             
         param = node                                                           
         call errmsg(16,param,dums,dumr,dumd)                                   
         go to 320                                                              
        end if                                                                  
c                                                                               
c                       check that the list node is not negative.               
c                                                                               
      if ( node .lt. 0 ) then                                                   
         param = node                                                           
         call errmsg(58,param,dums,dumr,dumd)                                   
         go to 320                                                              
      end if                                                                    
c                                                                               
c                       save the front node in the list                         
c                       if the user is inputing node sets, then                 
c                       save the value as a negative number.                    
c                                                                               
c                                                                               
      num_front_nodes = num_front_nodes + 1                                     
      if ( num_front_nodes .gt. max_front_nodes ) then                          
         call errmsg(331,max_front_nodes,dums,dumr,dumd)                        
         go to 320                                                              
      end if                                                                    
      if ( user_def_ct ) then                                                   
         front_nodes(num_front_nodes) = -node                                   
      else                                                                      
         front_nodes(num_front_nodes) = node                                    
      end if                                                                    
c                                                                               
 320  if ( iplist .ne. 0 ) go to 310                                            
      if ( true( dummy ) ) call splunj                                          
      call backsp( 1 )                                                          
c                                                                               
c                                                                               
c                       done with front node list. set the                      
c                       order of crack front interpolation.                     
c                       we now let the user set the element                     
c                       type. get the verify flag if                            
c                       present.                                                
c                                                                               
      front_order = 0                                                           
      if ( matchs('linear',4)     ) front_order = 1                             
      if ( matchs('quadratic',4)  ) front_order = 2                             
      if ( matchs('parabolic',4)  ) front_order = 2                             
      if ( matchs('l3disop',4)    ) front_order = 1                             
      if ( matchs('q3disop',4)    ) front_order = 2                             
      if ( matchs('verify',4)     ) verify_front = .true.                       
c                                                                               
c               set tolerance to identify coincident front nodes                
c               ================================================                
c                                                                               
      if( matchs('tolerance',3) ) then                                          
         if( matchs('=',1)         ) call splunj                                
         if( .not.numd(box_tol_relat) ) then                                    
            call errmsg(57,dum,dums,dumr,dumd)                                  
         end if                                                                 
      end if                                                                    
c                                                                               
      go to 10                                                                  
c                                                                               
c                                                                               
c               q-values by node or automatic ring definitions                  
c               ==============================================                  
c                                                                               
 400  continue                                                                  
      if ( matchs('-',1)      ) call splunj                                     
      if ( matchs('values',4) ) call splunj                                     
      if ( .not. matchs('automatic',4) ) go to 450                              
c                                                                               
c               list of rings given                                             
c                                                                               
      if ( matchs('rings',4) ) call splunj                                      
      call scan                                                                 
      call trlist(intlst,mxlsz,0,lenlst,errnum)                                 
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in the last three cases, command will be ig-            
c                       nored and a new domain parameter command                
c                       will be sought.                                         
c                                                                               
            if ( errnum .eq. 2 ) then                                           
               param = 1                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 3 ) then                                      
               param = 2                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 4 ) then                                      
               param = 4                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else                                                                
               if ( errnum .ne. 1 ) then                                        
                  param = 3                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  go to 10                                                      
               end if                                                           
            end if                                                              
c                                                                               
c                      we got a valid integer list. expand the                  
c                      list. set the ring number on or off in the               
c                      list of rings.                                           
c                                                                               
      icn    = 0                                                                
      iplist = 1                                                                
 410  call trxlst( intlst, lenlst, iplist, icn, ring_num )                      
c                                                                               
c                       check that the list ring does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
      if ( ring_num .gt. max_domain_rings) then                                 
         param = max_domain_rings                                               
         call errmsg2(76,param,dums,dumr,dumd)                                  
         go to 420                                                              
      end if                                                                    
c                                                                               
c                       check that the list ring is not negative.               
c                                                                               
      if ( ring_num .lt. 0 ) then                                               
         param = ring_num                                                       
         call errmsg2(77,param,dums,dumr,dumd)                                  
         go to 420                                                              
      end if                                                                    
c                                                                               
c                       save the ring number                                    
c                                                                               
      ring_list(ring_num) = 1                                                   
      num_auto_rings      = num_auto_rings + 1                                  
c                                                                               
 420  if ( iplist .ne. 0 ) go to 410                                            
      rings_given = .true.                                                      
      go to 10                                                                  
                                                                                
c                                                                               
c               list of nodes followed value.                                   
c               this command may be repeated as necessary within the            
c               domain definition.                                              
c                                                                               
 450  continue                                                                  
      call trlist(intlst,mxlsz,0,lenlst,errnum)                                 
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in the last three cases, command will be ig-            
c                       nored and a new domain parameter command                
c                       will be sought.                                         
c                                                                               
            if ( errnum .eq. 2 ) then                                           
               param = 1                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 3 ) then                                      
               param = 2                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 4 ) then                                      
               param = 4                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else                                                                
               if ( errnum .ne. 1 ) then                                        
                  param = 3                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  go to 10                                                      
               end if                                                           
            end if                                                              
c                                                                               
c                      we got a valid integer list. get the                     
c                      q-value to be assigned to nodes in list.                 
c                                                                               
      if ( true( dummy ) ) call splunj                                          
      call backsp( 1 )                                                          
      if ( .not. realn( qnow) ) then                                            
         call errmsg(97,param,dums,dumr,dumd)                                   
      end if                                                                    
c                                                                               
      icn    = 0                                                                
      iplist = 1                                                                
      if ( .not. allocated( compr_q_list ) ) then                               
            allocate( compr_q_list(nonode,2) )                                  
            compr_q_list(1:nonode,1) = 0                                        
            rword = rzero                                                       
            compr_q_list(1:nonode,2) = iword                                    
      end if                                                                    
 460  call trxlst(intlst,lenlst,iplist,icn,node)                                
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
      if ( node .gt. nonode) then                                               
         param = node                                                           
         call errmsg(16,param,dums,dumr,dumd)                                   
         go to 480                                                              
      end if                                                                    
c                                                                               
c                       check that the list node is not negative.               
c                                                                               
      if ( node .lt. 0 ) then                                                   
         param = node                                                           
         call errmsg(58,param,dums,dumr,dumd)                                   
         go to 480                                                              
      end if                                                                    
c                                                                               
c                       save the node in the list of nodes. column 1 is         
c                       the node, column two is q-value. q-values are           
c                       single precision reals.                                 
c                                                                               
      rword                       = qnow                                        
      last_compr                  = last_compr + 1                              
      compr_q_list(last_compr,1)  = node                                        
      compr_q_list(last_compr,2)  = iword                                       
c                                                                               
 480  if ( iplist .ne. 0 ) go to 460                                            
      qvals_given = .true.                                                      
c                                                                               
c                       done with q-values by node                              
c                                                                               
      go to 10                                                                  
c                                                                               
c                                                                               
c               print option flags                                              
c               ==================                                              
c                                                                               
 500  continue                                                                  
      if ( matchs( 'totals', 5  ) ) then                                        
         print_totals = .true.                                                  
         go to 10                                                               
      end if                                                                    
      if ( matchs( 'element', 4 ) ) then                                        
         if( matchs( 'values', 4 ) ) call splunj                                
         print_elem_values = .true.                                             
         go to 10                                                               
      end if                                                                    
c                                                                               
c                                                                               
c               domain type                                                     
c               ===========                                                     
c                                                                               
 600  continue                                                                  
      if ( matchs( 'type', 3 ) ) call splunj                                    
      if ( label(dummy) ) then                                                  
         call entits ( dums, nchar )                                            
      else if ( string(dummy) ) then                                            
         call entits ( dums, nchar )                                            
      else                                                                      
         goto 10                                                                
      endif                                                                     
      if ( dums .eq. 'a' ) domain_type = 1                                      
      if ( dums .eq. 'b' ) domain_type = 2                                      
      if ( dums .eq. 'c' ) domain_type = 3                                      
      if ( dums .eq. 'd' ) domain_type = 4                                      
      go to 10                                                                  
c                                                                               
c                                                                               
c               elements to be included                                         
c               =======================                                         
c               list of elements to be included in domain                       
c               when user specifies q-values explicitly.                        
c               this command can be repeated as necessary during                
c               definition of a domain                                          
c                                                                               
c                                                                               
 700  continue                                                                  
      call scan                                                                 
      call trlist(intlst,mxlsz,0,lenlst,errnum)                                 
            if ( errnum .eq. 2 ) then                                           
               param = 1                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 3 ) then                                      
               param = 2                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 4 ) then                                      
               param = 4                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else                                                                
               if ( errnum .ne. 1 ) then                                        
                  param = 3                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  go to 10                                                      
               end if                                                           
            end if                                                              
c                                                                               
c                      we got a valid integer list. expand the                  
c                      list of elements into the vector of element              
c                      numbers                                                  
c                                                                               
      icn    = 0                                                                
      iplist = 1                                                                
      if ( .not. allocated( q_element_maps ) ) then                             
           allocate( q_element_maps(noelem/30+1) )                              
           q_element_maps(1:noelem/30+1) = 0                                    
      end if                                                                    
 710  call trxlst(intlst,lenlst,iplist,icn,elemno)                              
c                                                                               
c                       check that the list element does not exceed             
c                       the number of elements in the structure.                
c                                                                               
      if ( elemno .gt. noelem) then                                             
         param = elemno                                                         
         call errmsg(35,param,dums,dumr,dumd)                                   
         go to 720                                                              
      end if                                                                    
c                                                                               
c                       check that the list element is not negative.            
c                                                                               
      if ( elemno .lt. 0 ) then                                                 
         param = elemno                                                         
         call errmsg(86,param,dums,dumr,dumd)                                   
         go to 720                                                              
      end if                                                                    
c                                                                               
c                       save the element number in the list bit map             
c                                                                               
      word = ( elemno-1 ) / 30 + 1                                              
      bit  = elemno - (word-1)*30                                               
      q_element_maps(word) = ior( q_element_maps(word), bits(bit) )             
c                                                                               
 720  if ( iplist .ne. 0 ) go to 710                                            
      go to 10                                                                  
c                                                                               
c                                                                               
c               debug flags for domain integral computations                    
c               ============================================                    
c                                                                               
 800  continue                                                                  
      if ( matchs('driver',4)   ) debug_driver   = .true.                       
      if ( matchs('elements',4) ) debug_elements = .true.                       
      if ( matchs('driver',4)   ) debug_driver   = .true.                       
      go to 10                                                                  
c                                                                               
c                                                                               
c               use 1 point integration for 8-node elements                     
c               ===========================================                     
c                                                                               
 900  continue                                                                  
      if ( matchs('1',1) ) then                                                 
        one_point_rule = .true.                                                 
      end if                                                                    
      go to 10                                                                  
c                                                                               
c                                                                               
c               ignore crack face loading during J computations                 
c               ===============================================                 
c                                                                               
 1000 continue                                                                  
      if ( matchs('crack',3) ) call splunj                                      
      if ( matchs('face',3) )  call splunj                                      
      if ( matchs('loading',3) )  then                                          
        ignore_face_loads = .true.                                              
      end if                                                                    
      go to 10                                                                  
c                                                                               
c                                                                               
c               omit crack-front elements during fgm J computations             
c               ===================================================             
c                                                                               
 1100 continue                                                                  
      if ( matchs('crack',3)    ) call splunj                                   
      if ( matchs('front',3)    ) call splunj                                   
      if ( matchs('elements',4) ) call splunj                                   
      if ( matchs('for',3)      ) call splunj                                   
      if ( matchs('fgms',3)     ) call splunj                                   
      if ( matchs('no',2)       ) omit_crack_front_elems = .false.              
      if ( matchs('yes',3)      ) omit_crack_front_elems = .true.               
      go to 10                                                                  
c                                                                               
c               set flag to output j and i values to packet.                    
c               ======================================                          
c                                                                               
 1200 continue                                                                  
      if ( matchs('packets',6) ) call splunj                                    
      if ( matchs('j',1)       ) then                                           
         output_packet_j = .true.                                               
         if( matchs('i',1) ) output_packet_i = .true.                           
      end if                                                                    
      if ( matchs('i',1)       ) then                                           
         output_packet_i = .true.                                               
         if( matchs('j',1) ) output_packet_j = .true.                           
      end if                                                                    
      go to 10                                                                  
c                                                                               
c                                                                               
c               enter uniform crack-face tractions                              
c               ==================================                              
c                                                                               
 1300 continue                                                                  
      if ( matchs('face',4)      ) call splunj                                  
      if (matchs('tractions',5)  ) call splunj                                  
c                                                                               
c                      parse through traction specifications                    
c                                                                               
 1310 continue                                                                  
      if(matchs('tx',2)) then                                                   
         dofn = 1                                                               
         go to 1320                                                             
      end if                                                                    
      if(matchs('ty',2)) then                                                   
         dofn = 2                                                               
         go to 1320                                                             
      end if                                                                    
      if(matchs('tz',2)) then                                                   
         dofn = 3                                                               
         go to 1320                                                             
      end if                                                                    
      go to 10                                                                  
c                                                                               
c                       store the face traction                                 
c                                                                               
 1320 continue                                                                  
      if( matchs('=',1) ) call splunj                                           
      if( .not.numd(forval) ) then                                              
         call errmsg(57,dum,dums,dumr,dumd)                                     
      else                                                                      
         cf_traction_flags(dofn) = .true.                                       
         cf_tractions(dofn)      = forval                                       
      end if                                                                    
      if( matchs(',',1) ) call splunj                                           
      go to 1310                                                                
c                                                                               
c                                                                               
c               set flag for output of plane stress                             
c               or plane strain results for interaction                         
c               integral                                                        
c               =======================================                         
c                                                                               
 1400 continue                                                                  
      stress_flag = .false.                                                     
      strain_flag = .false.                                                     
      if ( matchs('stress',6) ) then                                            
         stress_flag = .true.                                                   
         if ( matchs('plane',5)  ) call splunj                                  
         if ( matchs('strain',6) ) strain_flag = .true.                         
      end if                                                                    
      if ( matchs('strain',6) ) then                                            
         strain_flag = .true.                                                   
         if ( matchs('plane',5)  ) call splunj                                  
         if ( matchs('stress',6) ) stress_flag = .true.                         
      end if                                                                    
      if ( matchs('output',3)       ) call splunj                               
      if ( matchs('for',3)          ) call splunj                               
      if ( matchs('interaction',11) ) call splunj                               
      if ( matchs('integral',8)     ) call splunj                               
      if( .not.stress_flag ) out_pstress = .false.                              
      if( .not.strain_flag ) out_pstrain = .false.                              
      go to 10                                                                  
c                                                                               
c               compute stress intensity factors from J-values                  
c               ==============================================                  
c                                                                               
 1500 continue                                                                  
      if( matchs('to',2) ) call splunj                                          
      if( matchs('k',1) )  j_to_k = .true.                                      
      go to 10                                                                  
c                                                                               
c               set tag for crack-front curvature                               
c               =================================                               
c                                                                               
 1600 continue                                                                  
      if( matchs('circular',4) ) crack_curvature(2) = 1.0                       
      go to 10                                                                  
c                                                                               
c               dump domain definition                                          
c               ======================                                          
c                                                                               
 2000 continue                                                                  
      write(out,9000) domain_id(1:24)                                           
      if ( num_front_nodes .gt. 0 ) then                                        
        write(out,9100) num_front_nodes, front_nodes(1:num_front_nodes)         
      else                                                                      
        write(out,9110)                                                         
      end if                                                                    
      write(out,9120) crack_plane_normal                                        
      write(out,9720) tangent_vector_defined                                    
      if( tangent_vector_defined ) write(out,9710) crack_front_tangent          
      write(out,9130) symmetric_domain                                          
      write(out,9140) front_order, verify_front, domain_type                    
      write(out,9157) box_tol_relat                                             
      write(out,9150) print_elem_values, print_totals                           
      write(out,9155) one_point_rule                                            
      write(out,9156) ignore_face_loads                                         
      write(out,9160) qvals_given, rings_given                                  
      if ( qvals_given ) then                                                   
         write(out,9170)                                                        
         do list_entry = 1, last_compr                                          
           node  = compr_q_list(list_entry,1)                                   
           iword = compr_q_list(list_entry,2)                                   
           write(out,9180) node, rword                                          
         end do                                                                 
         write(out,9190)                                                        
         do elemno = 1, noelem                                                  
           word = ( elemno-1 ) / 30  + 1                                        
           bit  = elemno - ( word-1 ) * 30                                      
           if ( iand( bits(bit),q_element_maps(word) ) .ne. 0 )                 
     &         write(out,9200) elemno                                           
         end do                                                                 
      end if                                                                    
      if ( rings_given ) then                                                   
         write(out,9300)                                                        
         write(out,9310)                                                        
         do ring_num = 1, 300                                                   
           if ( ring_list(ring_num) .eq. 1 ) then                               
              write(out,9320) ring_num                                          
           end if                                                               
         end do                                                                 
      end if                                                                    
      write(out,9210) debug_driver, debug_elements                              
      write(out,9500)                                                           
c                                                                               
      go to 10                                                                  
c                                                                               
c        node sets to define the initially blunted crack tip                    
c        at a front position                                                    
c        =====================================================                  
c                                                                               
 2500 continue                                                                  
c                                                                               
      kount_set = 0                                                             
c                                                                               
      if ( matchs('set',3) ) call splunj                                        
      if ( .not. numi( node_set_id ) )                                          
     &           call errmsg(261,param,dums,dumr,dumd)                          
c                                                                               
c                                                                               
c               set-up list of nodes to be parsed.                              
c                                                                               
      call scan                                                                 
      call trlist(intlst,mxlsz,0,lenlst,errnum)                                 
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in the last three cases, command will be ig-            
c                       nored and a new domain parameter command                
c                       will be sought.                                         
c                                                                               
            if ( errnum .eq. 2 ) then                                           
               param = 1                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 3 ) then                                      
               param = 2                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else if ( errnum .eq. 4 ) then                                      
               param = 4                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
               go to 10                                                         
            else                                                                
               if ( errnum .ne. 1 ) then                                        
                  param = 3                                                     
                  call errmsg(24,param,dums,dumr,dumd)                          
                  go to 10                                                      
               end if                                                           
            end if                                                              
c                                                                               
c                      we got a valid integer list.                             
c                                                                               
      if ( true( dummy ) ) call splunj                                          
      call backsp( 1 )                                                          
c                                                                               
      icn    = 0                                                                
      iplist = 1                                                                
      if ( .not. allocated( node_set ) ) then                                   
         allocate( node_set(max_node_set,max_ct_node) )                         
         node_set(1:max_node_set,1:max_ct_node) = 0                             
      end if                                                                    
 2510 call trxlst(intlst,lenlst,iplist,icn,node)                                
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
      if ( node .gt. nonode) then                                               
         param = node                                                           
         call errmsg(259,param,dums,dumr,dumd)                                  
         go to 2520                                                             
      end if                                                                    
c                                                                               
c                       check that the list node is positive.                   
c                                                                               
      if ( node .le. 0 ) then                                                   
         param = node                                                           
         call errmsg(260,param,dums,dumr,dumd)                                  
         go to 2520                                                             
      end if                                                                    
c                                                                               
c                       check that the node set id is less than                 
c                       maximum allowable node sets.                            
c                                                                               
      if ( node_set_id .gt. max_node_set ) then                                 
         dumr = node_set_id                                                     
         call errmsg(263,param,dums,dumr,dumd)                                  
         go to 2520                                                             
      end if                                                                    
c                                                                               
c                       check that the number of ct nodes is less than          
c                       maximum allowable ct nodes.                             
c                                                                               
      if ( kount_set+1 .gt. max_ct_node ) then                                  
         dumr = kount_set+1                                                     
         call errmsg(262,param,dums,dumr,dumd)                                  
         go to 2520                                                             
      end if                                                                    
c                                                                               
      kount_set = kount_set + 1                                                 
      node_set(node_set_id,kount_set) = node                                    
      max_node_set_id = max( max_node_set_id, node_set_id )                     
c                                                                               
 2520 if ( iplist .ne. 0 ) go to 2510                                           
      call indom_chk_dup_nodes( node_set_id, out, dup_error )                   
      if( dup_error .gt. 0 ) input_ok = .false.                                 
      if ( true( dummy ) ) call splunj                                          
      call backsp( 1 )                                                          
c                                                                               
      go to 10                                                                  
c                                                                               
c               direction cosines for tangent vector at front.                  
c               =============================================                   
c                                                                               
 2600 continue                                                                  
      if ( matchs('vector',4) ) call splunj                                     
      tangent_vector_defined = .true.                                           
 2610 if ( endcrd(dummy) ) go to 10                                             
      if ( matchs('tx',2) ) then                                                
          if ( numd(crack_front_tangent(1)) ) go to 2610                        
          call errmsg(199,dum,dums,dumr,dumd)                                   
          go to 10                                                              
      end if                                                                    
      if ( matchs('ty',2) ) then                                                
          if ( numd(crack_front_tangent(2)) ) go to 2610                        
          call errmsg(199,dum,dums,dumr,dumd)                                   
          go to 10                                                              
      end if                                                                    
      if ( matchs('tz',2) ) then                                                
          if ( numd(crack_front_tangent(3)) ) go to 2610                        
          call errmsg(199,dum,dums,dumr,dumd)                                   
          go to 10                                                              
      end if                                                                    
      call errmsg(200,dum,dums,dumr,dumd)                                       
      go to 10                                                                  
c                                                                               
c                                                                               
c                                                                               
 9999 sbflg1 = .true.                                                           
      sbflg2 = .true.                                                           
c                                                                               
c                                                                               
 9000 format(//,5x,'Domain definition:',2x,a24,/)                               
 9100 format(8x,'Number of front nodes: ',i4,' Nodes on front:',                
     &     /,8x,10i7,/,8x,10i6,/,8x,10i6 )                                      
 9110 format(8x,'** no front nodes specified **')                               
 9120 format(8x,'Direction cosines: ',3f10.6 )                                  
 9130 format(8x,'Symmetric domain: ',l1 )                                       
 9140 format(8x,'Front order: ',i2,5x,'Verify front: ',l1,                      
     &          ' Domain type: ',i2 )                                           
 9150 format(8x,'Print element values: ',l1,5x,                                 
     &          'Print totals: ',l1 )                                           
 9155 format(8x,'Use 1 point rule for 8-node elements: ',l1)                    
 9156 format(8x,'Ignore crack face loading: ',l1)                               
 9157 format(8x,'box tolerance ratio: ',e13.6)                                  
 9160 format(8x,'Q-values specified: ',l1,' Rings specified: ',l1)              
 9170 format(8x,'Q-values by node: ')                                           
 9180 format(15x,i8,f10.3)                                                      
 9190 format(8x,'Elements included in domain:')                                 
 9200 format(15x,i8)                                                            
 9210 format(8x,'Debug driver: ',l1,5x,'Debug elements: ',l1)                   
 9300 format(8x,'Q-values by automatic rings: ')                                
 9310 format(8x,'Rings to be computed/printed:')                                
 9320 format(15x,i8)                                                            
 9500 format(/,5x,'** end of domain definition **' )                            
 9700 format ( /'Node set ',i4,' :',20i6)                                       
 9710 format(8x,'Tangent vector direction cosines: ',3f10.6 )                   
 9720 format(8x,'Tangent vector defined: ',l1 )                                 
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine indom_chk_dup_node           *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 2/17/10  rhdcw             *          
c     *                                                              *          
c     *     input of a node set completed. check for duplicates      *          
c     *     in the list                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine indom_chk_dup_nodes( set_id, out, dup_error )                  
      use j_data, only : node_set, max_ct_node                                  
      implicit integer (a-z)                                                    
c                                                                               
      logical found_duplicate                                                   
c                                                                               
      dup_error = 0                                                             
c                                                                               
c           count number of nodes in this node set on front                     
c                                                                               
      count = 0                                                                 
      do i = 1, max_ct_node                                                     
       if( node_set(set_id,i) .gt. 0 ) count = count + 1                        
      end do                                                                    
c                                                                               
c            traverse list looking for duplicates. simple (slow)                
c            algortihm but lists are very short (maybe 20-50 entries).          
c                                                                               
      found_duplicate = .false.                                                 
      do i = 1, count                                                           
        base_node = node_set(set_id,i)                                          
        do j = 1, count                                                         
          if ( i .eq. j ) cycle                                                 
          if( node_set(set_id,j) .ne. base_node ) cycle                         
          if( .not. found_duplicate ) then                                      
             write(out,9000) set_id                                             
             dup_error = 1                                                      
             found_duplicate = .true.                                           
          end if                                                                
          write(out,9010) base_node, j                                          
        end do                                                                  
       end do                                                                   
c                                                                               
      write(out,*) " "                                                          
      return                                                                    
c                                                                               
 9000 format(                                                                   
     & /1x,'>>>>> fatal error: 1 or more duplicate nodess found in',            
     & /,        '                    list for node set: ',i4 )                 
 9010 format(10x,'      node:',i9,                                              
     & ' has duplicate in list position: ',i4)                                  
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
