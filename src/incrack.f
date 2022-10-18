c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine incrack                      *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 10/18/22 rhd               *          
c     *                                                              *          
c     *                   input crack growth parameters              *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine incrack( sbflg1, sbflg2 )                                      
c
      use global_data, only : out, noelem, num_error, iprops, use_mpi,
     &                        dstmap
      use constants
      use allocated_integer_list, only : trlist_allocated
      use main_data, only : cnstrn_in                                           
      use damage_data
c                                                                               
      implicit none
c      
      logical :: sbflg1, sbflg2  
c 
c             locals
c           
      integer :: dof, elem, icn, iplist, lcktype, node, param,
     &           errnum, temp_int, dum, scan_order_list_size,
     &           scan_kill_order_length, scan_master_length,
     &           list_size
      logical, external :: matchs, numd, numi, endcrd,
     &                     true, matchs_exact, integr
      logical :: debug
      integer, allocatable :: scan_order_list(:), 
     &                        scan_kill_order_list(:),
     &                        scan_master_list(:)
c                                                                               
      double precision :: dumd, ctoa, temp_dble         
      double precision, parameter :: def_min_load_fact = 0.01d0                                                    
      character(len=1) :: dums                                                  
      real :: dumr                                                                 
c                                                                               
c         do some initialization                                                
c                                                                               
      debug = .false.                                                           
      if ( debug ) write (out,*) '>> inside incrack.'                           
c                                                                               
c         if we are back in the routine from just being here, then              
c         it is because there was a command line error.                         
c                                                                               
      if ( sbflg1 ) call errmsg( 201, dum, dums, dumr, dumd)                    
c                                                                               
c         branch on the crack growth parameter command                          
c                                                                               
 10   call readsc                                                               
c                                                                               
      if( matchs('type',4)    )        go to 200                               
      if( matchs('critical',4) )       go to 300                               
      if( matchs('release',4) )        go to 400                               
      if( matchs('print',5) )          go to 500                               
      if( matchs_exact('kill') )       go to 600                               
      if( matchs('sequential',4) )     go to 600                               
      if( matchs('force',4) )          go to 700                               
      if( matchs('crack',4) )          go to 800                               
      if( matchs('cell',3) )           go to 900                               
      if( matchs('angle',3) )          go to 1000                              
      if( matchs('dump',4) )           go to 1100                              
      if( matchs('characteristic',4) ) go to 1200                              
      if( matchs('alpha',4) )          go to 1300                              
      if( matchs('beta',4) )           go to 1400                              
      if( matchs_exact('gamma') )      go to 1450                              
      if( matchs('overshoot',4) )      go to 100                               
      if( matchs('automatic',4) )      go to 1500                              
      if( matchs('adaptive',4) )       go to 1500                              
      if( matchs('constant',5) )       go to 1600                              
      if( matchs('master',4) )         go to 1700                              
      if( matchs('number',4) )         go to 1800                              
      if( matchs('enforce',4) )        go to 1900                              
      if( matchs_exact('ppr') )        go to 2000
      if( matchs_exact('smcs') )       go to 2100
      if( matchs_exact('killed') )     go to 2200
      if( matchs_exact('stop') )       go to 2300
      if( matchs_exact('mesh') )       go to 2400
      if( matchs('regularization',7) ) go to 2400
      if( matchs_exact('use') )        go to 2500
      if( matchs('distortion',6) ) then
        call reset
        if( true() ) call splunj
        go to 2500
      end if
      if( matchs_exact('gurson') )     go to 2600
c                                                                               
      go to 9999                                                                
c                                                                               
c          ---------------------                                                
c          | Overshoot control |                                                
c          ---------------------                                                
c                                                                               
 100  continue                                                                  
c                                                                               
      if ( matchs('control',4)) call splunj                                     
c                                                                               
 110  continue                                                                  
      if ( matchs('off',3))     goto 120                                        
      if ( matchs('on',2))      goto 130                                        
      if ( matchs('percent',3)) goto 140                                        
      if ( matchs('minimum',3)) goto 150                                        
      if ( endcrd(dum) )        goto 10                                         
c                                                                               
      call errmsg( 201, dum, dums, dumr, dumd )                                 
c                                                                               
      goto 10                                                                   
c                                                                               
c          ========= turn on overshoot control =========                        
c                                                                               
 120  continue                                                                  
      overshoot_control_crk_grth = .false.                                      
      goto 110                                                                  
c                                                                               
c          ========= turn off overshoot control =========                       
c                                                                               
 130  continue                                                                  
      overshoot_control_crk_grth = .true.                                       
      goto 110                                                                  
c                                                                               
c          ========= set percent overshoot =========                            
c                                                                               
 140  continue                                                                  
      if ( matchs('overshoot',4)) call splunj                                   
      if ( .not. numd(overshoot_limit)) then                                    
         call errmsg ( 201, dum, dums, dumr, overshoot_limit )                  
      else                                                                      
         overshoot_limit = overshoot_limit / hundred                            
      endif                                                                     
      goto 110                                                                  
c                                                                               
c          ========= set minimum reduction =========                            
c                                                                               
 150  continue                                                                  
      if ( matchs('reduction',3)) call splunj                                   
      if ( .not. numd(min_load_fact)) then                                      
         call errmsg ( 201, dum, dums, dumr, min_load_fact )                    
      else                                                                      
         min_load_fact = min_load_fact / hundred                                
       if (min_load_fact .le. zero ) then                                       
            min_load_fact = def_min_load_fact                                   
          call errmsg ( 269, dum, dums, dumr, def_min_load_fact )               
         endif                                                                  
      endif                                                                     
      goto 110                                                                  
c                                                                               
c                                                                               
c          ---------------------                                                
c          | Crack Growth Type |                                                
c          ---------------------                                                
c                                                                               
c      type of crack growth <option>                                            
c                                                                               
c      <option> ::=   none | element_extinction | gurson |                      
c                     smcs | node_release | discrete |                          
c                     cohesive | user                                                 
c                                                                               
c      element_extinction maybe used as an optional phrase                      
c      ahead of gurson, smcs and user                                                
c                                                                               
c      element_extinction with nothing afterward implies gurson                 
c                                                                               
 200  continue                                                                  
      if ( matchs('of',2)     ) call splunj                                     
      if ( matchs('crack',3)  ) call splunj                                     
      if ( matchs('growth',4) ) call splunj                                     
c                                                                               
      if ( matchs('none',4) ) then                                              
         if ( crack_growth_type .ne. 0 ) call errmsg( 63, dum,                  
     &        dums, dumr, dumd )                                                
         crack_growth_type  = 0                                                 
         growth_by_release  = .false.                                           
         growth_by_kill     = .false.                                           
         growth_by_cohesive = .false.                                           
         go to 10                                                               
      end if                                                                    
      lcktype = 0                                                               
      if ( matchs('element_extinction',4) ) then                                
             lcktype = 1                                                        
             if ( matchs('gurson',4) ) call splunj                              
             if ( matchs('smcs',4)   ) lcktype = 3      
             if ( matchs('user',4)   ) lcktype = 5                        
      else if ( matchs('gurson',4) ) then                                       
             lcktype = 1                                                        
      else if ( matchs('smcs',4) ) then                                         
             lcktype = 3                                                        
      else if ( matchs('node_release',4) ) then                                 
             lcktype = 2                                                        
      else if ( matchs('discrete',4) ) then                                     
             lcktype = 2                                                        
      else if ( matchs('cohesive',4) ) then                                     
             lcktype = 4                                                        
      else if ( matchs('user',4) ) then                                     
             lcktype = 5    
      else                                                    
         call errmsg( 200, dum, dums, dumr, dumd )                              
         crack_growth_type  = 0                                                 
         growth_by_release  = .false.                                           
         growth_by_kill     = .false.                                           
         growth_by_cohesive = .false.                                           
         go to 10                                                               
      end if                                                                    
c              
      if ( lcktype .eq. 1 .or. lcktype .eq. 3                                   
     &      .or. lcktype .eq. 4 .or. lcktype .eq. 5 ) then                                          
c                                                                               
c                element extinction is chosen. if crack growth is               
c                off and elements were killed, then error. If node_release      
c                crack growth is on, also error.  Else,                         
c                extinction on.                                                 
c                                                                               
         if( crack_growth_type .eq. 0 ) then                                    
            if( .not. no_released_nodes ) then                                  
               call errmsg( 249, dum, dums, dumr, dumd )                        
            else if ( no_killed_elems ) then                                    
               crack_growth_type  = lcktype                                     
               growth_by_kill     = .true.                                      
               growth_by_release  = .false.                                     
               growth_by_cohesive = lcktype .eq. 4                              
               call incrack_dam_init ( debug )     
               if( lcktype .eq. 3 ) call allocate_damage ( 13 ) 
            else                                                                
               call errmsg( 214, dum, dums, dumr, dumd )                        
            endif                                                               
         else if ( crack_growth_type .eq. 2 ) then                              
            call errmsg( 249, dum, dums, dumr, dumd )                           
         end if                                                                 
c                                                                               
      elseif ( lcktype .eq. 2 ) then                                            
c                                                                               
c                 node_release crack growth is chosen.  If element              
c                 extinction is on, or if any elements are killed, then         
c                 error, otherwise turn node_release on.                        
c                                                                               
c                 Note: the data structures for node release are                
c                 allocated in the subroutine errchk.                           
c                                                                               
         if ( ( crack_growth_type.eq.0 .and. .not.no_killed_elems )             
     &        .or. growth_by_kill ) then                                        
            call errmsg( 249, dum, dums, dumr, dumd )                           
         else                                                                   
            crack_growth_type  = 2                                              
            growth_by_release  = .true.                                         
            growth_by_kill     = .false.                                        
            growth_by_cohesive = .false.                                        
         end if                                                                 
      end if                                                                    
c                                                                               
      go to 10                                                                  
c                                                                               
c          ------------------------------------------------------               
c          | Porosity limit -- element extinction option        |               
c          | Critical stress fraction to kill cohesive elements |               
c          ------------------------------------------------------               
c                                                                               
 300  continue                                                                  
      if ( matchs( 'effective', 4 ) ) call splunj                               
      if ( matchs( 'cohesive', 4 ) ) go to 350                                  
      if ( matchs( 'porosity', 4 ) ) call splunj                                
c                                                                               
c                make sure its a number                                         
c                                                                               
      if ( .not. numd ( porosity_limit ) ) then                                 
         call errmsg( 198, dum, dums, dumr, dumd )                              
       go to 10                                                                 
      end if                                                                    
c                                                                               
c                make sure its more than zero                                   
c                                                                               
      if ( porosity_limit .lt. zero ) then                                      
         call errmsg( 199, dum, dums, dumr, dumd )                              
       porosity_limit = zero                                                    
      end if                                                                    
c                                                                               
      go to 10                                                                  
c                                                                               
 350  continue                                                                  
      if ( matchs( 'displacement', 4 ) ) call splunj                            
      if ( matchs( 'multiplier', 4 ) ) call splunj                              
c                                                                               
c                make sure its a number                                         
c                                                                               
      if ( .not. numd ( critical_cohes_deff_fract ) ) then                      
         call errmsg( 329, dum, dums, dumr, dumd )                              
       go to 10                                                                 
      end if                                                                    
c                                                                               
c                make sure its more than 1                                      
c                                                                               
      if ( critical_cohes_deff_fract .le. one ) then                            
         call errmsg( 330, dum, dums, dumr, dumd )                              
       critical_cohes_deff_fract = 5.0                                          
      end if                                                                    
c                                                                               
      go to 10                                                                  
c                                                                               
c          -----------------------------------------------                      
c          | Number of release steps or release fraction |                      
c          -----------------------------------------------                      
c                                                                               
 400  continue                                                                  
c                                                                               
c                don't change number if either elements have been killed        
c                or nodes have been released.                                   
c                                                                               
      if ( .not. no_killed_elems .or. .not. no_released_nodes) then             
         call errmsg( 212, dum, dums, dumr, dumd )                              
         go to 10                                                               
      end if                                                                    
c                                                                               
      if ( matchs('steps',4) ) then                                             
        if ( .not. numi(temp_int) ) then                                        
          call errmsg( 203, max_dam_state, dums, dumr, dumd)                    
        else                                                                    
           if ( temp_int .le. 0 ) then                                          
              call errmsg( 202, max_dam_state, dums, dumr, dumd)                
           else                                                                 
              max_dam_state = temp_int                                          
           end if                                                               
        endif                                                                   
c                                                                               
      else if ( matchs('fraction',3) ) then                                     
         if ( .not. numd(release_fraction) ) then                               
            call errmsg( 239, dum, dums, dumr, dumd )                           
         end if                                                                 
c                                                                               
      else                                                                      
         call errmsg( 240, dum, dums, dumr, dumd )                              
      end if                                                                    
      go to 10                                                                  
c                                                                               
c          -----------------------------------------------------------          
c          | Print the status of the killed elements or crack front, |          
c          | or print the initial crack front or crack front nodes   |  
c          | or print top list of killable elements                  |        
c          -----------------------------------------------------------          
c                                                                               
 500  continue                                                                  
c                                                                               
c           branch on the proper command                                        
c                                                                               
 510  continue                                                                  
      if ( matchs('nodes',4) )   call splunj                                    
      if ( matchs('crack',5) )   call splunj                                    
      if ( matchs('status',4) )  call splunj                                    
c                                                                               
      if ( matchs('plane',5) )  goto 520                                        
      if ( matchs('front',5) )  goto 530                                        
      if ( matchs('off',3) )    goto 540                                        
      if ( matchs('on',2) )     goto 545  
      if ( matchs('top',3) )    goto 570                                      
      if ( endcrd(dum) )        goto 10                                         
c                                                                               
      call errmsg( 201, dum, dums, dumr, dumd )                                 
c                                                                               
      goto 10                                                                   
c                                                                               
c          ======= turn on printing of crack plane nodes =======                
c                                                                               
 520  continue                                                                  
      list_crkpln_nodes = .true.                                                
      goto 510                                                                  
c                                                                               
c          ======= turn on printing of crack front nodes =======                
c                                                                               
 530  continue                                                                  
      list_crkfrnt_nodes = .true.                                               
      goto 510                                                                  
c                                                                               
c          ======= turn off printing of crack status =======                    
c                                                                               
 540  continue                                                                  
      print_status = .false.   
      print_top_list = .false.                                                 
      goto 510                                                                  
c                                                                               
c          ======= turn on printing of crack status =======                     
c                                                                               
 545  continue                                                                  
      print_status = .true.        
      print_top_list = .false.                                             
c                                                                               
c            if we are using node release crack growth, we cant                 
c            have a list, so search for another command on the line.            
c                                                                               
      if ( growth_by_release  ) go to 510                                       
c                                                                               
c            input the order of the elements to be printed.  If not             
c            specified, then print all.                                         
c                                                                               
      if ( .not. matchs('order',5) ) go to 550                                  
      if ( matchs('elements',4) ) then                                          
        call splunj                                                             
      else                                                                      
        call backsp( 1 )                                                        
      end if                                                                    
      call scan
      if( allocated(scan_order_list) ) deallocate( scan_order_list )
      allocate( scan_order_list(20) )
      call trlist_allocated( scan_order_list, list_size, noelem,                              
     &                       scan_order_list_size, errnum )                               
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in these last 3 cases, the rest of the card             
c                       will be ignored and a new card will be sought.          
c                                                                               
      if ( errnum.eq.2 ) then                                                   
         param= 1                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 550                                                              
      else if ( errnum.eq.3 ) then                                              
         param= 2                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 550                                                              
      else if ( errnum.eq.4 ) then                                              
         param=4                                                                
         call errmsg ( 24, param, dums, dumr, dumd )                            
         go to 550                                                              
      else                                                                      
         if ( errnum.eq.1 ) then                                                
            call backsp ( 1 )                                                   
            go to 560                                                           
         end if                                                                 
         param= 3                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 550                                                              
      end if                                                                    
c                                                                               
c           if order was not specified, or list had an error, simulate all      
c                                                                               
 550  continue                                                                  
      call errmsg( 213, dum, dums, dumr, dumd )                                 
      scan_order_list(1) = 1                                                    
      scan_order_list(2) = -noelem                                              
      scan_order_list(3) = 1                                                    
      scan_order_list_size = 3                                                  
c                                                                               
c           translate list into a dynamically allocated array -- first          
c           must parse the list to check if all any elements are not            
c           killable                                                            
c                                                                               
 560  continue                                                                  
      iplist = 1                                                                
      icn = 0                                                                   
      num_print_list= 0                                                         
 565  continue                                                                  
      call trxlst( scan_order_list, scan_order_list_size, iplist,               
     &     icn, elem )  
      if( elem == 0 .or. elem > noelem ) then 
       write(out,9100) elem
       num_error = num_error + 1
       go to 10
      end if                                                       
      if( growth_by_kill ) then                                                 
         if ( iand (iprops(30,elem),2).ne.0 ) then                              
              num_print_list = num_print_list + 1                               
         end if                                                                 
      end if                                                                    
      if ( iplist .ne. 0 ) go to 565                                            
      if ( num_print_list .le. 0 ) then                                         
         call errmsg( 225, dum, dums, dumr, dumd )                              
         print_status = .false.                                                 
         go to 510                                                              
      end if                                                                    
c                                                                               
c           now allocate and fill the order array                               
c                                                                               
      call allocate_damage( 2 )                                                 
      call incrack_print_list_fill( scan_order_list,
     &                              scan_order_list_size, debug )                                                          
c                                                                               
      goto 510  
c
 570  continue
      if( .not. integr( num_top_list ) ) then
         call incrack_errmsg( 7 )
         print_top_list = .false.
         num_top_list = 0
         go to 10
      end if
      print_status = .false.
      print_top_list = .true.
      go to 10
                                                                
c                                                                               
c                                                                               
c          ------------------------------------------------------               
c          | Kill elements in order -- element extinction option|               
c          ------------------------------------------------------               
c                                                                               
c                                                                               
 600  continue                                                                  
      if ( matchs('sequentially',6)) call splunj                                
      if ( matchs('extinction',5)) call splunj                                  
c                                                                               
                                                                                
      if ( matchs('on',2) ) then                                                
         kill_order = .true.                                                    
      else if ( matchs('off',3) ) then                                          
         kill_order = .false.                                                   
         goto 10                                                                
      else                                                                      
         call errmsg( 201, dum, dums, dumr, dumd )                              
         kill_order = .false.                                                   
         goto 10                                                                
      endif                                                                     
c                                                                               
c                       now read the list of the order in which to              
c                       kill elements                                           
c                                                                               
      if ( .not.matchs('order',5) ) goto 650                                    
c                                                                               
      call scan    
      if( allocated( scan_kill_order_list ) ) 
     &          deallocate( scan_kill_order_list )
      allocate( scan_kill_order_list(20) )
      call trlist_allocated( scan_kill_order_list, list_size, noelem,                         
     &                       scan_kill_order_length, errnum )                                     
c                                                                               
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in these last 3 cases, the rest of the card             
c                       will be ignored and a new card will be sought.          
c                                                                               
      if ( errnum.eq.2 ) then                                                   
         param= 1                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 650                                                              
      else if ( errnum.eq.3 ) then                                              
         param= 2                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 650                                                              
      else if ( errnum.eq.4 ) then                                              
         param=4                                                                
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 650                                                              
      else                                                                      
         if ( errnum.eq.1 ) then                                                
            call backsp( 1 )                                                    
            go to 660                                                           
         end if                                                                 
         param= 3                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 650                                                              
      end if                                                                    
c                                                                               
c           if order was not specified, or list had an error, simulate all      
c                                                                               
 650  continue                                                                  
      call errmsg( 213, dum, dums, dumr, dumd )                                 
      scan_kill_order_list(1) = 1                                               
      scan_kill_order_list(2) = -noelem                                         
      scan_kill_order_list(3) = 1                                               
      scan_kill_order_length = 3                                                
c                                                                               
c           translate list into a dynamically allocated array -- first          
c           must parse the list to check if all any elements are not            
c           killable                                                            
c                                                                               
 660  continue                                                                  
      iplist = 1                                                                
      icn = 0                                                                   
      num_kill_order_list = 0                                                   
 665  continue                                                                  
      call trxlst( scan_kill_order_list, scan_kill_order_length,                
     &     iplist, icn, elem )                                                  
      if ( iand (iprops(30,elem),2).ne.0 ) then                                 
         num_kill_order_list = num_kill_order_list + 1                          
      endif                                                                     
      if ( iplist.ne.0 ) goto 665                                               
      if ( num_kill_order_list.le.0 ) then                                      
         call errmsg( 225, dum, dums, dumr, dumd )                              
         kill_order = .false.                                                   
         goto 10                                                                
      endif                                                                     
c                                                                               
c           now allocate and fill the order array                               
c                                                                               
      call allocate_damage( 3 )                                                 
      call incrack_kill_order_fill( scan_kill_order_list,                               
     &          scan_kill_order_length, debug )                                 
c                                                                               
      goto 10                                                                   
c                                                                               
c                                                                               
c          -----------------------------                                        
c          |      force release type   |                                        
c          -----------------------------                                        
c              1 = fixed number of load steps                                   
c              2 = traction-separation law                                      
c                                                                               
 700  continue                                                                  
c                                                                               
c                don't change number if any elements have been killed or        
c                nodes have been released.                                      
c                                                                               
      if ( .not. no_killed_elems .or. .not. no_released_nodes ) then            
         call errmsg( 212, dum, dums, dumr, dumd )                              
         go to 10                                                               
      end if                                                                    
c        
      use_estiff_at_death = .false.
      use_mesh_regularization = .false.
c                                                                             
      if ( matchs('release',3) )  call splunj                                   
      if ( matchs('type',3)    )  call splunj                                   
      if ( matchs('steps',4)   ) then                                           
           release_type = 1                                                     
      else if ( matchs('traction',4) ) then                                     
           release_type = 2                                                     
      else                                                                      
           call errmsg( 234, dum, dums, dumr, dumd )                            
      end if                                                                    
      go to 10                                                                  
c                                                                               
c          -------------------------------------------------------------        
c          | crack plane normal direction, coordinate, search tolerance |       
c          -------------------------------------------------------------        
c              used for traction separation law release                         
c              look two words deep onto line to make sure we have               
c              crack plane normal command                                       
c                                                                               
 800  continue                                                                  
      if ( matchs('plane',3)  )  call splunj                                    
      if ( .not. matchs('normal',4) ) then                                      
         call reset                                                             
         if ( true(dum) ) call splunj                                           
         go to 9999                                                             
      end if                                                                    
      call reset                                                                
      if ( true(dum) ) call splunj                                              
      if ( matchs('crack',5) ) call splunj                                      
c                                                                               
c                don't change number if either elements have been killed        
c                or nodes have been released                                    
c                                                                               
      if ( .not. no_killed_elems .or. .not. no_released_nodes ) then            
         call errmsg( 251, dum, dums, dumr, dumd )                              
         go to 10                                                               
      end if                                                                    
c                                                                               
      crkpln_srch_tol = 1.0e-05                                                 
      if ( matchs('plane',3)  ) call splunj                                     
      if ( matchs('normal',4) ) then                                            
        crk_pln_normal_idx = 0                                                  
        crack_plane_coord  = 0.0                                                
        if ( matchs('x',1) ) crk_pln_normal_idx = 1                             
        if ( matchs('y',1) ) crk_pln_normal_idx = 2                             
        if ( matchs('z',1) ) crk_pln_normal_idx = 3                             
        if ( crk_pln_normal_idx .eq. 0 ) then                                   
           call errmsg( 235, dum, dums, dumr, dumd )                            
           go to 10                                                             
        end if                                                                  
        if ( matchs('coord',4) ) then                                           
           if ( .not. numd(crack_plane_coord) ) then                            
              call errmsg( 236, dum, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
        if ( endcrd(dum) ) go to 10                                             
        if ( matchs('search',4) ) call splunj                                   
        if ( matchs('tolerance',3) ) call splunj                                
        if ( .not. numd(crkpln_srch_tol) ) then                                 
           call errmsg( 325, dum, dums, dumr, dumd )                            
        end if                                                                  
      else                                                                      
        call errmsg( 237, dum, dums, dumr, dumd )                               
      end if                                                                    
      go to 10                                                                  
c                                                                               
c          --------------------------------------------                         
c          | cell height -- element extinction option |                         
c          --------------------------------------------                         
c              used for traction separation law release                         
c                                                                               
 900  continue                                                                  
c                                                                               
c                don't change number if any elements have been killed.          
c                                                                               
      if ( .not. no_killed_elems ) then                                         
         call errmsg( 212, dum, dums, dumr, dumd )                              
         go to 10                                                               
      end if                                                                    
c                                                                               
      if ( matchs('height',3) .or. matchs('hieght',3) ) then                    
        if ( .not. numd(gurson_cell_size) ) then                                
           call errmsg( 238, dum, dums, dumr, dumd )                            
           go to 10                                                             
        end if                                                                  
      else                                                                      
        call errmsg( 238, dum, dums, dumr, dumd )                               
      end if                                                                    
      go to 10                                                                  
c                                                                               
c          --------------------------------------------------                   
c          | angle for CTOA release  -- node release option |                   
c          --------------------------------------------------                   
c                                                                               
c           also sets the distance back to measure ctoa when                    
c           using constant front growth                                         
c                                                                               
 1000  continue                                                                 
c                                                                               
       if ( matchs ('for',3)     ) call splunj                                  
       if ( matchs ('initiation',3) ) then                                      
          if ( .not. numd(ctoa) ) then                                          
             call errmsg ( 289, dum, dums, dumr, dumd )                         
          else                                                                  
             init_crit_ang = ctoa / two                                         
          endif                                                                 
          if ( matchs ('distance',4) ) then                                     
             if ( .not. numd(init_ctoa_dist) ) then                             
                call errmsg ( 287, dum, dums, dumr, dumd )                      
             endif                                                              
          endif                                                                 
          go to 10                                                              
       else                                                                     
          if ( matchs ('release',3) ) call splunj                               
          if ( .not. numd(ctoa) ) then                                          
             call errmsg ( 289, dum, dums, dumr, dumd )                         
          else                                                                  
             critical_angle = ctoa / two                                        
          endif                                                                 
          if ( matchs ('distance',4) ) then                                     
             if ( .not. numd(ctoa_dist) ) then                                  
                call errmsg ( 287, dum, dums, dumr, dumd )                      
             endif                                                              
          endif                                                                 
          go to 10                                                              
       endif                                                                    
c                                                                               
c          ----------------------------------------                             
c          | dump all crack growth data to screen |                             
c          ----------------------------------------                             
c                                                                               
c                                                                               
 1100  continue                                                                 
       call dam_debug                                                           
       go to 10                                                                 
c                                                                               
c                                                                               
c          --------------------------------------------                         
c          | characteristic length for release height |                         
c          |         - node release option            |                         
c          --------------------------------------------                         
c                                                                               
c                                                                               
 1200  continue                                                                 
       if ( matchs('length',3) ) call splunj                                    
       if ( .not. numd(char_length) ) then                                      
          call errmsg( 254, dum, dums, dumr, dumd )                             
       endif                                                                    
       go to 10                                                                 
c                                                                               
c                                                                               
c                                                                               
c          ---------------------------------------------                        
c          | alpha and beta parameters for smcs growth |                        
c          |         - element kill option             |                        
c          ---------------------------------------------                        
c                                                                               
c                                                                               
 1300  continue                                                                 
       if ( .not. numd(smcs_alpha) ) then                                       
          call errmsg( 254, dum, dums, dumr, dumd )                             
       end if                                                                   
       go to 10                                                                 
 1400  continue                                                                 
       if ( .not. numd(smcs_beta) ) then                                        
          call errmsg( 254, dum, dums, dumr, dumd )                             
       end if                                                                   
       go to 10                                                                 
 1450  continue                                                                 
       if ( .not. numd(smcs_gamma) ) then                                        
          call errmsg( 254, dum, dums, dumr, dumd )                             
       end if                                                                   
       go to 10                                                                 
c                                                                               
c          ----------------------------                                         
c          | automatic load reduction |                                         
c          ----------------------------                                         
c                                                                               
 1500 continue                                                                  
      if ( matchs ('load', 4    ) ) call splunj                                 
      if ( matchs ('reduction',5) ) call splunj                                 
      if ( matchs ('control',5) )   call splunj                                 
c                                                                               
 1510 continue                                                                  
      if ( matchs ('off', 3     ) ) goto 1520                                   
      if ( matchs ('on', 2      ) ) goto 1530                                   
      if ( matchs ('minimum', 3 ) ) goto 1510                                   
      if ( matchs ('maximum', 3 ) ) goto 1510                                   
      if ( matchs ('steps', 4   ) ) goto 1540                                   
      if ( matchs ('porosity', 4) ) goto 1550                                   
      if ( matchs ('plastic', 4 ) ) goto 1510                                   
      if ( matchs ('strain', 3  ) ) goto 1560                                   
      if ( matchs ('relative',4 ) ) goto 1570   
      if ( matchs_exact('alpha_min') ) goto 1575                               
      if ( matchs_exact('alpha_max') ) goto 1580   
      if ( matchs ( 'allowable_release', 4 ) ) goto 1585                        
      if ( endcrd ( dum         ) ) goto 10                                     
c                                                                               
      call errmsg( 201, dum, dums, dumr, dumd )                                 
c                                                                               
      goto 10                                                                   
c                                                                               
c          ========= turn off auto stepping =========                           
c                                                                               
 1520 continue                                                                  
      load_size_control_crk_grth = .false.                                      
      goto 1510                                                                 
c                                                                               
c          ========= turn on auto stepping =========                            
c                                                                               
 1530 continue                                                                  
      load_size_control_crk_grth = .true.                                       
      goto 1510                                                                 
c                                                                               
c          ========= set minimum steps between releases =========               
c                      (for CTOA crack growth)                                  
c                                                                               
 1540  continue                                                                 
      if ( .not. numi( temp_int ) ) then                                        
         call errmsg ( 201, dum, dums, dumr, dumd )                             
      else                                                                      
       if ( temp_int .le. 0 ) then                                              
          call errmsg ( 272, min_steps_for_release, dums, dumr,                 
     &           dumd)                                                          
         else                                                                   
            min_steps_for_release = temp_int                                    
         end if                                                                 
      end if                                                                    
      goto 1510                                                                 
c                                                                               
c                                                                               
c          ========= set % max porosity change between steps =========          
c                      (for gurson crack growth)                                
c                                                                               
 1550  continue                                                                 
      if ( matchs ('change', 6) ) call splunj                                   
c                                                                               
      if ( .not. numd( temp_dble ) ) then                                       
         call errmsg ( 201, dum, dums, dumr, dumd )                             
      else                                                                      
       if ( temp_dble .le. 0 ) then                                             
          call errmsg ( 275, dum, dums, dumr, dumd)                             
         else                                                                   
            max_porosity_change = temp_dble                                     
         end if                                                                 
      end if                                                                    
      goto 1510                                                                 
c                                                                               
c                                                                               
c          ========= set max plastic strain change between steps =========      
c                      (for smcs crack growth)                                  
c                                                                               
 1560 continue                                                                  
      if ( matchs ('change', 6) ) call splunj                                   
c                                                                               
      if ( .not. numd( temp_dble ) ) then                                       
         call errmsg ( 201, dum, dums, dumr, dumd )                             
      else                                                                      
       if ( temp_dble .le. 0 ) then                                             
          call errmsg ( 276, dum, dums, dumr, dumd)                             
         else                                                                   
            max_plast_strain_change = temp_dble                                 
         end if                                                                 
      end if                                                                    
      goto 1510                                                                 
c                                                                               
c          ========= set max change in effective relative                       
c                    displacement between steps =========                       
c                      (for cohesive crack growth)                              
c                  this is a fraction (e.g. = 0.20) means                       
c                  max change in effective displacement per                     
c                  step must be <= 20% of displacement at the                   
c                  peak stress value on the traction-separation                 
c                  curve.  Default value in startup is 0.2                      
c                                                                               
 1570 continue                                                                  
      if ( matchs( 'displacement', 6) ) call splunj                             
      if ( matchs('change', 3) ) call splunj                                    
c                                                                               
      if ( .not. numd( temp_dble ) ) then                                       
         call errmsg ( 201, dum, dums, dumr, dumd )                             
      else                                                                      
       if ( temp_dble .le. 0.0 ) then                                           
          call errmsg ( 328, dum, dums, dumr, dumd)                             
         else                                                                   
            max_deff_change = temp_dble                                         
         end if                                                                 
      end if                                                                    
      goto 1510 
c                                                                               
c          ========= additional adaptive parameters for SMCS
c
 1575 continue                                                                  
      if( .not. numd( temp_dble ) ) then                                       
         call incrack_errmsg( 23 )
      else                                                                      
       if( temp_dble .le. 0.0d0 .or. temp_dble .ge. 1.d0 ) then                                           
          call incrack_errmsg( 24 )
         else                                                                   
            smcs_adapt_alpha_min = temp_dble       
         end if                                                                 
      end if                                                                    
      goto 1510                                                                 
c
 1580 continue                                                                  
      if( .not. numd( temp_dble ) ) then                                       
         call incrack_errmsg( 25 )
      else                                                                      
       if( temp_dble .le. 1.0d0 ) then                                           
          call incrack_errmsg( 26 )
         else                                                                   
            smcs_adapt_alpha_max = temp_dble                                         
         end if                                                                 
      end if                                                                    
      goto 1510                                                                 
c
 1585 continue                                                                  
      if( .not. integr( temp_int ) ) then                                       
         call incrack_errmsg( 27 )
         go to 1510
      end if
      smcs_allowable_in_release = temp_int   
      goto 1510                                                                 
c                                                                               
c          ----------------------------                                         
c          | initialize const growth |                                          
c          ----------------------------                                         
c                                                                               
 1600 continue                                                                  
      if ( matchs('front',5) ) call splunj                                      
      if ( matchs('growth',4) ) call splunj                                     
      if ( matchs('on',2) ) then                                                
         const_front = .true.                                                   
      else if ( matchs('off',3) ) then                                          
         const_front = .false.                                                  
      else                                                                      
         call errmsg ( 201, dum, dums, dumr, dumd )                             
      endif                                                                     
c                                                                               
      goto 10                                                                   
c                                                                               
c          ------------------------                                             
c          | list of master nodes |                                             
c          ------------------------                                             
c                                                                               
 1700 continue                                                                  
c                                                                               
c             if crack plane normal has not yet been set, error                 
c                                                                               
      if ( crk_pln_normal_idx.lt.1 .or. crk_pln_normal_idx.gt.3) then           
         call errmsg (286, dum, dums, dumr, dumd)                               
         goto 10                                                                
      endif                                                                     
                                                                                
c                                                                               
      if ( matchs('node',4) ) call splunj                                       
      if ( matchs('list',4) ) call splunj                                       
c                                                                               
c             now read list of master nodes                                     
c     
      if( allocated( scan_master_list ) ) 
     &    deallocate( scan_master_list )
      allocate( scan_master_list(20) )                                                                         
      call scan                                                                 
      call trlist_allocated( scan_master_list, list_size, noelem,                             
     &     scan_master_length, errnum )                                         
c                                                                               
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       for errors 2 and 4, we will turn off const              
c                       growth.  For 3, just output a warning.                  
c                                                                               
      if ( errnum.eq.2 ) then                                                   
         param= 1                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 10                                                               
      else if ( errnum.eq.3 ) then                                              
         param= 2                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
      else if ( errnum.eq.4 ) then                                              
         param=4                                                                
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 10                                                               
      else                                                                      
         if ( errnum.eq.1 ) then                                                
            call backsp( 1 )                                                    
            go to 1760                                                          
         end if                                                                 
         param= 3                                                               
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 10                                                               
      end if                                                                    
c                                                                               
c           translate list into a dynamically allocated array -- first          
c           must parse the list to check if nodes are on crack plane            
c                                                                               
 1760 continue                                                                  
      iplist = 1                                                                
      icn = 0                                                                   
      num_crack_fronts = 0                                                      
 1765 continue                                                                  
      call trxlst( scan_master_list, scan_master_length,                        
     &     iplist, icn, node )                                                  
      dof = dstmap(node) + crk_pln_normal_idx - 1                               
      if (cnstrn_in(dof) .ne. d32460) then                                      
         num_crack_fronts = num_crack_fronts + 1                                
      else                                                                      
         call errmsg( 283, node, dums, dumr, dumd )                             
      endif                                                                     
      if ( iplist.ne.0 ) goto 1765                                              
c                                                                               
c           now allocate and fill the order array                               
c                                                                               
      if (num_crack_fronts .eq. 0) then                                         
         call errmsg( 282, dum, dums, dumr, dumd )                              
      else                                                                      
         call allocate_damage( 10 )                                             
         call incrack_master_list_fill( scan_master_list,                               
     &          scan_master_length, debug )                                     
      endif                                                                     
c                                                                               
      goto 10                                                                   
                                                                                
c                                                                               
c          ------------------------------------------                           
c          | number of nodes in thickness direction |                           
c          ------------------------------------------                           
c                                                                               
 1800 continue                                                                  
c                                                                               
      if (matchs ('of',2)) call splunj                                          
      if (matchs ('nodes',4)) call splunj                                       
      if (matchs ('along',3)) call splunj                                       
      if (matchs ('front',5)) call splunj                                       
      if (matchs ('through',5)) call splunj                                     
      if (matchs ('thickness',5)) call splunj                                   
c                                                                               
      if ( .not. numi(num_nodes_thick)) then                                    
         call errmsg( 103, dum, dums, dumr, dumd )                              
      endif                                                                     
c                                                                               
      goto 10                                                                   
c                                                                               
c          ------------------------------------------                           
c          | enforced node release on next step     |                           
c          ------------------------------------------                           
c                                                                               
 1900 continue                                                                  
c                                                                               
      if (matchs ('node',4))  call splunj                                       
      if (matchs ('release',4)) call splunj                                     
      if (matchs ('next',4)) call splunj                                        
      if (matchs ('step',4)) call splunj                                        
      enforce_node_release = .true.                                             
      go to 10                                                                  
c                                                                               
c          ------------------------------------------                           
c          | PPR option criterion/data              |                           
c          ------------------------------------------                           
c                                                                               
 2000 continue                                                                  
c                                                                               
c          ========= set fraction of limit displacement at                      
c                    element extinction. limit displacement                     
c                    occurs at zero tractions in PPR ==========                 
c                                                                               
      if( matchs('displacement',4))  call splunj                                
      if( matchs('fraction',4)) then                                            
        if( matchs('for',3)) call splunj                                        
        if( matchs('extinction',4)) call splunj                                 
        if( .not. numd( ppr_kill_displ_fraction ) ) then                        
           call incrack_errmsg( 1 )                                             
        end if                                                                  
        if( ppr_kill_displ_fraction .le. 0.0 .or.                               
     &      ppr_kill_displ_fraction .gt. 1.0 )                                  
     &      call incrack_errmsg( 3 )                                            
        go to 10                                                                
      end if                                                                    
      call incrack_errmsg( 2 )                                                  
      go to 10                                                                  
c                                                                               
c          ------------------------------------------                           
c          | SMCS options                           |                           
c          ------------------------------------------                           
c                                                                               
 2100 continue            
      if( matchs_exact('type') ) then                                 
        call incrack_smcs_parms
        go to 10
      end if
      if( matchs('states',5) ) then
        call incrack_smcs_states
        go to 10
      end if
      if( matchs_exact('deleted') ) then
        call incrack_smcs_deleted_list_file
        go to 10
      end if
      if( matchs_exact('removed') ) then
        call incrack_smcs_removed_list_file
        go to 10
      end if
      call incrack_errmsg( 8 )
      go to 10
c                                                                               
c          ------------------------------------------                           
c          | killed element limit                   |                           
c          ------------------------------------------                           
c                                                                               
 2200 continue                                                                  
c                                                                               
c          ========= killed element limit
c         
      if( matchs('element',3) ) call splunj
      if( matchs('limit',3) ) call splunj
      if( .not. numi( killed_element_limit ) ) then
            call incrack_errmsg( 6 )
      end if
      go to 10
c                                                                               
c          ------------------------------------------                           
c          | stop if elements killed <list>         |   
c          | stop of nodes released <list>          |                      
c          ------------------------------------------                           
c                                                                               
 2300 continue                                                                  
c                          
c              stop (if) killed (elements) < list >  
c              stop (if) nodes released <list>
c                 
      if( matchs_exact('if') ) call splunj
      if( matchs_exact('killed') ) then
          call incrack_stop_elist
          go to 10
      end if
      if( matchs('nodes',4) ) then
          call incrack_stop_nlist
          go to 10
      end if
      call incrack_errmsg( 21 )
      go to 10
c
c          ------------------------------------------                           
c          | (mesh) regularization                  |                           
c          ------------------------------------------                           
c                                                                               
 2400 continue   
c
c                not supported yet with MPI
c
      if( use_mpi ) then
         call incrack_errmsg( 65 )
         go to 10
      end if
                                                                     
c                                                                               
c                don't change number if elements have been killed        
c
      if ( .not. no_killed_elems ) then             
         call errmsg( 212, dum, dums, dumr, dumd )        
         call scan_flushline                      
         go to 10                                                               
      end if                                                                    
c
      release_type = -1
      call allocate_damage( 14 )
      if( matchs('regularization',7) ) call splunj
      call incrack_regularization                                                                       
      go to 10
c                                                                               
c
c          -----------------------------------------------------------                           
c          | use distortion metric on|off                            |
c          |      (write initial values file) 
c          |      Oddy critical ratio <value>                        | 
c          |      plastic strain limit <value>  display                     |
c          -----------------------------------------------------------                           
c                                                                               
 2500 continue   
c
c                not supported yet with MPI
c
      use_distortion_metric = .false.
      Oddy_print_initial = .false.
      distortion_plastic_limit = hundred  ! so never gets invoked 
c
      if( use_mpi ) then
         call incrack_errmsg( 66 )
         go to 10
      end if
c
      if( .not. matchs('distortion',6) ) then
         call incrack_errmsg( 67 )
         go to 10
      end if
c
      if( matchs('metrics',6) ) call splunj
      if( matchs_exact('off') ) go to 10
      if( .not. matchs_exact('on') ) then
         call incrack_errmsg( 68 )
         go to 10
      end if
c
      use_distortion_metric = .true.
      call allocate_damage( 15 )
      if( endcrd() ) go to 10  
c
      if( matchs_exact('display') ) write(out,9060) 
     &     use_distortion_metric, Oddy_print_initial, 
     &     Oddy_critical_ratio, distortion_plastic_limit 
c
      if( matchs_exact('write') ) then
        if( matchs('initial',4) ) call splunj
        if( matchs('values',3) ) call splunj
        if( .not. matchs_exact('file') ) then
         call incrack_errmsg( 73 )
         use_distortion_metric = .false.
         Oddy_print_initial = .false.
         call allocate_damage( 16 ) ! delete if allocated
         go to 10
        end if
        if( matchs_exact('off') ) go to 2510
        if( .not. matchs_exact('on') ) then
         call incrack_errmsg( 74 )
         use_distortion_metric = .false.
         Oddy_print_initial = .false.
         call allocate_damage( 16 ) ! delete if allocated
         go to 10
        end if
        Oddy_print_initial = .true.
        call allocate_damage( 17 ) 
      end if 
c
 2510 continue
      if( matchs_exact('display') ) write(out,9060) 
     &     use_distortion_metric, Oddy_print_initial, 
     &     Oddy_critical_ratio, distortion_plastic_limit 

c
      if( endcrd() ) go to 10
      if( .not. matchs_exact('Oddy') ) then
         call incrack_errmsg( 75 )
         use_distortion_metric = .false.
         Oddy_print_initial = .false.
         call allocate_damage( 16 ) ! delete if allocated
         go to 10
      end if
      if( matchs('critical',3) ) call splunj
      if( matchs('ratio',3) ) call splunj
      if( .not. numd(Oddy_critical_ratio) ) then
         call incrack_errmsg( 69 )
         use_distortion_metric = .false.
         Oddy_print_initial = .false.
         call allocate_damage( 16 ) ! delete if allocated
         go to 10
      end if
c
       if( matchs_exact('display') ) write(out,9060) 
     &     use_distortion_metric, Oddy_print_initial, 
     &     Oddy_critical_ratio, distortion_plastic_limit 
c
      if( endcrd() ) go to 10
c
      if( .not. matchs_exact('plastic') ) then
         call incrack_errmsg( 70 )
         use_distortion_metric = .false.
         Oddy_print_initial = .false.
         call allocate_damage( 16 ) ! delete if allocated
         go to 10 
      end if
      if( matchs('strain',3) ) call splunj
      if( matchs('limit',3) ) call splunj
      if( .not. numd(distortion_plastic_limit) ) then
         call incrack_errmsg( 71 )
         use_distortion_metric = .false.
         Oddy_print_initial = .false.
         call allocate_damage( 16 ) ! delete if allocated
         go to 10 
      end if
c
      if( matchs_exact('display') ) write(out,9060) 
     &     use_distortion_metric, Oddy_print_initial, 
     &     Oddy_critical_ratio, distortion_plastic_limit 
      go to 10
c                                                                               
c          ---------------------------------------------------------                           
c          | gurson deleted (elements) file on|off (name) <string> |                           
c          ---------------------------------------------------------                           
c                                                                               
 2600 continue            
      if( matchs_exact('deleted') ) then
        call incrack_gt_list_file
        go to 10
      end if
      call incrack_errmsg( 72 )
      go to 10
      
c
 9999 continue                                                                  
      if ( debug ) write (out,*) '>>>>>> leaving incrack.'           
      sbflg1 = .true.                                                           
      sbflg2 = .true.                                                           
      return  
c
 9000 format(/1x,'>>>>> error: no element list found',
     &   /,14x,'read new line ...',/)
 9010 format(/1x,'>>>>> error: integer list rules failed....',
     &   ' read new line ...',/)
 9020 format(/1x,'>>>>> error: too many elements in list....',
     &   ' read new line ...',/)
 9050 format(/1x,'>>>>> FATAL ERROR: delete_elements inconsistency',
     &   /,14x,'job terminated now ...',/)
 9060 format(/,5x,'... use distortion metrics: ',l1,
     &       /,5x,'... make file of initial Oddy values: ',l1,
     &       /,5x,'... critical Oddy ratio: ',f12.3,
     &       /,5x,'... plastic strain limit: ',f12.5, //)
 9100 format(/1x,'>>>>> error: invalid element number: ',i8)
c
      contains
c     ========
c
c          ------------------------------------------                           
c          | contains: incrack_regularization       |                           
c          ------------------------------------------                           
c                                                                                           
      subroutine incrack_regularization
c
      implicit none
c
c        on | off  type  <integer> L <value> up_max <value>,
c              num_points <integer> [ <up,d. pair] (display)
c    
c        up vlaues are normalized 0.0 <= up <= 1.0
c        d values are normalized  0.0 <= d <= 1.0
c        first and last points must be set to req'd limits.
c
      integer :: i
      logical :: ok, onflag, bad_type, bad
      double precision :: x,y 
      double precision, parameter :: local_tol = 0.00000001d0
c
      use_mesh_regularization = .false.
      use_estiff_at_death = .false.
      regular_type = -1
      regular_length = -1.0d0
      regular_npoints = 0
      regular_up_max = zero
c
      if( matchs_exact('off') )  then
        call scan_flushline
        return
      end if
      onflag = .false.
      if( matchs_exact('on') ) onflag = .true.
      if( .not. onflag ) then  ! invalid command word
        call incrack_errmsg( 28 )
        return
      end if
c
c              type <integer>
c
      if( endcrd(dum) ) then
          call incrack_errmsg( 50 )
          return
      end if        
c
      if( matchs_exact('type') ) then
        if( .not. integr( regular_type ) ) then
          call incrack_errmsg( 30 )
          return
        end if
      else
         call incrack_errmsg( 29 )
         return
      end if
c
      bad_type = regular_type <=0 .or. regular_type > 3
      if( bad_type ) then
         call incrack_errmsg( 31 )
      end if
c
c              length_scale or L <number>
c
      if( endcrd(dum) ) then
          call incrack_errmsg( 50 )
          return
      end if        
c
      ok = .false.
      if( matchs_exact('L') ) ok = .true.
      if( matchs('length_scale',6) ) ok = .true.
      if( .not. ok ) then 
         call incrack_errmsg( 32 )
         return
      end if
c
      if( .not. numd( regular_length ) ) then
         call incrack_errmsg( 33 )
         return
      end if
c
      if( regular_length <= zero ) then
         call incrack_errmsg( 34 )
         return
      end if 
c
      if( regular_type == 1 ) then
c
c              up_max <number>
c
        if( endcrd(dum) ) then
            call incrack_errmsg( 50 )
            return
        end if        
c
        if( .not. matchs_exact('up_max') ) then
           call incrack_errmsg( 35 )
           return
        end if
c
        if( .not. numd( regular_up_max ) ) then
           call incrack_errmsg( 36 )
           return
        end if
c
        if( regular_up_max <= zero ) then
           call incrack_errmsg( 37 )
           return
        end if 
c
      end if !  on regular_type
c
c              input data for the type
c
      if( endcrd(dum) ) then
          call incrack_errmsg( 50 )
          return
      end if        
c
      select case( regular_type )
c
      case( 1 )  ! piecewise linear. get #pts and data pairs
         if( .not. matchs_exact('num_points') ) then 
           call incrack_errmsg( 38 )
           return
         end if 
         if( .not. integr( regular_npoints ) ) then 
           call incrack_errmsg( 40 )
           return
         end if 
         ok = regular_npoints > 0 .and. regular_npoints <= 10
         if( .not. ok ) then 
           call incrack_errmsg( 39 )
           return
         end if
c
         do i = 1, regular_npoints
            x = zero
            y = zero
            if( matchs_exact(',') ) call splunj
            if( endcrd( dum) ) then 
              call incrack_errmsg( 50 )
              return
            end if
            if( matchs_exact(',') ) call splunj
            if( .not. numd( x ) ) then
              call incrack_errmsg( 41)  ! no up value
              return
            end if
            if( matchs_exact(',') ) call splunj
            if( .not. numd( y ) ) then
              call incrack_errmsg( 42 )  ! no d value
              return
            end if
            regular_points(i,1) = x
            regular_points(i,2) = y
         end do   ! over i. need to check all values
c
         bad = .false.
         if( regular_points(1,1) /= zero ) then
           call incrack_errmsg( 43 )  ! bad 1st up value
           bad = .true.
         end if
         if( regular_points(regular_npoints,2) /= one ) then
           call incrack_errmsg( 44 )  ! bad last up value
           bad = .true.
         end if
         if( regular_points(regular_npoints,1) /= one ) then
           call incrack_errmsg( 45 )  ! bad last d value
           bad = .true. 
         end if
c 
         do i = 2, regular_npoints-1  ! check values each pair
           x = regular_points(i,1)
           y = regular_points(i,2)
           ok = x > zero .and. x < one   .and.
     &          y > zero .and. y < one  .and.
     &          x > regular_points(i-1,1)
           if( .not. ok ) write(out,9040) i, x, y
           if( .not. ok ) bad = .true.
         end do
         if( bad )then
           call incrack_errmsg( 46)
           return
         end if
c
      case( 2 )  ! G_f. d = G/G_f
         if( .not. matchs_exact('Gf') ) then 
           call incrack_errmsg( 56 )
           return
         end if 
         if( .not. numd( regular_Gf ) ) then 
           call incrack_errmsg( 57 )
           return
         end if
         if( regular_Gf <= zero ) then
           call incrack_errmsg( 58 )
           return
         end if
         if( .not. matchs_exact('m') ) then 
           call incrack_errmsg( 60 )
           return
         end if 
         if( .not. numd( regular_m_power ) ) then 
           call incrack_errmsg( 61 )
           call scan_flushline
           return
         end if
         if( regular_m_power <= zero ) then
           call incrack_errmsg( 62 )
           return
         end if

      case( 3 )  ! G_f w/ d = exponential form of G/G_f
         if( .not. matchs_exact('Gf') ) then 
           call incrack_errmsg( 56 )
           return
         end if 
         if( .not. numd( regular_Gf ) ) then 
           call incrack_errmsg( 57 )
           return
         end if
         if( regular_Gf <= zero ) then
           call incrack_errmsg( 58 )
           return
         end if
         if( .not. matchs_exact('m') ) then 
           call incrack_errmsg( 60 )
           return
         end if 
         if( .not. numd( regular_m_power ) ) then 
           call incrack_errmsg( 61 )
           call scan_flushline
           return
         end if
         if( regular_m_power <= zero ) then
           call incrack_errmsg( 62 )
           return
         end if
         if( .not. matchs_exact('alpha') ) then 
           call incrack_errmsg( 63 )
           return
         end if 
         if( .not. numd( regular_alpha ) ) then 
           call incrack_errmsg( 64 )
           call scan_flushline
           return
         end if
c
      case default
         call incrack_errmsg( 35 )
         call die_abort
      end select 
c
c              input completed. no errors
c
      use_mesh_regularization = .true.
      use_estiff_at_death = .true.
c      
      if( endcrd(dum) ) return
      if( matchs_exact('display') ) then
         select case( regular_type )
         case( 1 ) 
           write(out,9000) regular_type, regular_length, regular_up_max
           write(out,9010) regular_npoints
           do i = 1, regular_npoints
             write(out,9020) i, regular_points(i,1),
     &                          regular_points(i,2)
           end do
         case( 2 ) 
           write(out,9001) regular_type, regular_length,
     &                     regular_GF, regular_m_power
         case( 3 ) 
           write(out,9002) regular_type, regular_length,
     &                     regular_GF, regular_m_power, regular_alpha
         end select
         write(out,9030) 
      end if
c
      call scan_flushline ! just in case anything left
      return
c
 9000 format(/,5x,'... mesh regularization parameter values: ',/,
     & 10x,'regular type:        ',i2, /,
     & 10x,'length-scale:      ',f10.5, /,
     & 10x,'up-max:            ',f10.5 )
 9001 format(/,5x,'... mesh regularization parameter values: ',/,
     & 10x,'regular type:        ',i2, /,
     & 10x,'length-scale:      ',f10.5, /,
     & 10x,'Gf:                ',f10.5, /,
     & 10x,'m exponent:        ',f10.5 )
 9002 format(/,5x,'... mesh regularization parameter values: ',/,
     & 10x,'regular type:        ',i2, /,
     & 10x,'length-scale:      ',f10.5, /,
     & 10x,'Gf:                ',f10.5, /,
     & 10x,'m exponent:        ',f10.5, /,
     & 10x,'alpha:             ',f10.5 )
 9010 format(/,10x,'... number of curve points, values (x,y):',i3)
 9020 format(15x,i2,2f10.5)
 9030 format(/,5x,'... end of display',/)
 9040 format(/1x,'>>>>> error: bad data pair values. ',i2,2f10.5)
! 9100 format(10x,'d0:                ',f10.5)
! 9110 format(10x,'alpha:             ',f10.5)
c
      end subroutine incrack_regularization
c 
c
c          ------------------------------------------                           
c          | contains: incrack_smcs_parms           |                           
c          ------------------------------------------                           
c                                                                                           
      subroutine incrack_smcs_parms
      implicit none
c
c
c        smcs type 1 gamma, alpha, beta, kappa, triaxiality_cutoff
c
c        smcs type 2 gamma, beta_1, beta_2, A_plus, A_minus, kappa,
c                    triaxiality_cutoff
c
c        smcs type 3 gamma, alpha_1, alpha_2, beta_1, beta_2,
c                    triaxiality_cutoff
c
c
      logical :: ok
c
      if( .not. integr(smcs_type) ) then
       call incrack_errmsg( 9 )
       return
      end if
c
      ok = smcs_type .ge. 1  .and.  smcs_type .le. 5
      if( .not. ok ) then
       call incrack_errmsg( 10 )
       return
      end if
c
      smcs_growth = .true.
c
      do
       if( endcrd(dum) ) return
       if( matchs_exact('display') ) then
         write(out,9000) smcs_type, smcs_gamma, smcs_alpha,
     &                   smcs_beta, smcs_kappa,
     &                   smcs_beta_1,
     &                   smcs_beta_2, smcs_a_plus, smcs_a_minus,
     &                   smcs_alpha_1, smcs_alpha_2,
     &                   smcs_cutoff_triaxiality,
     &                   max_eps_critical, smcs_type_4_A, smcs_type_4_n,
     &                   smcs_type_4_c1, smcs_type_4_c2, smcs_type_4_c3,
     &                   smcs_type_5_power, smcs_type_5_tp_critical 
         cycle
       end if

       if( matchs_exact('gamma') ) then
         if( numd(smcs_gamma) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('alpha') ) then
         if( numd(smcs_alpha) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('beta') ) then
         if( numd(smcs_beta) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('beta_1') ) then
         if( numd(smcs_beta_1) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('beta_2') ) then
         if( numd(smcs_beta_2) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('a_plus') ) then
         if( numd(smcs_a_plus) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('a_minus') ) then
         if( numd(smcs_a_minus) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('kappa') ) then
         if( numd(smcs_kappa) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('alpha_1') ) then
         if( numd(smcs_alpha_1) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('alpha_2') ) then
         if( numd(smcs_alpha_2) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs('cutoff_triaxiality',6) ) then
         if( numd(smcs_cutoff_triaxiality) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs('triaxiality_cutoff',6) ) then
         if( numd(smcs_cutoff_triaxiality) ) cycle
         call incrack_errmsg( 4 )
         return
       end if    
       if( matchs('max_eps_critical',12) ) then
         if( numd(max_eps_critical) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('A') ) then
         if( numd(smcs_type_4_A) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('n') ) then
         if( numd(smcs_type_4_n) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('c_1') ) then
         if( numd(smcs_type_4_c1) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('c_2') ) then
         if( numd(smcs_type_4_c2) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('c_3') ) then
         if( numd(smcs_type_4_c3) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('power') ) then
         if( numd(smcs_type_5_power) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       if( matchs_exact('tp_critical') ) then
         if( numd(smcs_type_5_tp_critical) ) cycle
         call incrack_errmsg( 4 )
         return
       end if
       call incrack_errmsg( 11 )
       return
      end do
c
      return
c
 9000 format(/,5x,'... smcs parameter values: ',/,
     & 10x,'smcs_type :          ',i2, /,
     & 10x,'gamma :           ',f10.4, /,
     & 10x,'alpha:            ',f10.4, /,
     & 10x,'beta:             ',f10.4, /,
     & 10x,'kappa:            ',f10.4, /,
     & 10x,'beta_1 :          ',f10.4, /,
     & 10x,'beta_2 :          ',f10.4, /,
     & 10x,'A+ :              ',f10.4, /,
     & 10x,'A- :              ',f10.4, /,
     & 10x,'alpha_1 :         ',f10.4, /,
     & 10x,'alpha_2 :         ',f10.4, /,
     & 10x,'cutoff_triaxiaiality:',f10.4, /,
     & 10x,'max_eps_critical :',2x,f9.5, /,
     & 10x,'type_4_A :        ',f10.4, /,
     & 10x,'type_4_n :        ',f10.4, /,
     & 10x,'type_4_c1 :       ',f10.4, /,
     & 10x,'type_4_c2 :       ',f10.4, /,
     & 10x,'type_4_c3 :       ',f10.4,/
     & 10x,'type_5_power:     ',f10.4, /,
     & 10x,'type_5_critical:  ',f10.4,/ )
c
      end subroutine incrack_smcs_parms
c
c          ------------------------------------------                           
c          | contains: incrack_smcs_states          |                           
c          ------------------------------------------                           
c                                                                                           
      subroutine incrack_smcs_states
      use allocated_integer_list
      implicit none
c
c        smcs states on | off text | stream steps <step list> 
c
      integer, allocatable :: temp_list(:)
      integer :: errnum, lenlst, list_size
c
      smcs_states  = .false.
      smcs_stream  = .false.
      smcs_text    = .false.
      if( allocated( smcs_states_intlst ) )
     &      deallocate( smcs_states_intlst )
c
      if( matchs('files',4) ) call splunj
      if( matchs_exact('on') ) then
         smcs_states = .true.
      elseif( matchs_exact('off') ) then 
         smcs_states = .false.
      else
         call incrack_errmsg( 15 )
         smcs_states = .false.
         return
      end if
      if( .not. smcs_states ) return
c
      if( matchs_exact('text') ) then
         smcs_stream = .false.
         smcs_text   = .true.
      elseif( matchs_exact('stream') ) then 
         smcs_stream = .true.
         smcs_text   = .false.
      else
         call incrack_errmsg( 16 )
         smcs_states = .false.
         return
      end if
c
      if( matchs('steps',4) ) call splunj
      allocate( temp_list(10) ) 
      temp_list = 0
      call scan                                                                 
      call trlist_allocated( temp_list, list_size, 0,   ! all not allowed                         
     &             lenlst, errnum )     
      if ( errnum == 2 ) then   ! parse rules failed                                                
         call incrack_errmsg( 12 )
         smcs_states  = .false.
         return
      else if ( errnum == 3 ) then  ! list overflow                                            
         call incrack_errmsg( 13 )
         smcs_states  = .false.
         return
      else if ( errnum == 4 ) then  ! no list                                             
         call incrack_errmsg( 14 )
         smcs_states  = .false.
         return
      endif
c
      allocate( smcs_states_intlst( lenlst ) )
      smcs_states_intlst(1:) = temp_list(1:lenlst)
c
c        not yet supported on MPI
c     
      if( use_mpi ) call incrack_errmsg( 17 )
c
      return
c
      end subroutine incrack_smcs_states
c
c          ------------------------------------------                           
c          | contains: incrack_gt_list_file         |                           
c          ------------------------------------------                           
c                                                                                           
      subroutine incrack_gt_list_file
c
      use damage_data, only : gt_list_file_flag, gt_list_file_name
      implicit none
c
      integer :: idummy, nchars, device
      logical, external :: label, string
      logical :: fexists
      character(len=80) :: temp_name
c
c        gurson deleted element file on | off text | name <string> 
c
      gt_list_file_flag = .false.
      gt_list_file_name = " "
c
      if( matchs('element',4) ) call splunj
      if( .not. matchs('file',4) ) then
        call incrack_errmsg( 19 )
        call scan_flushline
        return
      end if
c
      if( matchs_exact('on') ) then
         gt_list_file_flag = .true.
      elseif( matchs_exact('off') ) then 
         gt_list_file_flag = .false.
         call scan_flushline
         return
      else
         call incrack_errmsg( 15 )
         gt_list_file_flag = .false.
         call scan_flushline
         return
      end if
c
      if( matchs('name',4) ) call splunj
c
      if( label( idummy) ) then
        call entits( temp_name, nchars )
      elseif( string( idummy ) ) then
        call entits( temp_name, nchars )
      else
        call incrack_errmsg( 18 ) 
        call scan_flushline
        gt_list_file_flag = .false.
        return
      end if
c
      gt_list_file_name(1:) = temp_name(1:nchars)
      inquire( file=gt_list_file_name, exist=fexists )
      if( fexists ) then
         call incrack_errmsg( 20 ) 
         open(newunit=device, file=gt_list_file_name, status='old')
         close(device, status='delete')
      end if
c
      end subroutine incrack_gt_list_file
c
 
c          -------------------------------------------------                           
c          | contains: incrack_smcs_deleted_list_file       |                           
c          --------------------------------------------------                           
c                                                                                           
      subroutine incrack_smcs_deleted_list_file
c
      use damage_data, only : smcs_deleted_list_file_flag, 
     &                        smcs_deleted_list_file_name
      implicit none
c
      integer :: idummy, nchars, device
      logical, external :: label, string
      logical :: fexists
      character(len=80) :: temp_name
c
c        smcs deleted element file on | off text | name <string> 
c
      smcs_deleted_list_file_flag = .false.
      smcs_deleted_list_file_name = " "
c
      if( matchs('element',4) ) call splunj
      if( .not. matchs('file',4) ) then
        call incrack_errmsg( 19 )
        call scan_flushline
        return
      end if
c
      if( matchs_exact('on') ) then
         smcs_deleted_list_file_flag = .true.
      elseif( matchs_exact('off') ) then 
         smcs_deleted_list_file_flag = .false.
         call scan_flushline
         return
      else
         call incrack_errmsg( 15 )
         smcs_deleted_list_file_flag = .false.
         call scan_flushline
         return
      end if
c
      if( matchs_exact('name') ) call splunj
c
      temp_name = " " 
      if( label(idummy) ) then
        call entits( temp_name, nchars )
      elseif( string( idummy ) ) then
        call entits( temp_name, nchars )
      else
        call incrack_errmsg( 18 ) 
        call scan_flushline
        smcs_deleted_list_file_flag = .false.
        return
      end if
c
      smcs_deleted_list_file_name(1:) = " "
      smcs_deleted_list_file_name(1:) = temp_name(1:nchars)
      inquire( file=smcs_deleted_list_file_name, exist=fexists )
      if( fexists ) then
         call incrack_errmsg( 20 ) 
         open(newunit=device, file=smcs_deleted_list_file_name, 
     &        status='old')
         close(device, status='delete')
      end if
c
      end subroutine incrack_smcs_deleted_list_file
 
c          --------------------------------------------------                           
c          | contains: incrack_smcs_removed_list_file       |                           
c          ---------------------------------------------------                           
c                                                                                           
      subroutine incrack_smcs_removed_list_file
c
      use damage_data, only : smcs_removed_list_file_flag, 
     &                        smcs_removed_list_file_name
      implicit none
c
      integer :: idummy, nchars, device
      logical, external :: label, string
      logical :: fexists
      character(len=80) :: temp_name
c
c        smcs removed element file on | off text | name <string> 
c
      smcs_removed_list_file_flag = .false.
      smcs_removed_list_file_name = " "
c
      if( matchs('element',4) ) call splunj
      if( .not. matchs('file',4) ) then
        call incrack_errmsg( 19 )
        call scan_flushline
        return
      end if
c
      if( matchs_exact('on') ) then
         smcs_removed_list_file_flag = .true.
      elseif( matchs_exact('off') ) then 
         smcs_removed_list_file_flag = .false.
         call scan_flushline
         return
      else
         call incrack_errmsg( 15 )
         smcs_removed_list_file_flag = .false.
         call scan_flushline
         return
      end if
c
      if( matchs('name',4) ) call splunj
c
      if( label( idummy) ) then
        call entits( temp_name, nchars )
      elseif( string( idummy ) ) then
        call entits( temp_name, nchars )
      else
        call incrack_errmsg( 18 ) 
        call scan_flushline
        smcs_removed_list_file_flag = .false.
        return
      end if
c
      smcs_removed_list_file_name(1:) = temp_name(1:nchars)
      inquire( file=smcs_removed_list_file_name, exist=fexists )
      if( fexists ) then
         call incrack_errmsg( 20 ) 
         open(newunit=device, file=smcs_removed_list_file_name, 
     &        status='old')
         close(device, status='delete')
      end if
      inquire( file=smcs_removed_list_file_name, exist=fexists )
      if( fexists ) then
         open(newunit=device, file=smcs_removed_list_file_name,
     &         status='old')
         close(device, status='delete')
      end if
c
      end subroutine incrack_smcs_removed_list_file
c
      subroutine incrack_stop_elist 
      implicit none
      integer :: vector_length, errnum  
      logical, parameter :: ldebug = .false.                                                        
c                          
c              stop (if) killed (elements) < list >  
c                 
      if( matchs('elements',3) ) call splunj
      stop_killed_elist_length = 0
      if( allocated( deleted_elist_to_stop ) ) 
     &    deallocate( deleted_elist_to_stop )     
      if( matchs_exact('off') ) return
      if( .not. matchs_exact('on') ) then
         call incrack_errmsg( 22 )
         return
      end if     
      call scan
      call trlist_allocated( deleted_elist_to_stop, vector_length, 
     &                       0, stop_killed_elist_length, errnum )
c
c              1 = no error. list found
c              2 = parse rules failed in list
c              3 = list overflowed its maximum length of mxlsz
c              4 = no list found 
c
       if( errnum .eq. 4 ) then 
          write(out,9000) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
      if( errnum .eq. 2 ) then
          write(out,9010) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
      if( errnum .eq. 3 ) then
          write(out,9020) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
c
      if( errnum .ne. 1 ) then  ! something bad wrong
         write(out,9050)
         call die_gracefully
      end if
c     
      return
 9000 format(/1x,'>>>>> error: no element list found',
     &   /,14x,'read new line ...',/)
 9010 format(/1x,'>>>>> error: integer list rules failed....',
     &   ' read new line ...',/)
 9020 format(/1x,'>>>>> error: too many elements in list....',
     &   ' read new line ...',/)
 9050 format(/1x,'>>>>> FATAL ERROR: delete_elements inconsistency',
     &   /,14x,'job terminated now ...',/)
c
      end subroutine incrack_stop_elist
c
      subroutine incrack_stop_nlist 
      implicit none
      integer :: vector_length, errnum  
      logical, parameter :: ldebug = .false.                                                        
c                          
c              stop (if) nodes (released) < list >  
c                 
      if( matchs('released',3) ) call splunj
      stop_released_nlist_length = 0
      if( allocated( released_nlist_to_stop ) ) 
     &    deallocate( released_nlist_to_stop )     
      if( matchs_exact('off') ) return
      if( .not. matchs_exact('on') ) then
         call incrack_errmsg( 22 )
         return
      end if     
      call scan
      call trlist_allocated( released_nlist_to_stop, vector_length, 
     &                       0, stop_released_nlist_length, errnum )
c
c              1 = no error. list found
c              2 = parse rules failed in list
c              3 = list overflowed its maximum length of mxlsz
c              4 = no list found 
c
       if( errnum .eq. 4 ) then 
          write(out,9000) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
      if( errnum .eq. 2 ) then
          write(out,9010) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
      if( errnum .eq. 3 ) then
          write(out,9020) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
c
      if( errnum .ne. 1 ) then  ! something bad wrong
         write(out,9050)
         call die_gracefully
      end if
c     
      return
 9000 format(/1x,'>>>>> error: no node list found',
     &   /,14x,'read new line ...',/)
 9010 format(/1x,'>>>>> error: integer list rules failed....',
     &   ' read new line ...',/)
 9020 format(/1x,'>>>>> error: too many nodes in list....',
     &   ' read new line ...',/)
 9050 format(/1x,'>>>>> FATAL ERROR: released nodes inconsistency',
     &   /,14x,'job terminated now ...',/)
c
      end subroutine incrack_stop_nlist
c
      end subroutine incrack                                                                 
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine incrack_dam_init             *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 05/16/2019 rhd             *          
c     *                                                              *          
c     *     this subroutine allocates dam_ifv and dam_state.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_dam_init( debug )                                              
      use global_data ! old common.main
c                                                                               
      use damage_data                                                           
c                                                                               
      implicit none
c
      logical :: debug                                                    
c                                                                               
c             locally allocated arrays                                          
c                       
      integer :: elem, dum                                                        
      character(len=1) :: dums                                                  
      real :: dumr                                                                 
      double precision  :: dumd                                                                    
c                                                                               
      if ( debug ) write (out,*) '>>>> in incrack_dam_init'                             
c                                                                               
c                if the damage data structures have already been allocated,     
c                then leave                                                     
c                                         
      if ( num_kill_elem .ne. 0 ) then                                          
       call errmsg( 215, dum, 'dam_ifv', dumr, dumd )                           
         return                                                                 
      end if                                                                    
c                                                                               
c                must find out how many elements can be killed and              
c                fill dam_ptr                                                   
c                                                                               
      num_kill_elem = 0                                                         
      do elem = 1, noelem                                                       
        if ( iand( iprops(30,elem),2 ) .ne. 0 ) then                            
            num_kill_elem = num_kill_elem + 1                                   
            dam_ptr(elem) = num_kill_elem                                       
        end if                                                                  
      end do                                                                    
c                                                                               
c                error if no elements are killable.                             
c                                                                               
      if ( num_kill_elem .eq. 0 ) then                                          
         call errmsg( 208, dum, dums, dumr, dumd )                              
         crack_growth_type = 0                                                  
         growth_by_kill    = .false.                                            
         return                                                                 
      end if                                                                    
c                                                                               
c               now allocate needed varables                                    
c                                                                               
      call allocate_damage( 1 )                                                 
c                                                                               
c               set global counter for total elements killed at                 
c               point in analysis                                               
c                                                                               
      num_elements_killed = 0                                                   
c                                                                               
      if ( debug ) then                                                         
         call dam_debug                                                         
         write (out,*) '<<<< leaving incrack_dam_init'                                  
      end if                                                                    
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_dam_init_release             *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 09/05/95                   *          
c     *                                                              *          
c     *     this subroutine finds the nodes on the crack plane       *          
c     *     and allocates the data structures needed for the         *          
c     *     node_release crack growth algorithm.  Then it calls      *          
c     *     a routine to do the second half of the initialization    *          
c     *     process.                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_dam_init_release 
c                                              
      use global_data ! old common.main
      use main_data, only : crdmap                                              
      use damage_data
      use constants  
                                                               
c                                                                               
      implicit none                                                    
c                                                                               
c             local declarations
c           
      integer :: axis, node, start_of_last_row, i, j, dum                                                                   
      character(len=1) :: dums                                                  
      real :: dumr                                                                 
      double precision :: plane_tol, max_dist, dumd 
      integer, allocatable :: tmp_nodes(:)                                                  
      logical, parameter :: debug = .false.                                                             
c                                                                               
      if( debug ) write (out,*) '>>>> in incrack_dam_init_release'                      
c                                                                               
c                return if we have already initialized                          
c                                                                               
      if( num_crack_plane_nodes .ne. 0 ) return                                 
c                                                                               
      write(out,9000)                                                           
      write(out,9010)                                                           
                                                                                
c                                                                               
c                allocate a temporary array of flags for each node              
c                to tell if the node is on the crack plane or not.              
c                                                                               
      allocate( tmp_nodes(nonode) )                                             
c                                                                               
c                now loop over all of the nodes.  find all nodes                
c                that are within tolerance of the given crack                   
c                plane. exclude nodes that are on the plane but                 
c                which are also unconstrained in the crack plane                
c                normal direction.                                              
c                                                                               
      axis                  = crk_pln_normal_idx - 1                            
      num_crack_plane_nodes = 0                                                 
      plane_tol             = max( crkpln_srch_tol, crkpln_srch_tol *           
     &                             abs(crack_plane_coord) )                     
      max_dist = zero                                                           
c                                                                               
      do node = 1, nonode                                                       
         if( abs(c(crdmap(node)+axis) - crack_plane_coord)                      
     &        .le. plane_tol ) then                                             
           num_crack_plane_nodes            = num_crack_plane_nodes + 1         
           tmp_nodes(num_crack_plane_nodes) = node                              
           max_dist                         = max( max_dist,                    
     &                 abs(c(crdmap(node)+axis) - crack_plane_coord) )          
         end if                                                                 
      end do                                                                    
c                                                                               
c                if no crack plane nodes are found, do a fatal error            
c                message and skip allocation phase.                             
c                                                                               
      if( num_crack_plane_nodes .eq. 0 ) then                                   
         call errmsg ( 247, dum, dums, dumr, dumd )                             
         go to 9999                                                             
      end if                                                                    
c                                                                               
c                if list_crkpln_nodes is true, then list all the                
c                crack plane nodes                                              
c                                                                               
      write (out,9020) num_crack_plane_nodes, max_dist                          
      if( list_crkpln_nodes ) then                                              
         write(out,9030)                                                        
         start_of_last_row = 1                                                  
         if( num_crack_plane_nodes .gt. 5 ) then                                
            do i = 1, num_crack_plane_nodes/5                                   
               write(out,9040) (tmp_nodes(j),j=(i-1)*5+1,i*5)                   
            end do                                                              
            start_of_last_row = (int(num_crack_plane_nodes/5))*5 + 1            
         end if                                                                 
         if( start_of_last_row .le. num_crack_plane_nodes ) then                
            write(out,9040) (tmp_nodes(j),j=start_of_last_row,                  
     &                       num_crack_plane_nodes)                             
         end if                                                                 
      end if                                                                    
c                                                                               
c                                                                               
c                allocate the needed variables for the                          
c                node_release algorithm:                                        
c                           crack_plane_nodes -- nodes on crack plane           
c                           num_neighbors -- # of neighboring nodes             
c                                      for each crack plane node                
c                           neighbor_nodes -- the neighboring nodes             
c                           crack_front_nodes -- linked list structure          
c                                      for the crack front nodes                
c                           node_release_frac -- stores the current             
c                                      release fraction for each                
c                                      released node (only for                  
c                                      traction separation)                     
c                 then initialize those data structures (e.g find               
c                 the neighbor nodes, find the initial crack front, etc.)       
c                                                                               
      call allocate_damage ( 6 )                                                
      call incrack_dam_init_release2 ( tmp_nodes, debug )                               
c                                                                               
 9999 continue                                                                  
      deallocate( tmp_nodes )                                                   
c                                                                               
      return                                                                    
 9000 format (/1x,'>>>> Initializing crack growth by node release')             
 9010 format (/1x,'       * find crack plane nodes...')                         
 9020 format (/1x,'       * found ',i5,' crack plane nodes.  Maximum',          
     &        /1x,'          distance from crack plane to crack plane',         
     &        /1x,'          node: ',e13.6)                                     
 9030 format (/1x,'         list of crack plane nodes:')                        
 9040 format (13x,5(2x,i7))   
c                                                  
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *            subroutine incrack_dam_init_release2              *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 05/16/2019 rhd             *          
c     *                                                              *          
c     *     this subroutine sets up the permanent array to hold all  *          
c     *     the crack plane nodes.  It also sets up and fills the    *          
c     *     linked list to hold the crack front nodes, and finds the *          
c     *     neighbor nodes for all of the crack plane nodes.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_dam_init_release2( tmp_nodes, debug )                          
      use global_data ! old common.main
c                                                                               
      use node_release_data                                                     
      use main_data, only : incmap, incid, inverse_incidences                   
      use damage_data  
      use constants                                                         
c                                                                               
      implicit none                                                    
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: tmp_nodes(*)                                                    
      logical :: debug                                                             
c                                                                               
c                     local declarations                                        
c             
      integer :: node, entry, elem, i, elem_node, dum
      integer, external :: incrack_master
      logical :: crack_node, first_node                                            
      double precision :: dumd     
      real :: dumr                                                                 
      character(len=1) :: dums                                                  
c                                                                               
c                                                                               
      if ( debug ) write (out,*) '>>>> in incrack_dam_init_release'                     
c                                                                               
c                fill the crack_plane_nodes array                               
c                                                                               
      do node = 1, num_crack_plane_nodes                                        
         crack_plane_nodes(node) = tmp_nodes(node)                              
      end do                                                                    
c                                                                               
c                fill the inv_crkpln_nodes array.  This array                   
c                indicates if a node is not on the crack plane                  
c                (has a 0 as the entry), or, if the node IS on                  
c                the crack plane, gives the corresponding entry                 
c                in the node release data structures.                           
c                                                                               
      inv_crkpln_nodes(1:nonode) = 0                                            
c                                                                               
      do entry = 1, num_crack_plane_nodes                                       
         inv_crkpln_nodes(crack_plane_nodes(entry)) = entry                     
      end do                                                                    
c                                                                               
c                The crack face may move in the wrong direction                 
c                in certain circumstances, moving through the crack             
c                crack plane.  This can occur due to changes in the             
c                structure loading or through release anolmolies.  In           
c                order to detect this behaviour, we need to know whether        
c                the elements are on the positive or negative side              
c                of the crack plane.  Find out which is the case and            
c                store it in crack_plane_sign (= 1 if elements are on           
c                the positive side, = -1 if on the negative side).  Check       
c                this through an element connected to the first crack plane     
c                node.                                                          
c                                                                               
      node = crack_plane_nodes(1)                                               
      elem = inverse_incidences(node)%element_list(1)                           
      do i = 1, 8                                                               
         elem_node = incid(incmap(elem)+i-1)                                    
         if( inv_crkpln_nodes(elem_node) .ne. 0 ) cycle                         
         if( c(dstmap(elem_node)+crk_pln_normal_idx - 1)                        
     &        .gt. crack_plane_coord ) then                                     
              crack_plane_sign = one                                            
         else                                                                   
              crack_plane_sign = - one                                          
         end if                                                                 
         exit                                                                   
      end do                                                                    
      if ( debug ) write (out,'("CRACK PLANE SIGN:",f5.2)')                     
     &     crack_plane_sign                                                     
c                                                                               
c                initialize killed node information structures (f90 syntax)     
c                                                                               
      do i = 1, num_crack_plane_nodes                                           
         crkpln_nodes_react(i) = zero                                           
         crkpln_nodes_state(i) = 0                                              
      end do                                                                    
c                                                                               
c                To describe the current crack front and to provide a           
c                data structure the is easy to change as nodes are              
c                deleted and added to the crack front, use a singly             
c                linked list to hold the crack front nodes.  The                
c                list is stored in crack_front_nodes(nodes,2), where the        
c                crack_front_nodes(node,1) entry stores the global node         
c                number for crack front node "node", and the                    
c                crack_front_nodes(node,2) entry stores the pointer             
c                to the next crack front node in the list.  A value             
c                of -1 as the pointer means the end of the list.  Array         
c                entries that are not in the node list are kept in a            
c                garbage list, where the node number entry is -1 and the        
c                pointer entry stores the pointer to the next entry in the      
c                garbage list.  crack_front_start and crack_front_end are       
c                head and tail pointers to the list of crack front nodes,       
c                and crkfrnt_garbage_start and crkfrnt_garbage_end are          
c                the head and tail pointers to the garbage list.                
c                                                                               
c                To initialize the linked list, make the entire                 
c                list a garbage list, with crkfrnt_garbage_start                
c                pointing to the top of the array, all the entries              
c                pointing sequentially down the array, and                      
c                crkfrnt_garbage_end the array, pointing at the last            
c                array entry.  After this, adding and deleting from the         
c                list will just be a changing of the entry status               
c                from garbage to data and vice versa.                           
c                                                                               
      crkfrnt_garbage_start = 1                                                 
      crkfrnt_garbage_end = num_crack_plane_nodes                               
      do i = 1, num_crack_plane_nodes                                           
         crack_front_nodes (i,1) = -1                                           
         crack_front_nodes (i,2) = i+1                                          
      end do                                                                    
      crack_front_nodes (num_crack_plane_nodes,2) = -1                          
c                                                                               
      crack_front_start = -1                                                    
      crack_front_end = -1                                                      
c                                                                               
c                 find all of the neighbors for each crack plane node.          
c                 A neighbor is defined as a node on the crack plane            
c                 which shares an element edge with the node we are checking.   
c                 If one or more of the neighbors are unconstrained in          
c                 the direction normal to the crack plane, then the node        
c                 we are checking is an initial crack front node, and we        
c                 need to add it to the crack_front_nodes linked list.          
c                                                                               
c                 If we are using const_front growth, then the only             
c                 nodes put into the crack_front_nodes are the                  
c                 master nodes, one per crack front.                            
c                                                                               
      write (out,9000)                                                          
      first_node = .true.                                                       
      do i = 1, num_crack_plane_nodes                                           
         call find_neighbors ( crack_plane_nodes(i), num_neighbors(i),          
     &        neighbor_nodes(1,i), crack_node, inv_crkpln_nodes,                
     &        crk_pln_normal_idx )                                              
         if( .not. crack_node) cycle                                            
         if( const_front ) then                                                 
            if( incrack_master(crack_plane_nodes(i)) .eq. 0 ) cycle                     
         end if                                                                 
         call add_to_list( crack_plane_nodes(i) )                               
c                                                                               
c                     write crack front node if printing is requested           
c                                                                               
         if (.not. list_crkfrnt_nodes) cycle                                    
         if ( first_node ) then                                                 
            write (out,9030)                                                    
            first_node = .false.                                                
         endif                                                                  
         write (out,9060) crack_plane_nodes(i)                                  
      end do                                                                    
      if ( list_crkfrnt_nodes ) write(out,*)                                    
c                                                                               
c                 check the crack_front_nodes pointers to make sure that:       
c                     1) there is a crack on the crack plane, i.e.              
c                        a crack front was found, and                           
c                     2) the crack does not encompass the entire plane,         
c                        i.e. there is at least one constrained node            
c                 If either of these is false, print an error and               
c                 inhibit future computation of load steps.                     
c                                                                               
      if ( crack_front_start .eq. -1 ) then                                     
         call errmsg( 250, dum, dums, dumr, dumd )                              
      end if                                                                    
c                                                                               
      if ( debug ) then                                                         
         call dam_debug                                                         
         write (out,*) '<<<< leaving incrack_dam_init_release2'                         
      end if                                                                    
      return                                                                    
 9000 format (/1x,'       * find neighboring nodes...')                         
 9030 format (/1x,'         list of crack front nodes:')                        
 9060 format (15x,i7)                                                           
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_find_release_height          *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 5/16/2019 rhd              *          
c     *                                                              *          
c     *     this subroutine fills the dam_print_list array           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_find_release_height                                            
      use global_data ! old common.main
c                                                                               
      use damage_data
      use constants                                                           
c                                                                               
      implicit none                                                    
c                                                                               
c                     local declarations                                        
c      
      integer :: dum                                                                    
      double precision :: dumd                                                  
      real :: dumr                                                                 
      character(len=1) :: dums                                                  

c                                                                               
c                 the release height is calculated as follows:                  
c                 consider a line segment with length characteristic_length     
c                 (which is specified in the user input) at an angle to the     
c                 crack plane of release_fraction * critical_angle. One end     
c                 point of the line rests on the crack plane.  The release      
c                 height is the distance between the other endpoint and         
c                 the crack plane.                                              
c                                                                               
      if( no_released_nodes ) then                                              
         call allocate_damage ( 7 )                                             
         if( const_front ) then                                                 
            release_height = max(init_ctoa_dist,ctoa_dist) *                    
     &           tan( critical_angle * release_fraction *                       
     &           (pi/oneeighty) )                                               
         else                                                                   
            release_height = char_length * tan( critical_angle *                
     &           release_fraction * (pi/oneeighty) )                            
       end if                                                                   
         write(out,9010) release_height                                         
      else                                                                      
         call errmsg(256,dum,dums,dumr,dumd)                                    
      end if                                                                    
c                                                                               
      return                                                                    
 9010 format (/1x,'       * computed release height (for traction',             
     &        /1x,'           separation):',e13.6,/)                            
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_print_list_fill              *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 5/16/2019 rhd              *          
c     *                                                              *          
c     *     this subroutine fills the dam_print_list array           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_print_list_fill( scan_order_list,                              
     &         scan_order_list_size, debug )                                    
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_print_list                              
      use damage_data                                                           
c                                                                               
      implicit none
c                                                    
      integer :: scan_order_list_size 
      integer :: scan_order_list(scan_order_list_size)                                              
      logical :: debug   
c
      integer :: iplist, icn, list_entry, elem, i, dum                                                         
c                                                                               
      if ( debug ) write (out,*) '>>>> in incrack_print_list_fill'                      
c                                                                               
c           now fill the order array with the specified elements                
c                                                                               
      iplist     = 1                                                            
      icn        = 0                                                            
      list_entry = 0                                                            
 570  continue                                                                  
      call trxlst( scan_order_list, scan_order_list_size, iplist, icn,          
     &     elem )                                                               
      if ( iand (iprops(30,elem),2).ne.0 ) then                                 
         list_entry = list_entry + 1                                            
         dam_print_list(list_entry) = elem                                      
      end if                                                                    
c                                                                               
      if ( iplist.ne.0 ) goto 570                                               
c                                                                               
      if ( debug ) then                                                         
         write (out,*) '>>>>> dam_print_list:'                                  
       do i = 1, num_print_list                                                 
          write (out,'(2x,i3,1x,i7)') i,dam_print_list(i)                       
       end do                                                                   
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_kill_order_fill              *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 5/16/2019 rhd              *          
c     *                                                              *          
c     *     this subroutine fills the kill_order_list array.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_kill_order_fill( scan_kill_order_list,                         
     &      scan_kill_order_length, debug )                                     
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : kill_order_list                             
      use damage_data                                                           
c                                                                               
      implicit none
c                                                    
      integer :: scan_kill_order_length
      integer :: scan_kill_order_list(scan_kill_order_length)                                         
      logical :: debug                                                             
c   
      integer :: iplist, icn, list_entry, elem, i
c                                                                            
      if ( debug ) write (out,*) '>>>> in incrack_kill_order_fill'                      
c                                                                               
c           now fill the order array                                            
c                                                                               
      iplist     = 1                                                            
      icn        = 0                                                            
      list_entry = 0                                                            
c                                                                               
 670  continue                                                                  
      call trxlst( scan_kill_order_list, scan_kill_order_length,                
     &     iplist, icn, elem )                                                  
      if ( iand (iprops(30,elem),2).ne.0 ) then                                 
         list_entry = list_entry + 1                                            
         kill_order_list(list_entry) = elem                                     
         if ( debug ) write(out,'(" entry:",i3," elem:",i7)')list_entry,        
     &        elem                                                              
      end if                                                                    
c                                                                               
      if ( iplist.ne.0 ) goto 670                                               
c                                                                               
      if ( debug ) then                                                         
         write (out,*) '>>>>> kill_order_list:'                                 
       do i = 1, num_kill_order_list                                            
          write (out,'(2x,i3,1x,i7)') i,kill_order_list(i)                      
       end do                                                                   
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine incrack_master_list_fill              *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 10/15/94                   *          
c     *                                                              *          
c     *     this subroutine fills the master_nodes array.            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_master_list_fill( scan_master_list,                            
     &      scan_master_length, debug )                                         
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_nodes                                
      use main_data,         only : cnstrn_in                                   
      use damage_data 
      use constants                                                          
c                                                                               
      implicit none
c
      integer :: scan_master_length                                                    
      integer :: scan_master_list(scan_master_length)                                             
      logical :: debug  
c
      integer :: iplist, icn, list_entry, node, i, dof
c                                                                               
      if ( debug ) write (out,*) '>>>> in incrack_master_list_fill'                     
c                                                                               
c           now fill the order array                                            
c                                                                               
      iplist     = 1                                                            
      icn        = 0                                                            
      list_entry = 0                                                            
 670  continue                                                                  
      call trxlst( scan_master_list, scan_master_length,                        
     &     iplist, icn, node )                                                  
c                                                                               
      dof = dstmap(node) + crk_pln_normal_idx - 1                               
      if ( cnstrn_in(dof) .ne. d32460 ) then                                    
         list_entry = list_entry + 1                                            
         master_nodes(list_entry) = node                                        
         if ( debug ) write(out,'(" entry:",i3," node:",i7)')list_entry,        
     &        node                                                              
      endif                                                                     
c                                                                               
      if ( iplist.ne.0 ) goto 670                                               
c                                                                               
      if ( debug ) then                                                         
         write (out,*) '>>>>> master_nodes:'                                    
         do i = 1, num_crack_fronts                                             
            write (out,'(2x,i3,1x,i7)') i,master_nodes(i)                       
       end do                                                                   
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      function incrack_master                 *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 5/16/2019 rhd              *          
c     *                                                              *          
c     *                   is a node in master  list                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      function incrack_master( node ) result( iresult )                                              
c                                                                               
      use node_release_data, only : master_nodes                                
      use damage_data, only : const_front, num_crack_fronts                                         
c                                                                               
      implicit none
c
      integer :: node, iresult
c
      integer :: temp, i
c                                                                               
c                                                                               
      temp = 0                                                                  
      if (.not. const_front) goto 9999                                          
      do i=1, num_crack_fronts                                                  
         if (node .eq. master_nodes(i)) then                                    
            temp = i                                                            
            exit                                                                
         endif                                                                  
      enddo                                                                     
c                                                                               
 9999 continue                                                                  
      iresult = temp                                                             
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_init_ctoa_back               *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 5/16/2019 rhd              *          
c     *                                                              *          
c     *     this subroutine initializes the data structures required *          
c     *     to get the CTOA at an arbitrary distance behind the      *          
c     *     crack tip. This only works for const_front growth.       *          
c     *                                                              *          
c     *     Measuring the CTOA back from a master node requires      *          
c     *     the determination of a master line --  a line of nodes   *          
c     *     directly behind the master node, perpendicular to the    *          
c     *     crack front. A master line is created for each master    *          
c     *     node going back until a crack end is reached or          *          
c     *     num_nodes_back nodes, whichever is fewer.                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine incrack_init_ctoa_back                                                 
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_nodes, num_neighbors,                
     &     neighbor_nodes, inv_crkpln_nodes, master_lines                       
      use main_data, only : cnstrn_in                                   
      use damage_data  
      use constants                                                         
c                                                                               
      implicit none
c                            
      integer :: num_lines, i, j, node_data_entry, neighbor, 
     &           neighbor_node, dof, master_node  
      logical, parameter :: debug = .false.                                                
c                                                                               
      if (debug) write (*,*) '>>> in incrack_init_ctoa_back'                            
c                                                                               
c         if we've already been here, then get out                              
c                                                                               
      if ( master_lines_set ) return                                            
c                                                                               
c         loop over each of the master nodes, creating the master line          
c         for each node.                                                        
c                                                                               
      num_lines = 0                                                             
      do i=1, num_crack_fronts                                                  
         master_node = master_nodes(i)                                          
c                                                                               
c         put master node on list, then find its unconstrained                  
c         node and put it on the list.                                          
c                                                                               
         node_data_entry = inv_crkpln_nodes ( master_node )                     
         do neighbor = 1, num_neighbors(node_data_entry)                        
            neighbor_node = neighbor_nodes(neighbor,node_data_entry)            
            dof= dstmap(neighbor_node)+crk_pln_normal_idx-1                     
            if ( cnstrn_in(dof).ne.d32460 ) cycle                               
c                                                                               
c             now we have two nodes on the master line. Going in the            
c             direction from the master node to the other node we found, use    
c             the connectivity data to find the rest of the nodes on the master 
c             line in that direction.                                           
c                                                                               
            num_lines = num_lines + 1                                           
            call incrack_find_master_line( master_node, neighbor_node,                  
     &           num_nodes_back, num_lines)                                     
c                                                                               
         enddo                                                                  
      enddo                                                                     
c                                                                               
      master_lines_set = .true.                                                 
c                                                                               
      if (debug) then                                                           
         write (*,*) '>>> Here are the master_lists...'                         
         do i=1, num_lines                                                      
            write (*,'(6x,i7,6x,20i7)') i, (master_lines(i,j),j=1,              
     &           num_nodes_back + 1)                                            
         enddo                                                                  
         write (*,*) '<<< leaving incrack_init_ctoa_back'                               
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine incrck_find_master_line             *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 10/15/97                   *          
c     *                                                              *          
c     *     This subroutine, given a master node and its adjacent    *          
c     *     unconstrained node, finds the master line extending      *          
c     *     back from these nodes.  This enables the measurement     *          
c     *     of CTOA at a given distance behind the crack front.      *          
c     *     Stop adding nodes to line when the end of the crack      *          
c     *     plane section is reached, or when num_nodes_back is      *          
c     *     reached.                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine incrack_find_master_line( base_node, next_node,                        
     &                      num_nodes_back, num_line )                                            
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_lines       
      use constants                         
c                                                                               
      implicit none
c              
      integer ::  base_node, next_node, num_nodes_back, num_line
c                    
      integer :: idx, new_node, master_line_idx                  
      logical :: reached_end 
      logical, parameter ::  debug = .false.                                               
c                                                                               
      if (debug) write (*,*) '>>> in incrack_find_master_line'                          
c                                                                               
c         we are given two nodes on the master line. Using these two nodes      
c         as a starting point, find num_nodes_back nodes on the master line     
c         behind the nodes.                                                     
c                                                                               
      master_lines(num_line,1) = base_node                                      
      master_lines(num_line,2) = next_node                                      
      master_line_idx = 2                                                       
c                                                                               
      do idx = 3, num_nodes_back + 1                                            
         call incrack_find_master_line_node( base_node, next_node,
     &                                       new_node, reached_end )                                                     
         if (reached_end) exit                                                  
c                                                                               
         master_lines(num_line,idx) = new_node                                  
         base_node = next_node                                                  
         next_node = new_node                                                   
c                                                                               
      end do                                                                     
c                                                                               
      if (debug) write (*,*) '<<< leaving incrack_find_master_line'                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_find_master_line_node        *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 5/16/2019 rhd              *          
c     *                                                              *          
c     *     this subroutine, given two nodes on the master line,     *          
c     *     finds the next node in the master line.  If the node     *          
c     *     found is the last node on the line in the model, then    *          
c     *     send back reached_end = .true.                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine incrack_find_master_line_node ( base_node, 
     &         next_node, new_node, reached_end )                                                     
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : num_neighbors, neighbor_nodes,              
     &        inv_crkpln_nodes                                                  
      use main_data, only : inverse_incidences                                  
      use damage_data                                                           
c                                                                               
      implicit none
c
      integer :: base_node, next_node, new_node
      logical :: reached_end
c
      integer :: num_shared_elems, num_base_elems, num_next_elems,
     &           base_elem_list,  base_elem, next_elem_list, i, j,
     &           next_elem, node_data_entry, neighbor_node, 
     &           num_elems, elem_list, elem                                              
      logical :: found_new_node     
      logical, parameter :: debug = .false.                           
      integer :: shared_elems(2)                                                 
c                                                                               
      if (debug) then                                                           
         write (*,*) '   >>> in incrack_find_master_line_node'                          
         write (*,*) '    base_node:',base_node,' next_node:',next_node         
      endif                                                                     
c                                                                               
c           The two given nodes define an element edge.  Find which             
c           elements have this edge; for normal cases there should be           
c           two, unless the given nodes are on a symmetry plane, where          
c           there will be only one.                                             
c                                                                               
      num_shared_elems = 0                                                      
      num_base_elems = inverse_incidences(base_node)%element_count              
      num_next_elems = inverse_incidences(next_node)%element_count              
c                                                                               
      do base_elem_list = 1, num_base_elems                                     
c                                                                               
         base_elem =                                                            
     &     inverse_incidences(base_node)%element_list(base_elem_list)           
         do next_elem_list = 1, num_next_elems                                  
c                                                                               
            next_elem =                                                         
     &        inverse_incidences(next_node)%element_list(next_elem_list)        
            if ( next_elem .eq. base_elem ) then                                
               num_shared_elems = num_shared_elems + 1                          
               shared_elems( num_shared_elems ) = next_elem                     
               exit                                                             
            endif                                                               
         enddo                                                                  
      enddo                                                                     
c                                                                               
      if (debug) then                                                           
         write (*,*) '   shared elems:'                                         
         do i = 1, num_shared_elems                                             
            write (*,'(6x,i1,1x,i7)') i, shared_elems(i)                        
         enddo                                                                  
      endif                                                                     
c                                                                               
c          We now know which elements are connected to both of the given        
c          nodes.  Now, loop over the neighbors of next_node.  If any of        
c          the elements connected to the neighbor is one of the elements        
c          we just found, then the new node isn't on the master_line.           
c          Otherwise, it is.                                                    
c                                                                               
      reached_end = .false.                                                     
      new_node = 0                                                              
      node_data_entry = inv_crkpln_nodes ( next_node )                          
      do i = 1, num_neighbors(node_data_entry)                                  
c                                                                               
         neighbor_node = neighbor_nodes (i,node_data_entry)                     
         num_elems = inverse_incidences(neighbor_node)%element_count            
         if (neighbor_node .eq. base_node) cycle                                
         found_new_node = .true.                                                
c                                                                               
c               loop over elements connected to neighbor                        
c                                                                               
         do elem_list = 1, num_elems                                            
c                                                                               
            elem =                                                              
     &        inverse_incidences(neighbor_node)%element_list(elem_list)         
            do j = 1, num_shared_elems                                          
               if (shared_elems(j) .eq. elem)  found_new_node = .false.         
            enddo                                                               
c                                                                               
            if (.not. found_new_node) exit                                      
         enddo                                                                  
c                                                                               
         if (.not. found_new_node) cycle                                        
c                                                                               
         new_node = neighbor_node                                               
         exit                                                                   
c                                                                               
      end do                                                                    
c                                                                               
c         if we were not able to find another node for the master_list,         
c         then we have reached one of the ends of the master_line.              
c                                                                               
      if (new_node .eq. 0) reached_end = .true.                                 
      if (debug) then                                                           
         if (reached_end) then                                                  
            write (*,*) '   reached end.'                                       
         else                                                                   
            write (*,*) '   new_node:',new_node                                 
         endif                                                                  
         write (*,*) '   <<< leaving incrack_find_master_line_node'                     
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine incrack_errmsg               *          
c     *                                                              *          
c     *                       written by : RHD                       *          
c     *                                                              *          
c     *                   last modified : 8/25/21 rhd                 *          
c     *                                                              *          
c     *                service routine for error messages            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine incrack_errmsg( ierrno )                                       
      use global_data ! old common.main
c                                                                               
      implicit none
c
      integer :: ierrno   
c
      integer :: strlng
      character(len=50) :: string                                               
c                                                                               
      select case( ierrno )                                                     
c                                                                               
      case( 1 )                                                                 
         num_error = num_error + 1                                              
         call entits( string, strlng )                                          
         write(out,9001)  string(1:strlng)                                      
         input_ok = .false.                                                     
c                                                                               
      case( 2 )                                                                 
         num_error = num_error + 1                                              
         call entits( string, strlng )                                          
         write(out,9002)  string(1:strlng)                                      
         input_ok = .false.                                                     
c                                                                               
      case( 3 )                                                                 
         num_error = num_error + 1                                              
         write(out,9003)                                                        
         input_ok = .false.                                                     
c                                                                               
      case( 4 )                                                                 
         num_error = num_error + 1                                              
         write(out,9004)                                                        
         input_ok = .false. 
         call scan_flushline
c                                                                               
      case( 5 )                                                                 
         num_error = num_error + 1                                              
         write(out,9005)                                                        
         input_ok = .false. 
         call scan_flushline
c                                                                               
      case( 6 )                                                                 
         num_error = num_error + 1                                              
         write(out,9006)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 7 )                                                                 
         num_error = num_error + 1                                              
         write(out,9007)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 8 )                                                                 
         num_error = num_error + 1                                              
         write(out,9008)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 9 )                                                                 
         num_error = num_error + 1                                              
         write(out,9009)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 10 )                                                                 
         num_error = num_error + 1                                              
         write(out,9010)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 11 )                                                                 
         num_error = num_error + 1                                              
         write(out,9011)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 12 )                                                                 
         num_error = num_error + 1                                              
         write(out,9012)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 13 )                                                                 
         num_error = num_error + 1                                              
         write(out,9013)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 14 )                                                                 
         num_error = num_error + 1                                              
         write(out,9014)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 15 )                                                                 
         num_error = num_error + 1                                              
         write(out,9015)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 16 )                                                                 
         num_error = num_error + 1                                              
         write(out,9016)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 17 )                                                                 
         num_error = num_error + 1                                              
         write(out,9017)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 18 )                                                                 
         num_error = num_error + 1                                              
         write(out,9018)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 19 )                                                                 
         num_error = num_error + 1                                              
         write(out,9019)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 20 )                                                                 
         write(out,9020)                                                        
c
      case( 21 )                                                                 
         num_error = num_error + 1                                              
         write(out,9021)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 22 )                                                                 
         num_error = num_error + 1                                              
         write(out,9022)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 23 )                                                                 
         num_error = num_error + 1                                              
         write(out,9023)  'alpha_min'                                                      
         input_ok = .false. 
         call scan_flushline
c
      case( 24 )                                                                 
         num_error = num_error + 1                                              
         write(out,9024)                                                        
         input_ok = .false. 
         call scan_flushline
c
      case( 25 )                                                                 
         num_error = num_error + 1                                              
         write(out,9025)  'alpha_max'                                                      
         input_ok = .false. 
         call scan_flushline
c
      case( 26 )                                                                 
         num_error = num_error + 1                                              
         write(out,9026) 
         input_ok = .false. 
         call scan_flushline
c
      case( 27 )                                                                 
         num_error = num_error + 1                                              
         write(out,9027) 
         input_ok = .false. 
         call scan_flushline
c
      case( 28 )                                                                 
         num_error = num_error + 1                                              
         write(out,9028) 
         input_ok = .false. 
         call scan_flushline
c
      case( 29 )                                                                 
         num_error = num_error + 1                                              
         write(out,9029) 
         input_ok = .false. 
         call scan_flushline
c
      case( 30 )                                                                 
         num_error = num_error + 1                                              
         write(out,9030) 
         input_ok = .false. 
         call scan_flushline
c
      case( 31 )                                                                 
         num_error = num_error + 1                                              
         write(out,9031) 
         input_ok = .false. 
         call scan_flushline
c
      case( 32 )                                                                 
         num_error = num_error + 1                                              
         write(out,9032) 
         input_ok = .false. 
         call scan_flushline
c
      case( 33 )                                                                 
         num_error = num_error + 1                                              
         write(out,9033) 
         input_ok = .false. 
         call scan_flushline
c
      case( 34 )                                                                 
         num_error = num_error + 1                                              
         write(out,9034) 
         input_ok = .false. 
         call scan_flushline
c
      case( 35 )                                                                 
         num_error = num_error + 1                                              
         write(out,9035) 
         input_ok = .false. 
         call scan_flushline
c
      case( 36 )                                                                 
         num_error = num_error + 1                                              
         write(out,9036) 
         input_ok = .false. 
         call scan_flushline
c
      case( 37 )                                                                 
         num_error = num_error + 1                                              
         write(out,9037) 
         input_ok = .false. 
         call scan_flushline
c
      case( 38 )                                                                 
         num_error = num_error + 1                                              
         write(out,9038) 
         input_ok = .false. 
         call scan_flushline
c
      case( 39 )                                                                 
         num_error = num_error + 1                                              
         write(out,9039) 
         input_ok = .false. 
         call scan_flushline
c
      case( 40 )                                                                 
         num_error = num_error + 1                                              
         write(out,9040) 
         input_ok = .false. 
         call scan_flushline
c
      case( 41 )                                                                 
         num_error = num_error + 1                                              
         write(out,9041) 
         input_ok = .false. 
         call scan_flushline
c
      case( 42 )                                                                 
         num_error = num_error + 1                                              
         write(out,9042) 
         input_ok = .false. 
         call scan_flushline
c
      case( 43 )                                                                 
         num_error = num_error + 1                                              
         write(out,9043) 
         input_ok = .false. 
         call scan_flushline
c
      case( 44 )                                                                 
         num_error = num_error + 1                                              
         write(out,9044) 
         input_ok = .false. 
         call scan_flushline
c
      case( 45 )                                                                 
         num_error = num_error + 1                                              
         write(out,9045) 
         input_ok = .false. 
         call scan_flushline
c
      case( 46 )                                                                 
         num_error = num_error + 1                                              
         write(out,9046) 
         input_ok = .false. 
         call scan_flushline
c
      case( 50 )                                                                 
         num_error = num_error + 1                                              
         write(out,9050) 
         input_ok = .false. 
         call scan_flushline
c
      case( 56 )                                                                 
         num_error = num_error + 1                                              
         write(out,9056) 
         input_ok = .false. 
         call scan_flushline
c
      case( 57 )                                                                 
         num_error = num_error + 1                                              
         write(out,9057) 
         input_ok = .false. 
         call scan_flushline
c
      case( 58 )                                                                 
         num_error = num_error + 1                                              
         write(out,9058) 
         input_ok = .false. 
         call scan_flushline
c
      case( 60 )                                                                 
         num_error = num_error + 1                                              
         write(out,9060) 
         input_ok = .false. 
         call scan_flushline
c
      case( 61 )                                                                 
         num_error = num_error + 1                                              
         write(out,9061) 
         input_ok = .false. 
         call scan_flushline
c
      case( 62 )                                                                 
         num_error = num_error + 1                                              
         write(out,9062) 
         input_ok = .false. 
         call scan_flushline
c
      case( 63 )                                                                 
         num_error = num_error + 1                                              
         write(out,9063) 
         input_ok = .false. 
         call scan_flushline
c
      case( 64 )                                                                 
         num_error = num_error + 1                                              
         write(out,9064) 
         input_ok = .false. 
         call scan_flushline
c
      case( 65 )                                                                 
         num_error = num_error + 1                                              
         write(out,9065) 
         input_ok = .false. 
         call scan_flushline
c
      case( 66 )                                                                 
         num_error = num_error + 1                                              
         write(out,9066) 
         input_ok = .false. 
         call scan_flushline
c
      case( 67 )                                                                 
         num_error = num_error + 1                                              
         write(out,9067) 
         input_ok = .false. 
         call scan_flushline
c
      case( 68 )                                                                 
         num_error = num_error + 1                                              
         write(out,9068) 
         input_ok = .false. 
         call scan_flushline
c
      case( 69 )                                                                 
         num_error = num_error + 1                                              
         write(out,9069) 
         input_ok = .false. 
         call scan_flushline
c
      case( 70 )                                                                 
         num_error = num_error + 1                                              
         write(out,9070) 
         input_ok = .false. 
         call scan_flushline
c
      case( 71 )                                                                 
         num_error = num_error + 1                                              
         write(out,9071) 
         input_ok = .false. 
         call scan_flushline
c
      case( 72 )                                                                 
         num_error = num_error + 1                                              
         write(out,9072) 
         input_ok = .false. 
         call scan_flushline
c
      case( 73 )                                                                 
         num_error = num_error + 1                                              
         write(out,9073) 
         input_ok = .false. 
         call scan_flushline
c
      case( 74 )                                                                 
         num_error = num_error + 1                                              
         write(out,9074) 
         input_ok = .false. 
         call scan_flushline
c
      case( 75 )                                                                 
         num_error = num_error + 1                                              
         write(out,9075) 
         input_ok = .false. 
         call scan_flushline
c
      case default                                                              
        write(out,9999)                                                         
        call die_abort                                                                   
      end select                                                                
c                                                                               
      return                                                                    
c                                                                               
 9001 format(/1x,'>>>>> error: expecting fraction value for extinction,'        
     & /14x,'in interface element with PPR cohesive option.',                   
     & /14x,'scanning: ',a, '. command ignored...',/)                           
c                                                                               
 9002 format(/1x,'>>>>> error: unrecognized command for PPR cohesive',          
     & /14x,'crack growth. Scanning: ',a,'. command ignored',/)                 
c                                                                               
 9003 format(/1x,'>>>>> error: displacement fraction invalid value',            
     & /14x,'must be >0 and .le. 1.0. command ignored',/)                       
c                                                                               
 9004 format(/1x,'>>>>> error: expecting value for SMCS parameter',            
     & /14x,'line ignored',/)                       
c                                                                               
 9005 format(/1x,'>>>>> error: unrecognized SMCS parameter',            
     & /14x,'line ignored',/)                       
c                                                                               
 9006 format(/1x,'>>>>> error: expecting killed element limit value',            
     & /14x,'line ignored',/)                       
c                                                                               
 9007 format(/1x,'>>>>> error: expecting number of top elements ',            
     & /14x,'to print.  line ignored',/)                       
c                                                                               
 9008 format(/1x,'>>>>> error: expecting keyword type, states, ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9009 format(/1x,'>>>>> error: expecting integer value 1 or 2 ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9010 format(/1x,'>>>>> error: smcs type must be valid integer ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9011 format(/1x,'>>>>> error: unrecognized smcs parameter ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9012 format(/1x,'>>>>> error: invalid step list format... ',    
     & /14x,'(including -all- steps not allowed)',        
     & /14x,'line ignored',/)                       
c                                                                               
 9013 format(/1x,'>>>>> error: step list overflow ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9014 format(/1x,'>>>>> error: no steps list found ... ',            
     & /14x,'(including -all- steps not allowed)',        
     & /14x,'line ignored',/)                       
c                                                                               
 9015 format(/1x,'>>>>> error: expecting keyword on or off ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9016 format(/1x,'>>>>> error: expecting keyword text or stream ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9017 format(/1x,'>>>>> error: feature not supported yet with MPI... ',            
     & /14x,'line ignored',/)                       
c
 9018 format(/1x,'>>>>> error: expecting file name... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9019 format(/1x,'>>>>> error: expecting keyword file ... ',            
     & /14x,'line ignored',/)                  
c                                                                               
 9020 format(/1x,'>>>>> warning: existing file deleted...',/)                  
c                                                                               
 9021 format(/1x,'>>>>> error: unrecognized keyword ... ',            
     & /14x,'line ignored',/)                  
c                                                                               
 9022 format(/1x,'>>>>> error: expecting keyword on or off ... ',            
     & /14x,'line ignored',/)                  
c                                                                               
 9023 format(/1x,'>>>>> error: expecting number after keyword ',a,            
     & /14x,'line ignored',/)                  
c                                                                               
 9024 format(/1x,'>>>>> error: alpha_min must be > 0 & < 1',            
     & /14x,'line ignored',/)                  
c
 9025 format(/1x,'>>>>> error: expecting number after keyword ',a,            
     & /14x,'line ignored',/)                  
c                                                                               
 9026 format(/1x,'>>>>> error: alpha_max must be > 1',            
     & /14x,'line ignored',/)                  
c                                                                               
 9027 format(/1x,'>>>>> error: expecting # allowable elements',            
     & /14x,'line ignored',/)                  
c                                                                               
 9028 format(/1x,'>>>>> error: expecting keyword: on or off',            
     & /14x,'command ignored',/)                  
c                                                                               
 9029 format(/1x,'>>>>> error: expecting keyword: type',            
     & /14x,'command ignored',/)                  
c                                                                               
 9030 format(/1x,'>>>>> error: expecting <integer> type number',            
     & /14x,'command ignored',/)    
c              
 9031 format(/1x,'>>>>> error: invalid regularization type',            
     & /14x,'command ignored',/)        
c          
 9032 format(/1x,'>>>>> error: expecting keyword: L or length_scale',            
     & /14x,'command ignored',/)        
c          
 9033 format(/1x,'>>>>> error: expecting length value',            
     & /14x,'command ignored',/)        
c          
 9034 format(/1x,'>>>>> error: invalid value for length-scale',            
     & /14x,'command ignored',/)        
c          
 9035 format(/1x,'>>>>> error: expecting keyword: up_max',            
     & /14x,'command ignored',/)        
c          
 9036 format(/1x,'>>>>> error: expecting value for up_max',            
     & /14x,'command ignored',/)  
c      
 9037 format(/1x,'>>>>> error: invalid value for up_max',            
     & /14x,'command ignored',/)        
c      
 9038 format(/1x,'>>>>> error: expecting keyword: num_points',            
     & /14x,'command ignored',/)        
c      
 9039 format(/1x,'>>>>> error: invalid number of points. must be',            
     & /14x,'>0 and <=10. command ignored',/)        
c      
 9040 format(/1x,'>>>>> error: expecting number of points',            
     & /14x,'command ignored',/)        
c      
 9041 format(/1x,'>>>>> error: expecting up value',            
     & /14x,'command ignored',/)        
c      
 9042 format(/1x,'>>>>> error: expecting d value',            
     & /14x,'command ignored',/)        
c      
 9043 format(/1x,'>>>>> error: 1st up value must = 0.0')
c      
 9044 format(/1x,'>>>>> error: last up value must = 1.0')
c      
 9045 format(/1x,'>>>>> error: last d value must = 1.0')
c      
 9046 format(/1x,'>>>>> error: bad curve points',            
     & /14x,'command ignored',/)        
c
 9050 format(/1x,'>>>>> error: incomplete regularization data',            
     & /14x,'command ignored',/)                  
c
 9056 format(/1x,'>>>>> error: expecting keyword: Gf',            
     & /14x,'command ignored',/)                  
c
 9057 format(/1x,'>>>>> error: expecting Gf value',            
     & /14x,'command ignored',/)                  
c
 9058 format(/1x,'>>>>> error: Gf must be > 0',            
     & /14x,'command ignored',/)                  
c
 9060 format(/1x,'>>>>> error: expecting keyword: m',            
     & /14x,'command ignored',/)                  
c
 9061 format(/1x,'>>>>> error: expecting m exponent value',            
     & /14x,'command ignored',/)                  
c
 9062 format(/1x,'>>>>> error: m must be > 0',            
     & /14x,'command ignored',/)                  
c
 9063 format(/1x,'>>>>> error: expecting keyword: alpha',            
     & /14x,'command ignored',/)                  
c
 9064 format(/1x,'>>>>> error: expecting alpha value',            
     & /14x,'command ignored',/)                  
c     
 9065 format(/1x,'>>>>> error: regularization not yet supported with',            
     & ' MPI execution',
     & /14x,'No solution allowed.',/)                  
c     
 9066 format(/1x,'>>>>> error: distortion metric not yet supported with',            
     & ' MPI execution',
     & /14x,'No solution allowed.',/)     
c                                                                               
 9067 format(/1x,'>>>>> error: expecting keyword: distortion ... ',            
     & /14x,'line ignored',/)                  
c                                                                               
 9068 format(/1x,'>>>>> error: expecting keyword: on or off ... ',            
     & /14x,'line ignored',/)                  
c                                                                               
 9069 format(/1x,'>>>>> error: expecting value for critical ',
     &      ' distortion ... ', /14x,'line ignored',/)                  
c                                                                               
 9070 format(/1x,'>>>>> error: expecting keyword: plastic ... ',            
     & /14x,'line ignored',/)                  
c                                                                               
 9071 format(/1x,'>>>>> error: expecting value for plastic ',
     &      ' strain limit ... ', /14x,'line ignored',/)      
c                                                                               
 9072 format(/1x,'>>>>> error: expecting keyword: deleted ... ',            
     & /14x,'line ignored',/)                       
c                                                                               
 9073 format(/1x,'>>>>> error: expecting keyword: file for  ... ',            
     & /14x,'for write file option. line ignored',/)                       
c                                                                               
 9074 format(/1x,'>>>>> error: expecting keyword: on | off for  ... ',            
     & /14x,'for write file option. line ignored',/)                       
c                                                                               
 9075 format(/1x,'>>>>> error: expecting keyword: Oddy  ... ',            
     & /14x,'line ignored',/)                       
c                  
 9999 format(/1x,'>>>>> Fatal Error: routine incrck_errmsg.',                   
     &   /16x,   'should have not reach this point.')                           
c                                                                                
      end                                                                       
c         
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine incrack_delete_elements              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/4/2019 rhd              *          
c     *                                                              *          
c     *            user-directed element deletion for crack growth   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c
      subroutine incrack_delete_elements( sbflg1, sbflg2 )
      use global_data ! old common.main
      use allocated_integer_list
c
      use main_data
      use damage_data, only : num_user_kill_elems, crack_growth_type,
     &                        user_kill_list_now, max_dam_state
c
      implicit none
c
c              parameters
c
      logical :: sbflg1, sbflg2
c
c              locals
c
      integer :: list_size, lenlst, errnum, elem, icn,
     &           iplist, bad_count, temp_int
      integer, allocatable :: ele_intlst(:)
      logical :: debug, bad_list
      logical, external :: matchs, numi
      double precision :: dumd
      real :: dumr
      character(len=1) :: dums
c
c              if sub flag 1 is on, there is re-entry after an error
c              in input. if we have bad input somewhere here,
c              release the data structure on exit
c
c              delete (elements) <list (release) steps <integer>
c
      debug = .false.
      if( debug ) then
         write(out,*) '... entered delete elements (user-defined) ...'
         write(out,*) '    sbflg1, sbflg2: ', sbflg1, sbflg2
      end if
c 
      if( sbflg1 ) then
         write(out,9000) ; num_error = num_error + 1
         num_user_kill_elems = 0
         return
      end if
c
c      if( crack_growth_type .ne. 5 ) then
c         num_user_kill_elems = 0
c         write(out,9040)
c         write(*,*) '... crack_growth_type: ',crack_growth_type
c         num_error = num_error + 1
c         call scan_flushline
c      end if
c
      num_user_kill_elems = 0
      if( matchs('elem',4) ) call splunj   ! optional word
      allocate( ele_intlst(20) )
      call scan
      call trlist_allocated( ele_intlst, list_size, 
     &                       noelem, lenlst, errnum )
c
c              1 = no error. list found
c              2 = parse rules failed in list
c              3 = list overflowed its maximum length of mxlsz
c              4 = no list found 
c
      if( debug ) write(out,*)
     & '    ... trlist done, errnum, lenlst: ', errnum, lenlst
c
      if( errnum .eq. 4 ) then 
          write(out,9000) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
      if( errnum .eq. 2 ) then
          write(out,9010) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
      if( errnum .eq. 3 ) then
          write(out,9020) 
          call scan_flushline
          num_error = num_error + 1
          return
      end if
c
      if( errnum .ne. 1 ) then  ! something bad wrong
         write(out,9050)
         call die_gracefully
      end if
c

      icn = 0; iplist = 1; bad_list = .false.; bad_count = 0
      do while ( iplist .ne. 0 )
        call trxlst( ele_intlst, lenlst, iplist, icn, elem )
        if( elem <= 0 .or. elem > noelem ) then
            if( bad_count > 10 ) then 
              write(out,9060) 
              bad_list = .true.
              exit
            end if
            write(out,9030 ) elem
            bad_list = .true.
            bad_count = bad_count + 1
            num_error = num_error + 1
            cycle
        end if
        num_user_kill_elems = num_user_kill_elems + 1 
        if( num_user_kill_elems > 100 ) then
          write(out,9070) 100
          num_error = num_error + 1
          bad_list = .true.
          exit        
        end if
        user_kill_list_now(num_user_kill_elems) = elem
      end do  ! extracting elements from list
c                                                                               
c              read/processed. cleanup.
c
      crack_growth_type = 5
      if( bad_list ) then
          num_user_kill_elems = 0
           crack_growth_type = 0
      end if
      sbflg1 = .true.
      sbflg2 = .true.
c
      if( debug ) then
         write(out,*) '... leaving incrack_delete_elements ...'
         write(out,*) '  num_user_kill_elems: ',num_user_kill_elems
         if( num_user_kill_elems == 0 ) return
         do icn = 1, num_user_kill_elems
          write(out,*) '    ... delete element: ',
     &                 user_kill_list_now(icn)
         end do
      end if      
c
      return
c
 9000 format(/1x,'>>>>> error: no element list found',
     &   /,14x,'read new line ...',/)
 9010 format(/1x,'>>>>> error: integer list rules failed....',
     &   ' read new line ...',/)
 9020 format(/1x,'>>>>> error: too many elements in list....',
     &   ' read new line ...',/)
 9030 format(/1x,'>>>>> error: invalid element number ',i8)
 9040 format(/1x,'>>>>> error: crack growth type previously defined',
     &     /,    '             is not user-directed' )
 9050 format(/1x,'>>>>> FATAL ERROR: delete_elements inconsistency',
     &   /,14x,'job terminated now ...',/)
 9060 format(/1x,'>>>>> error: further messages suppressed ...' )
 9070 format(/1x,'>>>>> error: too many elements in list. limit: ',i4 )
c
      end
