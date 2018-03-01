c     ****************************************************************          
c     *                                                              *          
c     *                subroutine release_constraints                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 6/20/2016 rhd               *          
c     *                                                              *          
c     *     read/store data for the release constraints command      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine release_constraints( sbflg1, sbflg2 )                          
      use global_data ! old common.main
c                                                                               
      use main_data, only : cnstrn_in, release_cons_table,                      
     &                      release_cons_steps, mdiag, rload                    
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       parameters                                              
c                                                                               
      logical :: sbflg1, sbflg2                                                 
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: intlst(mxlsz), lenlst, param, errnum, dummy,                   
     &           iplist, icn, node, idof, dof, nsteps, i                        
      integer, save :: count_stored                                             
      logical :: release_flags(3), found_list, bad_list, bad                    
      logical, external :: match, matchs, endcrd, true, numd, integr,           
     &         realn, match_exact                                               
      logical :: debug                                                          
      logical, save :: delete_release_cons                                      
      real ::  dumr                                                             
      double precision :: dumd, zero, d32460, react                             
      character :: dums*1, curtyp *1, dof_names(3)*1                            
      data dof_names / 'u', 'v', 'w' /                                          
      data  d32460, zero /  32460.0d00, 0.0d00 /                                
c                                                                               
c                       if sub flag 1 is on, there is re-entry                  
c                       after an error in input. if we have bad                 
c                       input somewhere here,                                   
c                       release the data structure on exit                      
c                                                                               
      debug = .false.                                                           
      count_stored = 0                                                          
c                                                                               
      if( debug ) then                                                          
         write(out,*) '... entered release_constraints ...'                     
         write(out,*) '    sbflg1, sbflg2: ', sbflg1, sbflg2                    
      end if                                                                    
c                                                                               
      if( sbflg1 ) then                                                         
         write(out,9000) ; num_error = num_error + 1                            
         go to 100                                                              
      else                                                                      
         delete_release_cons = .false.                                          
      end if                                                                    
c                                                                               
c                       get number of steps used to release reaction            
c                       forces to zero. default = 1 step                        
c                                                                               
      release_cons_steps = 1                                                    
      if( matchs( 'constraints', 4 ) ) call splunj                              
      if( matchs( 'steps',4 ) ) call splunj                                     
      if( matchs( '=',1 ) ) call splunj                                         
      if( .not. integr( release_cons_steps ) ) then                             
         if( .not. endcrd(dummy) ) then                                         
           write(out,9120)                                                      
           num_error = num_error + 1                                            
         end if                                                                 
      end if                                                                    
c                                                                               
c                       allocate release constraints table if needed            
c                                                                               
      if( .not. allocated( release_cons_table ) ) then                          
        if( debug ) write(out,*) '... allocating release_cons_table'            
        allocate( release_cons_table(3,nonode) )                                
        do i = 1, nonode                                                        
           release_cons_table(1:3,i)%num_release_steps = 0                      
           release_cons_table(1:3,i)%remaining_steps_for_release = 0            
           release_cons_table(1:3,i)%reaction_force = zero                      
        end do                                                                  
      end if                                                                    
c                                                                               
 100  continue                                                                  
      do !   outer read loop over lines of node lists and dofs                  
c                                                                               
c                       get a list of nodes and the released                    
c                       directions (u,v,w)                                      
c                                                                               
      call readsc                                                               
      call release_cons_scan                                                    
      if( .not. found_list ) exit  ! out of read loop. back to main             
      if( bad_list ) then                                                       
        delete_release_cons = .true.                                            
        cycle                                                                   
      end if                                                                    
c                                                                               
c                       list of nodes and release directions appears to         
c                       be ok. Process node list to store reactions             
c                       in release constraints data structure                   
c                                                                               
      icn    = 0                                                                
      iplist = 1                                                                
      do while ( iplist .ne. 0 )                                                
         call trxlst( intlst, lenlst, iplist, icn, node )                       
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure                    
c                       and is positive.                                        
c                                                                               
         if( node .gt. nonode ) then                                            
            param = node                                                        
            call errmsg(16,param,dums,dumr,dumd)                                
            cycle                                                               
         end if                                                                 
         if( node .le. 0 )  then                                                
            param = node                                                        
            call errmsg(58,param,dums,dumr,dumd)                                
            cycle                                                               
         end if                                                                 
c                                                                               
c                       before storing make sure the dof for node has           
c                       (1) is not already in release, (2) has a                
c                       constraint imposed                                      
c                                                                               
         bad = .false.                                                          
         do idof = 1, 3                                                         
           curtyp = dof_names(idof)                                             
           if( .not. release_flags(idof) ) cycle                                
           nsteps = release_cons_table(idof,node)%num_release_steps             
           if( nsteps .ne. 0 ) then                                             
              write(out,9100) node, dof_names(1)                                
              bad = .true.                                                      
           end if                                                               
           dof = dstmap(node) + idof - 1                                        
           if( cnstrn_in(dof) .eq. d32460 ) then                                
              write(out,9110) node, curtyp                                      
              bad = .true.                                                      
           end if                                                               
         end do                                                                 
c                                                                               
         if( bad ) then                                                         
           num_error = num_error + 1                                            
           cycle                                                                
         end if                                                                 
c                                                                               
c                       pull reaction values from global load vector            
c                       and store. remove constraint from constraints           
c                       data structure.                                         
c                                                                               
         do idof = 1, 3                                                         
           if( .not. release_flags(idof) ) cycle                                
           dof = 3 * (node-1) + idof                                            
           release_cons_table(idof,node)%num_release_steps =                    
     &            release_cons_steps                                            
           release_cons_table(idof,node)%remaining_steps_for_release =          
     &            release_cons_steps                                            
           react = -rload(dof) + ifv(dof) + mdiag(dof)*a(dof)                   
           release_cons_table(idof,node)%reaction_force = react                 
           call release_cons_update_constraints( dof )                          
           count_stored = count_stored + 1                                      
         end do                                                                 
c                                                                               
      end do  ! while loop to process nodes in list                             
c                                                                               
      end do  ! get next input line                                             
c                                                                               
c                                                                               
      sbflg1         = .true.                                                   
      sbflg2         = .true.                                                   
      if( debug ) then                                                          
         write(out,*) '... count_stored: ',count_stored                         
         write(out,*) '... delete_release_cons: ', delete_release_cons          
      end if                                                                    
      if( delete_release_cons ) then                                            
        if( allocated( release_cons_table ) )                                   
     &      deallocate( release_cons_table )                                    
      end if                                                                    
      call release_empty_cons_table                                             
      if( debug ) write(out,*) '... leaving release_constraints ...'            
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/1x,'>>>>> error: syntax error in release constraints',            
     &   /,14x,'input. read new line ...')                                      
 9100 format(/1x,'>>>>> error: node: ',i7,' dof: ',a1,                          
     &  /14x,'already in release')                                              
 9110 format(/1x,'>>>>> error: node: ',i7,' dof: ',a1,                          
     &  /14x,'has no constraint to release')                                    
 9120 format(/1x,'>>>>> error: invalid number of release steps',                
     &   /,14x,'input. read new line ...')                                      
                                                                                
      contains                                                                  
c     ========                                                                  
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine release_cons_scan                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/222104 RHD              *          
c     *                                                              *          
c     *     scan/store list of nodes and constraint directions to    *          
c     *     be released                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine release_cons_scan                                              
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: strlng                                                         
      character :: string*80                                                    
      logical :: ok                                                             
c                                                                               
      if( debug ) write(out,*) '    ... entered release_cons_scan ...'          
c                                                                               
c                       input node list                                         
c                                                                               
      if( matchs( 'plane', 5 ) ) call splunj ! implement later                  
      call trlist( intlst, mxlsz, nonode, lenlst, errnum)                       
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       in these last two cases, the illegal list               
c                       will be ignored and a new node list will                
c                       be sought. a value of 4 indicates that no list          
c                       was found. in this case, release cons input             
c                       has ceased.                                             
c                                                                               
      if( debug ) write(out,*)                                                  
     & '    ... trlist done, errnum, lenlst: ',errnum,lenlst                    
c                                                                               
      bad_list = .false.                                                        
      found_list = .false.                                                      
c                                                                               
      if( errnum .eq. 4 ) return                                                
      if( errnum .eq. 2 ) then                                                  
         param = 1                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         bad_list = .true.                                                      
         call scan_flushline                                                    
         return                                                                 
      end if                                                                    
c                                                                               
      if( errnum .eq. 3 ) then                                                  
         param = 2                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         bad_list = .true.                                                      
         call scan_flushline                                                    
         return                                                                 
      end if                                                                    
c                                                                               
      if( errnum .eq. 1 ) then                                                  
          found_list = .true.                                                   
          call backsp(1)                                                        
          if( true(dummy) ) call splunj                                         
      else                                                                      
         write(out,9000)                                                        
         call die_abort                                                         
      end if                                                                    
c                                                                               
c                       get one or more directions to be released               
c                                                                               
      release_flags(1:3) = .false.                                              
      do                                                                        
         if( match_exact( ',',1 ) ) then                                        
           call splunj ! do nothing forces compiler to execute                  
           cycle                                                                
         end if                                                                 
         if( endcrd(dummy) ) exit                                               
         if( match_exact( 'u',1 ) ) then                                        
             release_flags(1) = .true.                                          
             cycle                                                              
         end if                                                                 
         if( match_exact( 'v',1 ) ) then                                        
             release_flags(2) = .true.                                          
             cycle                                                              
         end if                                                                 
         if( match_exact( 'w',1 ) ) then                                        
             release_flags(3) = .true.                                          
             cycle                                                              
         end if                                                                 
c                                                                               
         bad_list = .true.                                                      
         num_error = num_error + 1                                              
         call entits( string, strlng )                                          
         write(out,9010) string(1:strlng)                                       
         call scan_flushline                                                    
         return                                                                 
      end do                                                                    
c                                                                               
      ok = release_flags(1) .or. release_flags(2) .or. release_flags(3)         
      if( ok ) return                                                           
      write(out,9020)                                                           
      num_error = num_error + 1                                                 
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/1x,'>>>>> error: invalid return on trlist in ',                   
     &  /14x,'release_cons_scan. system error. job aborted.',//)                
 9010 format(/1x,'>>>>> error: unrecognized data in list of released',          
     &  /14x,'dof. scanning: ', a )                                             
 9020 format(/1x,'>>>>> error: no valid components to release found')           
c                                                                               
      end subroutine release_cons_scan                                          
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine release empty table               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/23/03                   *          
c     *                                                              *          
c     *     skips to the end of a logical line of the input file     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine release_empty_cons_table                                       
c                                                                               
c                       delete an empty release constraints table               
c                                                                               
      logical :: found_entry                                                    
      integer :: i, n1, n2, n3                                                  
c                                                                               
      if(  .not. allocated( release_cons_table ) ) return                       
c                                                                               
      found_entry = .false.                                                     
      do i = 1, nonode                                                          
          n1 = release_cons_table(1,i)%remaining_steps_for_release              
          n2 = release_cons_table(2,i)%remaining_steps_for_release              
          n3 = release_cons_table(3,i)%remaining_steps_for_release              
          if( n1 .ne. 0 ) found_entry = .true.                                  
          if( n2 .ne. 0 ) found_entry = .true.                                  
          if( n3 .ne. 0 ) found_entry = .true.                                  
          if( found_entry ) exit                                                
      end do                                                                    
      if( .not. found_entry ) then                                              
           deallocate( release_cons_table )                                     
           if( debug )  write(out,*) '... deleted empty table ...'              
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end subroutine release_empty_cons_table                                   
c                                                                               
      end subroutine release_constraints                                        
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *           subroutine release_cons_update_constraints         *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/24/03                   *          
c     *                                                              *          
c     * remove constraint on 1 dof from constraint data structure    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine release_cons_update_constraints( sdof )                        
      use global_data ! old common.main
c                                                                               
      use main_data, only : cnstrn_in, cnstrn                                   
      use damage_data, only : csttail                                           
c                                                                               
      implicit none                                                             
c                                                                               
      integer :: sdof                                                           
c                                                                               
      integer :: cst_ptr, above_ptr                                             
      double precision :: d32460                                                
      data d32460 / 32460.0d00 /                                                
c                                                                               
c             traverse the singly-linked list of constraints in the             
c             model. If the user has input new constraints, then some           
c             of the previously released nodes may be re-constrained.           
c             If a released node has been re-constrained, remove the            
c             constraint.                                                       
c                                                                               
c                set up indexes for constraint linked list                      
c                                                                               
      cst_ptr   = csthed                                                        
      above_ptr = -1                                                            
c                                                                               
c                enter top of constraint linked list loop. run till             
c                dof is found or end of constraints.                            
c                                                                               
      do while ( cst_ptr .ne. -1 )                                              
c                                                                               
      if( cst_ptr .eq. sdof ) then                                              
         cnstrn(cst_ptr)    = d32460                                            
         cnstrn_in(cst_ptr) = d32460                                            
c                                                                               
         if( above_ptr .eq. -1 ) then                                           
c                                                                               
c                           at top of list.  move head pointer.                 
c                                                                               
            csthed          = cstmap(cst_ptr)                                   
            cstmap(cst_ptr) = 0                                                 
            cst_ptr         = csthed                                            
         else                                                                   
c                                                                               
c                           in middle of list. correct link indexes             
c                                                                               
            cstmap(above_ptr) = cstmap(cst_ptr)                                 
            cstmap(cst_ptr)   = 0                                               
            if ( csttail .eq. cst_ptr ) csttail = above_ptr                     
            cst_ptr = cstmap(above_ptr)                                         
         end if                                                                 
         exit                                                                   
      end if                                                                    
c                                                                               
c                    examine next constraint in list                            
c                                                                               
      above_ptr = cst_ptr                                                       
      cst_ptr   = cstmap(cst_ptr)                                               
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine incon                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 06/22/2016 rhd             *          
c     *                                                              *          
c     *     input of nodal displacement constraints                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine incon( sbflg1, sbflg2, olddof )                                
      use global_data ! old common.main
c                                                                               
      use main_data, only : trn, trnmat, cnstrn, cnstrn_in,                     
     &                      inverse_incidences                                  
      use mod_mpc, only : mpcs_exist, num_user_mpc, user_mpc_table              
      use damage_data, only : csttail                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       locally allocated                                       
c                                                                               
      double precision                                                          
     &  convec(mxndof), trans(mxndof,mxndof), tval, cval, dumd, zero,           
     &  one, d32460, rlen1, rlen2, rlen3, rottol, t11, t12, t13,                
     &  t21, t22, t23, t31, t32, t33                                            
      real dumr                                                                 
      character dums, curtyp *1                                                 
      logical sbflg1, sbflg2, skew, inpflg(mxndof), defcon(3),                  
     &        cons_defined, rflag1, rflag2                                      
      logical matchs, endcrd, true, numd, integr, realn                         
      dimension intlst(mxlsz)                                                   
      data zero, one, d32460, cons_defined / 0.0, 1.0, 32460.0,                 
     &                                       .false. /                          
      data rottol / 0.0001 /                                                    
c                                                                               
c                       if sub flag 1 is on, there is reentry into              
c                       incon after an error in constraints input.              
c                                                                               
      if( sbflg1 ) then                                                         
         call errmsg(71,dum,dums,dumr,dumd)                                     
         go to 710                                                              
      end if                                                                    
c                                                                               
c                       make sure that any previous transformation              
c                       matrices created due to contact are re-applied          
c                       to rotate the corresponding degrees of freedom          
c                       to global coordinates.                                  
c                                                                               
      call contact_remove (.true.)                                              
c                                                                               
c                       warn user that all constraints are being                
c                       destroyed and that all constraints must be              
c                       be re-defined....  deallocate also deletes
c                       allocatable sub-objects (F2003)                                     
c                                                                               
      new_constraints = .true.                                                  
      if ( cons_defined ) then                                                  
         call errmsg(224,dum,dums,dumr,dumd)                                    
         call errmsg2(59,dum,dums,dumr,dumd)                                    
         if (mpcs_exist) then                                                   
            if (allocated(user_mpc_table)) deallocate(user_mpc_table)           
            num_user_mpc = 0                                                    
            mpcs_exist = .false.                                                
         end if                                                                 
      endif                                                                     
      cons_defined = .true.                                                     
c                                                                               
c                       initialize the constraints link list.                   
c                                                                               
      csthed = -1                                                               
      olddof = 0                                                                
      do i = 1, nodof                                                           
         cstmap(i)    = 0                                                       
         cnstrn_in(i) = d32460                                                  
      end do                                                                    
c                                                                               
c                       initialize transformation matrix indexes                
c                       and flags. deallocate any old transformation            
c                       matrices                                                
c                                                                               
      do i = 1, nonode                                                          
         if ( trn(i) ) call allo_trnmat(i,2,dum)                                
         trn(i) = .false.                                                       
      end do                                                                    
c                                                                               
c                      initialize defined constraint flags to catch             
c                      problems with one or more totally unconstrained          
c                      direction                                                
c                                                                               
      defcon(1) = .false.                                                       
      defcon(2) = .false.                                                       
      defcon(3) = .false.                                                       
c                                                                               
c                                                                               
                                                                                
 710  call readsc                                                               
 711  continue                                                                  
c                                                                               
c                       look for transformation matrix input. set               
c                       flag depending on whether or not it is found.           
c                                                                               
      if( matchs('transformation',5) ) then                                     
         if( matchs('matrix',5) ) call splunj                                   
         skew = .true.                                                          
      else                                                                      
         skew = .false.                                                         
      end if                                                                    
c                                                                               
c                       look for dump command to dump out all the               
c                       constraints                                             
c                                                                               
      if ( matchs('dump',4) ) then                                              
         call con_dump(olddof)                                                  
         go to 710                                                              
      end if                                                                    
c                                                                               
c                       look for a command to impose constraints                
c                       on all nodes on a plane                                 
c                                                                               
      if ( matchs('plane',4) ) then                                             
         call inconplane( olddof, defcon )                                      
         go to 710                                                              
      end if                                                                    
c                                                                               
c                       look for a command to start input of multi-point        
c                       constraint equations. returns when line does not        
c                       start with an integer. Did a reset, true before         
c                       return.                                                 
c                                                                               
      if ( matchs('multipoint',5) ) then                                        
         call incon_mpcs                                                        
         go to 711                                                              
      end if                                                                    
c                                                                               
c                       input node list                                         
c                                                                               
      call trlist(intlst,mxlsz,nonode,lenlst,errnum)                            
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       in these last two cases, the illegal list               
c                       will be ignored and a new node list will                
c                       be sought. a value of 4 indicates that no list          
c                       was found. in this case, constraints input              
c                       has ceased.                                             
c                                                                               
      if( errnum .eq. 2 ) then                                                  
         param= 1                                                               
         call errmsg(24,param,dums,dumr,dumd)                                   
         go to 710                                                              
      else if( errnum .eq. 3 ) then                                             
         param = 2                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         go to 710                                                              
      else if( errnum .eq. 4 ) then                                             
         go to 9999                                                             
      else                                                                      
         if( errnum .eq. 1 ) then                                               
            call backsp(1)                                                      
            if( true(dummy) ) go to 715                                         
         end if                                                                 
         param = 3                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         go to 710                                                              
      end if                                                                    
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     input of transformation matrix at list of      *        
c *                     nodes which defines constraint compatable      *        
c *                     global coordinates.                            *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 715  if( skew ) then                                                           
c                                                                               
c                       initialize nodal input transformation array             
c                       and the input flag vector.                              
c                                                                               
         do i = 1, mxndof                                                       
            inpflg(i) = .false.                                                 
            do j = 1, mxndof                                                    
               if( i .eq. j ) then                                              
                  trans(i,j) = one                                              
               else                                                             
                  trans(i,j) = zero                                             
               end if                                                           
            end do                                                              
         end do                                                                 
c                                                                               
c                       branch on row of matrix currently being input.          
c                                                                               
 717     if( matchs('row_1',5) ) then                                           
            row = 1                                                             
            go to 719                                                           
         end if                                                                 
c                                                                               
         if( matchs('row_2',5) ) then                                           
            row = 2                                                             
            go to 719                                                           
         end if                                                                 
c                                                                               
         if( matchs('row_3',5) ) then                                           
            row = 3                                                             
            go to 719                                                           
         end if                                                                 
c                                                                               
         if( matchs(',',1) ) go to 718                                          
c                                                                               
c                       if there is an end of card, matrix input has            
c                       ended. branch to store the matrix globally.             
c                       if not, ignore the current entity and search            
c                       for another row to input. check validity of             
c                       3x3 rotation matrix.                                    
c                                                                               
         if( endcrd(dum) ) then                                                 
          rlen1 = sqrt( trans(1,1)**2 + trans(1,2)**2 +trans(1,3)**2 )          
          rlen2 = sqrt( trans(2,1)**2 + trans(2,2)**2 +trans(2,3)**2 )          
          rlen3 = sqrt( trans(3,1)**2 + trans(3,2)**2 +trans(3,3)**2 )          
          rflag1 = abs(rlen1-one) .gt. rottol .or.                              
     &             abs(rlen2-one) .gt. rottol .or.                              
     &             abs(rlen3-one) .gt. rottol                                   
          t11 = trans(1,1)**2 + trans(2,1)**2 + trans(3,1)**2                   
          t12 = trans(1,1)*trans(1,2) + trans(2,1)*trans(2,2) +                 
     &          trans(3,1)*trans(3,2)                                           
          t13 = trans(1,1)*trans(1,3) + trans(2,1)*trans(2,3) +                 
     &          trans(3,1)*trans(3,3)                                           
          t21 = trans(1,1)*trans(1,2) + trans(2,1)*trans(2,2) +                 
     &          trans(3,1)*trans(3,2)                                           
          t22 = trans(1,2)**2 + trans(2,2)**2 + trans(3,2)**2                   
          t23 = trans(1,2)*trans(1,3) + trans(2,2)*trans(2,3) +                 
     &          trans(3,2)*trans(3,3)                                           
          t31 = trans(1,1)*trans(1,3) + trans(2,1)*trans(2,3) +                 
     &          trans(3,1)*trans(3,3)                                           
          t32 = trans(1,2)*trans(1,3) + trans(2,2)*trans(2,3) +                 
     &          trans(3,2)*trans(3,3)                                           
          t33 = trans(1,3)**2 + trans(2,3)**2 + trans(3,3)**2                   
          rflag2 = abs(t11-one) .gt. rottol .or.                                
     &             abs(t22-one) .gt. rottol .or.                                
     &             abs(t33-one) .gt. rottol .or.                                
     &             abs(t12) .gt. rottol     .or.                                
     &             abs(t13) .gt. rottol     .or.                                
     &             abs(t21) .gt. rottol     .or.                                
     &             abs(t23) .gt. rottol     .or.                                
     &             abs(t31) .gt. rottol     .or.                                
     &             abs(t32) .gt. rottol                                         
          if ( rflag1 .or. rflag2 ) then                                        
                  call errmsg( 266, param, dums, dumr, dumd )                   
          end if                                                                
          go to 730                                                             
         else                                                                   
            call errmsg(84,dum,dums,dumr,dumd)                                  
            if( true(dum) ) go to 717                                           
         end if                                                                 
c                                                                               
c                       if there is a comma at the end of a line, the           
c                       input line is continued.                                
c                                                                               
 718     continue                                                               
         if( endcrd(dum) ) then                                                 
            call readsc                                                         
         end if                                                                 
         go to 717                                                              
c                                                                               
c                       check to make sure that the current row has             
c                       not been input twice. set input flag if not.            
c                                                                               
 719     if( inpflg(row) ) then                                                 
            param = row                                                         
            call errmsg(135,param,dums,dumr,dumd)                               
            go to 717                                                           
         else                                                                   
            inpflg(row) = .true.                                                
         end if                                                                 
c                                                                               
c                       store the row input.                                    
c                                                                               
         col= 0                                                                 
 720     if( .not. numd(tval) ) then                                            
            go to 717                                                           
         else                                                                   
            col = col+1                                                         
            if( col .gt. mxndof ) then                                          
               col = mxndof                                                     
               write(out,9073) 'columns ', node, 'columns '                     
               num_error = num_error + 1                                        
               if( true(dum) ) go to 717                                        
            end if                                                              
            trans(row,col) = tval                                               
         end if                                                                 
         go to 720                                                              
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     input of nodal constraint values, in con-      *        
c *                     straint compatable global coordinates.         *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      else                                                                      
c                                                                               
c                       initialize the temporary constraint vector.             
c                                                                               
         do i = 1,mxndof                                                        
            convec(i) = d32460                                                  
            inpflg(i) = .false.                                                 
         end do                                                                 
c                                                                               
c                       branch on dof type input.                               
c                                                                               
 726     if( matchs('u',1) ) then                                               
            dofn= 1                                                             
            curtyp= 'u'                                                         
            go to 728                                                           
         end if                                                                 
c                                                                               
         if( matchs('v',1) ) then                                               
            dofn= 2                                                             
            curtyp= 'v'                                                         
            go to 728                                                           
         end if                                                                 
c                                                                               
         if( matchs('w',1) ) then                                               
            dofn= 3                                                             
            curtyp= 'w'                                                         
            go to 728                                                           
         end if                                                                 
c                                                                               
         if( matchs(',',1) ) go to 727                                          
c                                                                               
c                       if there is an end of card, constr. input has           
c                       ended. branch to store the temp. vec. globally.         
c                       if not, ignore the current entity and search            
c                       for another dof to input.                               
c                                                                               
         if(endcrd(dum)) then                                                   
            go to 730                                                           
         else                                                                   
            call errmsg(72,dum,dums,dumr,dumd)                                  
            if(true(dum)) go to 726                                             
         end if                                                                 
c                                                                               
c                       if there is a comma at the end of a line, the           
c                       input line is continued.                                
c                                                                               
 727     continue                                                               
         if(endcrd(dum)) then                                                   
            call readsc                                                         
         end if                                                                 
         go to 726                                                              
c                                                                               
c                       store the constraints in the temporary vector           
c                                                                               
 728     if(matchs('=',1)) call splunj                                          
         if(.not.numd(cval)) then                                               
            call errmsg(73,dum,dums,dumr,dumd)                                  
         else                                                                   
c                                                                               
c                       make sure that the current dof has not                  
c                       already been input on the same nodal con-               
c                       straint command.                                        
c                                                                               
            if( inpflg(dofn) ) then                                             
               call errmsg(109,dum,curtyp,dumr,dumd)                            
               go to 726                                                        
            end if                                                              
            defcon(dofn) = .true.                                               
            inpflg(dofn) = .true.                                               
            convec(dofn) = cval                                                 
         end if                                                                 
         go to 726                                                              
      end if                                                                    
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     store constraints or transformation matrices   *        
c *                     globally.                                      *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 730  icn    = 0                                                                
      iplist = 1                                                                
 731  call trxlst(intlst,lenlst,iplist,icn,node)                                
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
         if( node .gt. nonode ) then                                            
            param= node                                                         
            call errmsg(16,param,dums,dumr,dumd)                                
            go to 734                                                           
         end if                                                                 
c                                                                               
c                       check that the list node is not negative.               
c                                                                               
         if( node .lt. 0 )  then                                                
            param= node                                                         
            call errmsg(58,param,dums,dumr,dumd)                                
            go to 734                                                           
         end if                                                                 
c                                                                               
         felem = inverse_incidences(node)%element_list(1)                       
         type  = iprops(1,felem)                                                
         ndof  = iprops(4,felem)                                                
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     store transformation matrix globally           *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
         if( skew ) then                                                        
c                                                                               
c                       set flag that there is a trn. mat. for the node.        
c                       if needed, allocate matrix for this node.               
c                       if the node already has had a trn. mat. input,          
c                       then overwrite the previous matrix.                     
c                                                                               
            if( .not. trn(node) ) then                                          
               trn(node) = .true.                                               
               call allo_trnmat ( node, 1, dum )                                
            endif                                                               
c                                                                               
c                       store the transformation matrix globally                
c                       for the node.                                           
c                                                                               
            do row = 1, 3                                                       
               trnmat(node)%mat(row,1) = trans(row,1)                           
               trnmat(node)%mat(row,2) = trans(row,2)                           
               trnmat(node)%mat(row,3) = trans(row,3)                           
            end do                                                              
c                                                                               
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     store the constraint vector globally           *        
c *                                                                    *        
c *                     ndof  -   number of dof per node               *        
c *                                                                    *        
c *                                                                    *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
         else                                                                   
c                                                                               
            do i = 1, ndof                                                      
               if( convec(i).eq.d32460 ) cycle                                  
               dof = dstmap(node)+i-1                                           
c                                                                               
c                       make sure this dof hasn't been previously               
c                       constrained.                                            
c                                                                               
               if( cnstrn_in(dof) .ne. d32460 ) then                            
                  param= dof                                                    
                  call errmsg(136,param,dums,dumr,dumd)                         
                  cycle                                                         
               end if                                                           
c                                                                               
c                       place the dof in the constraint mapping data            
c                       structure.                                              
c                                                                               
               if( csthed .eq. -1 ) then                                        
                  csthed = dof                                                  
               else                                                             
                  cstmap(olddof) = dof                                          
               end if                                                           
               olddof         = dof                                             
               cnstrn_in(dof) = convec(i)                                       
            end do                                                              
c                                                                               
         end if                                                                 
c                                                                               
 734  if( iplist .ne. 0 ) go to 731                                             
c                                                                               
      go to 710                                                                 
c                                                                               
c                                                                               
c **********************************************************************        
c **********************************************************************        
c                                                                               
c                                                                               
 9999 sbflg1         = .true.                                                   
      sbflg2         = .true.                                                   
      cstmap(olddof) = -1                                                       
      csttail        = olddof                                                   
c                                                                               
c                  user must specify at least 1 u, v, w                         
c                  constraint, else warning message                             
c                                                                               
      do i = 1, 3                                                               
         if ( .not. defcon(i) ) then                                            
            call errmsg(186,dum,dums,dumr,dumd)                                 
            exit                                                                
         end if                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9073 format(/1x,'>>>>> error: the number of ',a8,' input for a row of',        
     &           ' the transformation'/14x,'matrix of node ',i7,                
     &           ' exceeds the maximum number of'/14x,'degrees of',             
     &           ' freedom allowed any node. the number of'/14x,a8,             
     &           ' stored will be the above maximum.'/)                         
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine inconplane                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/3/2014                  *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *     constraints imposed on a specified plane                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine inconplane( olddof, defcon )                                   
      use global_data ! old common.main
c                                                                               
      use main_data, only : cnstrn, cnstrn_in, crdmap                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       parameters                                              
c                                                                               
      logical defcon(*)                                                         
c                                                                               
c                                                                               
c                       locally allocated                                       
c                                                                               
      double precision                                                          
     & convec(mxndof), tval, cval, dumd, zero, proximity_distance,              
     & one, d32460, xmin, xmax, ymin, ymax, zmin, zmax, x, y, z,                
     & coordtol(3), plane_coord, ctol, coordvec(3)                              
      real dumr                                                                 
      character(len=1) :: dums                                                  
      logical inpflg(mxndof), verify                                            
      logical matchs, numd, local_debug, endcrd                                 
      data zero, d32460 / 0.0d00, 32460.0 /                                     
      data ctol, local_debug / 0.00001, .false. /                               
c                                                                               
c                  find min and max coordinates for the body                    
c                                                                               
      xmin =  1.0e30                                                            
      xmax = -1.0e30                                                            
      ymin =  1.0e30                                                            
      ymax = -1.0e30                                                            
      zmin =  1.0e30                                                            
      zmax = -1.0e30                                                            
c                                                                               
      do node = 1, nonode                                                       
        x    = c(crdmap(node))                                                  
        y    = c(crdmap(node)+1)                                                
        z    = c(crdmap(node)+2)                                                
        xmax = max( xmax, x )                                                   
        xmin = min( xmin, x )                                                   
        ymax = max( ymax, y )                                                   
        ymin = min( ymin, y )                                                   
        zmax = max( zmax, z )                                                   
        zmin = min( zmin, z )                                                   
      end do                                                                    
c                                                                               
      coordtol(1) = abs( xmax-xmin ) * ctol                                     
      coordtol(2) = abs( ymax-ymin ) * ctol                                     
      coordtol(3) = abs( zmax-zmin ) * ctol                                     
c                                                                               
      if ( local_debug ) then                                                   
        write(*,*) '>> inside inconplane:'                                      
        write(*,*) ' xmin,max: ',xmin,xmax                                      
        write(*,*) ' ymin,max: ',ymin,ymax                                      
        write(*,*) ' zmin,max: ',zmin,zmax                                      
        write(*,*) ' xtol, ytol, ztol: ', coordtol                              
      end if                                                                    
c                                                                               
c                                                                               
c                  which plane to constrain                                     
c                                                                               
      plane = 0                                                                 
      if ( matchs('x',1) ) plane = 1                                            
      if ( matchs('y',1) ) plane = 2                                            
      if ( matchs('z',1) ) plane = 3                                            
      if ( plane .eq. 0 ) then                                                  
        call errmsg2( 4, dum, dums, dumr, dumd )                                
        return                                                                  
      end if                                                                    
c                                                                               
c                  coordinate of plane is zero by default                       
c                                                                               
      plane_coord = zero                                                        
      if ( matchs('=',1) ) then                                                 
       if ( .not. numd(plane_coord) ) then                                      
         call errmsg2( 5, dum, dums, dumr, dumd )                               
         return                                                                 
       end if                                                                   
      end if                                                                    
c                                                                               
      if ( local_debug ) then                                                   
        write(*,*) '>> plane id, coord: ',plane,plane_coord                     
      end if                                                                    
c                                                                               
      if( matchs('proximity',4) ) then                                          
        if( numd( proximity_distance ) ) then                                   
          coordtol(1:3) =  proximity_distance                                   
        end if                                                                  
      end if                                                                    
      write(out,9100) coordtol(1:3)                                             
                                                                                
c                                                                               
c                  get the type of constraints to impose on nodes               
c                  that lie on the plane, store values, return                  
c                                                                               
c                  fixed, symmetry or list of u, v, w values                    
c                                                                               
      convec(1:3) = d32460                                                      
      inpflg(1:3) = .false.                                                     
      verify = .false.                                                          
      if ( matchs('verify',3) ) verify = .true.                                 
                                                                                
c                                                                               
      if ( matchs('fixed',3) ) then                                             
        if ( matchs('verify',3) ) verify = .true.                               
        convec(1:3) = zero                                                      
        inpflg(1:3) = .true.                                                    
        call inconplane_store                                                   
        return                                                                  
      end if                                                                    
c                                                                               
      if ( matchs('symmetry',3) ) then                                          
        if ( matchs('verify',3) ) verify = .true.                               
        convec(plane) = zero                                                    
        inpflg(plane) = .true.                                                  
        call inconplane_store                                                   
        return                                                                  
      end if                                                                    
c                                                                               
      do   !  u, v, w values                                                    
        if ( matchs('verify',3) ) then                                          
           verify = .true.                                                      
           cycle                                                                
        end if                                                                  
        if( endcrd( dummy ) ) then                                              
          call inconplane_store                                                 
          return                                                                
        end if                                                                  
        if ( matchs('u',1) ) then                                               
          dofn = 1                                                              
        elseif ( matchs('v',1) ) then                                           
          dofn = 2                                                              
        elseif ( matchs('w',1) ) then                                           
          dofn = 3                                                              
        else                                                                    
          call errmsg2( 6, dum, dums, dumr, dumd )                              
          return                                                                
        end if                                                                  
c                                                                               
        inpflg(dofn) = .true.                                                   
        if ( matchs('=',1) ) call splunj                                        
        if ( .not. numd(convec(dofn)) ) then                                    
          call errmsg2( 7, dum, dums, dumr, dumd )                              
          return                                                                
        end if                                                                  
      end do                                                                    
c                                                                               
 9100 format(10x,">>>> Distances (x,y,z) for proximity test: ",                 
     &    3f15.6)                                                               
                                                                                
      contains  ! avoids go tos in above code.                                  
c     ========                                                                  
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine inconplane_store                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/3/2014                  *          
c     *                                                              *          
c     *     store constraints imposed by plane command               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine inconplane_store                                               
c                                                                               
c                  search for nodes that lie on the specified plane.            
c                  impose constraints on nodes. output verification list        
c                  if requested. we don't warn here about constraints           
c                  imposed on nodes that already have constraints.              
c                                                                               
c                                                                               
      if ( local_debug ) then                                                   
        write(*,*) '>> constraints to impose on matching nodes'                 
        write(*,*) ' convec: ', convec                                          
        write(*,*) ' inpflg: ', inpflg                                          
      end if                                                                    
      if ( verify ) write(out,9000)                                             
      con_node_count = 0                                                        
c                                                                               
      do node = 1, nonode                                                       
        coordvec(1) = c(crdmap(node))                                           
        coordvec(2) = c(crdmap(node)+1)                                         
        coordvec(3) = c(crdmap(node)+2)                                         
        if ( abs( plane_coord -coordvec(plane) ) .gt.                           
     &        coordtol(plane) ) cycle                                           
        con_node_count = con_node_count + 1                                     
        do i = 1, 3                                                             
          if( convec(i) .eq. d32460 ) cycle                                     
          dof = dstmap(node) + i -1                                             
          if( csthed .eq. -1 ) then                                             
             csthed = dof                                                       
          else                                                                  
             cstmap(olddof) = dof                                               
          end if                                                                
          olddof         = dof                                                  
          cnstrn_in(dof) = convec(i)                                            
          defcon(i)      = .true.                                               
          if ( verify ) then                                                    
            write(out,9100) node,i,convec(i)                                    
          end if                                                                
        end do                                                                  
      end do                                                                    
      write(out,9200) con_node_count                                            
c                                                                               
      return                                                                    
c                                                                               
 9000 format(//,1x,                                                             
     &'**** Constraints imposed by just entered plane command',/,/,             
     & 5x,' node   dof   constraint value')                                     
 9100 format(5x,i5,3x,i1,8x,f12.7)                                              
 9200 format(/,'>> Plane command applied constraints to: ',i5,                  
     & ' nodes...',//)                                                          
c                                                                               
      end subroutine inconplane_store                                           
c                                                                               
      end subroutine inconplane                                                 
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine incon_mpcs                     *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 08/22/2017 rhd             *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *     multi-point constraint equations                         *          
c     *     - update to following fix so leading term can have a +,- *          
c     *     - major error fix 7/29/11 to make eqn arrangement match  *          
c     *       that for tied contact. Byron implemented enforcement   *          
c     *       of mpcs assuming they were in form of tied contact     *          
c     *       arrangement. see notes below                           *          
c     *     - fixed several errors when MPCs are changed during      *          
c     *       an analysis                                                       
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine incon_mpcs                                                     
      use global_data ! old common.main
      use mod_mpc, only : mpcs_exist, num_user_mpc, user_mpc_table              
      use main_data, only : modified_mpcs                                       
c                                                                               
c              locals                                                           
c                                                                               
      integer, parameter :: max_trm = 100 ! limit on terms per mpc eqn          
c                                                                               
c              scan support                                                     
c                                                                               
      integer :: dumi                                                           
      real    :: dumr                                                           
      double precision ::  dumd                                                 
      character(len=1) :: dums                                                  
      logical, external :: matchs, realn, numr, integr, endcrd, true            
c                                                                               
c              local storage for translating, checking                          
c                                                                               
      real :: pre_term_multiplier, now_term_coeff, eqn_rhs_constant,            
     &        rtemp                                                             
      integer, dimension(:) :: eqn_nodes(max_trm), eqn_dofs(max_trm)            
      integer :: now_dof, ierr                                                  
      logical found_mpc                                                         
      real, dimension(:) :: eqn_coeffs(max_trm)                                 
c                                                                               
c              Set a flag so that we know later on the constraints have         
c              been touched                                                     
c                                                                               
      modified_mpcs = .true.                                                    
c                                                                               
c              notify user we're discarding any existing mpcs                   
c              in the mpc module. allocate permanent data structure             
c              for user mpcs in the module (see comments in mod_mpc)
c              F2003 - also deallocates allocated sub-objects.             
c                                                                               
      if( mpcs_exist ) then                                                     
         call errmsg2( 59, dumi, dums, dumr, dumd )                             
         if( allocated(user_mpc_table) ) deallocate( user_mpc_table )           
         num_user_mpc = 0                                                       
      end if                                                                    
c                                                                               
      allocate( user_mpc_table(max_mpc), stat=ierr )                            
      if( ierr .ne. 0 ) then                                                    
         call errmsg2( 47, dumi, dums, dumr, dumd )                             
         call die_abort                                                         
      end if                                                                    
c                                                                               
      mpcs_exist = .true.                                                       
      nmpc  = 0                                                                 
c                                                                               
c               loop to process next mpc equation. line must start              
c               with an +, - or integer. if not, treat as end of                
c               mpc equations                                                   
c                                                                               
      do    ! over mpc equations                                                
         call readsc                                                            
         found_mpc = .false.                                                    
         if( matchs('+',1) ) then                                               
            found_mpc = .true.                                                  
         elseif( matchs('-',1) ) then                                           
            found_mpc = .true.                                                  
         elseif( integr(idummy) ) then                                          
            found_mpc = .true.                                                  
         else                                                                   
            found_mpc = .false.                                                 
         end if                                                                 
c                                                                               
         call reset                                                             
         if( true(dum) )  call splunj                                           
         if( .not. found_mpc ) return                                           
c                                                                               
c              loop across line to extract each node number, its                
c              dof and multipler, the (+,-) between terms and the               
c              rhs const. for now, the constant = 0.0. when done.               
c              store the equation terms in the module data structure.           
c              loop over multiple lines continued by commas.                    
c                                                                               
c              for curremt mpc implementation, the rhs                          
c              constant must = 0.0                                              
c                                                                               
c              since we have done a line reset, the leading term may            
c              now have a pre-multiplier ( '+' or '-' ) just like other         
c              terms                                                            
c                                                                               
         nterm          = 0                                                     
         now_node       = 0                                                     
         now_dof        = 0                                                     
         now_term_coeff = 0.0                                                   
         nmpc = nmpc + 1                                                        
         if( nmpc .gt. max_mpc ) call  incon_mpcs_resize                                              
c                                                                               
         do ! all terms for this mpc eqn                                        
c                                                                               
            if( matchs('=',1) ) then                                            
              eqn_rhs_constant = 0.0                                            
              if( numr(rtemp) ) then                                            
                 eqn_rhs_constant = rtemp                                       
                 if( rtemp .ne. 0.0 ) then                                      
                  call errmsg( 320, 1, dums, dumr, dumd )                       
                  call incon_flushline                                          
                  nmpc = nmpc - 1                                               
                  exit ! from processing this mpc eqn                           
                 end if                                                         
              end if                                                            
              call incon_mpcs_store( nterm, eqn_rhs_constant,                   
     &                     eqn_nodes, eqn_dofs, eqn_coeffs )                    
               exit  !  start at loop for new mpc eqn                           
            end if                                                              
c                                                                               
            nterm = nterm + 1                                                   
c                                                                               
            if( nterm .gt. max_trm ) then                                       
               call errmsg2( 43, max_trm, dums, dumr, dumd )                    
               call die_abort                                                   
            end if                                                              
c                                                                               
            pre_term_multiplier = 1.0                                           
            if( matchs('+',1) ) then                                            
                  pre_term_multiplier = 1.0                                     
               else if (matchs('-',1) ) then                                    
                  pre_term_multiplier = -1.0                                    
            end if                                                              
c                                                                               
c                     node number                                               
c                                                                               
            if( integr(now_node) ) then                                         
               if( now_node .gt. nonode .or. now_node .le. 0  ) then            
                  call errmsg2( 34, now_node, dums, dumr, dumd )                
                  call reset                                                    
                  if( true(dum) ) call splunj                                   
                  call incon_flushline                                          
                  nmpc = nmpc - 1                                               
                  exit ! from processing this mpc eqn                           
               end if                                                           
               eqn_nodes(nterm) = now_node                                      
            end if                                                              
c                                                                               
c                    coefficient after node number. must be                     
c                    present and a real number.                                 
c                                                                               
            if( realn(now_term_coeff) ) then                                    
             eqn_coeffs(nterm) = pre_term_multiplier * now_term_coeff           
            else                                                                
             call errmsg( 319, 1, dums, dumr, dumd )                            
             call reset                                                         
             if( true(dum) ) call splunj                                        
             call incon_flushline                                               
             nmpc = nmpc - 1                                                    
             exit ! from processing this mpc eqn                                
            end if                                                              
c                                                                               
c                    get dof: u, v, w                                           
c                                                                               
            now_dof = 0                                                         
            if( matchs('u',1) ) then                                            
                  now_dof = 1                                                   
               else if( matchs('v',1) ) then                                    
                  now_dof = 2                                                   
               else if( matchs('w',1) ) then                                    
                  now_dof = 3                                                   
            end if                                                              
            if( now_dof. eq. 0 ) then                                           
             call errmsg( 210, 1, dums, dumr, dumd )                            
             call reset                                                         
             if( true(dum) ) call splunj                                        
             call incon_flushline                                               
             nmpc = nmpc - 1                                                    
             exit ! from processing this mpc eqn                                
            end if                                                              
            eqn_dofs(nterm) = now_dof                                           
c                                                                               
c                comma followed by eol is continuation of current mpc           
c                                                                               
            if( matchs(',',1) ) then                                            
               if( endcrd(dum) ) then                                           
                  call readsc                                                   
                  cycle  ! top of loop looking for terms in this mpc            
               end if                                                           
            end if                                                              
c                                                                               
         end do ! loop looking for terms in an mpc                              
      end do  ! loop to process all mpc eqns                                    
c                                                                               
      return                                                                    
c                                                                               
c    Notes on storage and implementation of MPC equations.                      
c                                                                               
c      (a) the equations must be homogeneous (rhs = 0)                          
c                                                                               
c      (b) none of the mentioned dof can also have a zero                       
c          or non-zero (absolute) constraint applied (which                     
c          effectively makes the equation non-homogeneous.                      
c                                                                               
c      (c) equations must be normalized so that the multiplier                  
c          stored for the leading term is -1.0                                  
c                                                                               
c    The tied-contact capability was implemented before user                    
c    defined mpcs. The tied contact arrangement of equations has                
c    the form (for example):                                                    
c                                                                               
c       ux = 0.25 ub + 0.3 uc + 0.4uc + 0.05 ud                                 
c                                                                               
c    Here node 'x' is located on an element face which has corner nodes         
c    a, b, c, d. The u-displ at 'x' is then a linear combination of             
c    displacements of the 4 face nodes with coefficients determined by          
c    the geometric location of node 'x' on the face.                            
c                                                                               
c    The tied contact input system stores the above equation in the form:       
c                                                                               
c     -1.0 ux + 0.25 ub + 0.3 uc + 0.4uc + 0.05 ud= 0.0                         
c                                                                               
c    Note the leading coefficient is -1.0. Byron implemented the above          
c    mpc assuming the leading term is -1.0 and that ux is the dependent         
c    dof                                                                        
c                                                                               
c    Similarly, for two coincident nodes connected by tied contact say          
c    20 and 50, the tied contact input system will store:                       
c                                                                               
c      -1.0 u20 + 1.0 u50 = 0.0                                                 
c                                                                               
c    Our user-defined mpc's must be transformed to follow the exact             
c    storage for tied contact since Byron's code to enforce the mpcs            
c    is based on the tied contact data structure.                               
c                                                                               
c    A user mpc of the form:                                                    
c                                                                               
c        53 1.0 u - 43 1.0 u = 0                                                
c                                                                               
c    to make the two u displacements equal must be stored in the form:          
c                                                                               
c         -1.0 u53 + 1.0 u43 = 0                                                
c                                                                               
c    The user can write this same equation as                                   
c                                                                               
c      - 53 1.0 u + 43 1.0 u = 0                                                
c      - 53 -1.0 u - 43 1.0 u = 0                                               
c      - 53 1.0 u - 43 -1.0 u = 0                                               
c       53 -42.5 + 43 u 42.5 = 0                                                
c                                                                               
c    All forms must be stored as                                                
c                                                                               
c        -1.0 u53 + 1.0 u43 = 0.0                                               
c                                                                               
c    More complex example. The user wants the v-dispacement on all 4            
c    nodes of an element face to sum = 0. For face nodes 23,                    
c    99, 103, 42 the user might write as input:                                 
c                                                                               
c            23 1.0 v + 99 1.0 v + 103 1.0 v + 1.0 42 v = 0                     
c                                                                               
c    We must store this equation as:                                            
c                                                                               
c            -1.0 u23  - 1.0 v99 - 1.0 v103 - 1.0 v42 = 0.0                     
c                                                                               
c    The code scans for the pre-multiplier sign (+,-), the                      
c    node number, the coefficient and the dof.                                  
c    The node number, dof and the multiplier = pre-multiplier *                 
c    the coefficient are stored in simple vectors while scanning.               
c                                                                               
c    While storing into the global MPC data structure that Brian's              
c    code uses, scale the coefficients so that the leading term                 
c    multiplier is -1.0.                                                        
c                                                                               
      end  
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine incon_mpcs_resize                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 8/22/2017 rhd              *          
c     *                                                              *          
c     *     increase size of the user_mpc_table. make larger one,    *
c     *     deep copy old -> new, move allocation. release old table *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      
      subroutine incon_mpcs_resize
      use global_data ! old common.main
      use mod_mpc, only : user_mpc_table, mpc_eqn                
      implicit none                                                     
c
      integer :: old_mpc_size, i, j, nt
      logical, parameter :: local_debug = .false.
      type (mpc_eqn), allocatable, dimension (:) :: new_mpc_table
c     
      if( local_debug ) then
         write(*,*) '...  resizing mpc user table.....'
         write(*,*) '       old_max_mpc: ', max_mpc 
      end if
c          
      old_mpc_size = max_mpc
      max_mpc = 2 * old_mpc_size
      allocate( new_mpc_table(max_mpc) )
c     
c              deep copy as required
c 
      do i = 1, old_mpc_size
          nt = user_mpc_table(i)%num_terms
          new_mpc_table(i)%num_terms = nt
          new_mpc_table(i)%constant = user_mpc_table(i)%constant
          allocate( new_mpc_table(i)%node_list(nt),
     &              new_mpc_table(i)%dof_list(nt),
     &              new_mpc_table(i)%multiplier_list(nt)  )
          do j = 1, nt
           new_mpc_table(i)%node_list(j) = 
     &                   user_mpc_table(i)%node_list(j)
           new_mpc_table(i)%dof_list(j) = 
     &                   user_mpc_table(i)%dof_list(j)
           new_mpc_table(i)%multiplier_list(j) = 
     &                   user_mpc_table(i)%multiplier_list(j)
          end do ! on j
      end do
c      
c              initialize remainder of new table. probably not req'd
c 
      do i =  old_mpc_size+1, max_mpc
          new_mpc_table(i)%num_terms = 0
          new_mpc_table(i)%constant = 0.0
          new_mpc_table(i)%node_list => null()
          new_mpc_table(i)%dof_list  => null()
          new_mpc_table(i)%multiplier_list  => null()
      end do ! on i
c
c              release old table and move allocation. F2003 &
c              later does a deep release on derived types. move_alloc
c              does a deep release on new_mpc_table
c
      deallocate( user_mpc_table )
      call move_alloc( new_mpc_table, user_mpc_table )
c
      return
      end          
c
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine incon_mpcs_store                 *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 8/22/2017 rhd              *          
c     *                                                              *          
c     *     stores a user-definite MPC equation into global data     *          
c     *     data structure                                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine incon_mpcs_store( nterm, const, node, dof, multi )             
      use mod_mpc, only : num_user_mpc, user_mpc_table     
      use global_data, only : max_mpc                     
      include 'param_def'                                                       
c                                                                               
c              parameters                                                       
c                                                                               
      integer :: nterm, node(*), dof(*)                                         
      real ::    const, multi(*)                                                
c                                                                               
c              locals                                                           
c                                                                               
      integer ::  err, dumi                                                     
      real ::     dumr, factor                                                  
      double precision :: dumd                                                  
      character(len=1) :: dums                                                  
c                                                                               
      num_user_mpc = num_user_mpc + 1                                           
      if( num_user_mpc .gt. max_mpc )  call incon_mpcs_resize
c
      user_mpc_table(num_user_mpc)%num_terms = nterm                            
      user_mpc_table(num_user_mpc)%constant  = const                            
c                                                                               
      allocate( user_mpc_table(num_user_mpc)%node_list(nterm),                  
     &          user_mpc_table(num_user_mpc)%dof_list(nterm),                   
     &          user_mpc_table(num_user_mpc)%multiplier_list(nterm),            
     &          stat=err )                                                      
c                                                                               
      if( err .ne. 0 ) then                                                     
         call errmsg2( 47,dumi,dums,dumr,dumd )                                 
         call die_abort                                                         
      end if                                                                    
c                                                                               
c              normalize the multiplier on each dof in the                      
c              equations such that the leading term is -1.0                     
c                                                                               
      factor = 1.0 / abs( multi(1) )                                            
      if( multi(1) .gt. 0.0 ) factor = -1.0 * factor                            
c                                                                               
      do i = 1, nterm                                                           
         user_mpc_table(num_user_mpc)%node_list(i)       = node(i)              
         user_mpc_table(num_user_mpc)%dof_list(i)        = dof(i)               
         user_mpc_table(num_user_mpc)%multiplier_list(i) =                      
     &         multi(i)*factor                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine incon_flushline               *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                   last modified : 02/24/03                   *          
c     *                                                              *          
c     *     this subroutine skips a logical line of the input file   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine incon_flushline                                                
      integer  dum                                                              
      logical  matchs, realn, integr, endcrd, true                              
c                                                                               
      call readsc                                                               
      do                                                                        
         if (matchs(',',1)) then                                                
            if (endcrd(dum)) then                                               
               call readsc                                                      
               cycle                                                            
            end if                                                              
         end if                                                                 
         if (.not.matchs(',',1)) then                                           
            if (endcrd(dum)) return                                             
            if (true(dum)) then                                                 
               call splunj                                                      
               cycle                                                            
            end if                                                              
         end if                                                                 
      end do                                                                    
c                                                                               
      end                                                                       
