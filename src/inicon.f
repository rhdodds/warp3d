c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine inicon                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 02/25/00                   *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of the *          
c     *     desired initial conditions for the structure at time 0.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine inicon( sbflg1, sbflg2 )                                       
      use global_data ! old common.main
c                                                                               
      use main_data, only : trn, trnmat, temper_nodes,                          
     &                      temper_nodes_ref, temperatures_ref,                 
     &                      inverse_incidences                                  
      implicit integer (a-z)                                                    
      logical sbflg1, sbflg2                                                    
c                                                                               
c                       local declarations                                      
c                                                                               
      real dumr                                                                 
      double precision                                                          
     &   cval, mpfact, edva(mxvl,mxndof), dumd, zero,                           
     &   trnmte(mxvl,mxedof,mxndof)                                             
      character name*80, iclnam*8, dums, curtyp*1                               
      logical found, dvaflg(mxndof), trne(mxvl,mxndel)                          
      logical matchs, endcrd, true, label, numd, scanms                         
      dimension intlst(mxlsz)                                                   
      data zero /0.0/                                                           
c                                                                               
c                       if sub flag 1 is on, there is reentry into              
c                       inicon after an error in the input of type              
c                       of initial condition input.                             
c                                                                               
      if( sbflg1 ) then                                                         
         call errmsg(119,dum,dums,dumr,dumd)                                    
      end if                                                                    
c                                                                               
 1100 call readsc                                                               
c                                                                               
c                       branch on the type of initial condition to              
c                       be input.                                               
c                                                                               
 1105 if( matchs('nodal',4)        ) call splunj                                
      if( matchs('displacement',7) ) go to 1110                                 
      if( matchs('velocity',3)     ) go to 1111                                 
      if( matchs('acceleration',3) ) go to 1112                                 
      if( matchs('temperature',4)  ) go to 1113                                 
c                                                                               
c                       if the above key words are not encountered,             
c                       return to driver subroutine and look for a              
c                       high level command.                                     
c                                                                               
      go to 9999                                                                
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     input of initial displacements.                *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 1110 continue                                                                  
      cond = 1                                                                  
      go to 1120                                                                
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     input of initial velocities.                   *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 1111 continue                                                                  
      cond = 2                                                                  
      go to 1120                                                                
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     input of initial accelerations                 *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 1112 continue                                                                  
      cond = 3                                                                  
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     input of initial (nodal) temperatures          *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 1113 continue                                                                  
      cond = 4                                                                  
c                                                                               
c                                                                               
c                                                                               
c                       read the list of nodes whose initial conditions         
c                       are to be input.                                        
c                                                                               
 1120 call readsc                                                               
      if( matchs('nodes',4) ) call scan                                         
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
c                       was found. in this case, input of initial cond-         
c                       itions has ceased.                                      
c                                                                               
      if( errnum .eq. 2 ) then                                                  
         param = 1                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 1120                                                             
      else if( errnum .eq. 3 ) then                                             
         param = 2                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 1120                                                             
      else if( errnum .eq. 4 ) then                                             
         go to 1105                                                             
      else                                                                      
         if( errnum .eq. 1 ) then                                               
            call backsp( 1 )                                                    
            if( true(dummy) ) go to 1125                                        
         end if                                                                 
         param = 3                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         go to 1120                                                             
      end if                                                                    
c                                                                               
c                       there is a valid nodal list. read the dof               
c                       to be assigned initial conditions.                      
c                                                                               
c                       initialize the d/v/a/ and dof type arrays.              
c                       use dofn = 0 for temperature                            
c                                                                               
 1125 do i = 1, mxndof                                                          
         edva(1,i) = zero                                                       
         dvaflg(i) = .false.                                                    
      end do                                                                    
c                                                                               
c                       input the dof type. if found, branch to                 
c                       store its initial value.                                
c                                                                               
 1127 if( matchs('u',1) ) then                                                  
         dofn   = 1                                                             
         curtyp = 'U'                                                           
         go to 1129                                                             
      end if                                                                    
c                                                                               
      if( matchs('v',1) ) then                                                  
         dofn   = 2                                                             
         curtyp = 'V'                                                           
         go to 1129                                                             
      end if                                                                    
c                                                                               
      if( matchs('w',1) ) then                                                  
         dofn   = 3                                                             
         curtyp = 'W'                                                           
         go to 1129                                                             
      end if                                                                    
c                                                                               
      if( matchs('temperature',4) ) then                                        
         dofn   = 0                                                             
         curtyp = 'T'                                                           
         go to 1129                                                             
      end if                                                                    
c                                                                               
      if( matchs(',',1) ) go to 1128                                            
c                                                                               
c                       if there is an end of card, in. cond. input has         
c                       ended. branch to store the temp. vec. globally.         
c                       if not, ignore the current entity and search            
c                       for another dof to input.                               
c                                                                               
      if( endcrd(dum) ) then                                                    
         go to 1135                                                             
      else                                                                      
         call errmsg( 72, dum, dums, dumr, dumd )                               
         if( true(dum) ) go to 1127                                             
      end if                                                                    
c                                                                               
c                       if there is a comma at the end of a line, the           
c                       input line is continued.                                
c                                                                               
 1128 continue                                                                  
      if( endcrd(dum) ) then                                                    
         call readsc                                                            
      end if                                                                    
      go to 1127                                                                
c                                                                               
c                       store the initial conditions.                           
c                                                                               
 1129 if( matchs('=',1) ) call splunj                                           
      if( .not. numd(cval) ) then                                               
         call errmsg( 73, dum, dums, dumr, dumd)                                
      else                                                                      
c                                                                               
c                       make sure that the current dof has not                  
c                       already been input on the same nodal initial            
c                       conditions command.                                     
c                                                                               
         if ( dofn .eq. 0 ) go to 1127                                          
         if( dvaflg(dofn) ) then                                                
            call errmsg( 117, dum, curtyp, dumr, dumd )                         
            go to 1127                                                          
         end if                                                                 
c                                                                               
         dvaflg(dofn) = .true.                                                  
         edva(1,dofn) = cval                                                    
c                                                                               
      end if                                                                    
c                                                                               
      go to 1127                                                                
c                                                                               
c                       store initial conditions globally.                      
c                                                                               
 1135 icn    = 0                                                                
      iplist = 1                                                                
 1136 call trxlst( intlst, lenlst, iplist, icn, node )                          
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
         if( node .gt. nonode ) then                                            
            param = node                                                        
            call errmsg( 16, param, dums, dumr, dumd )                          
            go to 1140                                                          
         end if                                                                 
c                                                                               
c                       check that the list node is not negative.               
c                                                                               
         if( node .lt. 0 ) then                                                 
            param = node                                                        
            call errmsg( 58, param, dums, dumr, dumd )                          
            go to 1140                                                          
         end if                                                                 
c                                                                               
         if ( dofn .eq. 0 ) then                                                
            temper_nodes(node) = cval                                           
            temper_nodes_ref(node) = cval                                       
            if ( abs(cval) .ne. zero ) temperatures_ref = .true.                
            go to 1140                                                          
         end if                                                                 
c                                                                               
         elem = inverse_incidences(node)%element_list(1)                        
         type = iprops(1,elem)                                                  
         ndof = iprops(4,elem)                                                  
c                                                                               
c                       as the displacements, velocities, and                   
c                       accelerations are input in uniform global               
c                       coordinates, they must be transformed to                
c                       constraint compatable global coordinates                
c                       before global storage.                                  
c                                                                               
c                       extract transformation matrix for this node             
c                       and transform to ccg coordinates.                       
c                                                                               
         trne(1,1) = trn(node)                                                  
         if( trne(1,1) ) then                                                   
            do i = 1, ndof                                                      
               do j = 1, ndof                                                   
                  trnmte(1,i,j) = trnmat(node)%mat(i,j)                         
               end do                                                           
            end do                                                              
            call trnvec( edva, trnmte, trne, ndof, 1, 1, 1 )                    
         end if                                                                 
c                                                                               
c                       store the transformed temporary vector                  
c                       globally, which contains the dis/vel/acc com-           
c                       patable with the current node in consecutive            
c                       order.                                                  
c                                                                               
         do i = 1, ndof                                                         
            if( cond .eq. 1 ) then                                              
               u(dstmap(node)+i-1) = edva(1,i)                                  
            else if( cond .eq. 2 ) then                                         
               v(dstmap(node)+i-1) = edva(1,i)                                  
            else                                                                
               a(dstmap(node)+i-1) = edva(1,i)                                  
            end if                                                              
         end do                                                                 
c                                                                               
c                                                                               
 1140 if( iplist .ne. 0 ) go to 1136                                            
c                                                                               
c                       return to process more lists.                           
c                                                                               
      go to 1120                                                                
c                                                                               
c                                                                               
c **********************************************************************        
c **********************************************************************        
c                                                                               
c                                                                               
 9999 sbflg1= .true.                                                            
      sbflg2= .false.                                                           
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
