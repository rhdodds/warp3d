c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oudva                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 8/4/2014 rhd               *          
c     *                                                              *          
c     *   nodal output of displacements, velocities, accelerations,  *          
c     *   reactions and temperatures                                 *          
c     *                                                              *          
c     *   (1) patran file, (2) flat file, or (3) printed values for  *          
c     *   node numbers lists and options specified by the user.      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oudva( dva, ouflg, oupat, oubin, ouasc, wide, eform,           
     &                  prec, noheader, react_totals_only ,                     
     &                  out_packet_now, flat_file, stream_file,                 
     &                  text_file, compressed  )                                
      use global_data ! old common.main
c                                                                               
      implicit none                                                             
c                                                                               
      integer :: dva                                                            
      logical :: ouflg, oupat, oubin, ouasc, wide, eform, prec,                 
     &           noheader, react_totals_only,                                   
     &           out_packet_now, flat_file, stream_file,                        
     &           text_file, compressed                                          
c                                                                               
      integer :: all, lsttyp, idum, lenlst, errnum, param                       
      integer :: intlst(mxlsz)                                                  
      real :: dumr                                                              
      double precision :: dumd                                                  
      character(len=1) :: dums                                                  
      logical :: first                                                          
      logical, external :: matchs, endcrd, true                                 
c                                                                               
      ouflg = .false. ! set .true. if bad outcome                               
c                                                                               
c                       output to patran or a flat file for post                
c                       processing. no list of nodes/elements checked.          
c                       return                                                  
c                                                                               
c                       write the displ/vel/acc/                                
c                       reactions/temps for all model nodes.                    
c                       Patran option is either a binary or asci file.          
c                       Flat file can be text or stream (binary)                
c                                                                               
      if( oupat .or. flat_file) then                                            
         call oupdva( dva, oubin, ouasc, flat_file, stream_file,                
     &                text_file, compressed )                                   
         return                                                                 
      end if                                                                    
c                                                                               
c                       hardcopy and/or packet file output. ignore              
c                       words "for" or "forces"                                 
c                                                                               
      if( matchs('forces',5) ) call splunj                                      
      if( matchs('for',3)    ) call splunj                                      
c                                                                               
c                      possible scan token at this point                        
c                                                                               
c                      eol -> process 'all' nodes and leave                     
c                                                                               
c                      "nodes" -> followed by <list>. if no list treat          
c                                as "all" and loop back to top here             
c                                                                               
c                      "elements" -> followed by <list>. if no list             
c                                treat as "all" and loop back here.             
c                                output element number then nodes on            
c                                the element and displacements.                 
c                                                                               
c                      <list> -> parse list as node list, process,              
c                                loop back here for another possible            
c                                keyword "nodes" or "elements"                  
c                                and another list to process                    
c                                                                               
c                                                                               
c                      not one of the above, back space 1 token and             
c                      leave. probably another quantity type.                   
c                                                                               
      all = nonode      !  assume output for user list of nodes                 
      lsttyp = 1        !  output for user list of nodes                        
                                                                                
      if( endcrd(idum) ) then                                                   
        lenlst = 3                                                              
        intlst(1) = 1; intlst(2) = -all; intlst(3) = 1                          
        call ouhdva( lsttyp, dva, wide, eform, prec, intlst, lenlst,            
     &             noheader, react_totals_only, out_packet_now )                
        return                                                                  
      end if                                                                    
                                                                                
      do ! "nodes" <lists> and "elements" <lists>                               
                                                                                
      lsttyp = 1; all = nonode                                                  
      if( matchs('nodes',3) ) then                                              
         call scan   ! get possible <list>                                      
      elseif( matchs('elements',4) ) then                                       
         lsttyp = 2; all = noelem; call scan ! <get possible list>              
      end if                                                                    
c                                                                               
c                       examine current token for <list> of nodes or            
c                       elements. if not a list, assume "all" of nodes          
c                       or elements -- scanner does not move.                   
c                                                                               
      call trlist( intlst, mxlsz, all, lenlst, errnum )                         
c                                                                               
c                       =1 <list> or "all" found. token after list is in        
c                          scanner                                              
c                       =2 parse rules failed                                   
c                       =3 list space would overflow                            
c                       =4 no <list> or "all" found                             
c                                                                               
      select case ( errnum )                                                    
      case( 1 )  ! no error. <lis> could be "all"                               
      case( 2 )  ! parse rules for list failed                                  
         param = 1                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         ouflg = .true.                                                         
         exit   ! do over lists                                                 
      case( 3 )  ! overflow of space provided for list                          
         param = 2                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         ouflg = .true.                                                         
         exit   ! do over lists                                                 
      case( 4 )                                                                 
         lenlst = 3                                                             
         intlst(1) = 1; intlst(2) = -all; intlst(3) = 1                         
      case default                                                              
         write(out,9000) 1                                                      
         call die_gracefully                                                    
      end select                                                                
c                                                                               
c                       there is a <list> (could be implied all).               
c                       hard copy or packet                                     
c                       output the dis/vel/acc/react/temp for the list.         
c                                                                               
      call ouhdva( lsttyp, dva, wide, eform, prec, intlst, lenlst,              
     &             noheader, react_totals_only, out_packet_now )                
c                                                                               
c                       backup 1 token and set next test function to            
c                       advance and examine. this puts scanner into a           
c                       known state. token can be: (1) eol,                     
c                       (2) word nodes/elements again with another              
c                       <list>, (3) another output quantity.                    
c                                                                               
      call backsp( 1 )                                                          
      if( true(idum) ) call splunj                                              
      if( endcrd(idum) ) return                                                 
c                                                                               
      if( matchs( 'nodes', 4 ) ) then                                           
         call backsp(1)                                                         
         cycle  ! to top of loop here                                           
      elseif( matchs( 'elements', 4 ) ) then                                    
         call backsp(1)                                                         
         cycle ! to top of loop here                                            
      else   ! should be next type of output quantity. let                      
         call backsp( 1 )                                                       
         if( true(idum) ) call splunj                                           
         exit ! do over node/elem/lists                                         
      end if                                                                    
c                                                                               
c                       should never fall thru here with above logic            
c                                                                               
      write(out,9000) 2                                                         
      call die_gracefully                                                       
c                                                                               
      end do   !  over nodes/elements lists.                                    
c                                                                               
c                                                                               
9000  format(/1x,'>>>>> internal error:  routine oudva @ ',i1,/,                
     & 14x,'job terminated...',/)                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
