c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine compute                      *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 5/18/2017 rhd              *          
c     *                                                              *          
c     *     scan the compute command, make some checks and call      *          
c     *     various driver routines                                  *          
c     *     -> this code runs only on root process                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine compute                                                        
      use global_data ! old common.main
      use j_data, only: comput_j, comput_i                                      
      use erflgs                                                                
c                                                                               
      implicit none                                                             
c                                                                               
c                                                                               
      integer :: intlst(mxlsz), dum, param, dummy, nc, lodn, ldnum,             
     &           lenlst, errnum                                                 
      real :: dumr                                                              
      double precision :: dumd                                                  
      double precision, parameter :: zero = 0.0d0                               
      character :: dums*1                                                       
      character(len=80) :: name                                                 
      character(len=8) :: ldname                                                
      logical :: matchs, endcrd, true, label, found, scanms                     
      logical, save :: notes_msg = .true.                                       
c                                                                               
c                       branch on the type of computation to be                 
c                       performed.                                              
c                                                                               
      comput_j = .false.; comput_i = .false.                                    
      if( matchs('domain',5) ) comput_j = .true.                                
      if( matchs('interaction',8) ) comput_i = .true.                           
      if( matchs('domain',5) ) comput_j = .true.                                
      if( comput_j .or. comput_i ) then                                         
       call didriv; go to 9999                                                  
      end if                                                                    
c                                                                               
      if( .not. matchs('displacements',5) ) then                                
        call errmsg(121,dum,dums,dumr,dumd); go to 9999                         
      end if                                                                    
c                                                                               
c                       computation of displacements for a specified            
c                       series of time steps.                                   
c                                                                               
c                       make sure that the constraints to be applied            
c                       for the specified time steps exist.                     
c                                                                               
      if( .not. constr ) then                                                   
        param = 6; call errmsg(11,param,dums,dumr,dumd); go to 9999             
      end if                                                                    
c                                                                               
c                       input the loading to be used to define the              
c                       specified time steps.                                   
c                                                                               
      if( matchs('for',3) ) call splunj                                         
      if( .not.matchs('loading',4) ) then                                       
        call errmsg(54,param,dums,dumr,dumd);  go to 9999                       
      end if                                                                    
c                                                                               
      if( .not. label(dummy) ) then                                             
        call errmsg(121,param,dums,dumr,dumd);  go to 9999                      
      end if                                                                    
c                                                                               
c                       if the loading name input cannot be found,              
c                       error skip command. otherwise                           
c                       search for the step loading in the loading              
c                       library.                                                
c                                                                               
      found = .false.                                                           
      name = ' '                                                                
      ldname = ' '                                                              
      call entits(name,nc)                                                      
      if( nc > 8) nc = 8                                                        
      ldname(1:nc) = name(1:nc)                                                 
      lodn = lodhed/two16                                                       
 1222 continue                                                                  
      if( lodn == 32460 ) go to 1223                                            
       if( scanms(lodnam(lodn),ldname,8) ) then                                 
         if( .not. scanms(lodtyp(lodn),'TIMESTEP',8) ) then                     
             param = lodn; call errmsg(122,param,dums,dumr,dumd)                
             go to 9999                                                         
         end if                                                                 
         ldnum = lodn; found = .true.; go to 1223                               
      end if                                                                    
      lodn = lodlst(lodn)/two16                                                 
      go to 1222                                                                
c                                                                               
 1223 continue                                                                  
      if( .not. found ) then                                                    
        call errmsg(60,dum,dums,dumr,dumd); go to 9999                          
      end if                                                                    
c                                                                               
c                       there is a valid loading. read the list of              
c                       time/load steps to be computed.                         
c                                                                               
      if( matchs('for',3) ) call splunj                                         
      if( .not. matchs('steps',4) ) then                                        
        call errmsg(121,param,dums,dumr,dumd);  go to 9999                      
      end if                                                                    
c                                                                               
      call scan                                                                 
      call trlist( intlst,mxlsz,nonode,lenlst,errnum )                          
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       a value of 4 indicates that no list was found.          
c                       in these last three cases, the illegal list             
c                       will be ignored and a new compute command will          
c                       be sought.                                              
c                                                                               
      if( errnum .ne. 1 ) then                                                  
          if( errnum  == 2) then                                                
               param = 1                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
          else if( errnum == 3 ) then                                           
               param = 2                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
          else if( errnum == 4) then                                            
               param = 4                                                        
               call errmsg(24,param,dums,dumr,dumd)                             
          end if                                                                
          go to 9999                                                            
      end if                                                                    
c                                                                               
c                       the list of time steps is a valid one. compute          
c                       the displacements for the steps specified.              
c                       run cehcks/setups that cannot be done until             
c                       now - most often for just 1st time step or              
c                       first step in a restart.                                
c                                                                               
      notes_msg = .false.                                                       
      call compute_checks                                                       
      call stpdrv( intlst, lenlst, ldnum)                                       
      return                                                                    
c                                                                               
 9999 continue                                                                  
      call scan_flushline                                                       
      return                                                                    
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine compute_checks                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/12/2017 rhd              *          
c     *                                                              *          
c     *   perform additional, generally one-time set ups that must   *          
c     *   done after input processed but not before every step or    *          
c     *   output (e.g. just restarted)                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine compute_checks                                                 
      use global_data ! old common.main
      use  main_data, only : cp_matls_present, matprp, imatprp                  
      implicit integer (a-z)                                                    
c                                                                               
      integer :: matnum                                                         
      logical :: found_cp                                                       
c                                                                               
c                                                                               
c              check for presence of material(s) using crystal                  
c              plasticity.                                                      
c                                                                               
c              cp_matls_present = -1 we've not yet checked                      
c                               =  0 we checked and no CP materials             
c                               =  1 we check. there are CP materials           
c                                  and we need to call CP setup                 
c                                  routine if not already done.                 
c                               for MPI, workers need sizes from root           
c                                                                               
      if( cp_matls_present == -1 ) then                                         
         found_cp = .false.                                                     
         do matnum = 1, nummat  ! nummat in common.main                         
           if(  matprp(9,matnum) /= 10 ) cycle                                  
           found_cp = .true.                                                    
           exit                                                                 
         end do                                                                 
         cp_matls_present = 0                                                   
         if( found_cp ) then                                                    
            cp_matls_present = 1                                                
            call mm10_set_history_locs                                          
            call wmpi_compute_set_history_locs                                  
         end if                                                                 
      end if                                                                    
c                                                                               
c              add more setups here ... may need to add another                 
c              wmpi_... routine as well                                         
c                                                                               
      return                                                                    
      end                                                                       
