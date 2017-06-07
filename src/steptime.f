                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine steptime                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 1/17/2016 rhd              *          
c     *                                                              *          
c     *     determine if the then load step will likely exceel the   *          
c     *     user specified maximum wall clock time.if so, write a    *          
c     *     restart file and quit.                                   *          
c     *     if the user-did not specify a time limit in the          *          
c     *     nonlinear parameters, skip this processing               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine steptime( step, status )                                       
      use global_data ! old common.main
      implicit integer (a-z)                                                    
c                                                                               
      integer :: step, status                                                   
c                                                                               
      integer :: idummy, last                                                   
      real, save :: last_step_time, time_before_step                            
      logical, save :: ignore                                                   
      double precision :: dumd                                                  
      character(len=80) :: name                                                 
      real :: wcputime, wwalltime, dumr, percent, time_so_far                   
      external :: wcputime, wwalltime                                           
      logical ldum1, ldum2, debug, quit_now                                     
      data percent, debug / .90, .false./                                       
c                                                                               
      select case( status )                                                     
      case( 1 )                                                                 
c                                                                               
c               status = 1: no steps have been run yet.  set                    
c               last_step_time = 0 and ignore = true                            
c                                                                               
        ignore = .true.                                                         
        last_step_time = 0.0                                                    
        if( debug ) write(out,*) '>>> last step time set to zero'               
c                                                                               
      case( 2 )                                                                 
c                                                                               
c               status = 2: before a load step.  if ignore, then                
c               return.  else, check the last load step versus                  
c               total wall time. if it is within 'percent' of                   
c               time_limit, write a restart file and quit.                      
c                                                                               
        if( ignore ) return                                                     
        time_so_far = wwalltime(idummy)                                         
        quit_now = percent .lt. (time_so_far +                                  
     &             last_step_time)/(time_limit)                                 
        if( .not. quit_now ) then                                               
          time_before_step = wwalltime(idummy)                                  
        else                                                                    
           if( debug ) then                                                     
             write(out,'("time_limit     =",f10.3)') time_limit                 
             write(out,'("last_step_time =",f10.3)') last_step_time             
             write(out,'("time now       =",f10.3)') wwalltime(idummy)          
           endif                                                                
c                                                                               
           last = 8                                                             
           call name_strip( stname, last )                                      
           name = stname(1:last) // '_overtime_db'                              
           call errmsg( 195, step, name, dumr, dumd )                           
           call store ( 'itsblank', name, ldum1, ldum2 )                        
           call die_gracefully                                                  
        endif                                                                   
c                                                                               
      case( 3 )                                                                 
c                                                                               
c               status = 3: after a load step.  if ignore, then                 
c               ignore. calculate the time                                      
c               needed for a complete step.                                     
c                                                                               
        if( ignore ) return                                                     
        last_step_time = wwalltime(idummy) - time_before_step                   
        if( debug ) write(out,'("last_step_time =",f10.3)')                     
     &             last_step_time                                               
c                                                                               
      case( 4 )                                                                 
         ignore = time_limit .le. 0.0                                           
c                                                                               
      case default                                                              
         return                                                                 
      end select                                                                
c                                                                               
      return                                                                    
      end                                                                       
