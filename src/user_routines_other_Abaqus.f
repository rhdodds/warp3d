c                                                                               
c           user_routines_other_abaqus.f   Distribution version                 
c                                                                               
c           Updated:  6/13/2017 rhd                                                  
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine uexternaldb                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/14/2013                  *          
c     *                                                              *          
c     *             Abaqus compatible support routine                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine uexternaldb( lop, lrestart, time, dtime, kstep, kinc )         
      implicit double precision (a-h,o-z)                                       
c                                                                               
c      Recall that Abaqus 'step" is == 1 in WARP3D. Abaqus                      
c      "increment" is the WARP3D load (time) step number.                       
c                                                                               
c     lop  indicates that the subroutine is being called                        
c       = 0 at the start of the analysis.                                       
c       = 1 at the start of the current analysis increment.                     
c         the subroutine can be called multiple times at the beginning          
c         of an analysis increment if the increment fails to converge           
c         and a smaller time increment is required.                             
c       = 2 at the end of the current analysis increment.                       
c         when lop=2, all information that you need to restart the              
c         analysis should be written to external files.                         
c       = 3 at the end of the analysis.                                         
c       = 4 at the beginning of a restart analysis.                             
c         when lop=4, all necessary external files should be opened             
c         and properly positioned and all information required for              
c         the restart should be read from the external files.                   
c                                                                               
c     lrestart                                                                  
c       = 0 indicates that an analysis restart file is not being                
c         written for this increment.                                           
c       = 1 indicates that an analysis restart file is being written            
c         for this increment.                                                   
c       = 2 indicates that an analysis restart file is being written            
c         for this increment and that only one increment is being               
c         retained per step so that the current increment overwrites            
c         the previous increment in the restart file                            
c                                                                               
c     time(1) = value of current (simulation) step time                         
c     time(2) = value of current (simulation) total time                        
c               time(1) = time(2) in WARP3D since there is only 1               
c               Abaqus 'step"                                                   
c                                                                               
c     dtime = (simulation) time increment (step ime increment in                
c             WARP3D)                                                           
c                                                                               
c     kstep = 1 for WARP3D                                                      
c                                                                               
c     kinc = current increment number. when lop=4,                              
c            kinc gives the restart increment number.                           
c            kinc is the WARP3D step number. kinc = 1 is                        
c            1st simulation step                                                
c                                                                               
c                                                                               
      integer :: lop, lrestart, kstep, kinc, termout                            
      dimension time(2)                                                         
c                                                                               
      logical debug_local                                                       
                                                                                
                                                                                
      termout = 6                                                               
      debug_local = .false.                                                     
c                                                                               
      if( debug_local ) then                                                    
        write(termout,9000) lop, lrestart, time(1), time(2),                    
     &                      dtime, kstep, kinc                                  
      endif                                                                     
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/,3x,"Entered UEXTERNALDB:",                                       
     & /,10x,"lop, lrestart: ",2i3,                                             
     & /,10x,"time(1),(2): ",2e16.6,2x,"dtime: ",e16.6,                         
     & /,10x,"kstep, kinc: ",2i3,//)                                            
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine uexternaldb_store              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/23/2016 rhd              *          
c     *                                                              *          
c     *                   uexternaldb support routine                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine uexternaldb_store( fileno, iout )                              
c                                                                               
      use mod_user_routines                                                     
      implicit none                                                             
c                                                                               
      integer :: fileno, iout                                                   
c                                                                               
c              local variables                                                  
c                                                                               
      logical, parameter :: local_debug = .true.                                
                                                                                
c      write(fileno) active_profile_up                                          
c      if( local_debug ) write(iout,9000) active_profile_up                     
                                                                                
      return                                                                    
c                                                                               
 9000 format(".... uexternaldb_store ...",                                      
     & /,10x"active_profile_up: ",i6,// )                                       
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine uexternaldb_reopen             *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/23/2016 rhd              *          
c     *                                                              *          
c     *                   uexternaldb support routine                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine uexternaldb_reopen( fileno, iout )                             
c                                                                               
      use mod_user_routines                                                     
      implicit none                                                             
c                                                                               
      integer :: fileno, iout                                                   
c                                                                               
c              local variables                                                  
c                                                                               
      logical, parameter ::  local_debug = .true.                               
                                                                                
c      read(fileno) computed_profile_from_restart                               
c      if( local_debug ) write(iout,9000)                                       
c     &  computed_profile_from_restart                                          
                                                                                
      return                                                                    
c                                                                               
 9000 format(".... uexternaldb_reopen ...",                                     
     & /,10x"computed_profile_from_restart: ",i6,// )                           
      end                                                                       
