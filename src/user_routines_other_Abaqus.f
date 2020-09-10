c                                                                               
c           user_routines_other_abaqus.f   Distribution version                 
c                                                                               
c           Updated:  9/8/2020 rhd
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
                                                                   
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine sigini                            *          
c     *                                                              *          
c     *                  written by : rhd                            *          
c     *                                                              *          
c     *                last modified : 9/9/2020 rhd                  *          
c     *                                                              *          
c     *     Abaqus compatible routine to supply initial stresses     *          
c     *                                                              *          
c     ****************************************************************          
c
      module work_for_sigini
      implicit none
      integer :: variable !  whatewver is needed
      end module work_for_sigini
c
      subroutine sigini( sigma, coords, ntens, ncrds, noel, npt,
     &                   layer, kspt, lrebar, names, thread, out,
     &                   initial_stresses_file )
      use work_for_sigini
c
      implicit none
c
c         sigma, coords, ntens, ncrds, noel, npt
c         follow Abaqus documentation.
c         sigma ordering follows Abaqus scheme: 11, 22, 33, 12, 13, 23
c         always: ntens = 6, ncrds = 3, npt = 1
c         layer, kspt, lrebar are not defined on input,
c         names(1) = blank, names(2) = WARP3D element name,
c         e.g. q3disop
c         thread = OMP thread running this routine
c         out = Fortran I/O device number to write messages, 
c         debug output
c
c         routine is called with noel = 0 outside a threaded region
c         to allow setup of any read-only variables in the module.
c         the initial_stresses_file (if used) must be opened, read, 
c         closed for the call with noel = 0.
c
      integer, intent(in) :: ntens, ncrds, noel, npt,
     &                       layer, kspt, lrebar, thread, out
      double precision, intent(in) ::  coords(3)
      double precision, intent(out) :: sigma(6)
      character(len=80), intent(in) :: names(2)
      character(len=*), intent(in) :: initial_stresses_file
c
c             local variables
c
      double precision :: x, y, z
      logical :: file_exists
      logical :: t1, t2
c
      double precision, parameter :: zero = 0.0d0
      logical, parameter :: ldebug = .true.
c
      return  !     NOTE 
c
      if( ldebug ) then
        write(out,9000) noel, npt, thread, names(2), coords
      end if
c
c              do any setup. initial read-only module variables,
c              read data from file into module for read-only use
c              during (parallel) processing to return initial
c              stresses
c
      if( noel == 0 ) then
         inquire( file=initial_stresses_file, exist=file_exists )
         return ! no setup needed
      end if
c
      x = coords(1); y = coords(2); z = coords(3)
      t1 = x > 0.2501d0 .and. x <= 0.7495d0
      t2 = y > 0.0001d0 .and. y <= 0.2495d0
      if( t1 .and. t2 ) then
        sigma(1) = -57.6d0
        sigma(2) = -57.6d0
        sigma(3) = -45.6d0 
      end if

!      if( noel <= 200 ) then
!        sigma(1) = -57.6d0
!        sigma(2) = -57.6d0
!        sigma(3) = -45.6d0 
!      end if
c
      if( ldebug ) then
        write(out,9100) sigma
      end if
     
c
      return
c
 9000 format(/,'...... debug inside sigini .....',
     &       /,'       noel, npt, thread: ',i8, i2, i4,
     &       /,'       element type: ',a,
     &       /,'       (x,y,z): ',3f15.6,
     &       // )
 9100 format(  '       sigma: ',6f10.3)
c
      end


