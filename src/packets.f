c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine open_packets_file            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/25/26 rhd               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine open_packets_file( msg_flag )                                    
c
      use main_data, only : output_packets, packet_file_name,                   
     &                      packet_file_no                                      
c
      implicit none                                                    
c
      logical :: msg_flag                                                   
c
      integer :: dum
      real :: dumr                                                                 
      character :: dums                                                         
      double precision :: dumd                                                                   
      logical :: found                                                   
                                                                                
c                                                                               
c                        open the existing binary file of packet                
c                        results if name found in restart file.                 
c                                                                               
      if( output_packets ) then                                                
         inquire( file=packet_file_name, exist=found )                          
         if( found .and. msg_flag )                                             
     &                   call errmsg2( 26, dum, packet_file_name,               
     &                                 dumr, dumd )                             
         open(unit=packet_file_no, file=packet_file_name,                       
     &        access='sequential',form='unformatted',                           
     &        status='unknown',position='append')                               
      end if        
c                                                                  
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine close_packets_file           *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 10/14/00                   *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                                                                                                             
      subroutine close_packets_file( msg_flag )                                   
c
      use main_data, only : packet_file_name, packet_file_no  
                        
      implicit none
c      
      logical :: msg_flag                                               
c
      integer :: dum
      real :: dumr                                                                 
      character :: dums                                                         
      double precision :: dumd                                                                   
      logical :: connected
c                                                                               
c                        open the existing binary file of packet                
c                        results if name found in restart file. 
c                
      inquire( unit=packet_file_no, opened=connected )                          
      if( connected ) then                                                     
        close(unit=packet_file_no,status='keep')                                
        if( msg_flag )                                                         
     &    call errmsg2( 27, dum, packet_file_name, dumr, dumd )                 
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
