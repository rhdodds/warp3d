c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouhprn                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 9/17/2025 rhd              *          
c     *                                                              *          
c     *     prints the hardcopy output and/or                        *          
c     *     the record within the packet for a single element        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouhprn( elem, type, ngp, nnode, long, nodpts,                  
     &                   do_stress, wide, eform, prec, lnum, pgnum,                
     &                   lbltyp, strlbl, hedtyp, center_output, geonl,                
     &                   noheader, out_packet_now ) 
c                                 
      use output_value_indexes, only : num_short_strain,
     &        num_short_stress, num_long_strain, num_long_stress
      implicit none
c      
      integer :: elem, type, ngp, nnode, lnum, lbltyp
      character(len=8) :: strlbl(*)                                             
      character(len=*) :: hedtyp                                                
      logical :: long, nodpts, do_stress, wide, eform, prec,                          
     &           center_output, geonl, noheader, out_packet_now                    
c                                                                               
c                       local declartions                                       
c                                                                               
      logical :: newel
      integer :: nstrou, pgnum, bele                                                             
c                                                                               
c                                                                               
c                       set the number of strain or stress                      
c                       values for output.                                      
c                                                                               
      nstrou = num_short_strain                                                 
      if ( do_stress ) nstrou = num_short_stress                                   
      if( long .or. out_packet_now ) then                                       
         nstrou = num_long_strain
         if ( do_stress ) nstrou = num_long_stress
      end if                                                                    
c                                                                               
      bele  = 1                                                                 
      newel = .true.                                                            
c                                                                               
c                       set the stress-strain output labels for the 
c                       current element.          
c                                                                               
      call oulbst( do_stress, lbltyp, type, elem, strlbl, long, hedtyp,            
     &             geonl, 0  )                                                  
c                                                                               
c                       output the stress/strain data for the current           
c                       element.                                                
c                                                                               
      call ouhel( bele, elem, hedtyp, nnode, ngp, nodpts, nstrou, wide,         
     &            eform, prec, strlbl, pgnum, lnum, newel,                      
     &            center_output, noheader, out_packet_now )                     
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
