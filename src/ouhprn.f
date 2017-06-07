c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouhprn                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 3/11/04 (rhd)              *          
c     *                                                              *          
c     *     this subroutine prints the hardcopy output and/or        *          
c     *     the record within the packet for a single element        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouhprn( elem, type, ngp, nnode, long, nodpts,                  
     &                   stress, wide, eform, prec, lnum, pgnum,                
     &                   lbltyp, strlbl, hedtyp, num_short_stress,              
     &                   num_short_strain, center_output, geonl,                
     &                   noheader, out_packet_now )                             
      use global_data ! old common.main
      implicit integer (a-z)                                                    
      character(len=8) :: strlbl(*)                                             
      character(len=*) :: hedtyp                                                
      logical long, nodpts, stress, wide, eform, prec,                          
     &        center_output, geonl, noheader, out_packet_now                    
c                                                                               
c                       local declartions                                       
c                                                                               
      logical newel                                                             
c                                                                               
c                                                                               
c                       set the number of strain or stress                      
c                       values for output.                                      
c                                                                               
      nstrou = num_short_strain                                                 
      if ( stress ) nstrou = num_short_stress                                   
      if( long .or. out_packet_now ) then                                       
         nstrou = num_short_strain + 15                                         
         if ( stress ) nstrou = num_short_stress + 15                           
      end if                                                                    
c                                                                               
      bele  = 1                                                                 
      newel = .true.                                                            
c                                                                               
c                       set the stress labels for the current element.          
c                                                                               
      call oulbst( stress, lbltyp, type, elem, strlbl, long, hedtyp,            
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
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
