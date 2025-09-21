c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouhstr                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/17/2025 rhd              *          
c     *                                                              *          
c     *     drive output of element strains/stresses to printed      *          
c     *     hardcopy and/or to packets files                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouhstr( do_stress, wide, eform, prec,                          
     &                   noheader, out_packet_now, element_list,                
     &                   num_list_entries )   
c                                  
      use global_data, only : ltmstp, iprops, lprops, elblks
      use main_data, only: packet_file_no, cohesive_ele_types,                  
     &                     elems_to_blocks, bar_types, link_types                                      
c
      implicit none
c
      integer, intent(in) :: num_list_entries  
      integer, intent(in) :: element_list(num_list_entries)                                    
      logical, intent(in) :: wide, eform, prec, do_stress, noheader, 
     &                       out_packet_now                                                       
c                                                                               
c             local declarations                                                
c     
      integer :: local_count, index, elem, elem_type, num_enodes,
     &           num_int_points, output_loc, additional, 
     &           cohes_pkt_offset, pkt_type, lnum, pgnum, lbltyp,
     &           span, blk, rel_elem, int_order, mat_type, totdof,
     &           cohesive_type, felem, num_enode_dof
      character(len=8) :: strlbl(30), hedtyp*30                                 
      real :: dumr                                                                 
      double precision :: dumd                                                                 
      character :: dums                                                         
      logical :: bbar_flg, geo_non_flg, long_out_flg, solid_elem,                              
     &           nodpts_flg, center_output, cohesive_elem,                        
     &           at_intpts, at_enodes, at_center, bar_elem, link_elem                                  
c                                                                               
c             for packet output, loop over the user specified                   
c             element list to determine how many data records                   
c             will be in the packet. the packet contains                        
c             as many elements an in the list with data                         
c             at integration points, element nodes, or center.                  
c                                                                               
      if( out_packet_now ) then                                                 
c                                                                               
        local_count = 0                                                         
        do index = 1, num_list_entries                                          
          elem           = element_list(index)                                  
          elem_type      = iprops(1,elem)                                       
          num_enodes     = iprops(2,elem)                                       
          num_int_points = iprops(6,elem)                                       
          output_loc     = iprops(12,elem)                                      
          cohesive_elem  = cohesive_ele_types(elem_type)                        
          at_intpts      = output_loc .eq. 1                                    
          at_enodes      = output_loc .eq. 2                                    
          at_center      = output_loc .eq. 3                                    
          additional     = 0                                                    
          cohes_pkt_offset = 0                                                  
          if( cohesive_elem ) then                                              
             additional = num_int_points                                        
             if( at_center )  additional = 1                                    
             cohes_pkt_offset = 16                                              
          else   ! solid element                                                
             if( at_enodes ) additional = num_enodes                           
             if( at_center ) additional = 1                                     
             if( at_intpts ) additional = num_int_points                        
          end if                                                                
          local_count = local_count + additional                                
         end do                                                                 
         pkt_type = 15 + cohes_pkt_offset                                       
         if( .not. do_stress ) pkt_type = pkt_type + 1                          
         write(packet_file_no) pkt_type, local_count, ltmstp, 0                 
c                                                                               
      end if                                                                    
c                                                                               
c             initialize parameters controlling hardcopy output.                
c                                                                               
      lnum               = 56  ! forces printing of headers                     
      pgnum              = 0                                                    
      lbltyp             = 0   ! tracks last printed labels
      strlbl(1:30)       = ' '                                                  
c                                                                               
c             pass over the use list of elements for output.                    
c             call lower routines to compute secondary (derivable)              
c             output values, then output values to hardcopy and/or              
c             packet.                                                           
c                                                                               
      do index = 1, num_list_entries                                            
        elem              = element_list(index)                                 
        span              = 1                                                   
        blk               = elems_to_blocks(elem,1)                             
        rel_elem          = elems_to_blocks(elem,2)                             
        felem             = elblks(1,blk)                                       
        elem_type         = iprops(1,elem)                                      
        int_order         = iprops(5,elem)                                      
        mat_type          = iprops(25,elem)                                     
        num_enodes        = iprops(2,elem)                                      
        num_enode_dof     = iprops(4,elem)                                      
        totdof            = num_enodes * num_enode_dof                          
        geo_non_flg       = lprops(18,elem)                                     
        bbar_flg          = lprops(19,elem)                                     
        num_int_points    = iprops(6,elem)                                      
        long_out_flg      = lprops(16,elem)                                     
        output_loc        = iprops(12,elem)                                     
        nodpts_flg        = output_loc .eq. 2                                   
        center_output     = output_loc .eq. 3                                   
        if( out_packet_now ) long_out_flg = .true.                              
        cohesive_type     = iprops(27,elem)                                     
        cohesive_elem     = cohesive_ele_types(elem_type)   
        bar_elem          = bar_types(elem_type)
        link_elem         = link_types(elem_type)
c                                                                               
c             duplicate necessary element block data.                           
c  
        call oudups( span, elem, num_int_points, geo_non_flg,                   
     &               do_stress, cohesive_elem )                                 
c                                                                               
c             compute the element block stress/strain data to be                
c             output. set to get new header labels if needed.                   
c             solids and interface-cohesive treated separately.  
c    
c             lbltyp = 1 for solids.
c                    = 2 for bars
c                    = 3 for links
c                    = 5 +  cohesive_type for cohesive elements              
c                                                                               
        if( cohesive_elem ) then                                                
           if( lbltyp .ne. 5 + cohesive_type  ) then
            lbltyp = 5 + cohesive_type 
            lnum = 56 ! forces new page & label printing
          end if                                                                
          call ouhrks_cohesive( elem, blk, felem, elem_type,                    
     &       int_order, num_int_points, num_enodes, center_output,              
     &       do_stress, mat_type, cohesive_type, wide, eform,                   
     &       prec, lnum, pgnum, lbltyp, strlbl, hedtyp,  noheader,             
     &       out_packet_now, geo_non_flg, long_out_flg )                                      
          cycle                                                                 
        end if
c        
        if( bar_elem ) then
          if( lbltyp .ne. 2  ) then
            lbltyp = 2
            lnum = 56 ! forces new page & label printing
          end if
          call ouhrks_bar( elem, blk, felem, elem_type, num_enodes,    
     &        do_stress, mat_type, wide, eform, prec, lnum, pgnum,
     &        lbltyp, strlbl, hedtyp,  noheader, out_packet_now, 
     &        geo_non_flg )   
           cycle
        end if
c                                                                                    
        if( link_elem ) then
          if( lbltyp .ne. 3  ) then
            lbltyp = 3
            lnum = 56 ! forces new page & label printing
          end if
          call ouhrks_link( elem, blk, felem, elem_type, num_enodes,    
     &        do_stress, mat_type, wide, eform, prec, lnum, pgnum,
     &        lbltyp, strlbl, hedtyp,  noheader, out_packet_now, 
     &        geo_non_flg )   
           cycle
        end if                                                                            
c                                                                               
c             process solid element                                             
c                                                                               
        if( lbltyp .ne. 1 ) then
          lbltyp = 1
          lnum = 56 ! forces new page & label printing
        end if  
c
        call ouhrks( span, blk, felem, elem_type, int_order,                    
     &                num_int_points, num_enodes,                               
     &                geo_non_flg, long_out_flg,                                
     &                nodpts_flg, do_stress, mat_type,                          
     &                center_output )
c                                                                               
       call ouhprn(  elem, elem_type, num_int_points,                          
     &                num_enodes, long_out_flg,                                 
     &                nodpts_flg, do_stress, wide, eform, prec,                 
     &                lnum, pgnum, lbltyp, strlbl, hedtyp,                      
     &                center_output, geo_non_flg, noheader,                     
     &                out_packet_now )                                          
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ou_gastr                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/22/21 rhd               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c               moved here to stop compilers from complaining about
c               history_global is history....(blk)%ptr(location).
c               compilers don't like this type of old-style Fortran
c               array passing. it thinks ptr might be an array section.
c
c               by putting routine here and with separate compiles, the 
c               cont's see what we're doing (it is legal!) 
c     
c                                                                                    
      subroutine ou_gastr( history_local, history_global, ngp,                  
     &                     mxhist, mxngp, hist_size, span, mxvl )               
      implicit none                                                             
c                                                                               
c               parameter declarations                                          
c                                                                               
      integer, intent(in) :: ngp, mxhist, mxngp, hist_size, mxvl, span                         
      double precision, intent(in) :: history_global(hist_size,ngp,span)                                       
      double precision, intent(out) :: history_local(mxvl,mxhist,mxngp)
c                                                                               
c               local declarations                                              
c                                                                               
      integer i, j, k                                                           
c                                                                               
      if( ngp .ne. 8 ) then                                                     
        do k = 1, ngp                                                           
         do  j = 1, hist_size                                                   
!DIR$ VECTOR ALIGNED
            do  i = 1, span                                                     
               history_local(i,j,k) = history_global(j,k,i)                     
            end do                                                              
         end do                                                                 
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c                number of gauss points = 8, unroll.                            
c                                                                               
      do  j = 1, hist_size                                                      
!DIR$ VECTOR ALIGNED
        do  i = 1, span                                                         
            history_local(i,j,1) = history_global(j,1,i)                        
            history_local(i,j,2) = history_global(j,2,i)                        
            history_local(i,j,3) = history_global(j,3,i)                        
            history_local(i,j,4) = history_global(j,4,i)                        
            history_local(i,j,5) = history_global(j,5,i)                        
            history_local(i,j,6) = history_global(j,6,i)                        
            history_local(i,j,7) = history_global(j,7,i)                        
            history_local(i,j,8) = history_global(j,8,i)                        
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
