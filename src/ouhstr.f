c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouhstr                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
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
      use main_data, only: packet_file_no, cohesive_ele_types,
     &                     elems_to_blocks
      implicit integer (a-z)
      include 'common.main'
      logical wide, eform, prec, do_stress, noheader, out_packet_now
     &        solid_elem
      integer element_list(num_list_entries)
c
c             local declarations
c
      character*8 strlbl(30), hedtyp*30
      real dumr
      double precision
     &     dumd
      character dums
      logical  bbar_flg, geo_non_flg, long_out_flg,
     &         nodpts_flg, center_output, cohesive_elem,
     &         at_intpts, at_enodes, at_center
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
             if( at_enodes  ) additional = num_enodes
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
      lbltyp             = 0   ! not sure this is still needed
      strlbl(1:30)       = ' '
      last_solid_type    = 0   ! for tracking need to change headers
      last_cohesive_type = 0   !     "                  "
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
        num_short_stress  = 11
        num_short_strain  = 7
        cohesive_type     = iprops(27,elem)
        cohesive_elem     = cohesive_ele_types(elem_type)
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
        if( cohesive_elem ) then
          if( last_cohesive_type .ne. cohesive_type ) then
             lnum = 56; last_cohesive_type = cohesive_type
          end if
          last_solid_type = 0
          call ouhrks_cohesive( elem, blk, felem, elem_type,
     &       int_order, num_int_points, num_enodes, center_output,
     &       do_stress, mat_type, cohesive_type, wide, eform,
     &        prec, lnum, pgnum, lbltyp, strlbl, hedtyp,  noheader,
     &       out_packet_now, geo_non_flg )
          cycle
        end if
c
c             process solid element
c
        last_cohesive_type = 0
        if(  last_solid_type .eq. 0 ) then
           lnum = 56; last_solid_type = 1
        end if
        call ouhrks( span, blk, felem, elem_type, int_order,
     &                num_int_points, num_enodes,
     &                geo_non_flg, long_out_flg,
     &                nodpts_flg, do_stress, mat_type,
     &                center_output, num_short_stress,
     &                num_short_strain )
c
        call ouhprn(  elem, elem_type, num_int_points,
     &                num_enodes, long_out_flg,
     &                nodpts_flg, do_stress, wide, eform, prec,
     &                lnum, pgnum, lbltyp, strlbl, hedtyp,
     &                num_short_stress, num_short_strain,
     &                center_output, geo_non_flg, noheader,
     &                out_packet_now )
c
      end do
c
      return
      end


