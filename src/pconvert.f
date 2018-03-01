c
c **********************************************************************
c *                                                                    *
c *                          FILE: pconvert.f                          *
c *                                                                    *
c *  This file contains a module and various subroutines for the high  *
c *  level command "convert", used to convert binary packets to the    *
c *  more portable ascii format.                                       *
c *                                                                    *
c **********************************************************************
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine pconvdriv                    *
c     *                                                              *
c     *                       written by : agc                       *
c     *                                                              *
c     *                   last modified : 05/15/06                   *
c     *                                                              *
c     *    this subroutine drives the execution of the high level    *
c     *    command "convert" for converting binary packets to ascii. *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine pconvdriv
      use global_data ! old common.main
      use main_data, only: packet_file_no, packet_file_name,
     &                   ascii_packet_file_no, ascii_packet_file_name
      use pvars
c
      implicit integer (a-z)

      integer dumi, dot
      real dumr
      double precision
     &   dumd
      character :: dums
      logical duml, matchs, endcrd, numd, label, integr
      logical string, packets_requested, all_packets
      integer, dimension(1:max_packet_types,1:2) :: ou_packet_control
c
c             Initialize Packet Control Variable
c              1st col = (*,1) = was conversion for packet * requested?
c              2nd col = (*,2) = how many times was packet * in binary file?
c
      all_packets=.false.
      do dumi = 1,max_packet_types
         ou_packet_control(dumi,1)=0
         ou_packet_control(dumi,2)=0
      end do
c
c             Search input line for packet numbers to convert
c             and store them in array. With the exception of an array
c             overflow check, store every request regardless
c             if it is a valid packet or not. (Someday it could be
c             a valid packet.) If invalid, no packets with that number
c             will be found and the end status message will reflect that.
c
      outer: do
         if(integr(dumi)) then
             if( dumi .lt. max_packet_types .and. dumi.gt.0 ) then
                 ou_packet_control(dumi,1)=1
             end if
         else if(matchs('all',3)) then
c                      No need to set requests now. All
c                      "found" packets will become "requested".
             all_packets=.true.
             do dumi = 1,max_packet_types
                ou_packet_control(dumi,1)=0
             end do
         else
             exit outer
         end if
      end do outer
c
c                       Look for optional words in command
c
      if (matchs('to',2))    call splunj
      if (matchs('ascii',5)) call splunj
      if ( matchs('on',2) .or. matchs('to',2) ) call splunj
      if (matchs('file',4))  call splunj
c
c                       Find and Store ascii output filename
c
      if ( label(dummy) ) then
         ascii_packet_file_name = ' '
         call entits(ascii_packet_file_name, nc)
         call name_strip(ascii_packet_file_name, nc)
      else if( string(dummy) ) then
         ascii_packet_file_name = ' '
         call entits(ascii_packet_file_name, nc)
         call name_strip(ascii_packet_file_name, nc)
      else
c
c         set the default ascii packet file name based on the
c         binary packet file name. Find the "." for the file
c         extension, then append the "tpf" default extension
c         to the file name.
c
         dot = index(packet_file_name,'.',back=.true.) - 1
         if(dot.eq.0) then
           dot = len_trim(packet_file_name)
         end if
         ascii_packet_file_name(1:) = packet_file_name(1:dot)//'.tpf'
         call errmsg2(68, dumi, dums, dumr, dumd)
      end if
c
c                       Final check to make that there are
c                       actually requested packets if there are none,
c                       default to outputting all packets.
c
      packets_requested = .false.
      do dumi=1,max_packet_types
         if (ou_packet_control(dumi,1).eq.1) then
            packets_requested = .true.
            exit
         end if
      end do
c
      if (packets_requested .eqv. .false.) then
         if (all_packets .eqv. .false.) then
           call errmsg2(69, dumi, dums, dumr, dumd)
c                      No need to set requests now. All
c                      "found" packets will become "requested".
           all_packets=.true.
           do dumi = 1,max_packet_types
             ou_packet_control(dumi,1)=0
           end do
         end if
      end if
c
c                       Call subroutine to convert chosen packets
c
      call pconv(ou_packet_control,all_packets)
      return
c
      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine pconv                        *
c     *                                                              *
c     *                       written by : rhd/gt/agc                *
c     *                                                              *
c     *                   last modified : 05/15/06                   *
c     *                                                              *
c     *    this subroutine converts the specified packets from the   *
c     *    binary packet file to a separate ascii file. It's a       *
c     *    modified version of the WARP packet reader.               *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine pconv(ou_packet_control, all_packets)
      use global_data ! old common.main
c
      use main_data, only: packet_file_name, packet_file_no,
     &                  ascii_packet_file_name, ascii_packet_file_no
      use pvars
      implicit integer (a-z)

      integer dumi
      real dumr
      double precision
     &   dumd
      character :: dums
      logical duml, binary_file_was_open, all_packets
      integer ios,type,num_lines,step,iter,lines_read
      integer output_set,last_step,last_iter
      integer, dimension(1:max_packet_types,1:2) :: ou_packet_control
c
c             Initialize some variables.
c
      output_set = 0
      saved_step = 0
      debug = .false.
      type = 0
      iter = 0
      num_lines = 0
      step = 0

c
c             Binary packet file is probably open. Close it,
c             then reopen it to set the position to the start
c             for conversion.
c
        inquire( unit=packet_file_no, opened=binary_file_was_open)
        if (binary_file_was_open) then
          call close_packets_file(.false.)
          open( unit=packet_file_no, file=packet_file_name,
     &        status='old', access='sequential', form='unformatted',
     &        position='rewind', iostat=ios)
        else
          open( unit=packet_file_no, file=packet_file_name,
     &        status='old', access='sequential', form='unformatted',
     &        position='rewind', iostat=ios)
        end if
c
c             Print status message.
c
        if( ios .ne. 0) then
           call errmsg2(70, dumi, dums, dumr, dumd)
           return
        else if (binary_file_was_open .eqv. .false.) then
           write(out,9040)
        else
           write(out,9050)
        end if
c
c             Begin an infinite loop, reading packets until the end of
c             file is reached.
c
      infinite: do
c
c             Read binary packet header line. Get the type of packet,
c             the number of lines in the entry, the analysis step, and
c             the analysis iteration.
c
         read(unit=packet_file_no,iostat=ios) type,num_lines,step,iter
c
c              Check i/o status. (End of file, invalid file)
c
         if( ios .gt. 0) then
            call errmsg2(71, dumi, dums, dumr, dumd)
            exit infinite
         else if( ios .lt. 0 ) then
            write(out,9075)
            exit infinite
         end if
c
c             i/o is okay. check for an empty packet (num_lines = 0)
c
         if(num_lines.le.0)then
           cycle
         end if
c
c                   check if the step number has changed to start the
c                   next output set; save the output_set index
c                   to the module variable saved_step for use in the
c                   packet reading routines
c                   (8/25/05 GVT)
c
         if(step.ne.last_step)then
           output_set = output_set + 1
           saved_step = output_set
           last_step = step
           last_iter = iter
         end if      ! step.ne.last_step
c
c                   check to see if the current packet type from the
c                   binary packet file should be analyzed.
c                   ou_packet_control(type,1) will be set
c                   to .true. if the user has specified that type to be
c                   processed. Provide control if the "all_packets" flag
c                   is .true., in which case, every found packet is
c                   essentially requested.
c
         if ( all_packets .eqv. .true.) ou_packet_control(type,1)=1
c
         if( ou_packet_control(type,1).eq.1 ) then
c
c                   Keep track of the number of times each desired
c                   packet is found.  A message will be sent to list
c                   how many of each requested packet was found in the
c                   binary packet file.
c
            ou_packet_control(type,2) = ou_packet_control(type,2) + 1
c
c                   Jump to packet processing routine based on current
c                   packet number. New case numbers need to be added as
c                   packets are added.
c
            select case(type)
               case(1)
                  call displ_node_packet(type,num_lines,step,iter)
               case(2)
                  call veloc_node_packet(type,num_lines,step,iter)
               case(3)
                  call accel_node_packet(type,num_lines,step,iter)
               case(4)
                  call reacts_node_packet(type,num_lines,step,iter)
               case(5)
                  call gurson_elem_1(type,num_lines,step,iter)
               case(6)
                  call gurson_elem_2(type,num_lines,step,iter)
               case(7)
                  call cohes_packet(type,num_lines,step,iter)
               case(8)
                  call ctoa_type_3_packet(type,num_lines,step,iter)
               case(9)
                  call ctoa_type_4_packet(type,num_lines,step,iter)
               case(10)
                  call cohesive_elem_killed(type,num_lines,step,iter)
               case(11)
                  call displ_elem_packet(type,num_lines,step,iter)
               case(12)
                  call veloc_elem_packet(type,num_lines,step,iter)
               case(13)
                  call accel_elem_packet(type,num_lines,step,iter)
               case(14)
                  call reacts_elem_packet(type,num_lines,step,iter)
               case(15)
                  call stress_elem_packet(type,num_lines,step,iter)
               case(16)
                  call strain_elem_packet(type,num_lines,step,iter)
               case(17)
                  call j_integral_packet(type,num_lines,step,iter)
               case(18)
                  call ctoa_type_1_packet(type,num_lines,step,iter)
               case(19)
                  call ctoa_type_2_packet(type,num_lines,step,iter)
               case(20)
                  call gurson_elem_killed(type,num_lines,step,iter)
               case(21)
                  call smcs_elem_killed(type,num_lines,step,iter)
               case(22)
                  call smcs_elem_1(type,num_lines,step,iter)
               case(23)
                  call smcs_elem_2(type,num_lines,step,iter)
               case(24)
                  call step_patt_factors(type,num_lines,step,iter)
               case(27)
                  call i_integral_sif_packet(type,num_lines,step,iter)
               case(28)
                  call i_integral_t_packet(type,num_lines,step,iter)
            end select
c
         else
c                   This packet not wanted. Skip unwanted packets.
c
            do lines_read=1,num_lines
               read(unit=packet_file_no,iostat=ios)
     &                           type,num_lines,step,iter
               if( ios .gt. 0)then
                  call errmsg2(71, dumi, dums, dumr, dumd)
                  exit infinite
               else if( ios .lt. 0 ) then
                  write(out,9075)
                  exit infinite
               end if
            end do
c
         end if
c
      end do infinite
c
c
c     Close both packet files. Reopen binary file for more writing if
c     it was open to start with.
c
      call closer(packet_file_no)
      call closer(ascii_packet_file_no)
      if (binary_file_was_open) then
         call open_packets_file(.false.)
         write(out,9060)
      else
         write(out,9070)
      end if
c
c     Output message listing how many of each packet was found.
c
      write(out,9080)
      do dumi=1,max_packet_types
        if (ou_packet_control(dumi,1).eq.1) then
           write(out,9085) dumi, ou_packet_control(dumi,2)
        end if
      end do
      write(out,9090)
c
      return
c
c
 9040 format(/,15x,'> Binary packet file opened.')
 9050 format(/,15x,'> Binary packet file rewound and ready.')
 9060 format(  15x,'> Binary packet file ready for further writing.',/)
 9070 format(  15x,'> Binary packet file closed.,/')
 9075 format(  15x,'> End of binary packet file reached.')
 9080 format(  15x,'> Packet conversion summary:')
 9085 format(  20x,'> Packet number:',i5,'    Packets converted:',i5)
 9090 format(/)
      end
c
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: closer                                             *
c *                                                                    *
c *     can close both the binary packet file and the                  *
c *     ascii packet file.                                             *
c *                                                                    *
c **********************************************************************
c
      subroutine closer(file)
      use global_data ! old common.main

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit integer (a-z)

      integer, intent(in):: file
      logical connected
c
      inquire( unit=file, opened=connected)
      if( connected )then
        close(unit=file, status='keep')
c        if(file .eq. packet_file_no) write(out,1000)
        if(file .eq. ascii_packet_file_no) write(out,2000)
      end if
      return
c
 1000 format(15x,'> Binary packet file closed.')
 2000 format(15x,'> Ascii packet file closed.')
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: open_output_file                                   *
c *                                                                    *
c *     this subroutine is called from a specific packet types         *
c *     subroutine. it will open the output file where the results     *
c *     should be written.                                             *
c *                                                                    *
c **********************************************************************
c
c
      subroutine open_output_file
      use global_data ! old common.main

      use main_data, only: ascii_packet_file_name, ascii_packet_file_no
      implicit integer (a-z)

      integer ios,dot
      integer dumi
      real dumr
      double precision
     &   dumd
      character :: dums
c
c         open the text results file
c
      open(unit=ascii_packet_file_no, file=ascii_packet_file_name,
     &      access='sequential',form='formatted',
     &      status='replace', iostat=ios)
      if( ios .ne. 0)then
         call errmsg2(75, dumi, dums, dumr, dumd)
         return
      else
         write(out,9040)
      end if
c
      return
c
 9040 format(15x,'> Ascii packet file opened.')
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: ctoa_type_1_packet                                 *
c *                                                                    *
c *     reads packet file for released nodes during the step for       *
c *     non-constant front growth                                      *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine ctoa_type_1_packet( type, num_lines, step, iter )

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem, i, node
      logical connected
      double precision crit_angle, init_crit_angle, ctoa_range,
     &                 angle
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  need to read and put in file
c
      write(ascii_packet_file_no,1100) step
      read(packet_file_no,iostat=ios) crit_angle, init_crit_angle,
     &                                ctoa_range
      if( ios .gt. 0)then
           return
      else if( ios .lt. 0 )then
           return
      end if
c
      do lines_read = 1, num_lines-1
       read(packet_file_no,iostat=ios) node, angle
       if( ios .gt. 0)then
           return
       else if( ios .lt. 0 )then
           return
       end if
       write(ascii_packet_file_no,9010) node, angle
      end do
c
      write(ascii_packet_file_no,9020) crit_angle, init_crit_angle,
     &                                ctoa_range
      return
c
 9010 format(  '        node: ',i7,'  CTOA: ',f6.3,' (degrees)')
 9020 format('        * critical CTOA:',
     &       /,'            for growth    : ',f6.3,' (degrees)',
     &       /,'            for initiation: ',f6.3,' (degrees)',
     &       /,'            nodes released within ',f4.1,' % of ',
     &                        'critical angle.' )
 1100 format(
     & /,' >>> CTOA Node Release for Non-Constant Front:',
     & /,'        Before step: ',i7 )
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: ctoa_type_2_packet                                 *
c *                                                                    *
c *     reads packet file for released nodes during the step for       *
c *     constant front growth                                          *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine ctoa_type_2_packet( type, num_lines, step, iter )

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem, i, node
      logical connected
      double precision crit_angle, init_crit_angle, ctoa_range,
     &                 angle
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  need to read and put in file
c
      write(ascii_packet_file_no,1100) step
      read(packet_file_no,iostat=ios) crit_angle, init_crit_angle,
     &                                ctoa_range
      if( ios .gt. 0)then
           return
      else if( ios .lt. 0 )then
           return
      end if
c
      do lines_read = 1, num_lines-1
       read(packet_file_no,iostat=ios) node, angle
       if( ios .gt. 0)then
           return
       else if( ios .lt. 0 )then
           return
       end if
       if ( node .gt. 0 ) then
         write(ascii_packet_file_no,9010) node, angle
       else
         write(ascii_packet_file_no,9012) -node, angle
       end if
      end do
c
      write(ascii_packet_file_no,9020) crit_angle, init_crit_angle,
     &                                ctoa_range
      return
c
 9010 format(  '        master node:  ',i7,'  CTOA: ',f6.3,' (degrees)')
 9012 format(  '        interim node: ',i7,'  CTOA: ',f6.3,' (degrees)')
 9020 format('        * critical CTOA:',
     &       /,'            for growth    : ',f6.3,' (degrees)',
     &       /,'            for initiation: ',f6.3,' (degrees)',
     &       /,'            nodes released within ',f4.1,' % of ',
     &                        'critical angle.' )
 1100 format(
     &  ' >>> CTOA Node Release for Constant Front Procedure:',/,
     &  '        Before step: ',i7 )
c
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: ctoa_type_3_packet                                 *
c *                                                                    *
c *     reads packet file for ctoa constant growth. start of step      *
c *     summary for each front master node                             *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine ctoa_type_3_packet( type, num_lines, step, iter )

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in) :: type, num_lines, step, iter
c
      integer ios, lines_read, i, node, status
      logical connected
      double precision  angle, angle_fract
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  need to read and put in file
c
      write(ascii_packet_file_no,1100) step
      write(ascii_packet_file_no,*) 'master node       CTOA(degrees)',
     &     '  CTOA/critical CTOA'
      write(ascii_packet_file_no,*) '-----------       -------------',
     &     '  ------------------'
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) node, angle, status, angle_fract
       if( ios .gt. 0) then
           return
       else if( ios .lt. 0 ) then
           return
       end if
       if ( status .eq. 0 ) then
         write(ascii_packet_file_no,9010) node,angle,'(*)',angle_fract
       else
         write(ascii_packet_file_no,9010) node,angle,'   ',angle_fract
       end if
      end do
c
      write(ascii_packet_file_no,9200)
      return
c
 9010 format(6x,i7,7x,f7.4,a3,11x,f4.1,3x,'<=')
 1100 format(
     &  ' >>> CTOA Step Summary for Constant Front:',/,
     &  '        After step: ',i7 )
 9200 format('   NOTE: * indicates critical CTOA is for',
     &        ' initiation', / )
c
      end
c
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: ctoa_type_4_packet                                 *
c *                                                                    *
c *     reads packet file for ctoa non-constant growth. start of step  *
c *     summary for each node                                          *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine ctoa_type_4_packet( type, num_lines, step, iter )

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in) :: type, num_lines, step, iter
c
      integer ios, lines_read, i, node, status, neighbor
      logical connected
      double precision  angle, angle_fract
      character(len=1) :: marker
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  need to read and put in file
c
      write(ascii_packet_file_no,1100) step
      write(ascii_packet_file_no,*)
     &     'crack front node  neighbor  CTOA(degrees)',
     &     '  CTOA/critical CTOA'
      write(ascii_packet_file_no,*)
     &     '----------------  --------  -------------',
     &     '  ------------------'
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) node, neighbor, angle,
     &                                 status, angle_fract
       if( ios .gt. 0)then
           return
       else if( ios .lt. 0 )then
           return
       end if
       marker = ' '
       if ( status .eq. 0 ) marker = '*'
       write(ascii_packet_file_no,9010) node, neighbor, angle, marker,
     &                             angle_fract
      end do
c
      write(ascii_packet_file_no,9200)
      return
 9010 format(6x,i7,9x,i7,7x,f7.4,a1,11x,f4.1,3x,'<=')
 1100 format(
     & /,' >>> CTOA Step Summary for Non-Constant Front:',
     & /,'        After step: ',i7 )
 9200 format('   NOTE: * indicates critical CTOA is for',
     &        ' initiation', / )
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: gurson_elem_1                                      *
c *                                                                    *
c *     reads packet file for gurson element property packets          *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine gurson_elem_1(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision gurson(6)
      integer special_1, special_2
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected ) write(ascii_packet_file_no,1100)
c
c     need to read and put in file
c
      do lines_read=1,num_lines
         read(packet_file_no,iostat=ios)elem,(gurson(i), i=1,5),
     &        special_1,gurson(6),special_2
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         if ( special_1 .eq. 0 .and. special_2 .eq. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     '*',gurson(6),'*'
         if ( special_1 .ne. 0 .and. special_2 .eq. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     ' ',gurson(6),'*'
         if ( special_1 .eq. 0 .and. special_2 .ne. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     '*',gurson(6),' '
         if ( special_1 .ne. 0 .and. special_2 .ne. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     ' ',gurson(6),' '
      end do
c
      return
c
 1000 format(1x,i7,1x,i7,1x,f12.5,3x,f12.5,3(1x,e14.6),
     &                                   a1,(1x,e14.6),a1)

 1100 format('  step  element     inital f     current f       Ep    ',
     &       '      sigma bar      mean stress    ',
     &       'mises stress')
c
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: gurson_elem_2                                      *
c *                                                                    *
c *     reads packet file for gurson element property packets          *
c *             for automatic load reduction or adaptive load control  *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine gurson_elem_2(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision gurson(7)
      integer special_1, special_2
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      write(ascii_packet_file_no,1100)
c
c     need to read and put in file
c
      do lines_read=1,num_lines
         read(packet_file_no,iostat=ios)elem,(gurson(i), i=1,5),
     &        special_1,gurson(6),special_2,gurson(7)
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         if ( special_1 .eq. 0 .and. special_2 .eq. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     '*',gurson(6),'*',gurson(7)
         if ( special_1 .ne. 0 .and. special_2 .eq. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     ' ',gurson(6),'*',gurson(7)
         if ( special_1 .eq. 0 .and. special_2 .ne. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     '*',gurson(6),' ',gurson(7)
         if ( special_1 .ne. 0 .and. special_2 .ne. 0 )
     &     write(ascii_packet_file_no,1000)step,elem,(gurson(i), i=1,5),
     &     ' ',gurson(6),' ',gurson(7)
      end do
c
      return
c
 1000 format(1x,i7,1x,i7,1x,f12.5,3x,f12.5,3(1x,e14.6),
     &                a1,(1x,e14.6),a1,1x,f10.6 )

 1100 format('  step  element     inital f     current f       Ep    ',
     &       '      sigma bar      mean stress    ',
     &       'mises stress     delta f ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: gurson_elem_killed                                 *
c *                                                                    *
c *     reads packet file for gurson element property packets          *
c *             for automatic load reduction or adaptive load control  *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine gurson_elem_killed(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision gurson(4)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected ) write(ascii_packet_file_no,1100)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(gurson(i), i=1,4)
               if( ios .gt. 0)then
                  return
               else if( ios .lt. 0 )then
                  return
               end if
              write(ascii_packet_file_no,1000)step,elem,
     &                                 (gurson(i), i=1,4)
            end do
c
      return
c
 1000 format(1x,i7,1x,i7,1x,f12.5,1x,f12.5,3x,2(1x,e14.6))
 1100 format('  step  element   initial f    final f ',
     &       '    mean stress    ',
     &       'mises stress  ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: smcs_elem_1                                        *
c *                                                                    *
c *     reads packet file for gurson element property packets          *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine smcs_elem_1(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision smcs(4)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected ) write(ascii_packet_file_no,1100)
c
c                  need to read and put in file
c
      do lines_read=1,num_lines
         read(packet_file_no,iostat=ios)elem,(smcs(i), i=1,4)
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         write(ascii_packet_file_no,1000)step,elem,(smcs(i), i=1,4)
      end do
c
      return
c
 1000 format(1x,i7,1x,i7,4(1x,e14.6))

 1100 format('  step  element   eps-pls        eps-crit      ',
     &       'sig-mean       sig-mises')
c
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: smcs_elem_2                                        *
c *                                                                    *
c *     reads packet file for gurson element property packets          *
c *             for automatic load reduction or adaptive load control  *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine smcs_elem_2(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision smcs(5)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected ) write(ascii_packet_file_no,1100)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(smcs(i), i=1,5)
               if( ios .gt. 0)then
                  return
               else if( ios .lt. 0 )then
                  return
               end if
              write(ascii_packet_file_no,1000)step,elem,(smcs(i), i=1,5)
            end do
c
      return
c
 1000 format(1x,i7,1x,i7,5(1x,e14.6))

 1100 format('  step  element    eps-pls        eps-crit      ',
     &       'sig-mean       sig-mises       d(eps-pls)'  )
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: smcs_elem_killed                                   *
c *                                                                    *
c *     reads packet file for gurson element property packets          *
c *             for automatic load reduction or adaptive load control  *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine smcs_elem_killed(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision smcs(2)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected ) write(ascii_packet_file_no,1100)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(smcs(i), i=1,2)
               if( ios .gt. 0)then
                  return
               else if( ios .lt. 0 )then
                  return
               end if
              write(ascii_packet_file_no,1000)step,elem,
     &                     (smcs(i), i=1,2)
            end do
c
      return
c
 1000 format(1x,i7,2x,i7,5x,f8.5,4x,f8.5)
 1100 format('  step  element   pl-strain  pl-strain limit ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: cohes_packet                                       *
c *                                                                    *
c *     reads packet file for cohesive element property packets        *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine cohes_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision cohes(6)
      character(len=3) :: special_char
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      write(ascii_packet_file_no,1100)
      write(*,*) '.. code for cohes_packet requires updating'
      write(*,*) '.. execution terminated..'
      call die_gracefully
c
c     need to read and put in file
c
      do lines_read=1,num_lines
         read(packet_file_no,iostat=ios)elem,(cohes(i), i=1,6),
     &        special_char
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         write(ascii_packet_file_no,1000)step,elem,(cohes(i), i=1,6),
     &        special_char
      end do
c
      return
c
 1000 format(1x,i7,1x,i7,1x,6(e14.6),1x,a3)
 1100 format('  step  elem        Tn             Ts            ',
     &       'Dn           Ds          Teff          Deff   ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: cohesive_elem_killed                               *
c *                                                                    *
c *     reads packet file for cohesive element property packets        *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine cohesive_elem_killed(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i
      double precision cohes(2)

      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      write(ascii_packet_file_no,1100)
      write(*,*) '.. code for cohes_elem_killed requires updating'
      write(*,*) '.. execution terminated..'
      call die_gracefully
c
c     need to read and put in file
c
      do lines_read=1,num_lines
         read(packet_file_no,iostat=ios)elem,(cohes(i), i=1,2)
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         write(ascii_packet_file_no,1000)step,elem,(cohes(i), i=1,2)
      end do
c
      return
c
 1000 format(1x,i7,1x,i7,6x,
     &                      f5.2,8x,
     &                      f5.2 )
 1100 format('  step  elem    Deff/Dpeak   Teff/Tpeak ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: reacts_node_packet                                 *
c *                                                                    *
c *     reads packet file for nodal reaction packets                   *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine reacts_node_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision reacts(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

C      if( .not. connected )write(ascii_packet_file_no,1200)
      write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
               if( lines_read .lt. num_lines ) then
                  read(packet_file_no,iostat=ios)node,
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,node,
     &                                     (reacts(i), i=1,3)
               else
c
c                             the last entry on a reactions packet
c                             contains the TOTALS data
c
                  read(packet_file_no,iostat=ios)
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1100)step,
     &                                     (reacts(i), i=1,3)
               end if
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1100 format(1x,i7,1x,' TOT ',1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  Node Reactions',/,
     &       '  step  node   x-reaction    y-reaction    z-reaction   ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: reacts_elem_packet                                 *
c *                                                                    *
c *     reads packet file for element reaction packets                 *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine reacts_elem_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision reacts(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected )write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
               if( lines_read .lt. num_lines ) then
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,elem,node,
     &                                     (reacts(i), i=1,3)
               else
c
c                             the last entry on a reactions packet
c                             contains the TOTALS data
c
                  read(packet_file_no,iostat=ios)
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1100)step,
     &                                     (reacts(i), i=1,3)
               end if
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1100 format(1x,i7,1x,' TOT ELEM  ',1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-reaction   y-reaction    z-reaction   ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: displ_node_packet                                  *
c *                                                                    *
c *     reads packet file for nodal displacement packets               *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine displ_node_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision displ(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

c      if( .not. connected )write(ascii_packet_file_no,1200)
      write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)node,
     &                                    (displ(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,node,
     &                                     (displ(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  Node Displacements',/,
     &       '  step  node   x-disp        y-disp        z-disp   ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: displ_elem_packet                                  *
c *                                                                    *
c *     reads packet file for element displacement packets             *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine displ_elem_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision displ(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected )write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (displ(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,elem,node,
     &                                     (displ(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-disp      y-disp       z-disp     ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: veloc_node_packet                                  *
c *                                                                    *
c *     reads packet file for nodal velocity packets                   *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine veloc_node_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision veloc(3)
      logical connected
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected )write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)node,
     &                                    (veloc(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,node,
     &                                     (veloc(i), i=1,3)
            end do

c
      return
c
 1000 format(1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step   node   x-veloc        y-veloc       z-veloc   ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: veloc_elem_packet                                  *
c *                                                                    *
c *     reads packet file for element velocity packets                 *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine veloc_elem_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision veloc(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected )write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (veloc(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,elem,node,
     &                                     (veloc(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-veloc      y-veloc       z-veloc     ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: accel_node_packet                                  *
c *                                                                    *
c *     reads packet file for nodal acceleration packets               *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine accel_node_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision accel(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected )write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)node,
     &                                    (accel(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,node,
     &                                     (accel(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step   node   x-accel        y-accel       z-accel   ')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: accel_elem_packet                                  *
c *                                                                    *
c *     reads packet file for element acceleration packets             *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine accel_elem_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision accel(3)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

      if( .not. connected )write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (accel(i), i=1,3)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,elem,node,
     &                                     (accel(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,i7,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-accel      y-accel       z-accel     ')
c
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: stress_elem_packet                                 *
c *                                                                    *
c *     reads packet file for element stress packets, long output      *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine stress_elem_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision stress(26)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

c      if( .not. connected )write(ascii_packet_file_no,1200)
      write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
c
c
c                          'node' will contain one of the following
c
c                           - nodes number (nodpts warp3d elem logical)
c                           - gauss points (gausspts warp3d elem logical)
c                           - node = 0 (center_output warp3d elem logical)
c
c
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (stress(i), i=1,26)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,elem,node,
     &                                     (stress(i), i=1,11)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,i7,1x,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6 )
 1200 format('  Stress Results',/,'  step  elem node/gpt    ',
     &       'sx            sy            sz            ',
     &       'sxy           syz           sxz           Uo            ',
     &       'svm           c1            c2            c3')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: strain_elem_packet                                 *
c *                                                                    *
c *     reads packet file for element strain packets, long output      *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine strain_elem_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision strain(22)
      logical connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected )call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)

c      if( .not. connected )write(ascii_packet_file_no,1200)
      write(ascii_packet_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
c
c                     'node' will contain one of the following
c                       - nodes number (nodpts warp3d elem logical)
c                       - gauss points (gausspts warp3d elem logical)
c                       - node = 0 (center_output warp3d elem logical)
c
c
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (strain(i), i=1,22)
                  if( ios .gt. 0)then
                     return
                  else if( ios .lt. 0 )then
                     return
                  end if
                  write(ascii_packet_file_no,1000)step,elem,node,
     &                                     (strain(i), i=1,7)
            end do
c
c
c
      return
c
 1000 format(1x,i7,1x,i7,1x,i7,1x,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6 )
 1200 format('  Strain Results',/,'  step  elem node/gpt    ',
     &       'ex            ey            ez            ',
     &       'gxy           gyz           gxz           eff')
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: j_integral_packet                                  *
c *                                                                    *
c *     reads packet file for j-integral packets                       *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine j_integral_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c             dummy variables
c
      integer, intent(in) :: type,num_lines,step,iter
c
c             local variables
c
      integer :: ios, lines_read, i, ltmstp, nowring, skipped_killed,
     &           num_domains, k
      double precision  :: diterms(9), di_min, di_max, di_avg
      character(len=24) :: domain_id
      character(len=8)  :: stname, lsldnm
      logical           :: connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  read from packet file and write to output
c
      read(packet_file_no,iostat=ios) num_domains, domain_id, stname,
     &                                lsldnm, ltmstp
      if( ios .gt. 0)then
         return
      else if( ios .lt. 0 )then
         return
      end if
      domain_id = ADJUSTL( domain_id )
      stname    = ADJUSTL( stname )
      lsldnm    = ADJUSTL( lsldnm )
      write(ascii_packet_file_no,1100)
     &  TRIM(domain_id),num_domains,TRIM(stname),TRIM(lsldnm),ltmstp
      write(ascii_packet_file_no,1200)
c
c             read J-integral values for domain(s)
c
      do lines_read = 1, num_domains
         read(packet_file_no,iostat=ios) nowring,(diterms(i),i=1,9),
     &        skipped_killed
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         write(ascii_packet_file_no,1300) nowring, (diterms(i),i=1,9),
     &        skipped_killed
      end do
c
c             read max, min and average values
c
      if( num_domains .gt. 1 ) then
         read(packet_file_no,iostat=ios) di_avg, di_min, di_max
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         write(ascii_packet_file_no,1500)
         write(ascii_packet_file_no,1550) di_avg, di_min, di_max
      end if
c
c             read stress intensity factor values computed from
c             J-values.
c
      write(ascii_packet_file_no,1600)
      do i = 1,num_domains
         read(packet_file_no,iostat=ios) nowring,(diterms(k),k=1,5)
         if( ios .gt. 0)then
            return
         else if( ios .lt. 0 )then
            return
         end if
         write(ascii_packet_file_no,1700) nowring, (diterms(k),k=1,5)
      end do
c
      return
c
 1100 format('domain: ',a,2x,'num: ',i6,2x,'structure: ',a,2x,
     &       'loading: ',a,2x,'step_number: ',i12)
 1200 format(1x,'domain',5x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total J',
     &    2x,'killed ele' )
 1300 format(1x,i3,3x,9(1x,e11.4),3x,i3)
 1500 format(1x,'domain average',3x,'domain min',3x,'domain max')
 1550 format(1x,e11.4,5x,e11.4,2x,e11.4)
 1600 format('stress intensity factors from J (single-mode loading):',
     &       /,1x,'domain',2x,'KI pstrs',4x,'KI pstrn',
     &       4x,'KII pstrs',3x,'KII pstrn',3x,'KIII' )
 1700 format(1x,i3,3x,5(1x,e11.4))
c
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: i_integral_sif_packet                              *
c *                                                                    *
c *     reads packet file for interaction-integral                     *
c *     stress intensity factor results                                *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine i_integral_sif_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c             dummy variables
c
      integer, intent(in) :: type,num_lines,step,iter
c
c             local variables
c
      integer :: ios, lines_read, i, j, ltmstp, nowring, skipped_killed,
     &           num_domains
      double precision  :: iiterms(10), ii_min, ii_max, ii_avg,
     &                     j_pstrs, j_pstrn
      character(len=24) :: domain_id
      character(len=8)  :: stname, lsldnm
      logical           :: connected
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  read from packet file and write to output
c
      read(packet_file_no,iostat=ios) num_domains, domain_id, stname,
     &                                lsldnm, ltmstp
      if( ios .gt. 0)then
         return
      else if( ios .lt. 0 )then
         return
      end if
      domain_id = ADJUSTL( domain_id )
      stname    = ADJUSTL( stname )
      lsldnm    = ADJUSTL( lsldnm )

      write(ascii_packet_file_no,1100)
     &  TRIM(domain_id),num_domains,TRIM(stname),TRIM(lsldnm),ltmstp
c
      do j=1,5
         if( j.eq.1 ) write(ascii_packet_file_no,1201)
         if( j.eq.2 ) write(ascii_packet_file_no,1202)
         if( j.eq.3 ) write(ascii_packet_file_no,1203)
         if( j.eq.4 ) write(ascii_packet_file_no,1204)
         if( j.eq.5 ) write(ascii_packet_file_no,1205)
         do lines_read = 1,num_domains
            read(packet_file_no,iostat=ios) nowring,
     &           (iiterms(i),i=1,10), skipped_killed
            if( ios .gt. 0)then
               return
            else if( ios .lt. 0 )then
               return
            end if
            write(ascii_packet_file_no,1400) nowring,
     &        (iiterms(i),i=1,10), skipped_killed
         end do
c
         if( num_domains.gt.1 ) then
            read(packet_file_no,iostat=ios) ii_avg, ii_min, ii_max
            if( ios .gt. 0)then
               return
            else if( ios .lt. 0 )then
               return
            end if
            write(ascii_packet_file_no,1500)
            write(ascii_packet_file_no,1550) ii_avg, ii_min, ii_max
         end if
      end do
c
c             J-values computed from stress intensity factors
c
      write(ascii_packet_file_no,1600)
      do lines_read = 1,num_domains
         read(packet_file_no,iostat=ios) nowring, j_pstrs, j_pstrn
         write(ascii_packet_file_no,1700) nowring, j_pstrs, j_pstrn
      end do
c
c
      return
c
 1100 format('domain: ',a,2x,'num: ',i6,2x,'structure: ',a8,2x,
     &       'loading: ',a,2x,'step_number: ',i12)
 1201 format('KI= 1  KII= 0  KIII= 0  plane stress',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KI pstrs',3x,'killed ele' )
 1202 format('KI= 1  KII= 0  KIII= 0  plane strain',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KI pstrn',3x,'killed ele' )
 1203 format('KI= 0  KII= 1  KIII= 0  plane stress',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KII pstrs',2x,'killed ele' )
 1204 format('KI= 0  KII= 1  KIII= 0  plane strain',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KII pstrn',2x,'killed ele' )
 1205 format('KI= 0  KII= 0  KIII= 1',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',6x,
     &       'KIII',4x,'killed ele' )
 1400 format(1x,i3,2x,10(1x,e11.4),2x,i3)
 1500 format(1x,'domain average',3x,'domain min',3x,'domain max')
 1550 format(1x,e11.4,5x,e11.4,2x,e11.4)
 1600 format('J-values from KI, KII and KIII',
     &       /,1x,'domain',2x,'J pstrs',5x,'J pstrn')
 1700 format(1x,i3,3x,e11.4,1x,e11.4)
c
c
      end
c
c **********************************************************************
c *                                                                    *
c *     subroutine: i_integral_t_packet                                *
c *                                                                    *
c *     reads packet file for interaction-integral t-stress results    *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
      subroutine i_integral_t_packet(type,num_lines,step,iter)

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c             dummy variables
c
      integer, intent(in) :: type,num_lines,step,iter
c
c             local variables
c
      integer :: ios, lines_read, i, j, ltmstp, nowring, skipped_killed,
     &           num_domains
      double precision  :: iiterms(11), ii_min1, ii_max1, ii_avg1,
     &                     ii_min2, ii_max2, ii_avg2
      character(len=24) :: domain_id
      character(len=8)  :: stname, lsldnm
      logical           :: connected
c
c
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  read from packet file and write to output
c
      read(packet_file_no,iostat=ios) num_domains, domain_id, stname,
     &                                lsldnm, ltmstp
      if( ios .gt. 0)then
         return
      else if( ios .lt. 0 )then
         return
      end if
      domain_id = ADJUSTL( domain_id )
      stname    = ADJUSTL( stname )
      lsldnm    = ADJUSTL( lsldnm )
      write(ascii_packet_file_no,1100)
     &  TRIM(domain_id),num_domains,TRIM(stname),TRIM(lsldnm),ltmstp
c
      do j=1,2
         if( j.eq.1 ) write(ascii_packet_file_no,1201)
         if( j.eq.2 ) write(ascii_packet_file_no,1202)
         do lines_read=1,num_domains
            read(packet_file_no,iostat=ios) nowring,
     &           (iiterms(i),i=1,11), skipped_killed
            if( ios .gt. 0)then
               return
            else if( ios .lt. 0 )then
               return
            end if
            write(ascii_packet_file_no,1400) nowring,
     &           (iiterms(i),i=1,11), skipped_killed
         end do
c
         if( num_domains.gt.1 ) then
            read(packet_file_no,iostat=ios) ii_avg1, ii_min1, ii_max1,
     &           ii_avg2, ii_min2, ii_max2
            if( ios .gt. 0)then
                return
            else if( ios .lt. 0 )then
               return
            end if
            write(ascii_packet_file_no,1500)
            write(ascii_packet_file_no,1550) ii_avg1, ii_min1, ii_max1
            write(ascii_packet_file_no,1550) ii_avg2, ii_min2, ii_max2
         end if
      end do
c
c
      return
c
 1100 format('domain: ',a,2x,'num: ',i6,2x,'structure: ',a,2x,
     &       'loading: ',a,2x,'step_number: ',i12)
 1201 format('T11, T33: plane stress',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',4x,
     &       'T11 pstrs',5x,'T33',5x,'killed ele' )
 1202 format('T11, T33: plane strain',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',4x,
     &       'T11 pstrn',5x,'T33',5x,'killed ele' )
 1400 format(1x,i3,2x,11(1x,e11.4),3x,i3)
 1500 format(1x,'domain average',3x,'domain min',3x,'domain max')
 1550 format(1x,e11.4,5x,e11.4,2x,e11.4)
c
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: step_patt_factors                                  *
c *                                                                    *
c *     reads packet file for accumulated step pattern factors         *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine step_patt_factors( type, num_lines, step, iter )

      use main_data, only: packet_file_no, ascii_packet_file_no
      implicit none
c
c             dummy variables
c
      integer, intent(in) :: type, num_lines, step, iter
c
c             local variables
c
      integer :: ios, i
      double precision  :: factor
      character(len=12) :: pattern
      logical           :: connected
c
      inquire( unit=ascii_packet_file_no, opened=connected)
      if( .not. connected ) call open_output_file

      call packet_header_SRT(type,num_lines,step,iter)
c
c                  read from packet file and write to output
c
      write(ascii_packet_file_no,8000) step
      do i = 1, num_lines
         read(packet_file_no,iostat=ios) pattern, factor
         if( ios .gt. 0 ) then
           return
         else if( ios .lt. 0 ) then
           return
         end if
         pattern = adjustl( pattern )
         write(ascii_packet_file_no,8010) pattern, factor
      end do
      write(ascii_packet_file_no,*) ' '
c
      return
c
 8000 format('>>> Accumulated loading pattern factors for step: ',
     & i6 )
 8010 format(8x,a12,1x,f12.3)
c
      end
c
c **********************************************************************
c *                                                                    *
c *     subroutine: packet_header_SRT                                  *
c *                                                                    *
c *     write the packet header line to the text file;                 *
c *     the header line in the text contains the packet type ID        *
c *     number, the number of data lines, the current step and iter,   *
c *     and the current output set number; these packet header values  *
c *     can be used by the routine reading the text packet file.       *
c *                                                                    *
c *     created:  8/25/2005 GVT                                        *
c *                                                                    *
c **********************************************************************
c
      subroutine packet_header_SRT(packet_type,num_lines,step,iter)
c
c Declare Modules
c
      use main_data, only: ascii_packet_file_no
      use pvars, only: saved_step
c
c Declare Variables
c
      implicit none
      integer, intent(in)::packet_type,num_lines,step,iter
c
c Declare Local Variables
c
      character(LEN=80)word_packet,word_num,word_step,word_iter,
     &  word_saved_step
c
c Write the packet header line values to the text file
c
      write(word_packet,*)packet_type
      word_packet = ADJUSTL(word_packet)

      write(word_num,*)num_lines
      word_num = ADJUSTL(word_num)

      write(word_step,*)step
      word_step = ADJUSTL(word_step)

      write(word_iter,*)iter
      word_iter = ADJUSTL(word_iter)

      write(word_saved_step,*)saved_step
      word_saved_step = ADJUSTL(word_saved_step)

CC      write(ascii_packet_file_no,700)
CC     &  packet_type,num_lines,step,iter,saved_step

      write(ascii_packet_file_no,701)TRIM(word_packet),TRIM(word_num),
     &  TRIM(word_step),TRIM(word_iter),TRIM(word_saved_step)
c
c Format statements
c
c            try to keep the packet header line within 80 characters
c
CC  700 format(/,'packet:',i3,'  num_lines:',i10,'  step:',i8,
CC     &         '  iter:',i10,'  output_set:',i8)
  701 format(/,'packet: ',a,'  num_lines: ',a,'  step: ',a,
     &         '  iter: ',a,'  output_set: ',a)

      end

