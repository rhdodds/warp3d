c **********************************************************************
c *                                                                    *
c *     module: pvars                                                  *
c *                                                                    *
c *     contains global variables                                      *
c *                                                                    *
c *     Use saved_dp1, saved_dp2, saved_dp3, saved_step, saved_int1,   *
c *     saved_int2, and saved_int3, when data from one packet          *
c *     subroutine is needed in following packet subroutines.          *
c *     If more than 3 double precision or integer variables are       *
c *     needed, add additional variables to this module.               *
c *                                                                    *
c **********************************************************************
c
c
      module pvars
      implicit none
c
c
c
      double precision saved_dp1, saved_dp2, saved_dp3
      integer packet_file_no /1/, screen /6/, key /5/,
     &        results_file_no /2/, saved_step, saved_int1,
     &                             saved_int2, saved_int3
      character * 80 packet_file_name, results_file_name
      logical ou_packet_control(30), debug
c
c
c
c
      end module pvars
c
c
c
c **********************************************************************
c *                                                                    *
c *     program: packet_reader                                         *
c *                                                                    *
c *     Main program controlling packet reading.  Calls appropriate    *
c *     subroutines to operate on packets when needed, or skips over   *
c *     unwanted packets.                                              *
c *                                                                    *
c **********************************************************************
c
c
c
      program packet_reader
      use pvars
      implicit none
c
      integer ios,type,num_lines,step,iter,lines_read,counter
c
c
      counter = 0
c
c
c                   call subroutines to initialize variables,
c                   open packet file, and obtain the desired
c                   packet types from the user.
c
c
      call init
      call get_packet_type
c
c
c                   infinite loop, reading packets until the end of
c                   file is reached.
c
c
      debug = .false.
      do
c
c                   read packet header line.
c
c
         read(unit=packet_file_no,iostat=ios)type,num_lines,step,iter
c
c
c                   check for errors in reading the packet file,
c                   or for the end of file, terminate program when
c                   appropriate.
c
c
         if( debug ) write(screen,*)type,num_lines,step,iter,ios
c
         if( ios .gt. 0)then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9000)
            stop
         else if( ios .lt. 0 ) then  ! -1 means eof on prior read
            if( counter .eq. 0 ) write(screen,9200)
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9100)
            stop
         end if
c
c                   check to see if the current packet type should
c                   be analyzed.  ou_packet_control(type) will be set
c                   to true when the user specifies that type to be
c                   processed.
c
         if( ou_packet_control ( type ) ) then
c
c
c                   keep track of the number of times the desired
c                   packets are found, message will be sent if
c                   no desired packets are found on the packet file.
c
c
            counter = counter + 1
c
c                   case(#) need to be added as packets are added
c
c
            select case(type)
c
               case(1)
                  call displ_node_packet(type,num_lines,step,iter)
c
               case(2)
                  call veloc_node_packet(type,num_lines,step,iter)
c
               case(3)
                  call accel_node_packet(type,num_lines,step,iter)
c
               case(4)
                  call reacts_node_packet(type,num_lines,step,iter)
c
               case(5)
                  call gurson_elem_1(type,num_lines,step,iter)
c
               case(6)
                  call gurson_elem_2(type,num_lines,step,iter)
c
               case(7)
                  call cohes_packet(type,num_lines,step,iter)
c
               case(8)
                  call ctoa_type_3_packet(type,num_lines,step,iter)
c
               case(9)
                  call ctoa_type_4_packet(type,num_lines,step,iter)
c
               case(10)
                  call cohesive_elem_killed(type,num_lines,step,iter)
c
               case(11)
                  call displ_elem_packet(type,num_lines,step,iter)
c
               case(12)
                  call veloc_elem_packet(type,num_lines,step,iter)
c
               case(13)
                  call accel_elem_packet(type,num_lines,step,iter)
c
               case(14)
                  call reacts_elem_packet(type,num_lines,step,iter)
c
               case(15)
                  call stress_elem_packet(type,num_lines,step,iter)
c
               case(16)
                  call strain_elem_packet(type,num_lines,step,iter)
c
               case(17)
                  call j_integral_packet(type,num_lines,step,iter)
c
               case(18)
                  call ctoa_type_1_packet(type,num_lines,step,iter)
c
               case(19)
                  call ctoa_type_2_packet(type,num_lines,step,iter)
c
               case(20)
                  call gurson_elem_killed(type,num_lines,step,iter)
c
               case(21)
                  call smcs_elem_killed(type,num_lines,step,iter)
c
               case(22)
                  call smcs_elem_1(type,num_lines,step,iter)
c
               case(23)
                  call smcs_elem_2(type,num_lines,step,iter)
c
               case(24)
                  call step_patt_factors(type,num_lines,step,iter)
c
               case(27)
                  call i_integral_sif_packet(type,num_lines,step,iter)
c
               case(28)
                  call i_integral_t_packet(type,num_lines,step,iter)
c
               case(29)
                  call temp_nodes_packet(type,num_lines,step,iter)
c
               case(30)
                  call temp_elem_packet(type,num_lines,step,iter)
c
               case(31)
                  call cohesive_tractions_packet(type,num_lines,
     &                         step,iter)
c
               case(32)
                  call cohesive_disp_jumps_packet(type,
     &                   num_lines,step,iter)
            end select
c
         else
c
c
c                   skip unwanted packets.
c
c
            do lines_read=1,num_lines
               read(unit=packet_file_no,iostat=ios)
     &                           type,num_lines,step,iter
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  if( counter .eq. 0 ) write(screen,9200)
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
            end do
c
         end if
c
      end do
c
c
      call closer(packet_file_no)
      call closer(results_file_no)
      stop
c
 9000 format(/,1x,'>> ERROR: reading packet file...terminating m',/)
 9100 format(/,1x,'>> END OF FILE...program terminating in main',/)
 9200 format(/,1x,'>> ERROR: the desired packet type can not be',
     &       /,1x,'          found in the given packet file.')
c
c
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: init                                               *
c *                                                                    *
c *     initializes variables, gets packet file name from user, and    *
c *     opens the file.                                                *
c *                                                                    *
c **********************************************************************
c
c
      subroutine init
      use pvars, only: screen,key,packet_file_no,
     &        packet_file_name,ou_packet_control,debug
      implicit none
      integer ios
c
c
c                   ou_packet_control is an array of logical
c                   variables, which are set to true when the
c                   user specifies that a type of packet should
c                   be analyzed.  all components of the array
c                   are set to false initially.
c
      ou_packet_control = .false.
c
      debug = .false.
c
c
      write(screen,9000)
c
c
c                   loop until a valid packet file is given
c                   and opened.
c
c
      do
        write(screen,9010)
        read(key,9020) packet_file_name
        call stripf( packet_file_name )
        if ( packet_file_name(1:1) .eq. ' ' ) then
            packet_file_name(1:) = 'packet.out.1'
        end if
        open( unit=packet_file_no, file=packet_file_name, status='old',
     &        access='sequential', form='unformatted',
     &        position='rewind', iostat=ios)
        if( ios .ne. 0)then
           write(screen,9030)
        else
           write(screen,9040)
           exit
        end if
      end do
c
      return
c
c
 9000 format(/,
     &       ' ****************************************************',
     &  /,   ' *                                                  *',
     &  /,   ' *   WARP3D Packet Reader (example implementation)  *',
     &  /,   ' *                                                  *',
     &  /,   ' *             Last updated: 5-23-2013              *',
     &  /,   ' *                                                  *',
     &  /,   ' *                                                  *',
     &  /,   ' ****************************************************'
     &   )
 9010 format(/,1x, '>> WARP3D binary packet file',/,
     &               ' >>   (default: packet.out.1) ? ',$)
 9020 format(a80)
 9030 format(/,1x, '>> ERROR: Invalid packet file')
 9040 format(/,1x, '>> Packet file opened',/)
c
      end
c
c **********************************************************************
c *                                                                    *
c *     subroutine: closer                                             *
c *                                                                    *
c *     closes both the packet file and the output file.               *
c *                                                                    *
c **********************************************************************
c
      subroutine closer(file)
      use pvars, only: screen
      implicit none
      integer, intent(in):: file
c
      logical connected
c
c
      inquire( unit=file, opened=connected)
      if( connected )then
        close(unit=file, status='keep')
        if(file .eq. 1)write(screen,1000)
        if(file .eq. 2)write(screen,2000)
      end if
      return
 1000 format(/,1x,'>> Packets file closed')
 2000 format(/,1x,'>> Results file closed')
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- s t r i p f                      *
c *                                                               *
c *****************************************************************
c
      subroutine stripf( string )
      implicit integer (a-z)
      character *(*) string, copy*256
c
c             strip leading blanks from string and return it.
c
      last = len(string)
      copy(1:last) = string(1:last)
      do 100 cpos = 1, last
       if ( string(cpos:cpos) .ne. ' ' ) go to 200
 100  continue
      return
 200  continue
      string(1:) = copy(cpos:last)
      return
      end
c
c **********************************************************************
c *                                                                    *
c *     subroutine: get_packet_type                                    *
c *                                                                    *
c *     collects all of the desired packet types from the user.        *
c *                                                                    *
c **********************************************************************
c
      subroutine get_packet_type
      use pvars, only: screen,key,packet_file_no,
     &                 ou_packet_control
      implicit none
c
c
c
      character * 1 y_or_n
      integer output_packet_no
      logical additional_packets_needed, valid_packet_no
c
c
      additional_packets_needed = .true.
      valid_packet_no = .true.
c
c
c                   loop asks user what type of packet is
c                   desired, repeats until all deired packets
c                   have been recorded.
c
c
      do
         write(screen,9000)
         valid_packet_no=.true.
         do
           read(key,9100)output_packet_no
           if(output_packet_no .gt. 32 .or.
     &        output_packet_no .lt. 0)then
                write(screen,9200)
                write(screen,9300)
           else
                call not_a_packet_yet(output_packet_no,
     &                                valid_packet_no)
                if(valid_packet_no)then
                   write(screen,9300)
                else
                   exit
                end if
           end if
         end do
c
         ou_packet_control(output_packet_no)=.true.
c
         write(screen,9400)
         read(key,9500)y_or_n
         if(y_or_n .eq. 'N' .or.
     &      y_or_n .eq. 'n' .or.
     &      y_or_n .eq. ' ')exit
c
      end do
c
      return
c
 9000 format(/,1x,'        Packet Types/Numbers              ',/,
     &         1x,'                                          ',/,
     &         1x,'          Type                    Number  ',/,
     &         1x,'   -------------------            ------  ',/,
     &         1x,'   nodal displacements               1    ',/,
     &         1x,'   nodal velocities                  2    ',/,
     &         1x,'   nodal accelerations               3    ',/,
     &         1x,'   nodal reactions                   4    ',/,
     &         1x,'   gurson elm status                 5    ',/,
     &         1x,'   gurson elm status (adap ld cntrl) 6    ',/,
     &         1x,'   cohesive elem status              7    ',/,
     &         1x,'   ctoa constant front step summary  8    ',/,
     &         1x,'   ctoa front step summary           9    ',/,
     &         1x,'   killed cohesive element          10    ',/,
     &         1x,'   element displacements            11    ',/,
     &         1x,'   element velocities               12    ',/,
     &         1x,'   element acclerations             13    ',/,
     &         1x,'   element reactions                14    ',/,
     &         1x,'   element stresses                 15    ',/,
     &         1x,'   element strains                  16    ',/,
     &         1x,'   domain J-integral results        17    ',/,
     &         1x,'   ctoa non-constant front release  18    ',/,
     &         1x,'   ctoa constant front release      19    ',/,
     &         1x,'   killed gurson element            20    ',/,
     &         1x,'   killed stress-mod crt strain elm 21    ',/,
     &         1x,'   smcs elem status                 22    ',/,
     &         1x,'   smcs elem status (auto ld red)   23    ',/,
     &         1x,'   accumulated step pattern factors 24    ',/,
     &         1x,'   I-integral SIF results           27    ',/,
     &         1x,'   I-integral T-stress results      28    ',/,
     &         1x,'   nodal temperatures (node list)   29    ',/,
     &         1x,'   temperatures at nodes of element 30    ',/,
     &         1x,'   interface tractions of element   31    ',/,
     &         1x,'   interface displ jumps of element 32    ',/,
     &       /,1x,'>> Enter desired packet number: ',$)
 9100 format(i3)
 9200 format(/,1x,'>> ERROR: Invalid packet number',
     &       /,1x,'>>        Pick a valid packet number')
 9300 format(/,1x,'>> Enter desired packet number: ',$)
 9400 format(/,1x,'>> Are additional packets desired?',
     &            ' (y or n) (n is the default): ',$)
 9500 format(a1)
      end
c
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: not_a_packet_yet                                   *
c *                                                                    *
c *     checks to see if the desired packet has been added to this     *
c *     program.  as packet types are added, the packet number is      *
c *     removed from this list.                                        *
c *                                                                    *
c **********************************************************************
c
c
c
c
      subroutine not_a_packet_yet(number,control)
      use pvars, only: screen
      implicit none
c
      integer, intent(in)::number
      logical, intent(inout)::control
c
c
c
c             when a new packet is added, remove the number from
c             this check
c
c
      if(
     &    number .eq. 25                 .or.
     &    number .eq. 26                 .or.
     &    number .eq. 33                 .or.
     &    number .eq. 34                 .or.
     &    number .eq. 35
     &        )then
c
          write(screen,1000)number
c
      else
          control=.false.
      end if
c
      return
c
c
c
 1000 format(/,1x,'>> ERROR: Packet number ',i2,
     &          ' has not been implimented yet. Pick again')
c
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
      use pvars, only: screen, results_file_name, results_file_no, key
      implicit none
c
      integer ios
c
      do
        write(screen,9010)
        read(key,9020) results_file_name
        call stripf( results_file_name )
        if ( results_file_name(1:1) .eq. ' ' ) then
            results_file_name(1:) = 'results.out.1'
        end if
        open( unit=results_file_no, file=results_file_name,
     &        status='replace',iostat=ios)
        if( ios .ne. 0)then
           write(screen,9030)
        else
           write(screen,9040)
           exit
        end if
      end do
c
      return
c
c
c
 9010 format(/,1x, '>> Name of results file',/,
     &               ' >>   (default: results.out.1)? ',$)
 9020 format(a80)
 9030 format(/,1x, '>> ERROR: opening results file')
 9040 format(/,1x, '>> Results file opened'  )
c
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  need to read and put in file
c
      write(results_file_no,1100) step
      read(packet_file_no,iostat=ios) crit_angle, init_crit_angle,
     &                                ctoa_range
      if( ios .gt. 0)then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
      else if( ios .lt. 0 )then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
      end if
c
      do lines_read = 1, num_lines-1
       read(packet_file_no,iostat=ios) node, angle
       if( ios .gt. 0)then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
       else if( ios .lt. 0 )then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
       end if
       write(results_file_no,9010) node, angle
      end do
c
      write(results_file_no,9020) crit_angle, init_crit_angle,
     &                                ctoa_range
      return
c
 9010 format(  '        node: ',i6,'  CTOA: ',f6.3,' (degrees)')
 9020 format('        * critical CTOA:',
     &       /,'            for growth    : ',f6.3,' (degrees)',
     &       /,'            for initiation: ',f6.3,' (degrees)',
     &       /,'            nodes released within ',f4.1,' % of ',
     &                        'critical angle.' )
 1100 format(
     & /,' >>> CTOA Node Release for Non-Constant Front:',
     & /,'        Before step: ',i5 )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in cohes')
c
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  need to read and put in file
c
      write(results_file_no,1100) step
      read(packet_file_no,iostat=ios) crit_angle, init_crit_angle,
     &                                ctoa_range
      if( ios .gt. 0)then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
      else if( ios .lt. 0 )then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
      end if
c
      do lines_read = 1, num_lines-1
       read(packet_file_no,iostat=ios) node, angle
       if( ios .gt. 0)then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
       else if( ios .lt. 0 )then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
       end if
       if ( node .gt. 0 ) then
         write(results_file_no,9010) node, angle
       else
         write(results_file_no,9012) -node, angle
       end if
      end do
c
      write(results_file_no,9020) crit_angle, init_crit_angle,
     &                                ctoa_range
      return
c
 9010 format(  '        master node:  ',i6,'  CTOA: ',f6.3,' (degrees)')
 9012 format(  '        interim node: ',i6,'  CTOA: ',f6.3,' (degrees)')
 9020 format('        * critical CTOA:',
     &       /,'            for growth    : ',f6.3,' (degrees)',
     &       /,'            for initiation: ',f6.3,' (degrees)',
     &       /,'            nodes released within ',f4.1,' % of ',
     &                        'critical angle.' )
 1100 format(
     & /,' >>> CTOA Node Release for Constant Front Procedure:',
     & /,'        Before step: ',i5 )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in cohes')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
c
      integer, intent(in) :: type, num_lines, step, iter
c
      integer ios, lines_read, i, node, status
      logical connected
      double precision  angle, angle_fract
c
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  need to read and put in file
c
      write(results_file_no,1100) step
      write(results_file_no,*) 'master node       CTOA(degrees)',
     &     '  CTOA/critical CTOA'
      write(results_file_no,*) '-----------       -------------',
     &     '  ------------------'
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) node, angle, status, angle_fract
       if( ios .gt. 0)then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
       else if( ios .lt. 0 )then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
       end if
       if ( status .eq. 0 ) then
         write(results_file_no,9010) node, angle, '(*)', angle_fract
       else
         write(results_file_no,9010) node, angle, '   ', angle_fract
       end if
      end do
c
      write(results_file_no,9200)
      return
c
 9010 format(6x,i5,7x,f7.4,a3,11x,f4.2,3x,2h<=)
 1100 format(
     & /,' >>> CTOA Step Summary for Constant Front:',
     & /,'        After step: ',i5 )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in cohes')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
c
      integer, intent(in) :: type, num_lines, step, iter
c
      integer ios, lines_read, i, node, status, neighbor
      logical connected
      double precision  angle, angle_fract
      character * 1 marker
c
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  need to read and put in file
c
      write(results_file_no,1100) step
      write(results_file_no,*)
     &     'crack front node  neighbor  CTOA(degrees)',
     &     '  CTOA/critical CTOA'
      write(results_file_no,*)
     &     '----------------  --------  -------------',
     &     '  ------------------'
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) node, neighbor, angle,
     &                                 status, angle_fract
       if( ios .gt. 0)then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
       else if( ios .lt. 0 )then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
       end if
       marker = ' '
       if ( status .eq. 0 ) marker = '*'
       write(results_file_no,9010) node, neighbor, angle, marker,
     &                             angle_fract
      end do
c
      write(results_file_no,9200)
      return
 9010 format(6x,i5,9x,i5,7x,f7.4,a1,11x,f4.2,3x,2h<=)
 1100 format(
     & /,' >>> CTOA Step Summary for Non-Constant Front:',
     & /,'        After step: ',i5 )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in cohes')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
      if( .not. connected ) write(results_file_no,1100)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(gurson(i), i=1,5),
     &                   special_1,gurson(6),special_2
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
               if ( special_1 .eq. 0 .and. special_2 .eq. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    '*',gurson(6),'*'
               if ( special_1 .ne. 0 .and. special_2 .eq. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    ' ',gurson(6),'*'
               if ( special_1 .eq. 0 .and. special_2 .ne. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    '*',gurson(6),' '
               if ( special_1 .ne. 0 .and. special_2 .ne. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    ' ',gurson(6),' '
            end do
c
      return
c
 1000 format(1x,i5,1x,i6,1x,f12.5,3x,f12.5,3(1x,e14.6),
     &                                   a1,(1x,e14.6),a1)

 1100 format('  step  element     inital f     current f       Ep    ',
     &       '      sigma bar      mean stress    ',
     &       'mises stress')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in gur_1')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
      if( .not. connected ) write(results_file_no,1100)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(gurson(i), i=1,5),
     &              special_1,gurson(6),special_2,gurson(7)
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
               if ( special_1 .eq. 0 .and. special_2 .eq. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    '*',gurson(6),'*',gurson(7)
               if ( special_1 .ne. 0 .and. special_2 .eq. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    ' ',gurson(6),'*',gurson(7)
               if ( special_1 .eq. 0 .and. special_2 .ne. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    '*',gurson(6),' ',gurson(7)
               if ( special_1 .ne. 0 .and. special_2 .ne. 0 )
     &          write(results_file_no,1000)step,elem,(gurson(i), i=1,5),
     &                    ' ',gurson(6),' ',gurson(7)
            end do
c
      return
c
 1000 format(1x,i5,1x,i6,1x,f12.5,3x,f12.5,3(1x,e14.6),
     &                a1,(1x,e14.6),a1,1x,f10.6 )

 1100 format('  step  element     inital f     current f       Ep    ',
     &       '      sigma bar      mean stress    ',
     &       'mises stress     delta f ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in gur_2')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
      if( .not. connected ) write(results_file_no,1100)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(gurson(i), i=1,4)
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
              write(results_file_no,1000)step,elem,
     &                                 (gurson(i), i=1,4)
            end do
c
      return
c
 1000 format(1x,i5,1x,i6,1x,f12.5,1x,f12.5,3x,2(1x,e14.6))
 1100 format('  step  element   initial f    final f ',
     &       '    mean stress    ',
     &       'mises stress  ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating gur3')
 9100 format(/,1x,'>> END OF FILE...program terminating in gur_3')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
      if( .not. connected ) write(results_file_no,1100)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(smcs(i), i=1,4)
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
               write(results_file_no,1000)step,elem,(smcs(i), i=1,4)
            end do
c
      return
c
 1000 format(1x,i5,1x,i6,4(1x,e14.6))

 1100 format('  step  element   eps-pls        eps-crit      ',
     &       'sig-mean       sig-mises')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in smcs_1')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
      if( .not. connected ) write(results_file_no,1100)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(smcs(i), i=1,5)
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
              write(results_file_no,1000)step,elem,(smcs(i), i=1,5)
            end do
c
      return
c
 1000 format(1x,i5,1x,i6,5(1x,e14.6))

 1100 format('  step  element    eps-pls        eps-crit      ',
     &       'sig-mean       sig-mises       d(eps-pls)'  )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in smcs_2')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
      if( .not. connected ) write(results_file_no,1100)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               read(packet_file_no,iostat=ios)elem,(smcs(i), i=1,2)
               if( ios .gt. 0)then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9000)
                  stop
               else if( ios .lt. 0 )then
                  call closer(packet_file_no)
                  call closer(results_file_no)
                  write(screen,9100)
                  stop
               end if
              write(results_file_no,1000)step,elem,
     &                     (smcs(i), i=1,2)
            end do
c
      return
c
 1000 format(1x,i5,2x,i6,5x,f8.6,4x,f8.6)
 1100 format('  step  element   pl-strain  pl-strain limit ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in smcs_k')
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
      subroutine cohes_packet( type, num_lines, step, iter )
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
      integer, intent(in)::type, num_lines, step, iter
c
c                  local declarations
c                  ------------------
c
      integer ios, lines_read, elem, i, cohes_type,
     &        num_gp_post_peak_normal, num_gp_post_peak_shear
      double precision cohes(6)
      character * 3 special_char
      logical connected, exp1_intf, ppr
c
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  determine which cohesive option we have
c                  (exp1_intf or ppr - they cannot be mixed)
c                  write a header line above element list and
c                  status values
c
       read(packet_file_no,iostat=ios) elem, cohes_type
       exp1_intf = cohes_type .eq. 4
       ppr       = cohes_type .eq. 6
       if( exp1_intf ) write(results_file_no,1100)
       if( ppr       ) write(results_file_no,1110)
       backspace( packet_file_no )
c
c                  loop over the number of cohesive elements with
c                  results for this step.
c
       do lines_read = 1, num_lines
c
         if( exp1_intf ) then
            read(packet_file_no,iostat=ios) elem, cohes_type,
     &             cohes(1:6), special_char
            if( ios .gt. 0 .or. ios .lt. 0 ) go to 8000
            write(results_file_no,1000) step, elem, cohes(1:6),
     &                                  special_char
         end if
c
         if( ppr ) then
            read(packet_file_no,iostat=ios) elem, cohes_type,
     &             cohes(1:6), num_gp_post_peak_normal,
     &             num_gp_post_peak_shear
            if( ios .gt. 0 .or. ios .lt. 0 ) go to 8000
            write(results_file_no,1010) step, elem, cohes(1:6),
     &                                  num_gp_post_peak_normal,
     &                                  num_gp_post_peak_shear
         end if
      end do
c
      return
c
 8000 continue
      if( ios .gt. 0 )then
        call closer(packet_file_no)
        call closer(results_file_no)
        write(screen,9000)
        stop
      else if( ios .lt. 0 )then
       call closer(packet_file_no)
       call closer(results_file_no)
       write(screen,9100)
       stop
      end if

c
 1000 format(1x,i5,1x,i6,1x,6e14.6,1x,a3)
 1010 format(1x,i5,1x,i6,e14.6,1x,e14.6,e14.6,e14.6,f11.5,
     &       f13.5,3x,'(',i1,',',i1,')')
 1100 format('>>> interface-cohesive element status: exp1_intf',
     & /,'  step  elem        Tn             Ts            ',
     &       'Dn           Ds          Teff          Deff   ')
 1110 format('>>> interface-cohesive element status: ppr',
     & /, '  step  elem    T-norm         T-shear        D-norm   ',
     &       '    D-shear     Tn/Tn-peak   ',
     &       ' Ts/Ts-peak   ' )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in cohes')
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
      subroutine cohesive_elem_killed( type, num_lines, step, iter )
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
      integer, intent(in)::type,num_lines,step,iter
c
c                  local declarations
c                  ------------------
c
      integer ios, lines_read, elem,i,cohes_type
      double precision cohes(2)
      logical connected, write_header
c
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  there can be multiple killed cohesive elements during the
c                  just completed step (num_lines>1). loop to read each
c                  line, handle each coheice model formulation,
c                  print a header before start of element list.
c
      write_header = .true.
      do lines_read=1,num_lines
        read(packet_file_no,iostat=ios) elem, cohes(1), cohes(2),
     &                                  cohes_type
        if( ios .gt. 0 ) then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9000)
           stop
        else if( ios .lt. 0 ) then
           call closer(packet_file_no)
           call closer(results_file_no)
           write(screen,9100)
           stop
        end if
c
        if( cohes_type .eq. 4 ) then
           if( write_header ) then
            write(results_file_no,1100)
            write_header = .false.
           end if
           write(results_file_no,1000) step, elem, cohes(1), cohes(2)
        end if
        if( cohes_type .eq. 6 ) then
           if( write_header ) then
            write(results_file_no,1110)
            write_header = .false.
           end if
           write(results_file_no,1000) step, elem, cohes(1), cohes(2)
        end if
      end do
c
      return
c
 1000 format(1x,i5,1x,i6,6x,f5.3,8x,f5.3 )
 1100 format('>>> interface-cohesive newly killed elements: exp1_intf',
     & /,'  step  elem    Deff/Dpeak   Teff/Tpeak ')
 1110 format('>>> interface-cohesive newly killed elements: ppr',
     & /,'  step  elem    Tn/Tn-peak   Tt/Tt-peak ' )
 9000 format(/,1x,'>> ERROR: reading packet file...terminating')
 9100 format(/,1x,'>> END OF FILE...program terminating in cohes_kill')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               if( lines_read .lt. num_lines ) then
                  read(packet_file_no,iostat=ios)node,
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,node,
     &                                     (reacts(i), i=1,3)
               else
c
c                             the last entry on a reactions packet
c                             contains the TOTALS data
c
                  read(packet_file_no,iostat=ios)
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1100)step,
     &                                     (reacts(i), i=1,3)
               end if
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1100 format(1x,i5,1x,' TOT ',1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step   node  x-reaction    y-reaction    z-reaction   ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in reac')
 9100 format(/,1x,'>> END OF FILE...program terminating in reacts')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
               if( lines_read .lt. num_lines ) then
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,elem,node,
     &                                     (reacts(i), i=1,3)
               else
c
c                             the last entry on a reactions packet
c                             contains the TOTALS data
c
                  read(packet_file_no,iostat=ios)
     &                                    (reacts(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1100)step,
     &                                     (reacts(i), i=1,3)
               end if
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1100 format(1x,i5,1x,' TOT ELEM  ',1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-reaction   y-reaction    z-reaction   ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in reac')
 9100 format(/,1x,'>> END OF FILE...program terminating in reacts')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)node,
     &                                    (displ(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,node,
     &                                     (displ(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step   node   x-disp        y-disp       z-disp   ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in disp')
 9100 format(/,1x,'>> END OF FILE...program terminating in disp')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (displ(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,elem,node,
     &                                     (displ(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-disp      y-disp       z-disp     ')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in disp')
 9100 format(/,1x,'>> END OF FILE...program terminating in disp')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)node,
     &                                    (veloc(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,node,
     &                                     (veloc(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step   node   x-veloc        y-veloc       z-veloc   ')
 9000 format(/,1x,'>> ERROR: reading packet file..terminating in veloc')
 9100 format(/,1x,'>> END OF FILE...program terminating in veloc')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (veloc(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,elem,node,
     &                                     (veloc(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-veloc      y-veloc       z-veloc     ')
 9000 format(/,1x,'>> ERROR: reading packet file..terminating in veloc')
 9100 format(/,1x,'>> END OF FILE...program terminating in veloc')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)node,
     &                                    (accel(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,node,
     &                                     (accel(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step   node   x-accel        y-accel       z-accel   ')
 9000 format(/,1x,'>> ERROR: reading packet file..terminating in accel')
 9100 format(/,1x,'>> END OF FILE...program terminating in accel')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (accel(i), i=1,3)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,elem,node,
     &                                     (accel(i), i=1,3)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,
     &                      e14.6,
     &                      e14.6,
     &                      e14.6)
 1200 format('  step element  node',
     &       '  x-accel      y-accel       z-accel     ')
 9000 format(/,1x,'>> ERROR: reading packet file..terminating in accel')
 9100 format(/,1x,'>> END OF FILE...program terminating in accel')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
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
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,elem,node,
     &                                     (stress(i), i=1,26)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6 ,/, 19x, 15e14.6)
 1200 format('  step elem  node/gpt',
     &       '    sx            sy            sz          ',
     &       '    sxy           syz           sxz         Uo      ',
     &       '    svm            c1             c2           c3      ',
     & /,21x,
     &    '    I1            I2            I3           sig-1    ',
     &    '    sig-2         sig-3         L1            M1      ',
     &    '      N1            L2            M2            N2      ',
     &    '      L3            M3           N3')
 9000 format(/,1x,'>> ERROR: reading packet file..terminating in stres')
 9100 format(/,1x,'>> END OF FILE...program terminating in stress')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected )call open_output_file
      if( .not. connected )write(results_file_no,1200)
c
c                  need to read and put in file
c
c
c
            do lines_read=1,num_lines
c
c                          'node' will contain one of the following
c
c                           - nodes number (nodpts warp3d elem logical)
c                           - gauss points (gausspts warp3d elem logical)
c                           - node = 0 (center_output warp3d elem logical)
c
c
                  read(packet_file_no,iostat=ios)elem,node,
     &                                    (strain(i), i=1,22)
                  if( ios .gt. 0)then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9000)
                     stop
                  else if( ios .lt. 0 )then
                     call closer(packet_file_no)
                     call closer(results_file_no)
                     write(screen,9100)
                     stop
                  end if
                  write(results_file_no,1000)step,elem,node,
     &                                     (strain(i), i=1,7)
            end do
c
c
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6,e14.6,e14.6,
     &                      e14.6 )
 1200 format('  step   elem node/gpt',
     &       '    ex            ey            ez            ',
     &       'gxy           gyz           gxz         eff    ')
 9000 format(/,1x,'>> ERROR: reading packet file..terminating in strn')
 9100 format(/,1x,'>> END OF FILE...program terminating in strain')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  read from packet file and write to output
c
      read(packet_file_no,iostat=ios) num_domains, domain_id, stname,
     &                                lsldnm, ltmstp
      if( ios .gt. 0)then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9000)
         stop
      else if( ios .lt. 0 )then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9100)
         stop
      end if
      domain_id = ADJUSTL( domain_id )
      stname    = ADJUSTL( stname )
      lsldnm    = ADJUSTL( lsldnm )
      write(results_file_no,1100) domain_id(1:24), stname(1:8),
     &                            lsldnm(1:8), ltmstp
      write(results_file_no,1200)
c
c             read J-integral values for domain(s)
c
      do lines_read = 1, num_domains
         read(packet_file_no,iostat=ios) nowring,(diterms(i),i=1,9),
     &        skipped_killed
         if( ios .gt. 0)then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9200)
            stop
         else if( ios .lt. 0 )then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9100)
            stop
         end if
         write(results_file_no,1300) nowring, (diterms(i),i=1,9),
     &        skipped_killed
      end do
c
c             read max, min and average values
c
      if( num_domains .gt. 1 ) then
         read(packet_file_no,iostat=ios) di_avg, di_min, di_max
         if( ios .gt. 0)then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9200)
            stop
         else if( ios .lt. 0 )then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9100)
            stop
         end if
         write(results_file_no,1500)
         write(results_file_no,1550) di_avg, di_min, di_max
      end if
c
c             read stress intensity factor values computed from
c             J-values.
c
      write(results_file_no,1600)
      do i = 1,num_domains
         read(packet_file_no,iostat=ios) nowring,(diterms(k),k=1,5)
         if( ios .gt. 0)then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9200)
            stop
         else if( ios .lt. 0 )then
            call closer(packet_file_no)
            call closer(results_file_no)
            write(screen,9100)
            stop
         end if
         write(results_file_no,1700) nowring, (diterms(k),k=1,5)
      end do
c
      return
c
 1100 format(/,'domain: ',a24,2x,'structure: ',a8,2x,'loading: ',a8,
     &       2x,'load step number: ',i12)
 1200 format(1x,'domain',5x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total J',
     &    2x,'killed ele' )
 1300 format(1x,i3,3x,9(1x,e11.4),3x,i3)
 1500 format(1x,'domain average',3x,'domain min',3x,'domain max')
 1550 format(1x,e11.4,5x,e11.4,2x,e11.4)
 1600 format(/,'stress intensity factors from J (single-mode loading):',
     &       /,1x,'domain',2x,'KI pstrs',4x,'KI pstrn',
     &       4x,'KII pstrs',3x,'KII pstrn',3x,'KIII' )
 1700 format(1x,i3,3x,5(1x,e11.4))
 9000 format(/,1x,'>> ERROR: reading packet file domain integral',
     &       ' labels...terminating in domain integral')
 9100 format(/,1x,'>> END OF FILE...program terminating in domain',
     &       ' integral')
 9200 format(/,1x,'>> ERROR: reading packet file domain integral',
     &       ' values...terminating in domain integral')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  read from packet file and write to output
c
      read(packet_file_no,iostat=ios) num_domains, domain_id, stname,
     &                                lsldnm, ltmstp
      if( ios .gt. 0)then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9000)
         stop
      else if( ios .lt. 0 )then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9100)
         stop
      end if
      domain_id = ADJUSTL( domain_id )
      stname    = ADJUSTL( stname )
      lsldnm    = ADJUSTL( lsldnm )
      write(results_file_no,1100) domain_id(1:24), stname(1:8),
     &                            lsldnm(1:8), ltmstp
c
      do j=1,5
         if( j.eq.1 ) write(results_file_no,1201)
         if( j.eq.2 ) write(results_file_no,1202)
         if( j.eq.3 ) write(results_file_no,1203)
         if( j.eq.4 ) write(results_file_no,1204)
         if( j.eq.5 ) write(results_file_no,1205)
         do lines_read = 1,num_domains
            read(packet_file_no,iostat=ios) nowring,
     &           (iiterms(i),i=1,10), skipped_killed
            if( ios .gt. 0)then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9200)
               stop
            else if( ios .lt. 0 )then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9100)
               stop
            end if
            write(results_file_no,1400) nowring, (iiterms(i),i=1,10),
     &           skipped_killed
         end do
c
         if( num_domains.gt.1 ) then
            read(packet_file_no,iostat=ios) ii_avg, ii_min, ii_max
            if( ios .gt. 0)then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9200)
               stop
            else if( ios .lt. 0 )then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9100)
               stop
            end if
            write(results_file_no,1500)
            write(results_file_no,1550) ii_avg, ii_min, ii_max
         end if
      end do
c
c             J-values computed from stress intensity factors
c
      write(results_file_no,1600)
      do lines_read = 1,num_domains
         read(packet_file_no,iostat=ios) nowring, j_pstrs, j_pstrn
         write(results_file_no,1700) nowring, j_pstrs, j_pstrn
      end do
c
c
      return
c
 1100 format(/,'domain: ',a24,2x,'structure: ',a8,2x,'loading: ',a8,
     &       2x,'load step number: ',i12)
 1201 format(/,'KI= 1  KII= 0  KIII= 0  plane stress',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KI pstrs',3x,'killed ele' )
 1202 format(/,'KI= 1  KII= 0  KIII= 0  plane strain',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KI pstrn',3x,'killed ele' )
 1203 format(/,'KI= 0  KII= 1  KIII= 0  plane stress',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KII pstrs',2x,'killed ele' )
 1204 format(/,'KI= 0  KII= 1  KIII= 0  plane strain',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',3x,
     &       'KII pstrn',2x,'killed ele' )
 1205 format(/,'KI= 0  KII= 0  KIII= 1',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',6x,
     &       'KIII',4x,'killed ele' )
 1400 format(1x,i3,2x,10(1x,e11.4),2x,i3)
 1500 format(1x,'domain average',3x,'domain min',3x,'domain max')
 1550 format(1x,e11.4,5x,e11.4,2x,e11.4)
 1600 format(/,1x,'J-values from KI, KII and KIII',
     &       /,1x,'domain',2x,'J pstrs',5x,'J pstrn')
 1700 format(1x,i3,3x,e11.4,1x,e11.4)
 9000 format(/,1x,'>> ERROR: reading packet file interaction integral',
     &       ' labels...terminating in interaction integral')
 9100 format(/,1x,'>> END OF FILE...program terminating in interaction',
     &       ' integral')
 9200 format(/,1x,'>> ERROR: reading packet file interaction integral',
     &       ' values...terminating in interaction integral')
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
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
c             dummy variables
c
      integer, intent(in) :: type,num_lines,step,iter
c
c             local variables
c
      integer, parameter :: tot_iiterms = 11
      integer :: ios, lines_read, i, j, ltmstp, nowring, skipped_killed,
     &           num_domains
      double precision  :: iiterms(tot_iiterms), ii_min1, ii_max1,
     &                     ii_avg1, ii_min2, ii_max2, ii_avg2
      character(len=24) :: domain_id
      character(len=8)  :: stname, lsldnm
      logical           :: connected

c
c
c
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  read from packet file and write to output
c
      read(packet_file_no,iostat=ios) num_domains, domain_id, stname,
     &                                lsldnm, ltmstp
      if( ios .gt. 0)then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9000)
         stop
      else if( ios .lt. 0 )then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9100)
         stop
      end if
      domain_id = ADJUSTL( domain_id )
      stname    = ADJUSTL( stname )
      lsldnm    = ADJUSTL( lsldnm )
      write(results_file_no,1100) domain_id(1:24), stname(1:8),
     &                            lsldnm(1:8), ltmstp
c
      do j=1,2
         if( j.eq.1 ) write(results_file_no,1201)
         if( j.eq.2 ) write(results_file_no,1202)
         do lines_read=1,num_domains
            read(packet_file_no,iostat=ios) nowring,
     &           (iiterms(i),i=1,tot_iiterms), skipped_killed
            if( ios .gt. 0)then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9200)
               stop
            else if( ios .lt. 0 )then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9100)
               stop
            end if
            write(results_file_no,1400) nowring,
     &           (iiterms(i),i=1,tot_iiterms), skipped_killed
         end do
c
         if( num_domains.gt.1 ) then
            read(packet_file_no,iostat=ios) ii_avg1, ii_min1, ii_max1,
     &           ii_avg2, ii_min2, ii_max2
            if( ios .gt. 0)then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9200)
               stop
            else if( ios .lt. 0 )then
               call closer(packet_file_no)
               call closer(results_file_no)
               write(screen,9100)
               stop
            end if
            write(results_file_no,1500)
            write(results_file_no,1550) ii_avg1, ii_min1, ii_max1
            write(results_file_no,1550) ii_avg2, ii_min2, ii_max2
         end if
      end do
c
c
      return
c
 1100 format(/,'domain: ',a24,2x,'structure: ',a8,2x,'loading: ',a8,
     &       2x,'load step number: ',i12)
 1201 format(/,'T11, T33: plane stress',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',4x,
     &       'T11 pstrs',5x,'T33',5x,'killed ele' )
 1202 format(/,'T11, T33: plane strain',
     &       /,1x,'domain',4x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &       9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total I',4x,
     &       'T11 pstrn',5x,'T33',5x,'killed ele' )
 1400 format(1x,i3,2x,11(1x,e11.4),3x,i3)
 1500 format(1x,'domain average',3x,'domain min',3x,'domain max')
 1550 format(1x,e11.4,5x,e11.4,2x,e11.4)
 9000 format(/,1x,'>> ERROR: reading packet file interaction integral',
     &       ' T-stress labels...terminating in interaction integral')
 9100 format(/,1x,'>> END OF FILE...program terminating in interaction',
     &       ' integral T-stress')
 9200 format(/,1x,'>> ERROR: reading packet file interaction integral',
     &       ' T-stress values...terminating in interaction integral')
c
c
      end
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
      use pvars, only: screen,key,packet_file_no,results_file_no
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
      inquire( unit=results_file_no, opened=connected)
      if( .not. connected ) call open_output_file
c
c                  read from packet file and write to output
c
      write(results_file_no,8000) step
      do i = 1, num_lines
         read(packet_file_no,iostat=ios) pattern, factor
         if( ios .gt. 0 ) then
           call closer( packet_file_no )
           call closer( results_file_no )
           write(screen,9000)
           stop
         else if( ios .lt. 0 ) then
           call closer( packet_file_no )
           call closer( results_file_no )
           write(screen,9100)
           stop
         end if
         pattern = adjustl( pattern )
         write(results_file_no,8010) pattern, factor
      end do
      write(results_file_no,*) ' '
c
      return
c
 8000 format('>>> Accumulated loading pattern factors for step: ',
     & i6 )
 8010 format(8x,a12,1x,f12.3)
 9000 format(/,1x,'>> ERROR: reading packet file for accumulated',
     &       ' load step patterns...terminating in step patterns')
 9100 format(/,1x,'>> END OF FILE...program terminating in step',
     &       ' patterns')
c
      end
c
c **********************************************************************
c *                                                                    *
c *     subroutine: temp_nodes_packet                                  *
c *                                                                    *
c *     reads packet file for nodal temperature packets                *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine temp_nodes_packet(type,num_lines,step,iter)
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision temp, rel_temp
      logical connected
c
      inquire( unit=results_file_no, opened=connected )
      if( .not. connected ) then
        call open_output_file
        write(results_file_no,1200)
      end if
c
c                  read and write to output file
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) node, temp, rel_temp
       if( ios .gt. 0) then
          call closer(packet_file_no)
          call closer(results_file_no)
          write(screen,9000)
          stop
       else if( ios .lt. 0 )then
          call closer(packet_file_no)
          call closer(results_file_no)
          write(screen,9100)
          stop
       end if
       write(results_file_no,1000) step, node, temp, rel_temp
      end do
c
      return
c
 1000 format(1x,i5,1x,i8,1x,2f15.3)
 1200 format('  step   node',13x,'T',10x,'T - Tref')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in disp')
 9100 format(/,1x,'>> END OF FILE...program terminating in disp')
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: temp_elem_packet                                   *
c *                                                                    *
c *     reads packet file for temps at element nodes                   *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine temp_elem_packet(type,num_lines,step,iter)
      use pvars, only: screen,key,packet_file_no,results_file_no
      implicit none
c
      integer, intent(in)::type,num_lines,step,iter
c
      integer ios, lines_read, elem,i,node
      double precision temp, rel_temp
      logical connected
c
      inquire( unit=results_file_no, opened=connected )
      if( .not. connected ) then
         call open_output_file
         write(results_file_no,1200)
      end if
c
c                  read and write to output file
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) elem,node,
     &                        temp, rel_temp
       if( ios .gt. 0)then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9000)
         stop
       else if( ios .lt. 0 )then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9100)
         stop
       end if
       write(results_file_no,1000) step, elem,node,
     &                             temp, rel_temp
      end do
c
      return
c
 1000 format(1x,i5,1x,i5,1x,i5,1x,2f15.3)
 1200 format('  step element  node',10x,'T',11x,'T - Tref')
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in disp')
 9100 format(/,1x,'>> END OF FILE...program terminating in disp')
c
      end

c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: cohesive_tractions_packet                          *
c *                                                                    *
c *     reads packet file for tractions of interface-cohesive          *
c *     elements
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine cohesive_tractions_packet( type,
     &     num_lines, step, iter)
      use pvars, only: screen, key, packet_file_no, results_file_no
      implicit none
c
      integer, intent(in):: type, num_lines, step, iter
c
c                  local declarations
c
      integer nvalues, lines_read, ios, elem, point, cohesive_type
      double precision values(20)
      logical connected
c
      inquire( unit=results_file_no, opened=connected )
      if( .not. connected ) call open_output_file
c
c                  read and write to output file
c
c                  effective reading of tractions must be
c                  aware of the cohesive model type to
c                  know number & type of data values
c
c                  set type below. may wish to alter based on
c                  element number.
c
c                  cohesive types:
c                    1 = linear-elastic   2,3,5 not used
c                    4 = exponential model 6 = ppr
c                    7 = creep
c
c
      cohesive_type = 6
c
      if( cohesive_type == 1 ) nvalues = 6
      if( cohesive_type == 4 ) nvalues = 8
      if( cohesive_type == 6 ) nvalues = 8
      if( cohesive_type == 7 ) nvalues = 0

      write(results_file_no,1100) cohesive_type, nvalues
      write(results_file_no,1200)
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) elem, point, values(1:nvalues)
       if( ios .gt. 0 ) then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9000)
         stop
       else if( ios .lt. 0 )then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9100)
         stop
       end if
       write(results_file_no,1000) step, elem, point,
     &                             values(1:nvalues)
      end do
c
      return
c
 1000 format(1x,i5,i8,i6,8f15.6)
 1100 format(/,1x,'Interface-cohesive tractions packet (type=31)....',
     & /,1x,'cohesive model type: ',i2,2x,'nvalues: ',i2 )
 1200 format('  step element  point',10x,20('-'),' values ',20('-'))
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in disp')
 9100 format(/,1x,'>> END OF FILE...program terminating in disp')
c
      end
c
c
c **********************************************************************
c *                                                                    *
c *     subroutine: cohesive_disp_jumps_packet                         *
c *                                                                    *
c *     reads packet file for displacement jumps of interface-cohesive *
c *     elements                                                       *
c *                                                                    *
c *     THIS SUBROUTINE SHOULD BE MODIFIED BY THE USER TO OUTPUT       *
c *     THE DESIRED DATA TO THE RESULTS FILE.                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine cohesive_disp_jumps_packet( type,
     &     num_lines, step, iter)
      use pvars, only: screen, key, packet_file_no, results_file_no
      implicit none
c
      integer, intent(in):: type, num_lines, step, iter
c
c                  local declarations
c
      integer nvalues, lines_read, ios, elem, point, cohesive_type
      double precision values(20)
      logical connected
c
      inquire( unit=results_file_no, opened=connected )
      if( .not. connected ) call open_output_file
c
c                  read and write to output file
c
c                  effective reading of displacement jumps must be
c                  aware of the cohesive model type to
c                  know number & type of data values
c
c                  set type below. may wish to alter based on
c                  element number.
c
c                  cohesive types:
c                    1 = linear-elastic   2,3,5 not used
c                    4 = exponential model 6 = ppr
c                    7 = creep
c
c
      cohesive_type = 6
c
      if( cohesive_type == 1 ) nvalues = 4
      if( cohesive_type == 4 ) nvalues = 6
      if( cohesive_type == 6 ) nvalues = 6
      if( cohesive_type == 7 ) nvalues = 0

      write(results_file_no,1100) cohesive_type, nvalues
      write(results_file_no,1200)
c
      do lines_read = 1, num_lines
       read(packet_file_no,iostat=ios) elem, point, values(1:nvalues)
       if( ios .gt. 0 ) then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9000)
         stop
       else if( ios .lt. 0 )then
         call closer(packet_file_no)
         call closer(results_file_no)
         write(screen,9100)
         stop
       end if
       write(results_file_no,1000) step, elem, point,
     &                             values(1:nvalues)
      end do
c
      return
c
 1000 format(1x,i5,i8,i6,8f15.6)
 1100 format(/,1x,'Interface-cohesive displacement jumps ',
     &   'packet (type=32)....',
     & /,1x,'cohesive model type: ',i2,2x,'nvalues: ',i2 )
 1200 format('  step element  point',10x,20('-'),' values ',20('-'))
 9000 format(/,1x,'>> ERROR: reading packet file...terminating in disp')
 9100 format(/,1x,'>> END OF FILE...program terminating in disp')
c
      end
