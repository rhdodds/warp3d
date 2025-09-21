c     ****************************************************************          
c     *                                                              *          
c     *               subroutine oustr_pat_flat_file                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 8/6/2025 rhd               *          
c     *                                                              *          
c     *     output of stresses or strains to patran results file     *          
c     *     or to a flat file (text or stream)                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oustr_pat_flat_file( do_stress, ouflg, oubin, ouasc,              
     &                           ounod, flat_file, stream_file,                 
     &                           text_file, compressed  )                       
      use global_data, only : ltmstp, stname, lsldnm
      use elblk_data, only : elestr ! mxvl,mxoupr,mxoupt 
      use constants                       
      implicit none                                                             
c                                                                               
      logical :: ouflg, oubin, ouasc, do_stress, ounod,                            
     &           flat_file, stream_file, text_file, compressed                  

c                                                                               
c                       MPI:                                                    
c                         Results remain distributed onto ranks for.            
c                         Patran and flat files. Results are                    
c                         written in pieces,                                    
c                         with one piece (file) for each rank.                  
c                                                                               
      call wmpi_alert_slaves ( 35 )                                             
      call wmpi_bcast_log( do_stress )                                             
      call wmpi_bcast_log( ouflg )                                              
      call wmpi_bcast_log( oubin )                                              
      call wmpi_bcast_log( ouasc )                                              
      call wmpi_bcast_log( ounod )                                              
      call wmpi_bcast_int ( ltmstp ) ! load step number                         
      call wmpi_bcast_string( stname, 8 )                                       
      call wmpi_bcast_string( lsldnm, 8 )                                       
c                                                                               
      call wmpi_bcast_log( flat_file )                                          
      call wmpi_bcast_log( stream_file )                                        
      call wmpi_bcast_log( text_file )                                          
      call wmpi_bcast_log( compressed )                                         
c                                                                               
c                       output stress or strain result files to                 
c                       patran or flat files.                                   
c                                                                               
c                       These can be averaged values at model nodes             
c                       averaged values at element centers.                     
c                                                                               
c                       Patran files can be ascii (text) or sequential          
c                       binary (Fortran).                                       
c                                                                               
c                       Flat files can be text (and optionally                  
c                       compressed) or stream binary files.                     
c                       Compression valid only on Linux and MacOS.               
c                                                                               
c                       set up block arrays used in the module by the 
c                       lower level routines. release them when done                  
c                       with output      
c                       
c                       elestr is mxvl x max number output vales in code 
c                       need better design than to share this across
c                       lower level routines.                                      
c                                                                               
      call oustr_set_block_arrays  
c          
      elestr = zero  ! used as work in lower levels. prevent un-inits                                       
      if ( ounod ) then  ! nodal strain or stress values                        
          call oupstr_node( do_stress, oubin, ouasc, flat_file,                    
     &                      stream_file, text_file, compressed )                
      else  ! element strain or stress values                                   
          call oupstr_elem( do_stress, oubin, ouasc, flat_file,                    
     &                      stream_file, text_file, compressed )                
      end if      
c                                                                    
      call oustr_release_block_arrays                                           
c                                                                               
      ouflg = .false.                                                           
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oustr                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 12/3/2019 rhd              *          
c     *                                                              *          
c     *     drive printed output or packet file output of stresses   *          
c     *     or strains according to the lists and options specified  *          
c     *     by the user.                                             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oustr( do_stress, ouflg, oupat, oubin, ouasc, ounod,              
     &                  wide, eform, prec, noheader, out_packet_now )           
c
      use global_data, only : noelem, out 
      use elblk_data, only : elestr ! mxvl,mxoupr,mxoupt                        
      use allocated_integer_list
      use constants
      implicit none                                                             
c
      logical :: ouflg, oupat, oubin, ouasc, ounod, wide, eform, prec,          
     &           do_stress, noheader, out_packet_now                               
c                                                                               
c                       local declarations                                      
c                                                                               
      integer :: lenlst, errnum, icn, iplist, num_list_entries,                 
     &           bad_list, param, next, i, dummy, list_size                                
      integer, allocatable :: intlst(:), element_list(:)                                  
      logical, external :: matchs, true, endcrd
      logical :: eject_flag, do_strains                                                     
      real :: dumr                                                              
      double precision :: dumd 
      character(len=1) :: dums                                                  
c                                                                               
c                       MPI:                                                    
c                         we need to gather all the stresses or strains         
c                         back to the root processor for printing.              
c    
      do_strains = .not. do_stress                                                                           
      if( do_stress )  call wmpi_get_str ( 2 )                                                
      if( do_strains ) call wmpi_get_str ( 4 )                                                
c                                                                               
c                      read the integerlist of element numbers from             
c                      input line.                                              
c                                                                               
c                      run thru list to count the actual number of              
c                      elements to be output. check validity of element         
c                      numbers.                                                 
c                                                                               
c                      allocate a vector to hold the expanded                   
c                      list.                                                    
c                                                                               
c                      then load all element number from the                    
c                      interlist. Maintains the order of elements               
c                      implied by the integerlist.                              
c                                                                               
      eject_flag = .false.                                                      
      ouflg = .false.  ! = .true. if parse rules fail                           
c                                                                               
      if( matchs('for',3) )      call splunj                                    
      if( matchs('elements',4) ) call scan                                      
c   
      allocate( intlst(10) )                                                                            
      call trlist_allocated( intlst, list_size, noelem, 
     &                       lenlst, errnum )                      
c                                                                               
      select case( errnum )                                                     
      case( 1 )  ! list found ok                                                
         continue                                                               
      case( 2 )  ! parse rules failed                                           
         param = 1                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         ouflg = .true.                                                         
      case( 3 )  ! list would overflow space                                    
         param = 2                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         ouflg = .true.                                                         
      case( 4 )  ! no list or all. make for all elements                        
         call errmsg( 244, param, dums, dumr, dumd )                            
         lenlst = 3; intlst(1) = 1; intlst(2) = -noelem                         
         intlst(3) = 1                                                          
      case default                                                              
         param = 3                                                              
         call errmsg( 24, param, dums, dumr, dumd )                             
         ouflg = .true.                                                         
      end select                                                                
c                                                                               
      if( ouflg ) return                                                        
c                                                                               
      icn = 0; iplist = 1; num_list_entries = 0; bad_list = 0                   
      do                                                                        
       if( iplist .eq. 0 ) exit                                                 
       call trxlst( intlst, lenlst, iplist, icn, next)                          
       num_list_entries = num_list_entries + 1                                  
       if( next .eq. 0 .or. next .gt. noelem ) then                             
         if( bad_list .eq. 0 ) write(out,9010)                                  
         write(out,9020) next                                                   
         bad_list = bad_list + 1                                                
       end if                                                                   
      end do                                                                    
c                                                                               
      if( bad_list .gt. 0 ) then                                                
        write(out,9030);  call backsp( 1 )                                      
        ouflg = .true.; return                                                  
      end if                                                                    
c                                                                               
      allocate( element_list(num_list_entries) )                                
      icn = 0; iplist = 1; i = 1;                                               
      do                                                                        
       if( iplist .eq. 0 ) exit                                                 
       call trxlst( intlst, lenlst, iplist, icn, next)                          
       element_list(i) = next; i = i + 1                                        
      end do                                                                    
c                                                                               
c                       set up the "block"                                      
c                       arrays used in the module by the lower                  
c                       level routines. release them when done                  
c                       with output                                             
c                                                                               
      call oustr_set_block_arrays
      elestr = zero ! used as work array in lower routines                                               
      call ouhstr( do_stress, wide, eform, prec,                                   
     &             noheader, out_packet_now, element_list,                      
     &             num_list_entries )                                           
      call oustr_release_block_arrays                                           
c                                                                               
c                       done with output. set up scanner                        
c                       to continue processing same output                      
c                       command. could be an eol after we                       
c                       scanned the intergerlist above.                         
c                                                                               
      eject_flag = .true.                                                       
      call noscan   ! don't move on next test function                          
      if( .not. endcrd(dummy) ) then                                            
        call backsp( 1 )                                                        
        if( true(dummy) ) call splunj                                           
      end if                                                                    
c                                                                               
      if ( eject_flag ) write(out,fmt='(a1)') char(12)                          
c                                                                               
      deallocate( element_list, intlst )                                                
c                                                                               
      return                                                                    
c                                                                               
 9010 format(/1x,'>>>>> warning: these elements are invalid...'/)               
 9020 format(15x,i10)                                                           
 9030 format(/1x,'>>>>> warning: output request ignored...'/)                   
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine oustr_set_block_arrays               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/30/2019 rhd              *          
c     *                                                              *          
c     *     allocate block type arrays used by lower-level output    *          
c     *     routines                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine oustr_set_block_arrays  
c                                             
      use global_data, only : nelblk, out, mxvl
      use elem_block_data, only : history_blk_list, gausspts_blk_list           
      use elblk_data, only : elem_hist, elem_hist1, blk_size_gp,                
     &                       blk_size_hist                                     
c                                                                               
      implicit none                                                             
c                                                                               
      integer :: blk, ierror  
      logical, parameter :: local_debug = .false.                                                      
c                                                                               
c                  allocate a 3-D array block for the element histories         
c                  at n and n+1. find the maximum number of gauss points        
c                  for all blocks and the maximum history size per              
c                  gauss point for all blocks.                                  
c                                                                               
      blk_size_gp = 0                                                           
      blk_size_hist = 0                                                         
      do blk = 1, nelblk                                                        
         blk_size_gp   = max( blk_size_gp, gausspts_blk_list(blk) )             
         blk_size_hist = max( blk_size_hist, history_blk_list(blk) )            
      end do 
c
      if( local_debug ) then
        write(out,*) "... leaving  oustr_set_block_arrays:"
        write(out,*) "     nelblk, blk_size_gp, blk_size_hist",
     &           nelblk, blk_size_gp, blk_size_hist
      end if
c                                                                               
      allocate( elem_hist(mxvl,blk_size_hist,blk_size_gp),
     &          stat=ierror )
      if( ierror .ne. 0 ) then
        write(out,*) '.... oustr_set_block_arrays'
        write(out,*) '     allocate failure @ 1'
        call die_abort
      end if                                  
c                      
      allocate( elem_hist1(mxvl,blk_size_hist,blk_size_gp),
     &          stat=ierror )                    
      if( ierror .ne. 0 ) then
        write(out,*) '.... oustr_set_block_arrays'
        write(out,*) '     allocate failure @ 2'
        call die_abort
      end if                                  
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine oustr_release_block_arrays           *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/30/2019 rhd              *          
c     *                                                              *          
c     *     deallocate block type arrays used by lower-level output  *          
c     *     routines                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine oustr_release_block_arrays    
c                                       
      use global_data, only : out 
      use elblk_data, only : elem_hist, elem_hist1                              
c                                                                               
      implicit none  
      integer  :: ierror
c                                                                               
c                  deallocate a 3-D array block for the element histories       
c                  at n and n+1.                                                
c                                                                               
      deallocate( elem_hist, elem_hist1, stat=ierror )     
      if( ierror .ne. 0 ) then
        write(out,*) '.... oustr_release_block_arrays'
        write(out,*) '     deallocate failure'
        call die_abort
      end if                                  
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
