c     ****************************************************************          
c     *                                                              *          
c     *                      module local oupstr_node                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/1/2019 rhd               *          
c     *                                                              *          
c     ****************************************************************          
      module local_oupstr_node
c
      implicit none
c
      type :: node_entry      
         integer :: count                                                       
         double precision, dimension(:), allocatable :: node_values                 
      end type node_entry                                                       
c                                                                               
      type (node_entry), dimension (:), allocatable :: nodal_values
c
      end module local_oupstr_node
c
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupstr_node                  *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *     drives output of stress or strain nodal results to       *          
c     *     (1) patran files in either binary or formatted forms or  *          
c     *     (2) flat text or stream files                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupstr_node( do_stress, oubin, ouasc, flat_file,                  
     &                        stream_file, text_file, compressed )               
c
      use global_data, only : nonode, nelblk, myid, iprops, lprops,
     &                        outmap, use_mpi, numprocs, mxelmp,
     &                        out, elblks
      use local_oupstr_node, only : nodal_values
      use main_data, only : incmap, incid, cohesive_ele_types                   
      use output_value_indexes, only : num_short_strain, 
     &            num_short_stress, num_long_strain, num_long_stress
c
      implicit none
c      
      logical :: do_stress, oubin, ouasc, flat_file,                                  
     &           stream_file, text_file, compressed                                
c                                                                               
c                local declarations                                             
c                                                                               
      logical :: bbar_flg, geo_non_flg, long_out_flg, nodpts_flg,                 
     &           center_output, cohesive_elem, do_average, 
     &           threads_only, do_strains, is_cohesive                
      integer  :: iblk, node, blk, span, felem, elem_type,  
     &            int_order, mat_type, num_enodes, num_enode_dof,
     &            totdof, num_vals, ifelem, count,num_int_points,
     &            output_loc, ierror   
      integer :: elem_out_map(mxelmp) ! mxelmp in global_data
      double precision, parameter :: zero = 0.0d0
c                                                                               
c                data structure for averaged nodal results. vector              
c                of derived types. global vector allocated for all              
c                structure nodes. derived type allocated only for               
c                nodes touched by elements in this domain.                      
c                                                                               
c                patran nodal result file contains values only                  
c                for nodes touched by elements in this domain.                  
c                they are averaged by the number of elements that               
c                touch then in this domain.                                     
c                                                                               
c                       create structure size vector of derived types.          
c                       set count for each node to zero.                       
c               
      do_strains = .not. do_stress                                                                
      allocate( nodal_values(nonode) )                                          
      nodal_values(1:nonode)%count = 0                                             
c                                                                               
c                       compute the stress/strain nodal output                  
c                       vectors and assemble them into the global               
c                       nodal results for each element block                    
c                                                                               
      do blk = 1, nelblk                                                        
c                                                                               
         if( elblks(2,blk) .ne. myid ) cycle                                   
         span              = elblks(0,blk)                                      
         felem             = elblks(1,blk)                                      
         elem_type         = iprops(1,felem)                                    
         int_order         = iprops(5,felem)                                    
         mat_type          = iprops(25,felem)                                   
         num_enodes        = iprops(2,felem)                                    
         num_enode_dof     = iprops(4,felem)                                    
         totdof            = num_enodes * num_enode_dof                         
         geo_non_flg       = lprops(18,felem)                                   
         bbar_flg          = lprops(19,felem)                                   
         num_int_points    = iprops(6,felem)                                    
         long_out_flg      = lprops(16,felem)                                   
         output_loc        = iprops(12,felem)                                   
         nodpts_flg        = .true.                                             
         center_output     = output_loc .eq. 3                                  
         cohesive_elem     = cohesive_ele_types(elem_type)                      
c                                                                               
         if( do_strains) num_vals = num_long_strain
         if( do_stress ) num_vals = num_long_stress
c                                                                               
c                       skip cohesive elements for now. they do not             
c                       contribute to nodal results at structure level.         
c                                                                               
         if( cohesive_elem ) cycle                                             
c                                                                               
         elem_out_map(1:mxelmp) = outmap(elem_type,1:mxelmp)                             
c                                                                               
c                       duplicate necessary element block data.                 
c              
         is_cohesive = .false.                                                                 
         call oudups( span, felem, num_int_points, geo_non_flg, 
     &                do_stress, is_cohesive )                                               
c                                                                               
c                       compute the element block nodal stress/strain           
c                       data.                                                   
c                                                                               
         iblk   = blk                                                           
         ifelem = felem                                                         
         call ouprks( span, iblk, ifelem, elem_type, int_order,                 
     &                num_int_points, num_enodes, geo_non_flg,        
     &                do_stress, mat_type, center_output, .false. )                               
c                                                                               
c                       add the element block nodal stress/strain data          
c                       into the global nodal results data structure.           
c  
         call oupstr_node_1                  
c                                                                               
      end do 
c
c                       for threads only processing and/or mpi with just              
c                       one process, average the total results at               
c                       each node, using node count for the model.              
c                       then calculate the extended values of stress            
c                       and strain at each node in the model.                   
c                                                                               
c                       we can have nodes with only interface-cohesive          
c                       elements attached. make sure they have a vector         
c                       of zero values to write in the file.                    
c                                                                               
c                       for multi-process mpi jobs, we just write               
c                       the summed results + the node count to the              
c                       simplified result file (for the domain).                
c                                                                               
c                       the external program to combine mpi result              
c                       files handles nodes with missing values.                
c                                                                               
c                                                                               
      threads_only = .not. use_mpi                                                    
      do_average =  threads_only .or. (use_mpi .and. numprocs .eq. 1)    
      if( do_average ) call oupstr_node_do_average
c                                                                               
c                       output the averaged total nodal results                 
c                       to a (1) file compatable with patran for post           
c                       processing in binary or formatted file type             
c                       or (2) flat file in text or stream types.               
c                                                                               
c                       only results for nodes that appear in this              
c                       domain are written. see notes                           
c                       above for values actually written.                      
c                                                                                
      call oustpa( do_stress, oubin, ouasc, num_vals, 
     &             nonode, flat_file, stream_file, text_file,                   
     &             compressed )                                                 
c                                                                               
c                       deallocate all instances of the derived type            
c
      deallocate( nodal_values )
c                                                                               
      return  
c
      contains
c     ========                                                                  
c

      subroutine oupstr_node_1  ! for a block
c     ------------------------                      
c
      use elblk_data, only : elestr ! mxvl x max allowed output values                                             
      implicit none                                                    
c                    
      integer :: i, j, k, snode, ierror, map                                                                        
c                                                                               
c                       look at nodes present on each element of                
c                       this block. if the vector of node average               
c                       values does not exist for a node, create and            
c                       zero it. update the counter for each structural node    
c                       accessed in the processing of the element               
c                       block.                                                  
c    
      do j = 1, num_enodes                                                          
       do i = 1, span                                                           
         snode = incid(incmap(felem+i-1)+j-1) 
         if( snode <= 0  .or. snode > nonode ) then
            write(out,*) '... bad snode @ 1. snode: ',snode
            write(out,*) '    num_enodes, i, j:',num_enodes,i,j
            stop' @ 5'     
         end if                                              
         if( .not. allocated( nodal_values(snode)%node_values) ) then          
             allocate( nodal_values(snode)%node_values(num_vals) )                
             nodal_values(snode)%node_values(1:num_vals) = zero                   
         end if                                                                 
         nodal_values(snode)%count = nodal_values(snode)%count + 1              
       end do                                                                   
      end do      
c                                                                               
c                       add the results from each node of each                  
c                       element in this block into the global nodal results     
c                       data structure.                                         
c                                                                               
      do k = 1, num_vals                                                        
         map = elem_out_map(k)                                                        
         do j = 1, num_enodes                                                        
!$omp simd
          do i = 1, span                                                        
            snode = incid(incmap(felem+i-1)+j-1)   
            nodal_values(snode)%node_values(map) =  
     &        nodal_values(snode)%node_values(map) + elestr(i,k,j)   
          end do                                                                
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end subroutine oupstr_node_1                                                                  


      subroutine oupstr_node_do_average ! over all model nodes
c     ---------------------------------
c
      implicit none
c
      integer :: snode, count, j, map
      double precision :: dcount
c
c$OMP PARALLEL DO  PRIVATE( snode )
      do snode = 1, nonode  ! all model nodes                                                    
        if( nodal_values(snode)%count .ne. 0 ) cycle
        allocate( nodal_values(snode)%node_values(num_vals) )                 
        nodal_values(snode)%node_values(1:num_vals) = zero                    
        nodal_values(snode)%count = 1                                       
      end do 
c$OMP END PARALLEL DO           
c                                                                       
c$OMP PARALLEL DO  PRIVATE( snode, dcount, j, map )
      do snode = 1, nonode  ! all model nodes                                                    
        dcount = dble( nodal_values(snode)%count )      
        do j = 1, num_vals                                                     
          map = elem_out_map(j)                                                
          nodal_values(snode)%node_values(map) =
     &              nodal_values(snode)%node_values(map) / dcount
        end do      
        call ouext2( nodal_values(snode)%node_values(1), 1, 1, 
     &               do_stress )   ! get expanded results
      end do ! on snode       
c$OMP END PARALLEL DO    
c
      return                                                 

      end subroutine oupstr_node_do_average
c
      end subroutine oupstr_node

c     ****************************************************************
c     *                                                              *
c     *                      subroutine oustpa                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 8/5/25 rhd                 *
c     *                                                              *
c     *     output stress or strain nodal results to (1) patran file *
c     *     in either binary or formatted forms, or (2) flat file    *
c     *     in stream or text formats                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oustpa( stress, oubin, ouasc, num_vals,
     &                   num_struct_nodes, flat_file,
     &                   stream_file, text_file, compressed )
c
      use global_data, only : ltmstp, out, use_mpi, myid, numprocs,
     &                        stname, lsldnm, nonode 
      use local_oupstr_node, only : nodal_values
      implicit none
c
      integer :: num_vals, num_struct_nodes
      logical oubin, ouasc, stress, flat_file, stream_file, text_file,
     &       compressed
c
c                       local declarations
c
      integer :: dum1, dum3, bnfile, flat_file_number, fmfile,
     &           quantity, stepno, titl(80), titl1(80), i, snode
      logical :: patran_file
      real :: dum2
      character(len=4) :: title(80), title1(80)
      character(len=80) :: string, strng1, stepstring*7
      equivalence (title,titl), (title1,titl1)
      double precision, parameter :: small_tol=1.0d-80, zero=0.0d0
c
      patran_file = oubin .or. ouasc
      stepno = ltmstp ! more descriptive name for use here
c
c                       open Patran or flat files
c
      call ouocst( stress, stepno, oubin, ouasc, bnfile, fmfile, 1,
     &             out, use_mpi, myid, flat_file, stream_file,
     &             text_file, compressed, flat_file_number )
c
c                        set Patran title records.
c
      if( patran_file ) call oustpa_1
c
c                       set dummy values not needed since the
c                       results are not displacement results.
c
      dum1 = 0
      dum2 = 0.0
      dum3 = 0
c
c                       write to a Patran binary file
c
      if( patran_file .and. oubin ) call oustpa_2
c
c                       write to a Patran formatted file
c
      if( patran_file .and. ouasc ) call oustpa_3
c
      if( patran_file ) then  ! close and leave
         call ouocst( stress, stepno, oubin, ouasc, bnfile, fmfile,
     &                2, out, use_mpi, myid, flat_file, stream_file,
     &                text_file, compressed, flat_file_number )
         return
      end if
c
c                       write to a flat text or stream file
c
      if( .not. flat_file ) then ! big error somewhere
         write(out,9200)
         call die_abort
      end if
c
      if( text_file ) then
        quantity = 1
        if( stress ) quantity = 2
        call ouddpa_flat_header( 2, quantity, flat_file_number,
     &                           num_vals, '(30e15.6)' )
      end if
c
      if( use_mpi .and. numprocs > 1 ) then
         call oustpa_4
       else    !  just threaded or mpi w/ num_procs = 1
         call oustpa_5
      end if
c
c                       close flat file
c
      call ouocst( stress, stepno, oubin, ouasc, bnfile, fmfile, 2,
     &             out, use_mpi, myid, flat_file, stream_file, 
     &             text_file, compressed, flat_file_number )
c
      return
c
 9200 format(/1x,'>>>>> Fatal Error: routine oustpa.',
     &   /16x,   'should have not reach this point @ 1.',
     &   /16x,   'job terminated' )
c
      contains
c     --------
c
      subroutine oustpa_1
c     -------------------
c
      implicit none
c
c                        set Patran title records.
c
      string = ' '
      strng1 = ' '
      string = 'nodal strain results for structure '
      if( stress ) string = 'nodal stress results for structure '
      strng1 = 'loading '
c
      string(36:43) = stname
      strng1(9:16)  = lsldnm
      write(stepstring,'(i7)') stepno
      string(44:) = ', step '// stepstring
c
      do i = 1, 80
         title(i)  = string(i:i)
         title1(i) = strng1(i:i)
      end do
c
      return
c
      end subroutine oustpa_1
c
      subroutine oustpa_2
c     -------------------
c
      implicit none
c                       write to a Patran binary file
c
c
c                       records 1-3
c
      write(bnfile) titl, nonode, dum1, dum2, dum3, num_vals
      write(bnfile) titl1
      write(bnfile) titl1
c
c                    nodal results records. only those with non-zero
c                    counts are in this domain. we zero very small
c                    values to prevent 3 digit exponents in
c                    formatted output.
c
      do snode = 1, nonode
        if( nodal_values(snode)%count .eq. 0 ) cycle
        associate( x = > nodal_values(snode)%node_values )
        if( use_mpi .and. numprocs > 1 ) then
            write(bnfile) snode, nodal_values(snode)%count,
     &                    x(1:num_vals)
        else
            write(bnfile) snode, (sngl(x(i)),i=1,num_vals)
        end if
      end associate
      end do
c
      return
c
      end subroutine oustpa_2

      subroutine oustpa_3
c     -------------------
c
      implicit none
c
c                       records 1-4
c
      write(fmfile,900) titl
      write(fmfile,910) nonode, dum1, dum2, dum3, num_vals
      write(fmfile,900) titl1
      write(fmfile,900) titl1
c
c                    nodal results records. only those with non-zero
c                    counts are in this domain
c
      do snode = 1, nonode
       if( nodal_values(snode)%count .eq. 0 ) cycle
       associate( x = > nodal_values(snode)%node_values )
       where( abs(x) .lt. small_tol ) x = zero
       if( use_mpi .and. numprocs > 1 ) then
         write(fmfile,930) snode, nodal_values(snode)%count
         write(fmfile,940) x(1:num_vals)
       else
         write(fmfile,920) snode, x(1:num_vals)
       end if
       end associate
      end do
c
      return
 900  format(80a1)
 910  format(2i9,e15.6,2i9)
 920  format(i8,5(e13.6))
 930  format(2i8)
 940  format(e23.15)
c
      end subroutine oustpa_3

      subroutine oustpa_4
c     -------------------
c
      implicit none

       if( stream_file ) write(flat_file_number) num_vals
       if( text_file ) write(flat_file_number,*) "  ", num_vals
c
       do snode = 1, nonode
         if( nodal_values(snode)%count .eq. 0 ) cycle
         associate( x = > nodal_values(snode)%node_values )
         where( abs(x) .lt. small_tol ) x = zero
         if( stream_file ) write(flat_file_number) snode,
     &             nodal_values(snode)%count, x(1:num_vals)
         if( text_file ) write(flat_file_number,9100) snode,
     &             nodal_values(snode)%count, x(1:num_vals)
         end associate
       end do
c
      return
c
 9100 format(2i9,30d15.6)
c
      end subroutine oustpa_4

      subroutine oustpa_5
c     -------------------
c
      implicit none
c
      do snode = 1, nonode
         if( nodal_values(snode)%count .eq. 0 ) cycle
         associate( x => nodal_values(snode)%node_values )
         where( abs(x) .lt. small_tol ) x = zero
         if( stream_file) write(flat_file_number) x(1:num_vals)
         if( text_file) write(flat_file_number,9300) x(1:num_vals)
         end associate
      end do
c
      return
c
 9300 format(30e15.6)
c
      end subroutine oustpa_5
c
      end subroutine oustpa

  