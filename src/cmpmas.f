c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine cmpmas                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 8/11/2017 rhdd             *          
c     *                                                              *          
c     *     drives computation of the global, diagonal mass matrix   *          
c     *     (stored as a vector). the effective diagonal mass        *          
c     *     of each element is computed/saved (compact form)         *          
c     *     element diagonal mass formed from consistent element     *          
c     *     mass matrix.                                             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine cmpmas                                                         
      use global_data ! old common.main
c                                                                               
      use elem_block_data, only : mass_blocks, cdest_blocks,                    
     &                            edest_blocks                                  
      use main_data,       only : mdiag,                                        
     &                            fgm_node_values_defined,                      
     &                            fgm_node_values,                              
     &                            incid, incmap,                                
     &                            cohesive_ele_types,                           
     &                            linear_displ_ele_types,                       
     &                            implemented_ele_types,                        
     &                            axisymm_ele_types, bar_types,
     &                            link_types                             
c                                                                               
      implicit none                                                   
c                                                                               
c              locals
c
      integer :: iout, blk, felem, elem_type, def_int_order, 
     &           num_enodes, num_enode_dof, totdof, def_num_int_points,
     &           span, iok, i, int_order, num_int_points                                                         
      double precision,                                                         
     &     allocatable, dimension(:,:) :: mel                                   
c                                                                               
      double precision :: totvol, ce_block(mxvl,mxecor), 
     &                    rho_block(mxndel,mxvl)
      double precision, parameter :: zero=0.0d0, three=3.0d0          
      logical :: fgm_props, is_bar_elem, is_link_elem                                                     
c                                                                               
      call thyme( 4, 1 )                                                        
      iout = out
c                                                                               
c              MPI: tell the slave processors to join us here                   
c              non-MPI: a dummy routine.                                        
c                                                                               
c              the former local_mdiag                                           
c              data structure was deleted with dropping the EBE                 
c              solver                                                           
c                                                                               
      call wmpi_alert_slaves ( 8 )                                              
c                                                                               
c              intitialize the structure diagonal mass vector                   
c                                                                               
      totvol = zero                                                             
      mdiag(1:nodof) = zero                                                     
c                                                                               
c              allocate blocks of element mass matrices and                     
c              the pointer vector to the blocks. we only store                  
c              num_enodes terms for each element since the remaining            
c              2 * num_enodes terms are the same.                               
c                                                                               
      call rdmass                                                               
c                                                                               
c              compute and assemble the masses for each element                 
c              block. use a full integration order to get mass.                 
c                                                                               
c              MPI: elblks(2,blk) holds which processor owns the                
c                   block.                                                      
c              non-MPI: all blocks processed                                    
c                                                                               
      do blk = 1, nelblk                                                        
c                                                                               
         if( myid .ne. elblks(2,blk) ) cycle                                    
c                                                                               
c              'def_int_order' and 'def_num_int_points'                         
c              are the default values taken from the props table.               
c              routine mass_getint returns final values in                      
c              variables int_order and num_int_points.                          
c                                                                               
         felem          = elblks(1,blk)                                         
         elem_type      = iprops(1,felem)                                       
         def_int_order  = iprops(5,felem)                                       
         num_enodes     = iprops(2,felem)                                       
         num_enode_dof  = iprops(4,felem)                                       
         totdof         = num_enodes * num_enode_dof                            
         def_num_int_points = iprops(6,felem)                                   
         span           = elblks(0,blk) 
         is_bar_elem    = bar_types(elem_type)
         is_link_elem   = link_types(elem_type)                                        
c                                                                               
         call mass_getint( elem_type, def_int_order, def_num_int_points,        
     &                     int_order, num_int_points, iout)                     
c                                                                               
         call dupmas( span, num_enodes, cdest_blocks(blk)%ptr(1,1),             
     &                c, totdof, ce_block, mxvl )                               
c                                                                               
c              mass density is defined at model nodes for the fgm               
c              capability. gather values for element nodes in block.            
c                                                                               
         if( fgm_node_values_defined )                                          
     &      call dupmas_fgm( span, num_enodes, incid(incmap(felem)),            
     &                       fgm_node_values, nonode, rho_block,                
     &                       mxndel  )                                          
c                                                                               
         fgm_props = fgm_node_values_defined                                    
c                                                                               
         allocate( mel(totdof,span),stat=iok ) ! full size block                
         if( iok .ne. 0 ) then                                                  
           write(iout,9100) iok, 1                                              
           call die_abort                                                       
         end if                                                                 
c                              
         call rknmas( span, felem, iprops(1,felem), int_order,                  
     &                num_int_points, num_enodes, totdof,                       
     &                props(1,felem),                                           
     &                mel, beta_fact,                                           
     &                totvol, fgm_props,                                        
     &                cohesive_ele_types(elem_type),                            
     &                linear_displ_ele_types(elem_type),                        
     &                axisymm_ele_types(elem_type),                             
     &                implemented_ele_types(elem_type),                         
     &                ce_block, rho_block, out )                                     
c                                                                               
         call addmas( span, edest_blocks(blk)%ptr(1,1), totdof,                 
     &                mdiag, mel )                                              
c                                                                               
         mass_blocks(blk)%ptr(1:num_enodes,1:span) =                            
     &                    mel(1:num_enodes,1:span)                              
         deallocate( mel, stat=iok )                                            
         if( iok .ne. 0 ) then                                                  
           write(iout,9100) iok, 2                                              
           call die_abort                                                       
         end if                                                                 
c                                                                               
      end do ! over blk                                                         
c                                                                               
c              MPI: combine all the processor local contributions               
c                   of mdiag and totvol into a full version on the              
c                   root processor.                                             
c              non-MPI: everything here                                         
c                                                                               
      call wmpi_reduce_vec ( mdiag, nodof )                                     
      call wmpi_reduce_vec ( totvol, 1 )                                        
c                                                                               
c              slaves leave now                                                 
c                                                                               
      if ( slave_processor ) return                                             
c                                                                               
c              set a global variable defining total mass of model               
c                                                                               
      total_mass = sum( mdiag(1:nodof) )
c
      write(iout,9000) total_mass / three                                       
      write(iout,9010) totvol                                                   
c                                                                               
      call thyme(4,2)                                                           
c                                                                               
      return                                                                    
c                                                                               
 9000 format(10x,'>>> total model mass (computed): ',e14.6 )                    
 9010 format(10x,'>>> total model volume:          ',e14.6 )                    
 9100 format('>> FATAL ERROR: cmpmas, memory allocate/deallocate',              
     & ' failure',                                                              
     &  /,   '                status= ',i5,' @ ',i5,                            
     &  /,   '                job terminated' )                                 
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine mass_getint                  *          
c     *                                                              *          
c     *                       written by : rau                       *          
c     *                                                              *          
c     *                    last modified : 08/11/2017 rhd            *          
c     *                                                              *          
c     *                                                              *          
c     *     this subroutine returns the order of integration         *          
c     *     and number of integration points given the element       *          
c     *     type. these values may be different from that            *          
c     *     given in the props table due to the higher order         *          
c     *     needed to integrate and determine the lump mass.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine mass_getint( elem_type, def_int_ord, def_num_pts,               
     &                       int_ord, num_int_pts, iout )                             
      implicit integer (a-z)                                                    
c                                                                               
      select case( elem_type )                                                  
c                                                                               
c        element type is 'q3disop'                                              
c                                                                               
         case(1)                                                                
            int_ord = 9                                                         
            num_int_pts = 14                                                    
c                                                                               
c        element type is 'l3disop'                                              
c                                                                               
         case(2)                                                                
            int_ord     = 1                                                     
            num_int_pts = 8                                                     
c                                                                               
c        element type is 'ts12isop'                                             
c                                                                               
         case(3)                                                                
            int_ord     = 9                                                     
            num_int_pts = 14                                                    
c                                                                               
c        element type is 'ts15isop'                                             
c                                                                               
         case(4)                                                                
            int_ord     = 9                                                     
            num_int_pts = 14                                                    
c                                                                               
c        element type is 'ts9isop'                                              
c                                                                               
         case(5)                                                                
            int_ord     = 9                                                     
            num_int_pts = 14                                                    
c                                                                               
c        element type is 'tet10'                                                
c                                                                               
         case(6)                                                                
            int_ord     = 15                                                    
            num_int_pts = 15                                                    
c                                                                               
c        element type is 'wedge15'                                              
c                                                                               
         case(7)                                                                
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                               
c        element type is 'tri6'                                                 
c                                                                               
         case(8)                                                                
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                               
c        element type is 'quad8'                                                
c                                                                               
         case(9)                                                                
            int_ord     = 1                                                     
            num_int_pts = 8                                                     
c                                                                               
c        element type is 'axiquad8'                                             
c                                                                               
         case(10)                                                               
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                               
c        element type is 'axitri6'                                              
c                                                                               
         case(11)                                                               
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                               
c        element type is 'inter8'                                               
c                                                                               
         case(12)                                                               
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                               
c        element type is 'tet4'                                                 
c                                                                               
         case(13)                                                               
            int_ord     = 4                                                     
            num_int_pts = 4                                                     
c                                                                               
c        element type is 'trint6'                                               
c                                                                               
         case(14)                                                               
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                               
c        element type is 'trint12'                                              
c                                                                               
         case(15)                                                               
            int_ord     = def_int_ord                                           
            num_int_pts = def_num_pts                                           
c                                                                                                                                                            
c        element type is 'bar2'                                              
c                                                                               
         case(18)                                                               
            int_ord     = 1                                           
            num_int_pts = 1                                           
c                                                                                                                                                            
c        element type is 'link2'                                              
c                                                                               
         case(19)                                                               
            int_ord     = 1                                           
            num_int_pts = 1                                           
c                                                                               
         case default                                                           
            write(iout, 1000)                                                   
            call die_abort                                                      
      end select                                                                
c                                                                               
      return                                                                    
c                                                                               
 1000 format(5x, ///'>>>> FATAL ERROR: cmpmas.f received an',                   
     &  /,       '                  invalid element type.')                     
c                                                                               
c                                                                               
      end subroutine mass_getint                                                
