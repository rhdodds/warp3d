c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine inelbk                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 11/26/2018 rhd             *          
c     *                                                              *          
c     *     input of element blocking information. computations for  *          
c     *     element blocking and optionally domains if requested     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine inelbk( sbflg1, sbflg2 )                                       
      use global_data ! old common.main
      implicit none                                                             
c                                                                               
      logical :: sbflg1,sbflg2                                                  
c                                                                               
      integer :: i, auto_size, ival, blk, param, proc, span, felem,             
     &           elem, msgcnt, dum, nc                                          
      integer, allocatable :: checkblk(:)                                       
      real :: dumr                                                              
      double precision :: dumd                                                  
      logical :: auto_blking, ok, display, auto_domains                         
      logical, external :: matchs, integr, endcrd                               
      logical, save :: first_proc = .true.                                      
      character :: dums*1, item*80                                              
c                                                                               
c                       if the subroutine has been previously                   
c                       entered and exited, then there was an                   
c                       error in the data of the last processed                 
c                       card. specifically, that error was the                  
c                       absence of the block number. under this                 
c                       circumstance, print an error message and                
c                       continue with input.                                    
c                                                                               
      if( sbflg1 ) call errmsg(43,dum,dums,dumr,dumd)                           
c                                                                               
c                       initialize blocking data structures.                    
c                                                                               
      nelblk = 0                  
      mxnmbl = noelem/128 + 1   ! starting size
      call mem_allocate( 17 )
c                                              
      do i = 1, mxnmbl                                                          
         elblks(0,i) = 0   ! number elements in block                           
         elblks(1,i) = 0   ! first element in block                             
         elblks(2,i) = 0   ! optional domain  number for element                
         elblks(3,i) = -1  ! doesn't seem to be used elsewhere                  
      end do                                                                    
c                                                                               
      if( numprocs == 1 ) then ! threads only, or mpi 1 domain                                
         do i = 1, mxnmbl                                                       
            elblks(2,i) = 0                                                     
         end do                                                                 
      else                                                       
         do i = 1, mxnmbl                                                       
            elblks(2,i) = -1                                                    
         end do                                                                 
      end if                                                                    
c                                                                               
      auto_size    = mxvl  !  from param_def                                    
      auto_blking  = .false.                                                    
      display      = .false.                                                    
      auto_domains = .false.                                                    
c                                                                               
      if( endcrd(dumr) ) then ! only word blocking on command                   
         call inelbk_user                                                       
      else                                                                      
         call inelbk_automatic ! probably automatic blocking/domains            
         call inelbk_user ! just in case use  still gives blocking              
      end if                                                                    
c                                                                               
c                      end of blocking input. run checks to make                
c                      sure that all elements have been assigned to a block     
c                      and that elements are assigned to only one block.        
c                      tries to find various errors that users can make !       
c                                                                               
      allocate( checkblk(noelem) )                                              
      checkblk(1:noelem) = 0                                                    
c                                                                               
      do blk = 1, nelblk                                                        
       felem = elblks(1,blk)                                                    
       span  = elblks(0,blk)                                                    
       do i = 1, span                                                           
         elem = felem + i - 1                                                   
         checkblk(elem) = checkblk(elem) + 1                                    
       end do                                                                   
      end do                                                                    
c                                                                               
      msgcnt = 0                                                                
      do elem = 1, noelem                                                       
        if( checkblk(elem) .ne. 1 ) then                                        
          if( msgcnt == 0 ) write(out,9000)                                     
          if( checkblk(elem) > 1 ) write(out,9010) elem                         
          if( checkblk(elem) == 0 ) write(out,9020) elem                        
          msgcnt = msgcnt + 1                                                   
        end if                                                                  
      end do                                                                    
c                                                                               
      deallocate( checkblk )                                                    
c                                                                               
      if( msgcnt > 0 ) then                                                     
         write(out,9030)                                                        
         call die_abort                                                         
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(//,'>>>> Fatal Error: inconsistent blocking for elements' )        
 9010 format(10x,'Element: ',i7,' assigned to more than 1 block' )              
 9020 format(10x,'Element: ',i7,' is not assigned to a block' )                 
 9030 format(//,'>>>> Job terminated due to fatal blocking errors' )            
c                                                                               
      contains                                                                  
c     ========                                                                  
c     ****************************************************************          
c     *                                                              *          
c     *      contains:   inelblk_user                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine inelbk_user                                                    
      implicit none                                                             
c                                                                               
      do while ( .true. )                                                       
      call readsc                                                               
c                                                                               
c                      read in element block. if no block is                    
c                      given, then either the block has been                    
c                      forgotten or blocking input has ended.                   
c                      in either case, set sbflg1 to true and                   
c                      exit the subroutine. if the block has                    
c                      been forgotten, then the rest of the                     
c                      line will be skipped.                                    
c                                                                               
      if( .not.integr(blk) ) exit                                               
c                                                                               
c                       valid block number ? resize                                    
c                                                                               
      if( blk < 0 ) then                                                        
         param = blk                                                            
         call errmsg(76,param,dums,dumr,dumd)                                   
         cycle                                                                  
      end if                                                                    
      if( blk > mxnmbl ) then  
         call inelbk_resize                                                 
c         param = blk                                                            
c         call errmsg(74,param,dums,dumr,dumd)                                   
c         cycle                                                                  
      end if                                                                    
c                                                                               
c                      read in the block span and the first                     
c                      element in the block.                                    
c                                                                               
      if( .not. integr(span) ) then                                             
         call errmsg(77,dum,dums,dumr,dumd)                                     
         cycle                                                                  
      end if                                                                    
c                                                                               
      if( (span < 1) .or. (span > mxvl) ) then                                  
         call errmsg(185,span,dums,dumr,dumd)                                   
         cycle                                                                  
      end if                                                                    
c                                                                               
      if( .not. integr(felem) ) then                                            
         call errmsg(159,dum,dums,dumr,dumd)                                    
         cycle                                                                  
      end if                                                                    
c                                                                               
c                      if present, read in processor assignment                 
c                      for given block.                                         
c                                                                               
      if( integr(proc) ) then                                                   
         if( proc .lt. 0 ) then                                                 
            call errmsg( 294, dum, dums, dumr, dumd )                           
            proc = 0                                                            
         end if                                                                 
         if( use_mpi ) then                                                     
            elblks(2,blk)= proc                                                 
         else                                                                   
            elblks(2,blk)= 0                                                    
         end if                                                                 
      end if                                                                    
c                                                                               
c                      set the appropriate block data structures.               
c                                                                               
      elblks(0,blk) = span                                                      
      elblks(1,blk) = felem                                                     
      nelblk        = max(nelblk,blk)                                           
c                                                                               
      end do                                                                    
c                                                                               
      sbflg1 = .true.                                                           
      sbflg2 = .true.                                                           
c                                                                               
      return                                                                    
      end subroutine inelbk_user                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *      contains:   inelblk_automatic                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine inelbk_automatic                                               
      implicit none                                                             
c                                                                               
      if( matchs( 'automatic',4) ) then                                         
         auto_blking = .true.                                                   
         if( use_mpi ) auto_domains = .true.                                    
      else                                                                      
          write(out,9040)  ! unknown blocking option                            
          return                                                                
      end if                                                                    
c                                                                               
      if( matchs('size',4) ) then                                               
        if( matchs('=',1) ) call splunj                                         
        if( integr( ival ) ) then                                               
          ok = ival >= 1  .and. ival <= mxvl                                    
          if( ok ) then                                                         
            auto_size = ival                                                    
          else                                                                  
            write(out,9050) ival, auto_size                                     
          end if                                                                
        else                                                                    
          call entits( item, nc )                                               
          write(out,9080) item(1:nc)                                            
          write(out,9070)                                                       
        end if                                                                  
      end if                                                                    
c                                                                               
c              scan optional words to maintain compatibility                    
c                                                                               
      if( matchs('domains',4) ) call splunj                                     
      if( matchs('automatic',4) ) call splunj                                   
c                                                                               
      if( matchs('display',4) ) then                                            
        display = .true.                                                        
      elseif( endcrd() ) then                                                   
        call splunj                                                             
      else ! unknown keyword                                                    
        call entits( item, nc )                                                 
        write(out,9060) item(1:nc)                                              
        write(out,9070)                                                         
      end if                                                                    
c                                                                               
      if( matchs('display',4) ) display = .true.                                
c                                                                               
c                      do auto blocking and optionally domains                  
c                                                                               
      call inelbk_simple_blocking( auto_size, auto_domains, display )           
c                                                                               
      return                                                                    
c                                                                               
 9040 format(//,1x                                                              
     &'>>>> Unknown option on blocking command - command ignore')               
 9050 format(//,1x,'>>>> Invalid blocksize: ',i4,' using: ',i4,                 
     &        /,1x,'     cannot exceed mxvl in file param_def')                 
 9060 format(//,1x,'>>>> Unrecognized blocking option: ',a20)                   
 9070 format(1x,   '     option ignored')                                       
 9080 format(//,1x,'>>>> Unrecognized size option: ',a20)                       
c                                                                               
      end subroutine inelbk_automatic                                           
      end subroutine inelbk                                                     
c     ****************************************************************          
c     *                                                              *          
c     *          subroutine inelbk_simple_blocking                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *             last modified : 2/27/2017 rhd                    *          
c     *                                                              *          
c     *     automatic generation of element blocking: simple version *          
c     *     for threaded w/o vectorized blocking. optional           *          
c     *     assignment of blocks to domains w/ simple algorithm      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine inelbk_simple_blocking( auto_size, auto_domains,               
     &                                   display )                              
      use global_data ! old common.main
      implicit none                                                             
      integer :: auto_size, auto_num_domains                                    
      logical :: auto_domains, display                                          
c                                                                               
c                     locals also visible by contains                           
c                                                                               
      logical :: geonon, bbar, cohesive, cyclic_plasticity,                     
     &           umat,                                                          
     &           blk_geonon, blk_bbar, blk_cohesive,                            
     &           blk_cyclic_plasticity, blk_umat,                               
     &           newblk, compatible, blk_dmg, dmg                               
c                                                                               
      integer ::  blk_matmodel, matmodel, blk_eletype, eletype,                 
     &            blk_intord, intord, blk_intnum, intnum, i,                    
     &            blk_matnum, matnum, current_size, felem,                      
     &            element, param, idomain, hblk, lblk, blks_per_domain,         
     &            blk, domain                                                   
      integer ::  counts_blks_ea_domain(0:max_procs-1,2)                        
      character(len=1) :: dums                                                  
      double precision :: dumd                                                  
      real :: dumr                                                              
c                                                                               
c                     first generation of automatic assignment of               
c                     elements to blocks.                                       
c                                                                               
c                     a) sequential pass thru elements                          
c                     b) assign to current block or open new                    
c                        block if full or element/material                      
c                        combinations are not compatible with                   
c                        current block.                                         
c                     c) no element renumbering or vectorized                   
c                        blocking in this first set of                          
c                        features.                                              
c                                                                               
c                                                                               
      nelblk       = 1   ! in global_data                                 
      current_size = 1                                                          
      felem        = 1                                                          
      elblks(1,1)  = 1  ! first element in block                                
      elblks(0,1)  = 1  ! number elements in block                              
      if( noelem == 1 ) return                                                  
c                                                                               
      call inelbk_load_block_props( 1 ) ! more than 1 elements                  
c                                                                               
c                     1. element fits in current block?                         
c                     2. if yes, is it compatible with elements now             
c                        in the block?                                          
c                     3. if yes, update current number of elements              
c                        in the block, next element                             
c                     4. otherwise start a new block. set first                 
c                        element in the block, init block size,                 
c                        load props for first element in block for              
c                        subsequent comparisons                                 
c                                                                               
      do element = 2, noelem                                                    
         newblk       = .false.                                                 
         current_size = current_size + 1                                        
         if( current_size .gt. auto_size ) newblk = .true.                      
         if( .not. newblk ) then                                                
           call inelbk_load_elem_props( element )                               
           call inelbk_chk_match( compatible )                                  
           if( .not. compatible ) newblk = .true.                               
         end if                                                                 
         if( .not. newblk ) then                                                
           elblks(0,nelblk)  = current_size                                     
           cycle                                                                
         endif                                                                  
         nelblk = nelblk + 1                                                    
         if( nelblk .gt. mxnmbl ) then                                          
            call inelbk_resize                                                      
c            call errmsg(74,param,dums,dumr,dumd)                                
c            call die_abort                                                      
         end if                                                                 
         felem             = element                                            
         elblks(1,nelblk)  = felem                                              
         elblks(0,nelblk)  = 1                                                  
         current_size       = 1                                                 
         call inelbk_load_block_props( felem )                                  
      end do ! on element                                                       
c                                                                               
c                     simplest assignment of blocks to domains.                 
c                     probably works well if user optimized element             
c                     number in code like Patran                                
c                                                                               
      if( auto_domains ) then                                                   
        auto_num_domains = 1                                                    
        if( use_mpi ) auto_num_domains = numprocs                               
        lblk = 1                                                                
        blks_per_domain = nelblk / auto_num_domains                             
        do idomain = 1, auto_num_domains                                        
          hblk = min( nelblk, idomain * blks_per_domain )                       
          if( idomain .eq. auto_num_domains ) hblk = nelblk                     
          elblks(2,lblk:hblk) = idomain - 1                                     
          lblk = hblk + 1                                                       
        end do                                                                  
      end if                                                                    
c                                                                               
c                     display blocking table if requested                       
c                                                                               
      if( .not. display ) return                                                
      write(out,9000) nelblk, auto_size                                         
      write(out,9010)                                                           
      do i = 1, nelblk                                                          
        write(out,9020) i, elblks(1,i),  elblks(0,i), elblks(2,i)               
      end do                                                                    
      write(out,*) ' '                                                          
c                                                                               
      if( .not. use_mpi ) return                                                
      counts_blks_ea_domain(0:max_procs-1,1:2) = 0                              
      do blk = 1, nelblk                                                        
       domain = elblks(2,blk)                                                   
       counts_blks_ea_domain(domain,1) =                                        
     &        counts_blks_ea_domain(domain,1) + 1                               
       counts_blks_ea_domain(domain,2) =                                        
     &        counts_blks_ea_domain(domain,2) + elblks(0,domain)                
      end do                                                                    
      write(out,9030)                                                           
      write(out,9040)                                                           
      do i = 0, numprocs-1                                                      
        write(out,9050)  i, counts_blks_ea_domain(i,1),                         
     &                   counts_blks_ea_domain(i,2)                             
      end do                                                                    
      write(out,*) ' '                                                          
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/,'>> Generated element blocking table:',                          
     & /,      1x,'  number of blocks, target size: ',i7,i5)                    
 9010 format(/,5x,                                                              
     &'block        1st element in blk        # elements in block',             
     &5x,'assigned domain' )                                                    
 9020 format(1x,i8, 10x,i9,25x,i4,10x,i6)                                       
 9030 format(/,'>> Domain statistics:')                                         
 9040 format(/,5x,'domain     # element blks    total # elements' )             
 9050 format(5x,i5,8x,i4,13x,i7)                                                
                                                                                
                                                                                
      contains                                                                  
c     ========                                                                  
                                                                                
c ********************************************************************          
c *                                                                  *          
c *    contains: inelbk_load_block_props                             *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine inelbk_load_block_props( now_elem )                            
      implicit none                                                             
      integer :: now_elem                                                       
c                                                                               
      blk_matmodel  = iprops(25,now_elem)                                       
      blk_eletype   = iprops(1,now_elem)                                        
      blk_intord    = iprops(5,now_elem)                                        
      blk_intnum    = iprops(6,now_elem)                                        
      blk_geonon    = lprops(18,now_elem)                                       
      blk_bbar      = lprops(19,now_elem)                                       
      blk_cohesive  = lprops(10,now_elem)                                       
      blk_matnum = iprops(38,now_elem)                                          
      blk_cyclic_plasticity = matmodel .eq. 5                                   
      blk_umat = matmodel .eq. 8                                                
      blk_dmg = iprops(42,now_elem) .ne. -1                                     
c                                                                               
      return                                                                    
      end subroutine  inelbk_load_block_props                                   
                                                                                
c ********************************************************************          
c *                                                                  *          
c *    contains: routine inelbk_load_elem_props                      *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine inelbk_load_elem_props( now_elem )                             
      implicit none                                                             
      integer :: now_elem                                                       
c                                                                               
      matmodel  = iprops(25,now_elem)                                           
      eletype   = iprops(1,now_elem)                                            
      intord    = iprops(5,now_elem)                                            
      intnum    = iprops(6,now_elem)                                            
      geonon    = lprops(18,now_elem)                                           
      bbar      = lprops(19,now_elem)                                           
      cohesive  = lprops(10,now_elem)                                           
      matnum = iprops(38,now_elem)                                              
      cyclic_plasticity = matmodel .eq. 5                                       
      umat = matmodel .eq. 8                                                    
      dmg = iprops(42, now_elem) .ne. -1                                        
c                                                                               
      return                                                                    
      end subroutine  inelbk_load_elem_props                                    
                                                                                
c ********************************************************************          
c *                                                                  *          
c *    contains: routine inelbk_chk_match                            *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine inelbk_chk_match( match )                                      
      implicit none                                                             
      logical match                                                             
c                                                                               
      match = .true.                                                            
c                                                                               
c                     in our simple blocking, the material                      
c                     associated with the element in the input must             
c                     be the same as for the block to be                        
c                     compatible. this is a more restrictive                    
c                     requirement than is rigidly required in                   
c                     WARP3D for the earliest developed material models.        
c                                                                               
c                     Example: elements in a block all use                      
c                     the "mises" model. The user defines multiple              
c                     materials with model mises where the materials            
c                     have different properties (E, nu, etc.)                   
c                     WARP3D allows elements in the same block to               
c                     be associated with the different "mises"                  
c                     materials. Later, more complex material models            
c                     did not support this capability (cyclic, umat,..)         
c                     as the data structures became too inflexible.             
c                                                                               
c                     Here, for simplicity, we require that elements            
c                     in a block must have the same user material               
c                     defined in the input.                                     
c                                                                               
c                                                                               
      if( blk_eletype .ne. eletype ) match = .false.                            
      if( .not. match ) return                                                  
c                                                                               
      if( blk_matnum .ne. matnum ) match = .false.                              
      if( .not. match ) return                                                  
c                                                                               
      if( blk_geonon .neqv. geonon ) match = .false.                            
      if( .not. match ) return                                                  
c                                                                               
      if( blk_intord .ne. intord ) match = .false.                              
      if( blk_intnum .ne. intnum ) match = .false.                              
      if( .not. match ) return                                                  
c                                                                               
      if( blk_bbar .neqv. bbar ) match = .false.                                
      if( .not. match ) return                                                  
c                                                                               
c                       Mark modification -- also need to check if              
c                       you are spoofing this material type for a               
c                       damage calculation or not.                              
c                                                                               
      if ( blk_dmg .neqv. dmg) match = .false.                                  
      if( .not. match ) return                                                  
      return                                                                    
      end subroutine inelbk_chk_match                                           
c                                                                               
      end subroutine inelbk_simple_blocking                                     
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine inelbk_resize                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 11/26/2018 rhd             *          
c     *                                                              *          
c     *               increase size of blocking table                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine inelbk_resize
      use global_data ! old common.main
      implicit none                                                             
c
c               current size is mxnmbl. reallocate to 1.x * current
c               size. copy in existing table, added rows.
c   
      integer :: new_mxnmbl, i
      integer, allocatable :: new_elblks(:,:)
      logical, parameter :: local_debug = .false.
c
      new_mxnmbl = int( 1.2 * mxnmbl ) + 5
      allocate( new_elblks(0:3,1:new_mxnmbl) )
      if( local_debug ) then
        write(out,*) "... resizing elblks. old, new: ", mxnmbl, 
     &                new_mxnmbl
      end if
c
      new_elblks(0:3,1:mxnmbl) =  elblks(0:3,1:mxnmbl)
c
      do i = mxnmbl + 1, new_mxnmbl
         new_elblks(0,i) = 0   ! number elements in block                           
         new_elblks(1,i) = 0   ! first element in block                             
         new_elblks(2,i) = 0   ! optional domain  number for element                
         new_elblks(3,i) = -1  ! doesn't seem to be used elsewhere     
      end do 

      if( numprocs > 1 ) then  ! mpi w/ domains used                               
         do i = 1, mxnmbl + 1, new_mxnmbl                                                       
            new_elblks(2,i) = -1                                                    
         end do                                                                 
      end if                                                                    
c
      if( local_debug ) write(out,*) '... moving allocation'
      call move_alloc( new_elblks, elblks )
      mxnmbl = new_mxnmbl
c
      return
      end

             
       
