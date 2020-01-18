c     ****************************************************************
c     *                                                              *
c     *                module allocated_integer_list                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/24/2019 rhd             *
c     *                                                              *
c     *     supplement for trlist. this version requires passing an  *
c     *     allocated vector to store the internal form of the list. *
c     *     the trlist_allocated routine resizes the passed          *
c     *     as needed to store the input list. this eliminates the   *
c     *     long standing challenge of defining a list vector        *
c     *     of some maximum possible length and then having the      *
c     *     input exceed that size.                                  *
c     *                                                              *
c     *     putting routines in a module eliminates need for         *
c     *     interface blocks where ever this routine is used         *
c     *                                                              *
c     ****************************************************************
c
      module allocated_integer_list
      use main_data, only : user_lists ! gfortran wants this here
!                                        not inside contains
      implicit none
c
      public  :: trlist_allocated
      private :: scan_list, list_resize
c
      contains
c     ****************************************************************
c     *                                                              *
c     *                subroutine trlist_allocated                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/2/2019 rhd              *
c     *                                                              *
c     *     alternative for trlist to support user-defined lists     *
c     *     the list storage is dynamically allocated and resized    *
c     *     as needed
c     *                                                              *
c     ****************************************************************
c
      subroutine trlist_allocated( list, mlist, iall, nlist, ierr )
c
c          scan action to input a list of integer terms.
c          the input can be a conventional integerlist or
c          a string containing the name of a previously defined
c          user list.
c
c          see trscan_list for all details on a conventional
c          integerlist.
c
c          dummy arguments
c              list    (output) - list of parsed input - as described
c                                 size increased here as needed so
c                                 no overflow ever occurs.
c              mlist   (out)    - size of list vector on return. 
c                                 incoming value not used. 
c              iall    (input)  - value of 'all'
c                                 = 0 - 'all' is not acceptable
c              nlist   (output) - number of terms stored in list
c              ierr    (output) - error code
c                                   = 1 - no error
c                                   = 2 - parse rules failded
c                                   = 3 - list overflow
c                                   = 4 - list not found
c
c
c          parsing rules:
c           - on entry, the calling routine has made the current scan entity
c             either the start of a conventional integer list or a string
c           - on exit we must have scan entity be the next item on line
c             after the string. this is compatible with the processing of
c             conventional integerlists.
c           - this routine does not touch the internal scan logical
c             flag tracking tru/false tests. here we do not know
c             what the user code expects so we leave exactly like simple
c             intergerlist
c
      implicit none
      include 'param_def'
c
      integer, allocatable :: list(:)
      integer ::  mlist, iall, nlist, ierr 
c
      integer :: i, idummy, nchars, list_col, stored_length
      character :: lname*24, name*80
      logical, external :: isstring, scanms
      logical, parameter :: debug = .false.
c
c          make sure the vector to store list has minimum size of 10
c          to start. could make this 20, 30, 50, ..
c          
      if( .not. allocated( list ) ) then
         mlist = 10
         allocate( list(mlist) ) 
      else
         mlist = size( list ) 
         if( mlist < 10 ) then
           deallocate( list )
           mlist = 10
           allocate( list(mlist) )
         end if
      end if
      list = 0  ! all entries
c
c          if we have a regular <integerlist> just process as before and
c          return. list can be increased in size. mlist has new size.
c          nlist has used number of entries in list.
c
      call scan_list( list, mlist, iall, nlist, ierr )
!      write(*,*) '.. back fronm scan_list. mlist, nlist:',mlist,nlist
!      if( nlist > 80 ) then
!!        do i = 1, nlist
!         write(*,*) i, list(i)
!         end do
!      end if
      if( ierr .ne. 4 ) return ! regular integer list
c
c          <integerlist> not found. check for user defined list name
c          in a string. scan entity is
c          already the string if it is there. use a scan function that
c          does not touch the internal scanner flag (next)
c
      if( .not. isstring(idummy) ) return ! don't know what we have
c
c          user defined list. find list in table, insert into list
c          space passed in. advance scanner to next entity on return to
c          match behavior of conventional integerlist.

      lname(1:24) = ' '; call entits( name, nchars )
      if( nchars > 24 ) nchars = 24; lname(1:nchars) = name(1:nchars)
      if( debug )  write(*,*) "... list id: ", lname
c
      list_col  = 0
      do i = 1, max_user_lists
       if( scanms( user_lists(i)%name, lname, 24 ) ) then
         list_col = i
         exit
       end if
      end do
c     
      ierr = 4  ! user list not found.
      if( list_col == 0 ) return
c
c          list found. check for overflow of space provided for
c          list.reallocate size as needed. extract values from 
c          stored lists and return.
c          put the next line entity into the scanner.
c
      stored_length = user_lists(i)%length_list
      if( stored_length == 0 ) then
         ierr = 2
         call ulist_error( 27 )
         call scan
         return
      end if
      if( mlist < stored_length ) then
        deallocate( list )
        allocate( list(stored_length) )
        mlist = stored_length
      end if
      nlist = stored_length
      list(1:nlist) = user_lists(list_col)%list(1:nlist)
      ierr = 1
      if( debug ) then
        write(*,*) '.. list_col, stored_length ', list_col, nlist
        write(*,*) 'list: ', list(1:nlist)
      end if
      call scan
      return
c
      end subroutine trlist_allocated

c *******************************************************************           
c *                                                                 *           
c * scan_list. same as trscan_list but with allocated list          *                                                       *           
c *                                                                 *           
c *******************************************************************           
c
      subroutine scan_list( list, mlist, iall, nlist, ierr ) 
      use scaner
c
      implicit none                
c                                                                               
c          scan action to input a list of integer terms                         
c              each action may be delimitted by a comma                         
c              a comma preceeding an eol indicates continuation                 
c              terms may be:                                                    
c                   a) <integer>                                                
c                   b) <integer1> - <integer2>                                  
c                   c) <integer1> to <integer2>                                 
c                   d) <integer1> - <integer2> by <integer3>                    
c                   e) <integer1> to <integer2> by <integer3>                   
c                   f) all                                                      
c              type a stores <integer>  in list                                 
c              type b and c stores <integer1>, -<integer2>, 1 in list           
c              type d and e stores <integer1>, -<integer2>, <integer3>          
c              in list                                                          
c              type f stores 1, -iall, 1 in list                                
c                                                                               
c         dummy arguments                                                       
c              list      (output) - list of parsed input - as described         
c              mlist     (input)  - allowable size of list                      
c              iall      (input)  - value of 'all'                              
c                                   = 0 - 'all' is not acceptable               
c              nlist     (output) - number of terms stored in list              
c              ierr      (output) - error code                                  
c                                   = 1 - no error                              
c                                   = 2 - parse rules failded                   
c                                   = 3 - list overflow 
c                                      no longer possible w/
c                                      list resizing here.                        
c                                   = 4 - list not found                        
c                                                                               
c         called subprograms                                                    
c              scan                                                             
c              rdline                                                           
c              scanmc                                                           
c                                                                               
c         local variables                                                       
c              istate            - fsa states                                   
c              nstate            - next state table                             
c              iclass            - class of input                               
c              ifsa              - action table                                 
c              iact              - action to do                                 
c              iby               - hollerth 'by'                                
c              ito               - hollerth 'to'                                
c              jall              - hollerth 'all'                               
c              istart            - flag set to true on first scan               
c              iovfl             - flag set to true on overflow                 
c                                                                               
c         algorithm terms                                                       
c              states   - 1 = start                                             
c                         2 = item (integer)                                    
c                         3 = delimiter (,)                                     
c                         4 = iteration (-, to)                                 
c                         5 = increment (by)                                    
c                         6 = delta (increment integer)                         
c              classes  - 1 = integer                                           
c                         2 = alpha 'to' or seperator '-'                       
c                         3 = alpha 'by'                                        
c                         4 = seperator ','                                     
c                         5 = end of line                                       
c                         6 = else                                              
c              actions  - 0 = switch state                                      
c                         1 = done                                              
c                         2 = error                                             
c                         3 = save inger, - value implies iteration te          
c                         4 = read line                                         
c                         5 = save -1*integer                                   
c                         6 = save 1, save integer                              
c                         7 = save 1                                            
c                         8 = save 1, done                                      
c                         9 = test overflow                                     
c                        10 = save iteration increment                          
c                                                                               
c      
      integer, allocatable :: list(:)
      integer :: mlist, iall, nlist, ierr                                                                         
      integer :: ifsa(6,6), nstate(6,6), istate, iclass, iact                             
      logical :: istart, iovfl                                                     
      logical, external :: scanmc                                                            
      real :: rby(1)  = real( Z'20207962', kind=kind(rby) ) ! by
      real :: rto(1)  = real( Z'20206F74', kind=kind(rto) ) ! to
      real :: rall(1) = real( Z'206C6C61', kind=kind(rall) ) ! all
c                                                                               
c         transfer table                                                        
c                                                                               
      data ifsa/                                                                
     1            3, 1, 1, 0, 1, 1,                                             
     2            3, 9, 2, 0, 1, 1,                                             
     3            3, 2, 2, 2, 4, 2,                                             
     4            5, 2, 2, 2, 2, 2,                                             
     5            6, 2, 0, 7, 8, 8,                                             
     6           10, 2, 2, 0, 4, 2/                                             
c                                                                               
c         next state table                                                      
c                                                                               
      data nstate/                                                              
     1            2, 0, 0, 1, 0, 0,                                             
     2            2, 4, 0, 3, 0, 0,                                             
     3            2, 0, 0, 0, 1, 0,                                             
     4            5, 0, 0, 0, 0, 0,                                             
     5            2, 0, 6, 3, 0, 0,                                             
     6            2, 0, 0, 6, 6, 0/                                             
c                                                                               
c         initialy in start state                                               
c                                                                               
      istart = .true.                                                           
      iovfl = .false.                                                           
      nlist = 1                                                                 
      istate = 1                                                                
c                                                                               
c         scan and determine class                                              
c                                                                               
 100  iclass = 6                                                                
      if(.not.istart)call scan                                                  
      if(mode.ne.9)go to 110                                                    
      iclass = 5                                                                
      go to 140                                                                 
 110  if(mode.ne.6)go to 120                                                    
      if(ivalue.eq.3)iclass = 2                                                 
      if(ivalue.eq.16)iclass = 4                                                
      go to 140                                                                 
 120  if( mode.ne.3 )go to 130                                                  
      if( scanmc(entity,rby,2) .and. nchar.eq.2 ) iclass = 3                 
      if( scanmc(entity,rto,2) .and. nchar.eq.2 ) iclass = 2                 
      if( scanmc(entity,rall,3) .and. istart .and. iall.ne.0                 
     &    .and. nchar.eq.3 ) go to 350                                          
      go to 140                                                                 
 130  if(mode.ne.1)go to 140                                                    
      iclass = 1                                                                
c                                                                               
c         determine next state and action                                       
c                                                                               
 140  iact = ifsa(iclass,istate)+1                                              
      istate = nstate(iclass,istate)                                            
      if(iact.ne.2)istart=.false.                                               
c                                                                               
c         action transfer, skip store on overflow                               
c                                                                               
      if(iovfl.and.(iact.eq.4.or.iact.eq.6.or.iact.eq.7.or.                     
     1              iact.eq.8.or.iact.eq.10))go to 100                          
      go to (100, 210, 220, 230, 250, 260, 270, 280, 290, 310, 320              
     1       ), iact                                                            
c                                                                               
c         done                                                                  
c                                                                               
 210  ierr = 1                                                                  
      nlist = nlist-1                                                           
      if(istart)ierr = 4                                                        
      if(iovfl)ierr = 3                                                         
      if(list(1).lt.0.and.nlist.gt.1)ierr = 2                                   
      return                                                                    
c                                                                               
c         error - parse                                                         
c                                                                               
 220  ierr = 2                                                                  
      nlist = nlist-1                                                           
      return                                                                    
c                                                                               
c         save integer, -value becomes iteration bound                          
c                                                                               
 230  if(nlist.le.mlist)go to 240                                               
      call list_resize( list, mlist )
 240  list(nlist) = ivalue                                                      
      nlist = nlist+1                                                           
      if(ivalue.lt.0.and.nlist.eq.2)go to 100                                   
      if(ivalue.lt.0)istate = 5                                                 
      go to 100                                                                 
c                                                                               
c         read line                                                             
c                                                                               
 250  call rdline                                                               
      go to 100                                                                 
c                                                                               
c         save -1*integer                                                       
c                                                                               
 260  list(nlist) = -ivalue                                                     
      if(ivalue.lt.0)go to 220                                                  
      nlist = nlist+1                                                           
      go to 100                                                                 
c                                                                               
c         save 1, save integer                                                  
c                                                                               
 270  if( nlist > mlist ) call list_resize( list, mlist ) 
      list(nlist) = 1      
      if( nlist + 1 > mlist ) call list_resize( list, mlist )                                                     
      list(nlist+1) = ivalue                                                    
      nlist = nlist+2                                                           
      go to 100                                                                 
c                                                                               
c         save 1                                                                
c                                                                               
 280  list(nlist) = 1                                                           
      nlist = nlist+1                                                           
      go to 100                                                                 
c                                                                               
c         save 1, done                                                          
c                                                                               
 290  if(iovfl)go to 300                                                        
      list(nlist) = 1                                                           
 300  ierr = 1                                                                  
      if(iovfl)ierr = 3                                                         
      if(list(1).lt.0.and.nlist.ne.1)ierr = 2                                   
      return                                                                    
c                                                                               
c         test for overflow in iteration term                                   
c                                                                               
 310  if(nlist+1.gt.mlist) call list_resize( list, mlist )
      go to 100                                                                 
c                                                                               
c         save integer for increment                                            
c                                                                               
 320  if(nlist.le.mlist)go to 330                                               
      call list_resize( list, mlist ) 
 330  list(nlist) = ivalue                                                      
      nlist = nlist+1                                                           
      go to 100                                                                 
c                                                                               
c        all, done                                                              
c                                                                               
 350  if(mlist.lt.3)go to 360                                                   
      list(1) = 1                                                               
      list(2) = -iall                                                           
      list(3) = 1                                                               
      nlist = 3                                                                 
      ierr = 1                                                                  
      call scan                                                                 
      if(mode.ne.6)return                                                       
      if(ivalue.ne.16)return                                                    
      call scan                                                                 
      if(mode.ne.9)return                                                       
      call rdline                                                               
      call scan                                                                 
      return                                                                    
 360  ierr = 3                                                                  
      nlist = 0                                                                 
      return                                                                    
c                                                                               
      end subroutine scan_list                                                                    

c *******************************************************************           
c *                                                                 *           
c * list_resize                                                     *           
c *                                                                 *           
c *******************************************************************           
c
      subroutine list_resize( list, size )
      implicit none
c
      integer, allocatable :: list(:)
      integer :: size
c
      integer, allocatable :: bigger_list(:)
      integer :: new_size
c
      new_size = size*2    !  max( size * 2, 100 )
      allocate( bigger_list(new_size) )
      bigger_list(1:size) = list(1:size)!!
      call move_alloc( bigger_list, list )
!      write(*,*) ' <<<< list_resize. old, new:',size,new_size
      size = new_size
c
      return
      end subroutine list_resize

      end module

      