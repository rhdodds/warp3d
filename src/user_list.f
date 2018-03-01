c     ****************************************************************
c     *                                                              *
c     *                      subroutine trwlist                      *
c     *                                                              *
c     *                 written by : dw and rhd                      *
c     *                                                              *
c     *          process a user-defined/named integer list           *
c     *                                                              *
c     ****************************************************************
c
      subroutine trwlist( sbflg1, sbflg2 )
      use global_data ! old common.main
      use main_data, only : user_lists, crdmap
      implicit integer (a-z)
      logical  :: debug, string, scanms, matchs, do_display,
     &        found, match_exact, matchs_exact, true, display_coords,
     &        endcrd
      character :: lname*24, name*80
c
c                 on entry, last scanner test was true (word
c                 list). get string id for user-defined list name.
c                 if no string found, error message and return.
c
      debug = .false.
      if( matchs_exact( "debug" ) )  debug = .true.
      if( debug ) write(out,*) "...inside trwlist ..."
      if( .not. string( idummy ) ) then
        call ulist_error( 1 ); call scan_flushline; return
      end if
c
c                 check if list already exists and/or find an available column
c                 for new list in table. if just display request, show list and
c                 leave. lists are never deleted from the table. just defined an
c                 possibly over written.
c
      lname(1:24) = ' '
      call entits( name, nchars )
      if( nchars .gt. 24 ) nchars = 24
      lname(1:nchars) = name(1:nchars)
      if( debug )  write(out,*) "... list id: ", lname
c
      list_col  = 0; avail_col = 0
      do i = 1, max_user_lists
       if( user_lists(i)%length_list .eq. 0 ) then
          if( avail_col .eq. 0 ) avail_col= i
       end if
       if( scanms( user_lists(i)%name, lname, 24 ) ) then
         list_col = i; exit
       end if
      end do
c
      if( matchs( 'display', 4 ) ) then
        if(  list_col .eq. 0 ) then
          call ulist_error( 2 )
        else
          display_coords = .false.
          if( matchs( 'coordinates', 4 ) ) display_coords = .true.
          call ulist_display( list_col, c(1), display_coords,
     &                        debug, out )
        end if
        return
      end if
c
c                 new list definition in process. If already exists,
c                 delete and re-define
c
      if( list_col .ne. 0 ) then
        call ulist_error( 3 )
        user_lists(list_col)%length_list = 0
        if( allocated(user_lists(list_col)%list) )
     &       deallocate(user_lists(list_col)%list)
      else
        list_col = avail_col; user_lists(list_col)%name = lname
      end if
c
c                 process remainder of command to construct
c                 the integer list. if nothing left on command,
c                 error and leave. check for different forms
c                 of list command, call lower routine to proess.
c
      if( endcrd() ) then
        call ulist_error( 4 ); call scan_flushline; return
      end if
c
      do_display = .false.; found = .false.; display_coords = .false.
c
      if( matchs_exact( 'x' ) ) then
        found = .true.
      elseif( matchs_exact( 'y' ) ) then
        found = .true.
      elseif( matchs_exact( 'z' ) )  then
        found = .true.
      end if
      if( found ) then
        call backsp( 1 )
        if( true( idummy ) ) call splunj
        call ulist_general( list_col, debug, do_display,
     &                      display_coords, nonode )
        go to 100
       end if
c
      if( matchs( 'cylinder', 6 ) ) then
        call ulist_cylinder( list_col, debug, do_display,
     &                       display_coords, nonode )
        go to 100
      end if
c
      if( matchs( 'plane', 5 ) ) then
        call ulist_plane( list_col, debug, do_display,
     &                    display_coords, nonode )
        go to 100
      end if
c
      if( matchs( 'sphere', 5 ) ) then
        call ulist_sphere( list_col, debug, do_display,
     &                     display_coords, nonode )
        go to 100
      end if
c
      call ulist_simple( list_col, debug, out,
     &                   4*max(nonode,noelem,100) )
      if( matchs( 'display',6 ) ) do_display = .true.
      if( matchs( 'coordinates', 4 ) ) display_coords = .true.
c
 100  continue
      if( do_display ) call ulist_display( list_col, c(1),
     &                            display_coords, debug, out )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ulist_simple                 *
c     *                                                              *
c     *                       written by : dw and rhd                *
c     *                                                              *
c     *          scan and store a conventional integer list          *
c     *                                                              *
c     ****************************************************************
c
      subroutine ulist_simple( list_col, debug, out, max_list )
      use main_data, only : user_lists
      implicit integer (a-z)
      logical debug
      integer, allocatable :: local_list(:)
c
c                 the default option for user defined integer list
c                 when no geometric selection option is specified.
c                 input must be a regular interger list at this point.
c                 value of 'all' is not acceptable since no context
c                 for what all means (nodes, elements, etc) is
c                 known.
c
      if( debug ) write(out,*) "..ulist_simple called..."
      allocate( local_list(max_list) )
c
c                 get an integer list and check for errors. "all"
c                 not allowed since we would not know all of
c                 what at this point?
c
      call trscan_list( local_list, max_list, 0, nlist, ierr )
      if( ierr .ne. 1 )  then
        if( ierr .eq. 2 )  call ulist_error( 5 )
        if( ierr .eq. 3 )  call ulist_error( 6 )
        if( ierr .eq. 4 )  call ulist_error( 7 )
        call scan_flushline
        return
      end if
c
c                 store the list into permanent data structure
c
      user_lists(list_col)%length_list = nlist
      allocate( user_lists(list_col)%list(nlist) )
      user_lists(list_col)%list(1:nlist) = local_list(1:nlist)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ulist_display                *
c     *                                                              *
c     *                       written by : dw and rhd                *
c     *                                                              *
c     *          display entries in a user-defined integer list      *
c     *                                                              *
c     ****************************************************************
c
      subroutine ulist_display( list_col, coords, display_coords,
     &                          debug, out )
      use main_data, only : user_lists, crdmap
      implicit integer(a-z)
      logical debug, display_coords
      double precision :: coords(*)
c
      integer line(9)
      double precision ::  xcoord, ycoord, zcoord
c
      if( debug ) write(out,*) "...ulist_display is called..."
c
      length = user_lists(list_col)%length_list
      if( length .eq. 0 ) then
         call ulist_error( 8 ); return
      end if
c
c                 loop to extract each entry from the compressed
c                 integer list format. write 9 values per line
c
      icn = 0; iplist = 1;count = 0; line_no = 0
      write(out,9000)
c
      do
       if( iplist .eq. 0 ) exit
       call trxlst( user_lists(list_col)%list(1), length, iplist,
     &              icn, next)
       count = count + 1; line(count) = next
       if( count .eq. 9 ) then
           line_no = line_no + 1; write(out,9100) line_no, line
           count = 0
        end if
      end do
c
      line_no = line_no + 1
      if( count .ne. 0 ) write(out,9100) line_no, line(1:count)
c
      if( .not. display_coords ) then
         write(out,9200)
         return
      end if
c
      icn = 0; iplist = 1
      write(out,9220)
c
      do
       if( iplist .eq. 0 ) exit
       call trxlst( user_lists(list_col)%list(1), length, iplist,
     &              icn, next)
       node = next
       position = crdmap(node)
       xcoord = coords(position+0); ycoord = coords(position+1);
       zcoord = coords(position+2)
       write(out,9210) node, xcoord, ycoord, zcoord
      end do
c
      write(out,9200)

      return
c
 9000 format(1x,">>>>> Entries in list ....")
 9100 format(1x,i6,":",9i9)
 9200 format(1x,"..... End of list ....",//)
 9210 format(10x,i7,3f20.9)
 9220 format(/,1x,"      Coordinates of nodes in list ....")
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ulist_error                  *
c     *                                                              *
c     *                       written by : dw and rhd                *
c     *                                                              *
c     *      error message routine for user-define integer lists     *
c     *                                                              *
c     ****************************************************************
c
      subroutine ulist_error( message )
      use global_data ! old common.main
      use main_data, only : user_lists
      implicit integer(a-z)
      character(len=80) :: string
c
      select case( message )
      case( 1 )
         call entits( string, strlng )
         write(out,9001) string(1:strlng)
         num_error = num_error + 1
      case( 2 )
         write(out,9002)
      case( 3 )
         write(out,9003)
      case( 4 )
         write(out,9004)
         num_error = num_error + 1
      case( 5 )
         write(out,9005)
         num_error = num_error + 1
      case( 6 )
         write(out,9006)
         num_error = num_error + 1
      case( 7 )
         call entits( string, strlng )
         write(out,9007) string(1:strlng)
         num_error = num_error + 1
      case( 8 )
         write(out,9008)
      case( 9 )
         call entits( string, strlng )
         write(out,9009) string(1:strlng)
         num_error = num_error + 1
      case( 10 )
         call entits( string, strlng )
         write(out,9010) string(1:strlng)
         num_error = num_error + 1
      case( 11 )
         call entits( string, strlng )
         write(out,9011) string(1:strlng)
         num_error = num_error + 1
      case( 12 )
         call entits( string, strlng )
         write(out,9012) string(1:strlng)
         num_error = num_error + 1
      case( 13 )
         call entits( string, strlng )
         write(out,9013) string(1:strlng)
         num_error = num_error + 1
      case( 14 )
         call entits( string, strlng )
         write(out,9014) string(1:strlng)
         num_error = num_error + 1
      case( 15 )
         call entits( string, strlng )
         write(out,9015) string(1:strlng)
         num_error = num_error + 1
      case( 16 )
         call entits( string, strlng )
         write(out,9016) string(1:strlng)
         num_error = num_error + 1
      case( 17 )
         write(out,9017)
         num_error = num_error + 1
      case( 18 )
         write(out,9018)
         num_error = num_error + 1
      case( 19 )
         call entits( string, strlng )
         write(out,9019) string(1:strlng)
         num_error = num_error + 1
      case( 20 )
         call entits( string, strlng )
         write(out,9020) string(1:strlng)
         num_error = num_error + 1
      case( 21 )
         call entits( string, strlng )
         write(out,9021) string(1:strlng)
         num_error = num_error + 1
      case( 22 )
         call entits( string, strlng )
         write(out,9022) string(1:strlng)
         num_error = num_error + 1
      case( 23 )
         write(out,9023)
         num_error = num_error + 1
      case( 24 )
         write(out,9024)
         num_error = num_error + 1
      case( 25 )
         call entits( string, strlng )
         write(out,9025) string(1:strlng)
         num_error = num_error + 1
      case( 26 )
         write(out,9026)
         num_error = num_error + 1
      case( 27 )
         call entits( string, strlng )
         write(out,9027) string(1:strlng)
      case default
        write(out,9999)
        stop
      end select
c
      return
c
 9001 format(/1x,'>>>>> error: list name expected. must be a string',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9002 format(/1x,'>>>>> warning: list name not found',
     & /14x,'command ignored...')
 9003 format(/1x,'>>>>> warning: existing list over-written')
 9004 format(/1x,'>>>>> error: unexpected end of command',
     & /14x,'command ignored...')
 9005 format(/1x,'>>>>> error: integer list rules failed',
     & /14x,'command ignored...')
 9006 format(/1x,'>>>>> error: list too long ',
     & /14x,'command ignored...')
 9007 format(/1x,'>>>>> error: expecting integer list or list option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9008 format(/1x,'>>>>> warning: specified list is empty')
 9009 format(/1x,'>>>>> error: expecting coordinate value',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9010 format(/1x,'>>>>> error: expecting tolerance value',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9011 format(/1x,'>>>>> error: unknown option in list command',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9012 format(/1x,'>>>>> error: expecting point 1 for cylinder option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9013 format(/1x,'>>>>> error: expecting keyword second to start ',
     & 'point 2 for cylinder option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9014 format(/1x,'>>>>> error: expecting radius value for cylinder',
     & ' or sphere option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9015 format(/1x,'>>>>> error: expecting an end-of-line ',
     & 'for cylinder or sphere option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9016 format(/1x,'>>>>> error: expecting keyword radius ',
     & 'for cylinder or sphere option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9017 format(/1x,'>>>>> error: all coordinates not given for 2 points',
     & ' on the cylinder axis',
     & /14x,'or radius not positive.... command ignored...')
 9018 format(/1x,'>>>>> error: cylinder axis length < 1.0d-10',
     & /14x,'command ignored...')
 9019 format(/1x,'>>>>> error: expecting point on plane',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9020 format(/1x,'>>>>> error: expecting definition of normal to plane',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9021 format(/1x,'>>>>> error: expecting component of normal vector',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9022 format(/1x,'>>>>> error: expecting an end-of-line ',
     & 'for plane option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9023 format(/1x,'>>>>> error: all coordinates not given for point',
     & ' on plane or components of normal vector',
     & /14x,'command ignored...')
 9024 format(/1x,'>>>>> error: normal vector length < 1.0d-10',
     & /14x,'command ignored...')
 9025 format(/1x,'>>>>> error: expecting origin value for sphere',
     & ' option',
     & /14x,'scanning: ',a,
     & /14x,'command ignored...')
 9026 format(/1x,'>>>>> error: coordinates of origin not specified',
     & /14x,'or radius not positivecommand ignored...')
 9027 format(/1x,'>>>>> warning: the user-list named: ',a,
     & /16x,'has no entries...')
 9999 format(/1x,'>>>>> Fatal Error: routine ulist_error.',
     &   /16x,   'should have not reach this point.')
c
      end

c     ****************************************************************
c     *                                                              *
c     *               subroutine ulist_general                       *
c     *                                                              *
c     *                   written by : dw and rhd                    *
c     *                                                              *
c     *        make a list of nodes that lie on a specified          *
c     *        cartesian plane, along a straight line or at a point  *
c     *                                                              *
c     ****************************************************************

c
      subroutine ulist_general( list_col, debug, do_display,
     &                          display_coords, max_list_size )
      use global_data ! old common.main
c
c       list "list name" [x/y/z (=) _] (tolerance(=)_) (display)
c
      use main_data, only : user_lists, crdmap

      implicit integer(a-z)
c
      double precision
     &   value, xvalue, yvalue, zvalue, tolerance, temp_t, zero,
     &   xcoord, ycoord, zcoord, x_toler, y_toler, z_toler
      logical :: do_display, xdefined, ydefined, zdefined, matchs,
     &           numd, endcrd, matchs_exact, debug, display_coords
      dimension matchlist(max_list_size)
      data zero /0.0d00/
c
      if( debug ) write(out,*) "...in ulist_line_point..."
c
      xdefined = .false.; ydefined = .false.; zdefined = .false.
      xvalue = zero; yvalue = zero; zvalue = zero
      do_display = .false.; tolerance = 0.000001
c
c                 loop to extract values for line/point. display
c                 request at end will be handled by driver.
c
      do
c
        if( endcrd())  exit  ! scan loop
        if( matchs_exact( "," ) ) call splunj
c
        if( matchs_exact( "x" ) ) then
          if( matchs_exact( "=" ) ) call splunj
          if( numd( xvalue ) ) then
            xdefined = .true.
            cycle
          else
            call ulist_error( 9 ); call scan_flushline; return
          end if
        end if
c
        if( matchs_exact( "y" ) ) then
         if( matchs_exact( "=" ) ) call splunj
         if( numd( yvalue ) ) then
           ydefined = .true.
           cycle
         else
           call ulist_error( 9 ); call scan_flushline; return
         end if
        end if
c
        if( matchs_exact( "z" ) ) then
         if( matchs_exact( "=" ) ) call splunj
         if( numd( zvalue ) ) then
           zdefined = .true.
           cycle
         else
           call ulist_error( 9 ); call scan_flushline; return
         end if
        end if
c
        if( matchs( "tolerance", 5 ) ) then
          if( matchs_exact( "=" ) ) call splunj
          if( numd( temp_t ) ) then
            tolerance = temp_t
          else
            call ulist_error( 10 ); call scan_flushline; return
          end if
          cycle
        end if
c
        if( matchs( "display", 6 ) ) then
          do_display = .true.
          if( matchs( "coordinates", 4 ) ) display_coords = .true.
          cycle
        end if
c
c                 no match to option of list command. flush and return
c
        call ulist_error( 11 ); call scan_flushline; return
c
      end do
c
      if( debug ) then
         write(out,9000) xdefined, ydefined, zdefined,
     &                   xvalue, yvalue, zvalue, tolerance,
     &                   do_display
      end if
c
      call ulist_toler( tolerance, x_toler, y_toler, z_toler,
     &                  nonode, crdmap, c, debug, out )
      write(out,9220)  x_toler, y_toler, z_toler
c
c                 examine all model nodes to compare against the
c                 point or line coordinate specification
c                 provided by user. compare with relative
c                 tolerances. if flag = 3 for node, it matched
c                 the tests, then store node number in list. note
c                 that matchlist to store nodes that match is
c                 an automatic vector.
c
      count = 0
      do node = 1, nonode
        flag = 0; position = crdmap(node)
        xcoord = c(position+0); ycoord = c(position+1);
        zcoord = c(position+2)
c
        if( xdefined ) then
          if( abs(xcoord - xvalue) <= x_toler ) then
             flag = flag + 1
          end if
        else
             flag = flag + 1
        end if
c
        if( ydefined ) then
          if( abs( ycoord - yvalue) <= y_toler ) then
             flag = flag + 1
          end if
        else
             flag = flag + 1
        end if
c
        if( zdefined ) then
          if( abs( zcoord - zvalue) <= z_toler ) then
             flag = flag + 1
          end if
        else
             flag = flag + 1
        end if
c
        if( flag == 3 ) then
          count = count + 1
          matchlist(count) = node
        end if
c
      end do
c
c                 if we found nodes that matched the specification
c                 store vlaues as an integer list in user defined lists
c                 table.
c
      if( count > 0 ) then
         user_lists(list_col)%length_list = count
         allocate( user_lists(list_col)%list(count) )
         do i = 1, count
          user_lists(list_col)%list(i) = matchlist(i)
         end do
         write(out,9200) count
       else
         write(out,*) " "
         write(out,9210)
         write(out,*) " "
         do_display = .false.
       end if
c
       return
c
 9000 format(2x,".... scanned values for lines option ....",
     & /,5x,"xdefined, ydefined, zdefined:",3l2,
     & /,5x,'x, y, z values: ',3f15.5,
     & /,5x,'tolerance, do_display: ',f15.8, l5,
     & //)
 9100 format(2x,".... tolerances for comparing coordinates ....",
     & /,5x,"nonode: ",i8,
     & /,5x,"x, y, z tolerance: ", 3e15.8, //)
 9200 format(1x,">>>>> number of nodes placed inlist: ", i8 )
 9210 format(1x,">>>>> WARNING: no nodes match the line/point ",
     & " specification")
 9220 format(1x,">>>>> (x,y,z) distances for proximity checks: ",
     &      3f15.8 )
c
      end

c     ****************************************************************
c     *                                                              *
c     *               subroutine ulist_toler                         *
c     *                                                              *
c     *              written by : dw and rhd                         *
c     *                updated:  8/12/2014 rhd. change to imp none   *s
c     *                                                              *
c     *        make a list of nodes that lie on a specified          *
c     *        or a plane  straight lines                            *
c     *                                                              *
c     ****************************************************************

c
      subroutine ulist_toler( tolerance, x_toler, y_toler, z_toler,
     &                        nonode, crdmap, coords, debug, out )
c
      implicit none
c
c                 arguments
c
      integer :: crdmap(*), out, nonode
      double precision ::
     & tolerance, x_toler, y_toler, z_toler, coords(*)
      logical :: debug
c
c                 locals
c
      double precision ::
     &  xcoord, ycoord, zcoord, xmin, xmax, ymin, ymax,
     &  zmin, zmax
      integer :: position, node
c
      if( debug ) write(out,*) "...ulist_toler..."
c
c                 user provides an absolute tolerance value.
c                 find size of box enclosing structure. construct
c                 relative tolerance for testing from box size
c                 * user absolute tolerance
c
      xmin =  1.0d20; xmax = -1.0d20; ymin =  1.0d20; ymax = -1.0d20
      zmin =  1.0d20; zmax = -1.0d20
c
      do node = 1, nonode
        position = crdmap(node)
        xcoord   = coords(position+0)
        ycoord   = coords(position+1)
        zcoord   = coords(position+2)
        xmin     = min( xmin, xcoord )
        xmax     = max( xmax, xcoord )
        ymin     = min( ymin, ycoord )
        ymax     = max( ymax, ycoord )
        zmin     = min( zmin, zcoord )
        zmax     = max( zmax, zcoord )
      end do
c
      x_toler = max( abs( xmax-xmin ) * tolerance, tolerance )
      y_toler = max( abs( ymax-ymin ) * tolerance, tolerance )
      z_toler = max( abs( zmax-zmin ) * tolerance, tolerance )
c
      if( debug ) write(out,9100) nonode, x_toler, y_toler, z_toler,
     &                            tolerance
      return

 9100 format(2x,".... tolerances for comparing coordinates ....",
     & /,5x,"nonode: ",i8,
     & /,5x,"x, y, z tolerance: ", 4e15.8, //)
c
      end

c     ****************************************************************
c     *                                                              *
c     *               subroutine ulist_cylinder                      *
c     *                                                              *
c     *                   written by : dw and rhd                    *
c     *                   last modified: 8/12/2014 rhd               *
c     *                                                              *
c     *        make a list of nodes that lie on a specified          *
c     *        cylinder (w/ open ends)                               *
c     *                                                              *
c     ****************************************************************

c
      subroutine ulist_cylinder( list_col, debug, do_display,
     &                           display_coords, max_list_size )
      use global_data ! old common.main
c
c      list "list name" cylinder first x(=)_ y(=)_ z(=)_ second
c             x(=)_ y(=)_ z(=)_ radius(=)_ (tolerance(=)_) (display)
c
      use main_data, only : user_lists, crdmap

      implicit integer(a-z)
c
      double precision
     &   tolerance, temp_t, zero,
     &   pt1(3), pt2(3), pt3(3),
     &   radius, flag, axis(3), axis_length, radius_tol,
     &   dx, dy, dz, d13, nx13, ny13, nz13, cos_theta, perp_dist,
     &   one, sin_theta
      logical :: do_display, matchs, check,
     &           numd, endcrd, matchs_exact, debug, display_coords
      dimension matchlist(max_list_size)
      data zero /0.0d00/, flag /-999999.0/, one /1.0d00/
c
      if( debug ) write(out,*) "...in ulist_cylinder..."
c
c                 loop to extract values for line/point. display
c                 request at end will be handled by driver.
c
      do_display = .false.; tolerance = 0.000001
      pt1(1:3) = flag; pt2(1:3) = flag
c
c                  -- must be input to define pt 1 on generator axis
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs_exact( "first" ) ) then
        call ulist_error( 12 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "x" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt1(1) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "y" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt1(2) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "z" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt1(3) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
      end if
c
c                  -- must be input to define pt 2 on generator axis
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs( "second", 5 ) ) then
        call ulist_error( 13 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "x" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt2(1) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "y" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt2(2) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "z" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt2(3) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
      end if
c
c                  -- must be radius of cyclinder
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs( "radius", 6 ) ) then
        call ulist_error( 16 ); call scan_flushline; return
      end if
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( radius ) ) then
        call ulist_error( 14 ); call scan_flushline; return
      end if
c
c                  -- optional tolerance
c
      if( matchs_exact( "," ) ) call splunj
      if( matchs( "tolerance", 5 ) ) then
        if( matchs_exact( "=" ) ) call splunj
        if( numd( temp_t ) ) then
          tolerance = temp_t
        else
          call ulist_error( 10 ); call scan_flushline; return
        end if
      end if
c
c                  -- optional display request
c
      if( matchs_exact( "," ) ) call splunj
      if( matchs( "display", 6 ) ) do_display = .true.
      if( matchs( 'coordinates', 4 ) ) display_coords = .true.
c
c                  -- if not an end of input at this point,
c                     we have an un-recognized item
c
      if( .not. endcrd() ) then
        call ulist_error( 15 ); call scan_flushline
          do_display = .false.; return
      end if
c
      if( debug ) write(out,9000) pt1, pt2, radius, tolerance,
     &                            do_display
c
c                  get a relative tolerance to use in comparisons
c                  just absolute tol  x radius for cyclinder.
c
c                  all coordinates for pts 1 and 2 on the cylinder
c                  axis must be specifed - no defaults. Axis length
c                  cannot be very near zero.
c
c                  normalize vector defining direction of cylinder
c                  axis.
c
      radius_tol = tolerance * radius
      write(out,9220) radius_tol
      check = pt1(1) .eq. flag .or. pt1(2) .eq. flag .or.
     &        pt1(3) .eq. flag .or. pt2(1) .eq. flag .or.
     &        pt2(2) .eq. flag .or. pt2(3) .eq. flag
     &        .or. radius <= zero
      if( check ) then
         call ulist_error( 17 ); do_display = .false.; return
      end if
c
      axis(1) = pt2(1) - pt1(1)
      axis(2) = pt2(2) - pt1(2)
      axis(3) = pt2(3) - pt1(3)
      axis_length = axis(1)*axis(1) + axis(2)*axis(2) + axis(3)*axis(3)
      axis_length = sqrt(axis_length)
      if( axis_length < 1.0d-10 ) then
         call ulist_error( 18 ); do_display = .false.; return
      end if
      axis(1) = axis(1) / axis_length
      axis(2) = axis(2) / axis_length
      axis(3) = axis(3) / axis_length
      if( debug ) write(out,9010) axis, radius_tol
c
c                 examine all model nodes.
c                  - compute components of vector from pt1 on
c                    cyl axis to the node (ie, pt3) and length
c                  - skip node if this distance < radius_tol (node
c                    is v. close to cyl axis pt1)
c                  - unit vector from pt1 to node
c                  - cosine of angle between cyl axis and line from
c                    pt1 on axis to node. axis direction and line both
c                    normalized just prior to unit length
c                  - sine of angle between pt1 and node
c                  - perpendicular distance from node to cyl axis =
c                    sin_angle * d13
c                  - check distance against radius tolerance
c
      count = 0
      do node = 1, nonode
        position = crdmap(node); pt3(1)  = c(position+0)
        pt3(2)  = c(position+1); pt3(3)  = c(position+2)
        dx = pt3(1) - pt1(1); dy = pt3(2) - pt1(2); dz = pt3(3) - pt1(3)
        d13 = sqrt( dx*dx + dy*dy + dz*dz )
        if( d13 < radius_tol ) cycle
        nx13 = dx / d13; ny13 = dy / d13; nz13 = dz / d13
        cos_theta = axis(1)*nx13 + axis(2)*ny13 + axis(3)*nz13
        sin_theta = sqrt( one - cos_theta*cos_theta )
        perp_dist = d13*sin_theta
        if( abs( perp_dist-radius ) <= radius_tol ) then
          count = count + 1; matchlist(count) = node
        end if
      end do
c
c                 if we found nodes that matched the specification
c                 store vlaues as an integer list in user defined lists
c                 table.
c
      if( count > 0 ) then
         user_lists(list_col)%length_list = count
         allocate( user_lists(list_col)%list(count) )
         do i = 1, count
          user_lists(list_col)%list(i) = matchlist(i)
         end do
         write(out,9200) count
       else
         write(out,9210)
       end if
c
       return
c
 9000 format(2x,".... scanned values for cylinder option ....",
     & /,5x,"pt1(x,y,z): ",3f15.6,
     & /,5x,"pt2(x,y,z): ",3f15.6,
     & /,5x,"radius, tolerance: ",2f15.6,
     & /,5x,"do_display: ",l5,
     & //)
 9010 format(2x,".... unit vector for cylinder axis ....",
     & /,5x,"nx, ny, nz: ",3f15.6,
     & /,5x,"radius_tol: ",f15.6,
     & //)
 9100 format(2x,".... tolerances for comparing coordinates ....",
     & /,5x,"nonode: ",i8,
     & /,5x,"x, y, z tolerance: ", 3e15.8, //)
 9200 format(1x,">>>>> number of nodes placed in list: ", i8 )
 9210 format(1x,">>>>> no nodes match the cylinder specification")
 9220 format(1x,">>>>> physical distance for proximity checks: ",
     &       f15.8)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                     subroutine trlist                        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/26/2011                 *
c     *                                                              *
c     *     replacement for trlist to support user-defined lists     *
c     *                                                              *
c     ****************************************************************
c
      subroutine trlist( list, mlist, iall, nlist, ierr )
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
c              list      (output) - list of parsed input - as described
c              mlist     (input)  - allowable size of list
c              iall      (input)  - value of 'all'
c                                   = 0 - 'all' is not acceptable
c              nlist     (output) - number of terms stored in list
c              ierr      (output) - error code
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
      use main_data, only : user_lists
      implicit integer(a-z)
      include 'param_def'
      dimension list(*)
      character lname*24, name*80
      logical isstring, scanms, debug
      data debug / .false. /
c
c          if we have a regular <integerlist> just process as before and
c          return
c
      call trscan_list( list, mlist, iall, nlist, ierr )
      if( ierr .ne. 4 ) return
c
c          <integerlist> not found. check for user defined list name
c          in a string. scan entity is
c          already the string if it is there. use a scan function that
c          does not touch the internal scanner flag (next)
c
      if( .not. isstring(idummy) ) return
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
         list_col = i; go to 100
       end if
      end do
c
      ierr = 4
      return
c
c          list found. check for overflow of space provided for
c          list. extract values from stored lists and return.
c          put the next line entity into the scanner.
c
 100  continue
      stored_length = user_lists(i)%length_list
      if( stored_length == 0 ) then
         ierr = 2
         call ulist_error( 27 )
         call scan
         return
      end if
      if( mlist < stored_length ) then
        ierr = 3; return
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
      end
c     ****************************************************************
c     *                                                              *
c     *               subroutine ulist_plane                         *
c     *                                                              *
c     *              written by : dw and rhd                         *
c     *              last modified: 8/12/2014 rhd                    *
c     *                                                              *
c     *    make a list of nodes that lie on a specified plane        *
c     *                                                              *
c     ****************************************************************

c
      subroutine ulist_plane( list_col, debug, do_display,
     &                        display_coords, max_list_size )
      use global_data ! old common.main
c
c     list "name" plane point x(=)_ y(=)_ z(=)_ normal
c             nx(=)_ ny(=)_ nz(=)_ (tolerance(=)_) (display)
c
      use main_data, only : user_lists, crdmap
c
      implicit integer(a-z)
c
      double precision
     &   tolerance, temp_t, zero, dot,
     &   x_toler, y_toler, z_toler, pt(3),
     &   flag, normvec(3), length, rel_tol,
     &   dx, dy, dz, d12, norm_tol
      logical :: do_display, matchs, check,
     &           numd, endcrd, matchs_exact, debug, display_coords
      dimension matchlist(max_list_size)
      data zero /0.0d00/, flag /-999999.0/, norm_tol / 1.0d-10 /
c
      if( debug ) write(out,*) "...in ulist_plane..."
c
      do_display   = .false.; tolerance = 0.000001
      pt(1:3) = flag; normvec(1:3) = flag
c
c                  -- must be input to define pt on the plane
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs_exact( "point" ) ) then
        call ulist_error( 19 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "x" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt(1) ) ) then
        call ulist_error( 9 )
        call scan_flushline
        return
      end if
c
      if( matchs_exact( "y" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt(2) ) ) then
        call ulist_error( 19 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "z" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt(3) ) ) then
        call ulist_error( 19 ); call scan_flushline; return
      end if
c
c                  -- must be input to define normal to plane
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs_exact( "normal" ) ) then
        call ulist_error( 20 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "nx" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( normvec(1) ) ) then
        call ulist_error( 21 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "ny" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( normvec(2) ) ) then
        call ulist_error( 21 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "nz" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( normvec(3) ) ) then
        call ulist_error( 21 ); call scan_flushline; return
      end if
c
c                  -- optional tolerance
c
      if( matchs( "tolerance", 5 ) ) then
        if( matchs_exact( "=" ) ) call splunj
        if( numd( temp_t ) ) then
          tolerance = temp_t
        else
          call ulist_error( 10 ); call scan_flushline
        end if
      end if
c
c                  -- optional display request
c
      if( matchs( "display", 6 ) ) do_display = .true.
      if( matchs( 'coordinates', 4 ) ) display_coords = .true.
c
c                  -- if not an end of input at this point,
c                     we have an un-recognized item
c
      if( .not. endcrd() ) then
       call ulist_error( 22 ); call scan_flushline; do_display = .false.
       return
      end if
c
      if( debug ) write(out,9000) pt, normvec, tolerance, do_display
c
c                  get a distance in user coordinates as the
c                  tolerance to use in checking if node lies on plane
c                  see definition in ulist_toler
c
      call ulist_toler( tolerance, x_toler, y_toler, z_toler,
     &                  nonode, crdmap, c, debug, out )
      rel_tol = sqrt( x_toler*x_toler + y_toler*y_toler +
     &                z_toler*z_toler )
c
c                  verify all that all coordinates for pt and
c                  components of normal vector are defined.
c                  normalize vector perpendicular to plane.
c
      check = pt(1) .eq. flag .or. pt(2) .eq. flag .or.
     &        pt(3) .eq. flag .or. normvec(1) .eq. flag .or.
     &        normvec(2) .eq. flag .or. normvec(3) .eq. flag
      if( check ) then
         call ulist_error( 23 ); do_display = .false.; return
      end if
c
      length = sqrt( normvec(1)*normvec(1) + normvec(2)*normvec(2) +
     &               normvec(3)*normvec(3) )
      if( length <= norm_tol ) then
         call ulist_error( 24 ); do_display = .false.; return
      end if
      normvec(1) = normvec(1) / length
      normvec(2) = normvec(2) / length
      normvec(3) = normvec(3) / length
      if( debug ) write(out,9010) normvec, rel_tol
      write(out,*) " "
      write(out,9220) normvec
      write(out,9222) rel_tol
c
c                  check all nodes in model.
c                    - get a vector from specified point on plane to
c                      the node
c                    - get dot product of that vector onto unit
c                      normal to plane. that is distance in user
c                      coordinates of node to the plane.
c                    - compare user point-plane distance to
c                      rel_tol computed above
c
      count = 0
      do node = 1, nonode
        position = crdmap(node)
        dx   = c(position+0) - pt(1) ! x, y, z for node
        dy   = c(position+1) - pt(2)! components of vector pt -> pt2
        dz   = c(position+2) - pt(3)
        dot  = abs( normvec(1)*dx + normvec(2)*dy + normvec(3)*dz )
        if( dot > rel_tol ) cycle
        count = count + 1
        matchlist(count) = node
      end do
c
c                 if we found nodes that matched the specification
c                 store vlaues as an integer list in user defined lists
c                 table.
c
      if( count > 0 ) then
         user_lists(list_col)%length_list = count
         allocate( user_lists(list_col)%list(count) )
         do i = 1, count
          user_lists(list_col)%list(i) = matchlist(i)
         end do
         write(out,9200) count
       else
         write(out,9210)
       end if
c
       return
c
 9000 format(2x,".... scanned values for plane option ....",
     & /,5x,"pt:     ",3f15.5,
     & /,5x,"nrmvec: ",3f15.5,
     & /,5x,'tolerance, do_display: ',f15.8, l5,
     & //)
 9010 format(2x,".... unit vector normal to plane ....",
     & /,5x,"nx, ny, nz: ",3f15.6,
     & /,5x,"rel_tol: ",f15.6,
     & //)
 9200 format(1x,">>>>> number of nodes placed in list: ", i8 )
 9220 format(1x,">>>>> unit normal vector to plane: ",3f10.6)
 9222 format(1x,">>>>> physical distance for proximity checks: ",
     &          f15.8  )
 9210 format(1x,">>>>> no nodes match the plane specification")
 9300 format( 5x, i7, 4f15.6)
 9310 format( 15x, 4f15.6)
c
      end

c     ****************************************************************
c     *                                                              *
c     *               subroutine ulist_sphere                        *
c     *                                                              *
c     *                   written by : dw and rhd                    *
c     *                                                              *
c     *    make a list of nodes that lie on a specified sphere       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ulist_sphere( list_col, debug, do_display,
     &                         display_coords, max_list_size )
      use global_data ! old common.main
c
c      list "list name" sphere origin x(=)_ y(=)_ z(=)_
c            radius =   tolerance  display
c
      use main_data, only : user_lists, crdmap

      implicit integer(a-z)
c
      double precision
     &   tolerance, temp_t, zero,
     &   pt(3), radius, flag, radius_tol, dx, dy, dz, d12
      logical :: do_display, matchs, check,
     &           numd, endcrd, matchs_exact, debug, display_coords
      dimension matchlist(max_list_size)
      data zero /0.0d00/, flag /-999999.0/
c
      if( debug ) write(out,*) "...in ulist_sphere..."
c
c                 loop to extract values. display
c                 request at end will be handled by driver.
c
      do_display = .false.; tolerance = 0.000001; pt(1:3) = flag
c
c                  -- must be input to define origin
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs_exact( "origin" ) ) then
        call ulist_error( 25 ); call scan_flushline; return
      end if
c
      if( matchs_exact( "x" ) ) call splunj
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( pt(1) ) ) then
        call ulist_error( 9 ); call scan_flushline; return
       end if
c
       if( matchs_exact( "y" ) ) call splunj
       if( matchs_exact( "=" ) ) call splunj
       if( .not. numd( pt(2) ) ) then
         call ulist_error( 9 ); call scan_flushline; return
       end if
c
       if( matchs_exact( "z" ) ) call splunj
       if( matchs_exact( "=" ) ) call splunj
       if( .not. numd( pt(3) ) ) then
         call ulist_error( 9 ); call scan_flushline; return
       end if
c
c                  -- must be radius of sphere
c
      if( matchs_exact( "," ) ) call splunj
      if( .not. matchs_exact( "radius" ) ) then
        call ulist_error( 16 ); call scan_flushline; return
      end if
      if( matchs_exact( "=" ) ) call splunj
      if( .not. numd( radius ) ) then
          call ulist_error( 14 ); call scan_flushline; return
      end if
c
c                  -- optional tolerance
c
      if( matchs_exact( "," ) ) call splunj
      if( matchs( "tolerance", 5 ) ) then
        if( matchs_exact( "=" ) ) call splunj
        if( numd( temp_t ) ) then
          tolerance = temp_t
        else
          call ulist_error( 10 ); call scan_flushline; return
        end if
      end if
c
c                  -- optional display request
c
      if( matchs_exact( "," ) ) call splunj
      if( matchs( "display", 6 ) )  do_display = .true.
      if( matchs( 'coordinates', 4 ) ) display_coords = .true.
c
c                  -- if not an end of input at this point,
c                     we have an un-recognized item
c
      if( .not. endcrd() ) then
       call ulist_error( 15 ); call scan_flushline; do_display = .false.
       return
      end if
c
      if( debug ) write(out,9000) pt, radius, tolerance, do_display
c
c                  get a relative tolerance to use in comparisons
c                  based on (x,y,z) dimensions of the model.
c
      radius_tol = tolerance * radius
      write(out,9220) radius_tol
c
      check = pt(1) .eq. flag .or. pt(2) .eq. flag .or.
     &        pt(3) .eq. flag .or. radius <= zero
      if( check ) then
         call ulist_error( 26 ); do_display = .false.; return
      end if
c
c                 examine all model nodes.
c                  - compute components of vector from origin
c                    to node and length
c                  - check length against radius tolerance
c
      count = 0
      do node = 1, nonode
        position = crdmap(node)
        dx  = c(position+0) - pt(1)
        dy  = c(position+1) - pt(2)
        dz  = c(position+2) - pt(3)
        d12 = sqrt( dx*dx + dy*dy + dz*dz )
        if( abs(d12-radius) > radius_tol ) cycle
        count = count + 1
        matchlist(count) = node
      end do
c
c                 if we found nodes that matched the specification
c                 store vlaues as an integer list in user defined lists
c                 table.
c
      if( count > 0 ) then
         user_lists(list_col)%length_list = count
         allocate( user_lists(list_col)%list(count) )
         do i = 1, count
          user_lists(list_col)%list(i) = matchlist(i)
         end do
         write(out,9200) count
       else
         write(out,9210)
       end if
c
       return
c
 9000 format(2x,".... scanned values for sphere option ....",
     & /,5x,"pt(x,y,z): ",3f15.6,
     & /,5x,"radius, tolerance: ",2f15.6,
     & /,5x,"do_display: ",l5,
     & //)
 9100 format(2x,".... tolerances for comparing coordinates ....",
     & /,5x,"nonode: ",i8,
     & /,5x,"x, y, z tolerance: ", 3e15.8, //)
 9200 format(1x,">>>>> number of nodes in list: ", i8 )
 9210 format(1x,">>>>> no nodes match the sphere specification")
 9220 format(1x,">>>>> physical distance for proximity checks: ",
     &        f15.8 )
c
      end
