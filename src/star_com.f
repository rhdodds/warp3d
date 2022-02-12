                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine star_com                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 02/10/22 rhd               *          
c     *                                                              *          
c     *     interprets the special star commands:                    *          
c     *     commands that are preceeded by an astrick and are        *          
c     *     trapped directly by scan.                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine star_com                                                       
      use global_data ! old common.main
      use scan_macros
      implicit none          
c        
      integer :: dum, nblank, reclen, endchr                                  
      double precision :: dumd                                                                    
      real :: t1, dumr                                                   
      real, external :: wcputime                                                         
      character(len=1) :: dums                                                  
      logical :: promsw, echosw, comsw, atrdsw, eolsw, eofsw, menusw,
     &           ptsw, signsw  
      logical, external :: matchs, label, string, endcrd, matchs_exact                                                             
      logical, parameter :: debug = .false.
c
      if( debug ) write(out,9000)
c                                                                               
c               command: *input -- getting input from file                       
c                                                                               
      if( matchs('input',5) ) then                                                
         call infile                                                            
         if( debug ) write (out,9010)    
         return
      end if
      if( matchs('inout',5) ) then  ! because Bob types wrong many times                                              
         call infile                                                            
         if( debug ) write (out,9010) 
         return
      end if
c
c               command: *output -- writing output to a file                     
c                                                                               
      if( matchs('output',6) ) then                                         
         call outfil    
         if( debug ) write (out,9010)
         return
      end if
c                                                                               
c               command: *time -- writes out total wall time    
c                         call WARP3D routine                  
c                                                                               
      if( matchs('time',4) ) then                                            
         t1 = wcputime(1)                                                       
         call errmsg (182,dum,dums,t1,dumd)
         if( debug ) write (out,9010)
         return
      end if
c                                                                               
c               command: *reset -- resets after a fatal error 
c                         call WARP3D routine                  
c                                                                               
      if( matchs('reset',5) ) then                                           
         input_ok = .true.  
         num_error = 0                                                    
         call errmsg(184,dum,dums,dumr,dumd)  
         if( debug ) write (out,9010)
         return
      end if
c                                                                               
c               command: *echo -- sets the scan echo on or off                   
c                                                                               
      if( matchs('echo',4) ) then                                            
         nblank= 20                                                             
         reclen= 80                                                             
         endchr= 1h$                                                            
         promsw= .false.                                                        
         comsw= .false.                                                         
         atrdsw= .false.                                                        
         eolsw= .true.                                                          
         eofsw= .true.                                                          
         menusw= .false.                                                        
         ptsw= .false.                                                          
         signsw= .false.                                                        
         call scinit(nblank,reclen,endchr,promsw,echosw,comsw,atrdsw,           
     &        eolsw,eofsw,menusw,ptsw,signsw)                                   
c                                                                               
         if( matchs('off',3) ) then                                              
            echosw= .false.                                                     
            call scinit(nblank,reclen,endchr,promsw,echosw,comsw,               
     &           atrdsw,eolsw,eofsw,menusw,ptsw,signsw)                         
         else                                                                   
            echosw= .true.                                                      
            call scinit(nblank,reclen,endchr,promsw,echosw,comsw,               
     &           atrdsw,eolsw,eofsw,menusw,ptsw,signsw)                         
         endif                                                                  
         if( debug ) write (out,*) '<<<<<< leaving star_com'      
         return
c                                                                               
      end if
c                                                                               
c               command: *parameter               
c
      if( matchs('parameter',5) ) then 
         if( .not. allocated(macros) ) allocate( macros(max_macros) )                                          
         call star_com_parameter
         if( debug ) write(out,9010)
         return
      end if
c                                                                               
c               unknown command: call WARP3D routine              
c                                                                               
      call errmsg (206,dum,dums,dumr,dumd)                                   
      if( debug ) write (out,9010) 
      return
c
 9000 format(1x,'<<<<<< entering star_com')   
 9010 format(1x,'<<<<<< leaving star_com')       
c
      contains
c     ========
c
      subroutine star_com_parameter
      implicit none
c
      integer :: id_nchars, param_nchars, n, n1, found_col, i
      character(len=80) :: s, s1, parameter_id, parameter_strng, 
     &                     bad_string
      logical, parameter :: ldebug = .false.
      logical :: found
c
c               commands:
c                 *parameter display
c                 *parameter [ <label> (=) <string> (,) ]
c 
      if( matchs("display",4) ) then
        if( num_macros == 0 ) then
           write(out,9200) 
           return
        end if
        write(out,9210) num_macros
        do i = 1, num_macros
          n = macros(i)%nchars_id 
          s = macros(i)%id   
          n1 = macros(i)%nchars_value
          s1 = macros(i)%value
          write(out,9220) i, s(1:n), s1(1:n1)
        end do  
        write(out,*) " "
        return
      end if
c     
      do while( .true. ) ! get all param definitions on line
c
        if( matchs_exact(",") ) call splunj()
        if( endcrd(dum) ) return
        if( .not. label(dum)  ) then
          call entits( bad_string, n )
          write(out,9100) bad_string(1:n)
          input_ok = .false.
          return
        end if
c
        call entits( parameter_id, id_nchars )  
        if( matchs("=",1) ) call splunj()
        if( .not. string(dum) ) then
          call entits( bad_string, n )
          write(out,9110) bad_string(1:n)
          input_ok = .false.
          return
        end if
        call entits( parameter_strng, param_nchars ) 
c
        if( num_macros == 0 ) then
          num_macros = 1
          macros(num_macros)%nchars_id    = id_nchars
          macros(num_macros)%nchars_value = param_nchars
          macros(num_macros)%id           = parameter_id
          macros(num_macros)%value        = parameter_strng
          if( ldebug ) then
           write(out,9000) num_macros, id_nchars, param_nchars
           write(out,9010) parameter_id(1:id_nchars)
           write(out,9020) parameter_strng(1:param_nchars)
          end if 
          cycle
        end if
c
        found = .false.
        do i = 1, num_macros
          n = macros(i)%nchars_id 
          s = macros(i)%id   
          if( n .ne. id_nchars ) cycle
          if( s(1:n) .ne. parameter_id(1:n) ) cycle
          found = .true.
          found_col = i
          exit 
        end do
        if( found ) then
          macros(found_col)%nchars_value = param_nchars
          macros(found_col)%value = parameter_strng
          if( ldebug ) then
            write(out,9000) num_macros, id_nchars, param_nchars
            write(out,9030) found_col
            write(out,9010) parameter_id(1:id_nchars)
            write(out,9020) parameter_strng(1:param_nchars)
          end if
          cycle
        end if 
c
        num_macros = num_macros + 1
        if( num_macros > max_macros ) then
          input_ok = .false.
          write(out,9120) max_macros
          call die_abort
        end if
c
        macros(num_macros)%nchars_id    = id_nchars
        macros(num_macros)%nchars_value = param_nchars
        macros(num_macros)%id           = parameter_id
        macros(num_macros)%value        = parameter_strng
        if( ldebug ) then
         write(out,9000) num_macros, id_nchars, param_nchars
         write(out,9010) parameter_id(1:id_nchars)
         write(out,9020) parameter_strng(1:param_nchars)
        end if
c
      end do 
c
      return
c
 9000 format(2x,".... num_macros, id_nchars, param_nchars: ",3i5)
 9010 format(2x,".... parameter id: ",a)
 9020 format(2x,".... value       : ",a)
 9030 format(2x,".... found_col: ",i3)
 9100 format(/1x,'>>>>> ERROR: expecting name of parameter.',
     &       /1x,'             scanning: ',a,/)
 9110 format(/1x,'>>>>> ERROR: expecting string for parameter value.',
     &       /1x,'             scanning: ',a,/)
 9120 format(/1x,'>>>>> ERROR: max # parameters exceeded. limit: ',i5,
     &       /1x,'             job terminated',///)
 9200 format(/1x,'>>>>> No parameters defined ....',/)
 9210 format(/1x,'>>>>> Currently defined parameters: ',i4)
 9220 format(10x,i4,2x,a,2x,a)
c
      end subroutine star_com_parameter
      end subroutine star_com
                                                                               
                                                                                
