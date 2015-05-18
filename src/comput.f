c     ****************************************************************
c     *                                                              *
c     *                      subroutine comput                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/24/2015 rhd              *
c     *                                                              *
c     *     this subroutine supervises the computation of the quan-  *
c     *     tities requested by the user.                            * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine comput
      use j_data, only: comput_j, comput_i
      implicit integer (a-z)
$add common.main
      real dumr
#dbl      double precision ::
#sgl      real ::
     &   dumd, zero
      character :: dums
      character*80 :: name
      character*8 :: ldname
      logical :: fatal, numnod, numel, coor, elprop, elinc, constr,
     &           matchs, endcrd, true, label, found, block,
     &           scanms, notes_msg
      dimension :: intlst(mxlsz)
      common/erflgs/ numnod,numel,fatal,coor,elprop,elinc,constr,block 
      data zero /0.0d00/, notes_msg / .true. /
c
c                       branch on the type of computation to be
c                       performed.
c         
      comput_j = .false.; comput_i = .false.
      if( matchs('domain',5) ) comput_j = .true.
      if( matchs('interaction',8) ) comput_i = .true.
      if( matchs('domain',5) ) comput_j = .true.
      if( comput_j .or. comput_i ) then
       call didriv; go to 9999
      end if
c
      if( .not. matchs('displacements',5) ) then
        call errmsg(121,dum,dums,dumr,dumd); go to 9999
      end if
c
c                       computation of displacements for a specified
c                       series of time steps.
c
c                       make sure that the constraints to be applied
c                       for the specified time steps exist.
c                       
      if( .not. constr ) then
        param = 6; call errmsg(11,param,dums,dumr,dumd); go to 9999
      end if
c
c                       input the loading to be used to define the
c                       specified time steps.
c                       
      if( matchs('for',3) ) call splunj
      if( .not.matchs('loading',4) ) then
        call errmsg(54,param,dums,dumr,dumd);  go to 9999
      end if
c
      if( .not. label(dummy) ) then
        call errmsg(121,param,dums,dumr,dumd);  go to 9999
      end if
c
c                       if the loading name input cannot be found,
c                       error skip command. otherwise
c                       search for the step loading in the loading
c                       library.
c
      found = .false.
      name = ' '    
      ldname = ' '
      call entits(name,nc)
      if( nc > 8) nc = 8
      ldname(1:nc) = name(1:nc)        
      lodn = lodhed/two16
 1222 continue
      if( lodn == 32460 ) go to 1223
       if( scanms(lodnam(lodn),ldname,8) ) then
         if( .not. scanms(lodtyp(lodn),'TIMESTEP',8) ) then
             param = lodn; call errmsg(122,param,dums,dumr,dumd)
             go to 9999
         end if
         ldnum = lodn; found = .true.; go to 1223
      end if
      lodn = lodlst(lodn)/two16
      go to 1222
c       
 1223 continue
      if( .not. found ) then
        call errmsg(60,dum,dums,dumr,dumd); go to 9999
      end if
c
c                       there is a valid loading. read the list of
c                       time/load steps to be computed.
c                     
      if( matchs('for',3) ) call splunj
      if( .not. matchs('steps',4) ) then
        call errmsg(121,param,dums,dumr,dumd);  go to 9999
      end if
c
      call scan
      call trlist( intlst,mxlsz,nonode,lenlst,errnum )
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       a value of 4 indicates that no list was found.
c                       in these last three cases, the illegal list
c                       will be ignored and a new compute command will
c                       be sought. 
c
      if( errnum .ne. 1 ) then
          if( errnum  == 2) then
               param = 1
               call errmsg(24,param,dums,dumr,dumd)
          else if( errnum == 3 ) then
               param = 2
               call errmsg(24,param,dums,dumr,dumd)
          else if( errnum == 4) then
               param = 4
               call errmsg(24,param,dums,dumr,dumd)
          end if
          go to 9999
      end if
c
c                       the list of time steps is a valid one. compute
c                       the displacements for the steps specified.
c  
      notes_msg = .false.
      call stpdrv( intlst, lenlst, ldnum)
      return
c
 9999 continue
      call scan_flushline; return
      end
