c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 02/17/2017 rhd             *
c     *                                                              *
c     *     this subroutine checks various program variables and     *
c     *     arrays for errors after the input of data pertaining to  *
c     *     the variables and arrays in question.                    *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errchk( lsn, chkprm, debug1 )
      use global_data ! old common.main
      implicit integer (a-z)
      logical debug1
c
      if (lsn .eq. 32) then
            call chk_crystal(chkprm)
            go to 9999
      end if
c                       branch on subroutine number
c
      if (debug1) write (out,*) '   branching on lsn  (errchk) ', lsn
      go to (100,200,300,400,500,600,700,800,900,1000,1100,1200,
     &       1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,
     &       2300,2400,2500) lsn
c
c         ******** Structure name was input ********
c
 100  continue
      go to 9999
c
c         ******** Material Descriptions ********
c
 200  continue
      call errchk_2( chkprm )
      go to 9999
c
c         ******** Number of nodes, elements was input ********
c
 300  continue
      go to 9999
c
c         ******** Coordinates were input ********
c
 400  continue
      call errchk_4( debug1 )
      goto 9999
c
c         ******** Element Description was input ********
c
 500  continue
      call errchk_5( debug1 )
      goto 9999
c
c         ******** Incidences were input ********
c
 600  continue
      call errchk_6
      go to 9999
c
c         ******** Constraints were given ********
c
 700  continue
      call errchk_7
      go to 9999
c
c         ******** Loading was specified********
c
 800  continue
      call errchk_8( chkprm, debug1 )
      goto 9999
c
c         ******** Display command ********
c
 900  continue
      go to 9999
c
c         ******** Solution analysis parameters were input ********
c
 1000 continue
      call errchk_10
      go to 9999
c
c         ******** Inital conditions command ********
c
 1100 continue
      go to 9999
c
c         ******** Compute command ********
c
 1200 continue
      go to 9999
c
c         ******** Output command ********
c
 1300 continue
      go to 9999
c
c         ******** dummy ********
c
 1400 continue
      go to 9999
c
c         ******** Save database was requested ********
c
 1500 continue
      go to 9999
c
c         ******** Retrieve database was requested ********
c
 1600 continue
      go to 9999
c
c         ******** dummy********
c
 1700 continue
      go to 9999
c
c         ******** Blocking was input ********
c
 1800 continue
      call errchk_18( debug1 )
      go to 9999
c
c         ******** Exit Warp ********
c
 1900 continue
      go to 9999
c
c         ******** Domain Integral was input ********
c
 2000 continue
      go to 9999
c
c         ******** Crack growth parameters were input ********
c

 2100 continue
      call errchk_21
      goto 9999
c
c         ******** Debug was set (on or off) ********
c
 2200 continue
      go to 9999
c
c         ******** stress-strain curve has been input ********
c
 2300 continue
      call errchk_23
      goto 9999
c
c         ******** license agreement has been requested ********
c
 2400 continue
      go to 9999
c
c         ******** contact planes have been defined ********
c
 2500 continue
      call errchk_25
      goto 9999
c
c
c
 9999 continue
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_2                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 5/18/2017                  *
c     *                                                              *
c     *    this subroutine checks the values input for the           *
c     *    material definitions.  If any are incorrect, WARP3D       *
c     *    will stop at the first compute command.                   *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errchk_2( matnum )
      use global_data ! old common.main
      use main_data, only : matprp, lmtprp
      use erflgs
      implicit none
c
c
      integer :: dum, matnum
      real :: dumr
      real, parameter :: fgm_mark = -99.0
      character :: dums
      double precision :: dumd
      logical :: is_matl_cohesive, is_umat, is_matl_cp
      character(len=50) erprms
c
      erprms = ' '
c
c          this is an older but still used routine. It does a quick check on
c          come key property values for material models for solid elements.
c          it does not check correctly for cohesive materials or CP.
c          just return.
c
      is_matl_cohesive = ( matprp(9,matnum) .eq. 4.0 )
      is_matl_cp       = ( matprp(9,matnum) .eq. 10.0 )
      if( is_matl_cp ) call chk_cp(matnum)
      if( is_matl_cohesive .or. is_matl_cp ) return
c
c           check E
c
      if( matprp(1,matnum) .ne. fgm_mark ) then
        if( matprp(1,matnum) .le. 0.0 ) then
           erprms = 'greater than zero.'
           call errmsg4( matnum, 'e', erprms,
     &                   num_fatal, out, matnam, input_ok )
        end if
      end if
c
c           check nu
c
      if( matprp(2,matnum) .ne. fgm_mark ) then
        if( matprp(2,matnum) .lt. -1.0 .or.
     &       matprp(2,matnum) .gt. 0.5 ) then
             erprms = 'between -1 and .5.'
             call errmsg4( matnum, 'nu', erprms,
     &                     num_fatal, out, matnam, input_ok )
        end if
      end if
c
c           check tan_e and n_power bases on which material model we are
c           using. if model is type:
c
c              1 -- only have linear hardening.  check tan_e.
c              2 -- only have power law hardening. check n_power
c              3 -- gurson/mises.  either one -- check both. If
c                    stress strain curve is given, skip the check.
c
c
      if( matprp(9,matnum) .eq. 1 ) then
c
         if( .not.lmtprp(8,matnum) .and.
     &        matprp(4,matnum) .lt. 0.0 ) then
            if(  matprp(4,matnum) .ne. fgm_mark ) then
               erprms = 'greater or equal to zero.'
               call errmsg4( matnum, 'tan_e', erprms,
     &                   num_fatal, out, matnam, input_ok )
            end if
         end if
c
      else if( matprp(9,matnum) .eq. 2 ) then
c
         if( matprp(11,matnum) .le. 0.0 ) then
            if(  matprp(11,matnum) .ne. fgm_mark ) then
              erprms = 'greater than zero.'
           call errmsg4( matnum, 'n_power', erprms,
     &                   num_fatal, out, matnam, input_ok )
            end if
         end if
c
      else if( matprp(9,matnum) .eq. 3 .and.
     &          .not. lmtprp(24,matnum) ) then
c
         if( matprp(11,matnum) .le. 0.0 .and.
     &        matprp(4,matnum) .lt. 0.0 ) then
            if( matprp(4,matnum) .ne. fgm_mark )
     &        call errmsg ( 258, dum, dums, dumr, dumd)
         end if
c
      endif
c
c           check beta
c
      if( matprp(3,matnum).lt.0.0 .or. matprp(3,matnum).gt.1.0) then
           erprms = 'between 0.0 and 1.0.'
           call errmsg4( matnum, 'beta', erprms,
     &                   num_fatal, out, matnam, input_ok )
      endif
c
c           check rho
c
      if( matprp(7,matnum) .ne. fgm_mark ) then
        if( matprp(7,matnum) .lt. 0.0 ) then
           erprms = 'greater or equal to zero.'
           call errmsg4( matnum, 'rho', erprms,
     &                   num_fatal, out, matnam, input_ok )
        end if
      end if
c
c           check m_power
c
      if( matprp(10,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'm_power', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check ref_eps
c
      if( matprp(12,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'ref_eps', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check f_0
c
      if( matprp(14,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'f_0', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check q1
c
      if( matprp(15,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'q1', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check q2
c
      if( matprp(16,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'q2', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check q3
c
      if( matprp(17,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'q3', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check s_n
c
      if( matprp(19,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 's_n', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check e_n
c
      if( matprp(20,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'e_n', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check f_n
c
      if( matprp(21,matnum) .lt. 0.0 ) then
         erprms = 'greater or equal to zero.'
         call errmsg4( matnum, 'f_n', erprms,
     &                 num_fatal, out, matnam, input_ok )
      end if
c
c           check for conflicting sets of thermal expansion
c           properties. we do not allow both isotropic and
c           anisotropic thermal properties to be specified.
c
      if( lmtprp(25,matnum) ) then
        if( matprp(6,matnum) .ne. 0.0 )
     &       call errmsg ( 274, matnum, 'f_n', dumr, dumd )
      end if
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_4                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 5/18/2017 rhd              *
c     *                                                              *
c     *    this subroutine checks the coordinates for constistency.  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errchk_4( debug1 )
      use global_data ! old common.main
c
      use main_data, only : crdmap
      use erflgs
c
      implicit none
      integer :: i
      real :: dumr
      character :: dums*1
      double precision :: dumd
      logical debug1
c
c                       check to make sure that all structural
c                       nodes have been given coordinates.
c
      if( debug1 ) write (out,*) '  checking all nodes (errchk) '
      coor = .true.
      do i = 1, nonode
         if( crdmap(i) == 0 ) then
            coor = .false.
            call errmsg( 15, i, dums, dumr, dumd )
         end if
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_5                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 5/18/2017 rhd              *
c     *                                                              *
c     *     this subroutine checks the element description for       *
c     *     consistency.                                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errchk_5( debug1 )
      use global_data ! old common.main
      use main_data, only : elstor
      use erflgs
c
      implicit none
      integer :: i
      real :: dumr
      character :: dums*1
      double precision :: dumd
      logical :: debug1
c
c
c                       check to make sure that all elements have
c                       been given the necessary properties and its
c                       stress/strain vectors have been initialized.
c
      if( debug1 ) write (out,*) '  checking all elements (errchk) '
      elprop = .true.
      do i = 1, noelem
         if( elstor(1,i) == 0 ) then
            elprop = .false.
            call errmsg(34,i,dums,dumr,dumd)
         end if
      end do
c
c                       data for all elements has been temporarily
c                       stored. store this information permanently.
c
      if( debug1 ) write(out,*)'    calling prcsel from errchk ',elprop
c
      if( elprop ) call prcsel
c
      if( debug1 ) write(out,*)'     returned from prcsel '
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_6                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 5/18/2017 rhd              *
c     *                                                              *
c     *     this subroutine checks the incidences for consistency    *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errchk_6
      use global_data ! old common.main
      use main_data, only : incmap, incid
      use erflgs
c
      implicit none
c
      integer :: dup_cnt, i, j, kk, incptr, nnode, curinc, num, rem,
     &           err, dumi
      integer, allocatable, dimension (:) :: dup_ele
      integer :: inc(mxndel)
      real :: dumr, mlt
      character :: dums*1
      double precision :: dumd
      logical :: dupnod_msg
c
c                       check the element incidences and the mapping
c                       vector for errors for every element.
c
      dup_cnt = 0
      elinc = .true.
      do i = 1, noelem
c
c                       check for holes in the mapping vector.
c
         if ( incmap(i) .eq. 0 ) then
            elinc = .false.
            call errmsg( 40, i, dums, dumr, dumd )
         else
            incptr = incmap(i)
            nnode = iprops(2,i)
            do j = 1, nnode
c
c                       check for holes in the incidences
c
               if ( incid(incptr+j-1).eq.0 ) then
                  elinc = .false.
                  call errmsg( 41, i, dums, dumr, dumd )
               end if
               inc(j) = incid(incptr+j-1)
            end do
            dupnod_msg = .true.
            do j = 1, nnode
               curinc = inc(j)
               do kk = j+1, nnode
c
c                       check for repeated incidences. warning only.
c
                  if ( curinc .eq. inc(kk) .and. dupnod_msg ) then

                      if (dup_cnt .eq. 0) then
                         allocate( dup_ele(noelem), stat=err)
                         if (err .ne. 0) then
                            call errmsg(42,dumi,dums,dumr,dumd)
                            call die_abort
                         end if
                      end if
                      dup_cnt = dup_cnt + 1
                      dup_ele(dup_cnt) = i

c                     call errmsg( 42, i, dums, dumr, dumd )

                     dupnod_msg = .false.
                  end if
               end do
            end do
         end if
      end do
c
      if (dup_cnt .gt. 0) then
         write(out,9000)
         mlt = real(dup_cnt)/8.0
         num = int(mlt)*8
         rem = dup_cnt - num
         do i = 1, num, 8
            write(out,9001) dup_ele( i ),dup_ele(i+1),dup_ele(i+2),
     &                      dup_ele(i+3),dup_ele(i+4),dup_ele(i+5),
     &                      dup_ele(i+6),dup_ele(i+7)
         end do
         if (rem .eq. 1)  write(out,9002) dup_ele(num+1)
         if (rem .eq. 2)  write(out,9003) dup_ele(num+1),dup_ele(num+2)
         if (rem .eq. 3)  write(out,9004) dup_ele(num+1),dup_ele(num+2),
     &                                    dup_ele(num+3)
         if (rem .eq. 4)  write(out,9005) dup_ele(num+1),dup_ele(num+2),
     &                                    dup_ele(num+3),dup_ele(num+4)
         if (rem .eq. 5)  write(out,9006) dup_ele(num+1),dup_ele(num+2),
     &                                    dup_ele(num+3),dup_ele(num+4),
     &                                    dup_ele(num+5)
         if (rem .eq. 6)  write(out,9007) dup_ele(num+1),dup_ele(num+2),
     &                                    dup_ele(num+3),dup_ele(num+4),
     &                                    dup_ele(num+5),dup_ele(num+6)
         if (rem .eq. 7)  write(out,9008) dup_ele(num+1),dup_ele(num+2),
     &                                    dup_ele(num+3),dup_ele(num+4),
     &                                    dup_ele(num+5),dup_ele(num+6),
     &                                    dup_ele(num+7)
         deallocate(dup_ele)
      end if
c
c                       call subroutine to set necessary mapping vectors and
c                       arrays and various other necessary data in preparation
c                       for the solution of the current structure.
c
      if ( elinc ) call setup( fatal )
c
 9000 format(/1x,'>> warning: there are duplicate entries in',
     &           ' the incidences for',/,
     &           '            the following elements:'/)
 9001 format(8(i9))
 9002 format(i9)
 9003 format(2(i9))
 9004 format(3(i9))
 9005 format(4(i9))
 9006 format(5(i9))
 9007 format(6(i9))
 9008 format(7(i9))
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_7                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 10/09/95                   *
c     *                                                              *
c     *   this subroutine checks constraints input for consistency.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine errchk_7
      use global_data ! old common.main
c
      use main_data, only : trn, trnmat, cnstrn_in, inverse_incidences
      use erflgs
c
      implicit none
c
      integer :: node, ndof, col, row, dum, dof
      real :: dumr
      character :: dums*1
      double precision :: dumd
      double precision, parameter :: zero=0.d0, d32460=32460.0d0
c
c                       check transformation matrix for each pertinent node.
c
      constr = .true.
      newtrn = .true.
      zrocon = .true.
c
      do node = 1, nonode
         if ( trn(node) ) then
c
            ndof   = iprops(4,inverse_incidences(node)%element_list(1))
            do row = 1, ndof
               do col = 1, ndof
c
c                       check the transformation matrix for holes
c
                  if ( trnmat(node)%mat(row,col).eq.d32460 ) then
                     constr = .false.
                     newtrn = .false.
                     call errmsg( 81, node, dums, dumr, dumd )
                  end if
               end do
            end do
         end if
      end do
c
c                       check to make sure that at least one con-
c                       straint has been input and that the con-
c                       straint values for each pertinent dof are
c                       feasible. set the flag indicating all zero
c                       constraints.
c
      if ( csthed .eq. -1 ) then
         call errmsg( 82, dum, dums, dumr, dumd )
         constr = .false.
      else
         dof = csthed
      end if
c
 720  if ( dof .eq. -1 ) go to 725
c
      if ( cnstrn_in(dof) .eq. d32460 ) then
         call errmsg( 83, dof, dums, dumr, dumd )
         constr = .false.
      else if ( cnstrn_in(dof) .ne. zero ) then
         zrocon = .false.
      end if
      dof = cstmap(dof)
c
      go to 720
c
 725  continue
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_8                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 2/8/2018 rhd               *
c     *                                                              *
c     *  this subroutine checks the loading input for constistency.  *
c     *                                                              *
c     ****************************************************************
c
      subroutine errchk_8( chkprm, debug1 )
      use global_data ! old common.main
      use main_data, only : stpchk
      implicit none
c
      integer :: i, chkprm, lodnum, step
      logical, external :: scanms
      logical :: debug1
c
c                       branch on time step versus regular loadings.
c                       chkprm is the loading number.
c
      if( scanms( lodtyp(chkprm), 'TIMESTEP', 8 ) ) then
c
c                       timestep definition:
c                         check to make sure all of the time steps
c                         existing for the loading number chkprm
c                         are defined.
c
         do i = 1, histep
            if( stpchk(i) ) cycle
            num_error = num_error + 1
            lodnum = chkprm
            step = i
            write(out,9084) step, lodnam(lodnum)
         end do
c
c                           set lowstp and histep, the range of time
c                           steps defined for this loading.
c
         stprng(chkprm,1) = lowstp
         stprng(chkprm,2) = histep
c
c                       regular loading definition:
c                        if face and body loadings have been defined,
c                        store the temporary data in the permanent data
c                        structure. Then deallocate temp data structures.
c
      else
         call allocate_perm_load( 1, chkprm )
         call allocate_temp_load( 2 )
      endif
c
      return

9084  format(/1x,'>>>>> error: step number ',i7,' of loading ',a8,
     &           ' has not been'/14x,'defined. an attempt to use this',
     &           ' time step will not be'/14x,'allowed.'/)
c

      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_10                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/12/2015                  *
c     *                                                              *
c     *  check that a convergence test is defined                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine errchk_10
      use global_data ! old common.main
      implicit integer (a-z)
      logical :: found
      real dumr
      character :: dums
      double precision
     &   dumd
c
c                       at least one convergence test must be defined.
c                       otherwise stop.
c
      found = .false.
      do i = 1, mxcvtests  ! parameter varariable
        if( convrg(i) ) found = .true.
      end do
      if ( .not. found ) call errmsg3( out, 13 )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_18                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/5/2017 rhd               *
c     *                                                              *
c     *    this subroutine checks the blocking input                 *
c     *    for constistency.                                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine errchk_18 (debug1)
      use global_data ! old common.main
      use main_data, only : elstor, incmap, incid
      use erflgs
c
      implicit none
c
      integer :: spnsum, i, k, num_bad_flags_used, blk_no, elem_loop,
     &           element,  matmodel, eletype, intord, intnum, matnum,
     &           matmodel2, eletype2, intord2, intnum2, curset,
     &           curset2, dum, node_loop, num_nodes, node, matnum2
      integer, allocatable :: node_refs(:)
      real :: dumr, beta_cohesive
      character :: dums*1
      double precision :: dumd
      logical :: cohesive
      logical geonon, geonon2, bbar, bbar2, debug1, segcur, segcur2,
     &        bad_block, bad_flags(20), cyclic_plasticity
c
c                       check the element blocking data structure.
c
      block  = .true.
      spnsum = 0
      allocate( node_refs(nonode) )
c
      do i = 1, nelblk
         if( (elblks(0,i).eq.0).or.(elblks(1,i).eq.0) ) then
            block = .false.
            call errmsg( 161, dum, dums, dumr, dumd )
            go to 9999
         else
            spnsum = spnsum + elblks(0,i)
         end if
      end do
c
      if ( spnsum .ne. noelem ) then
         block = .false.
         call errmsg( 173, dum, dums, dumr, dumd )
         go to 9999
      end if
c
c                      now check to see if blocking is scalar or vector.
c                      Also check if blocking is correct -- all
c                      elements in block are same material model type,
c                      element type, same integration order, same status of
c                      geometric nonlinearity and bbar. if material
c                      stress-strain curve is segmental, the curve set
c                      for all elements in the block must be identical.
c
      scalar_blocking = .false.
      num_bad_flags_used = 10
c
      write (out,9000)
      do blk_no = 1, nelblk
         if ( debug1 ) write (out,*) '> block:', blk_no
         node_refs(1:nonode) = 0
c
c                                  check properties
c
         do elem_loop = 1, elblks(0,blk_no)
            element = elblks(1,blk_no) + elem_loop - 1
            if ( debug1 ) write (out,*) '>> element:',element
            if ( elem_loop .eq. 1 ) then
               matmodel  = iprops(25,element)
               eletype   = iprops(1,element)
               intord    = iprops(5,element)
               intnum    = iprops(6,element)
               geonon    = lprops(18,element)
               bbar      = lprops(19,element)
               segcur    = iand( iprops(24,element), 4 ) .ne. 0
               curset    = -1
               if ( segcur) curset = iprops(21,element)
               cohesive  = lprops(10,element)
               beta_cohesive = props(23,element)
               matnum = iprops(38,element)
               cyclic_plasticity = matmodel .eq. 5
            else
               matmodel2 = iprops(25,element)
               eletype2  = iprops(1,element)
               intord2   = iprops(5,element)
               intnum2   = iprops(6,element)
               geonon2   = lprops(18,element)
               bbar2     = lprops(19,element)
               segcur2   = iand( iprops(24,element), 4 ) .ne. 0
               curset2   = -2
               if ( segcur2 ) curset2  = iprops(21,element)
               matnum2 = iprops(38,element)
               bad_flags(1) = matmodel2 .ne. matmodel
               bad_flags(2) = eletype2 .ne. eletype
               bad_flags(3) = intord2 .ne. intord
               bad_flags(4) = intnum2 .ne. intnum
               bad_flags(5) = .not. ( geonon2 .eqv. geonon )
               bad_flags(6) = .not. ( bbar2 .eqv. bbar )
               bad_flags(7) = .not. ( segcur .eqv. segcur2 )
               bad_flags(8) = .false.
               if ( segcur .and. ( segcur .eqv. segcur2) ) then
                  if ( curset .ne. curset2 ) bad_flags(8) = .true.
               end if
               bad_flags(9) = .false.
               if( cohesive ) then
                 if( props(23,element) .ne. beta_cohesive )
     &              bad_flags(9) = .true.
               end if
               bad_flags(10) = .false.
               if( cyclic_plasticity ) then
                  if( matnum .ne. matnum2 ) bad_flags(10) = .true.
               end if
               bad_block = .false.
               do k = 1, num_bad_flags_used
                 if ( bad_flags(k) ) bad_block = .true.
               end do
               if ( bad_block ) then
                 call errmsg( 216, blk_no, dums, dumr, dumd )
                 do k = 1, num_bad_flags_used
                   if ( bad_flags(k) ) then
                      call errmsg2( 19, k, dums, dumr, dumd )
                   end if
                 end do
                 call errmsg2( 20, 1, dums, dumr, dumd )
                 call die_gracefully
                 stop
               end if
            end if
         end do
c
c                                  check nodes for scalar or
c                                  vectorized blocking.
c
         do elem_loop = 1, elblks(0,blk_no)
            element   = elblks(1,blk_no) + elem_loop - 1
            num_nodes = iprops(2,element)
            if ( incmap(element) .eq. 0 ) then
               write(out,9900)
               call die_abort
            end if
c
            do node_loop = 1, num_nodes
               node = incid(incmap(element)+node_loop-1)
               if ( node_refs(node) .eq. 0 ) then
                  node_refs(node) = 1
               else if ( node_refs(node) .eq. 1 ) then
                  scalar_blocking = .true.
                  write(out,9100)
                  goto 9999
               endif
            end do
c
         end do
c
      end do
c
      write (out,9200)
c
 9999 continue
      return
c
 9000 format('>> Check element blocking ...')
 9100 format('     ... Elements have scalar blocking',/)
 9200 format('     ... Elements have vectorized blocking. No elements',
     &     /,'         in a block share a common node',/)
 9900 format('>>> ERROR: blocking input must follow element',
     &    /, '           incidences to enable consistency',
     &    /, '           checking. analysis terminated...',//)
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_21                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 5/28/02 rhd                *
c     *                                                              *
c     *    this subroutine checks the parameters for crack growth    *
c     *    for consistency.                                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine errchk_21
      use global_data ! old common.main
c
      use damage_data
      implicit integer (a-z)
c
      real dumr
      character :: dums
      double precision
     &   dumd, zero, tenth
      data zero, tenth / 0.0, 0.1 /
c
c            If the crack growth type is element_extinction, then:
c             - if load step size control has been activated and the
c               accompanying data structures are not yet allocated,
c               allocate them.
c             - if the release type is 2 (traction_separation), then
c               make sure crack plane normal, the gurson cell size, and
c               the release fraction were all set properly.  If not, the
c               error routines mark the input bad to prevent any compute
c               actions.
c
      if( growth_by_kill ) then
c
         if( load_size_control_crk_grth ) then
c
          if ( .not. g_stp_cntrl_allocated )
     &             call allocate_damage( 9 )
            if ( crack_growth_type .eq. 1 .and.
     &             max_porosity_change .eq. zero) then
               max_porosity_change = tenth *  porosity_limit
            else
               max_porosity_change = max_porosity_change *
     &              porosity_limit
            end if
c
         end if
c
         if( release_type .eq. 2 ) then
c
            if( crk_pln_normal_idx .le. 0 ) then
               call errmsg(241,dum,dums,dumr,dumd)
            end if
c
            if( gurson_cell_size .le. zero ) then
               call errmsg(242,dum,dums,dumr,dumd)
            end if
c
            if( release_fraction .le. zero ) then
               call errmsg(243,1,dums,dumr,dumd)
            end if
c
         end if
c
c            If the crack growth type is node_release, then check to make
c            sure that the critical angle and the crack plane normal have
c            been input before initializing the node_release data
c            structures.  Also, if the release type is traction separation,
c            then check the characteristic length and the release fraction
c            before calculating the release height.  If any of these checks
c            fail, the error routines mark the input bad to prevent compute
c            actions. if the user is requesting an enforced node
c            release, we can skip all these cheks - they were
c            done previously.
c
      else if( growth_by_release ) then
c
         if ( enforce_node_release ) return
c
         if( crk_pln_normal_idx .le. 0 ) then
            call errmsg( 241, dum, dums, dumr, dumd )
         end if
c
         if( critical_angle .le. zero ) then
            call errmsg( 246, 1, dums, dumr, dumd )
         end if
c
         if( init_crit_ang .le. zero ) then
            call errmsg( 246, 2, dums, dumr, dumd )
         end if
c
         if( const_front ) then
c
            if( num_nodes_thick .le. 0 ) then
               call errmsg( 284, dum, dums, dumr, dumd )
               call die_gracefully
               stop
            else if( num_crack_fronts .le. 0 ) then
               call errmsg( 281, dum, dums, dumr, dumd )
               call die_gracefully
               stop
            else if( init_ctoa_dist .le. zero ) then
               call errmsg( 288, dum, dums, dumr, dumd )
               call die_gracefully
               stop
            else if( ctoa_dist .le. zero ) then
               call errmsg( 288, dum, dums, dumr, dumd )
               call die_gracefully
               stop
            else if( char_length .le. zero ) then
               call errmsg( 243, 2, dums, dumr, dumd )
               call die_gracefully
               stop
            else
c
               if( master_lines_set ) then
                  call errmsg( 292, dum, dums, dumr, dumd )
               else
                  call allocate_damage( 11 )
               end if
c
            end if
c
c
 10         continue
c
         end if
c
         if( input_ok ) call dam_init_release
c
         if(const_front) call init_ctoa_back
c
         if( release_type .eq. 2 ) then
c
          if( release_fraction .le. zero ) then
               call errmsg( 243, 1, dums, dumr, dumd )
            else if( char_length .le. zero ) then
               call errmsg( 243, 2, dums, dumr, dumd )
            else
               call find_release_height
            end if
c
         end if
c
         if( overshoot_control_crk_grth .and.
     &        .not. overshoot_allocated ) call allocate_damage( 8 )
c
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_23                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 10/09/95                   *
c     *                                                              *
c     *   this subroutine checks the segmental stress-strain curves  *
c     *   for constistency.                                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine errchk_23
      use global_data ! old common.main
      use segmental_curves
      implicit integer (a-z)
      real dumr
      character :: dums
      double precision
     &   dumd
c
c
c
      if ( num_points .eq. 0 ) then
         call errmsg( 221, dum, dums, dumr, dumd )
         return
      endif
      num_seg_points(num_curve) =  num_points
      seg_curve_def(num_curve)  = .true.
c
c            store the current maximum number of points in a curve and the
c            maximum number of curves. this information is used in store/reopen
c            routines to work with only the necessary amount of data
c
      max_current_pts    = max(num_points,max_current_pts)
      max_current_curves = max(max_current_curves,num_curve)
      num_points         = 0
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine errchk_25                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 5/10/04                    *
c     *                                                              *
c     *   this subroutine checks the contact plane definitions       *
c     *   for constistency.                                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine errchk_25
      use contact, only : use_contact, maxcontact, contact_shape,
     &                    num_contact
      implicit integer (a-z)
      include 'param_def'
      real dumr
      character :: dums
      double precision
     &   dumd, zero
      data zero /0.0/
c
      use_contact = .false.
c
c              Here we check each plane
c
      do i = 1, maxcontact
         if( contact_shape(i) .ne. 0 ) then
            use_contact = .true.
            write (*,*) '   -> contact surface ',i,' is defined.'
            num_contact = i
         end if
      end do
c
      if ( .not. use_contact ) then
         write (*,*) '>>> No contact planes defined.'
      end if
c
c             IF we are using MPI:
c                send all processors copies of the contact information
c             If we are using the serial version:
c                return
c
      call wmpi_send_contact (.false.)
c
 9999 continue
      return
      end


c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine chk_cp                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/21/12                    *
c     *                                                              *
c     *   checks crystal plasticity material definitions for errors  *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_cp(matnum)
      use global_data ! old common.main
      use main_data, only: matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only: c_array
      implicit integer (a-z)
      integer, intent(in) :: matnum
      integer, dimension(2) :: chksz
      logical :: valid, exists
      character :: matname*24
c
      matname = matnam(matnum)
c           I'm just going to handle errors in this subroutine, the
c           error subroutine is getting way too long.
c
c           List of things to check:
c                 n_crystals >= 1                     c
c                 tolerance > 0                       c
c                 rho >= 0                            c
c                 alpha >= 0                          c
c                 if list, list is valid
c                 if offset, list and offset_el valid
c
      if (imatprp(101,matnum) .lt. 1) then
            write (out,1001) matname
 1001       format(/,1x,'>>>> Fatal error in material ', a24,
     &             ': the number of crystals',
     &             ' must be greater than zero!'/)
            call die_gracefully
      end if
      if (imatprp(6,matnum) .lt. 0) then
            write (out,1002) matname
 1002       format(/,1x,'>>>> Fatal error in material ', a24,
     &             ': thermal expansion coef.',
     &             ' must be greater or equal to zero!'/)
            call die_gracefully
      end if
      if (imatprp(7,matnum) .lt. 0) then
            write (out,1003) matname
 1003       format(/,1x,'>>>> Fatal error in material ', a24,
     &             ': mass density',
     &             ' must be greater or equal to zero!'/)
            call die_gracefully
      end if
      if (dmatprp(100,matnum) .le. 0.0) then
            write (out,1004) matname
 1004       format(/,1x,'>>>> Fatal error in material ', a24,
     &             ': tolerance',
     &             ' must be greater than zero!'/)
            call die_gracefully
      end if
c
      if (imatprp(104,matnum) .eq. 1) then
c           Check valid crystal
            if ((imatprp(105,matnum) .lt. 1) .or.
     &          (imatprp(105,matnum) .gt. max_crystals) .or.
     &          (.not. c_array(imatprp(105,matnum))%valid)) then
                  write (out,1007) matname, imatprp(105,matnum)
            end if
      elseif (imatprp(104,matnum) .eq. 2) then
c           Defer crystal check until later, but check filename
            inquire(FILE=smatprp(112,matnum), EXIST=exists)
            if ( .not. exists) then
                  write (out,1005) matname, smatprp(112,matnum)
                  call die_gracefully
            end if
      end if
c
      if (imatprp(107,matnum) .eq. 1) then
c           nearly anything is valid here
      elseif (imatprp(107,matnum) .eq. 2) then
c           Check for file
            inquire(FILE=smatprp(112,matnum), EXIST=exists)
            if (.not. exists) then
                  write (out,1005) matname,smatprp(112,matnum)
                  call die_gracefully
            end if
      end if


 1005 format(/1x,'>>>> Fatal error in material ', a24,
     &            ': file ', a24, ' not found.'/)
 1006 format(/1x,'>>>> Fatal error in material ', a24,
     &            ': offset element ', i7,
     &            ' is not valid'/)
 1007 format(/1x,'>>>> Fatal error in material ', a24,
     &            ': invalid crystal #',i3,/)
      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine chk_list                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/21/12                    *
c     *                                                              *
c     *   helper method check is list is defined, and compares it's  *
c     *     length to a provided value (if value is zero, don't      *
c     *     compare                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_list(list,length,valid)
      use main_data, only : tables
      implicit none
      character, intent(in) :: list*24
      integer, intent(in), dimension(2) :: length
      logical, intent(out) :: valid
      integer :: i
      logical :: found,sz
      write (*,*) list
c           First check existence
      found = .false.
      do i=1,size(tables)
            if (tables(i)%table_name .eq. list) then
                  found = .true.
                  exit
            end if
      end do
c           Optionally check size
      if (length(1) .eq. 0) then
            sz = .true.
      else
            if ((tables(i)%num_rows .eq. length(1)) .and.
     &            (tables(i)%num_cols .eq. length(2)) ) then
                  sz = .true.
            else
                  sz = .false.
            end if
      end if

      valid = found .and. sz

      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine chk_crystal                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/12/2017 rhd              *
c     *                                                              *
c     *   checks a crystal definition for errors                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_crystal(cnum)
      use global_data ! old common.main
      use crystal_data, only: c_array
      implicit integer (a-z)
      integer, intent(in) :: cnum

      if (c_array(cnum)%elastic_type .ne. 3) then ! ti6242

        if (c_array(cnum)%C11 .le. 0.0) then
            write(out,9001) cnum, 'C11'
            go to 1001
        end if

        if (c_array(cnum)%C12 .le. 0.0) then
            write(out,9001) cnum, 'C12'
            go to 1001
        end if

        if (c_array(cnum)%C13 .le. 0.0) then
            write(out,9001) cnum, 'C13'
            go to 1001
        end if

        if (c_array(cnum)%C33 .le. 0.0) then
            write(out,9001) cnum, 'C33'
            go to 1001
        end if

        if (c_array(cnum)%C44 .le. 0.0) then
            write(out,9001) cnum, 'C44'
            go to 1001
        end if

        if (c_array(cnum)%C55 .le. 0.0) then
            write(out,9001) cnum, 'C55'
            go to 1001
        end if

      end if ! elasticity = ti6242

      if (c_array(cnum)%harden_n .le. 0.0) then
            write(out,9001) cnum, 'harden_n'
            go to 1001
      end if

      if (c_array(cnum)%iD_v .lt. 0.0) then
            write (out,9003) cnum, 'iD_v'
            go to 1001
      end if

      if (c_array(cnum)%h_type .eq. 1) then
c           simple voche

      if (c_array(cnum)%k_o .lt. 0.0) then
            write (out,9003) cnum, 'k_o'
            go to 1001
      end if

c     ! Now softening is allowed
c      if (c_array(cnum)%theta_o .le. 0.0) then
c            write(out,9001) cnum, 'theta_o'
c            go to 1001
c      end if

      if (c_array(cnum)%tau_y .le. 0.0) then
            write(out,9001) cnum, 'tau_y'
            go to 1001
      end if

c     ! Now softening is allowed
c      if (c_array(cnum)%tau_v .le. 0.0) then
c            write(out,9001) cnum, 'tau_v'
c            go to 1001
c      end if

      if (c_array(cnum)%voche_m .le. 0.0) then
            write(out,9001) cnum, 'voche_m'
            go to 1001
      end if

      elseif (c_array(cnum)%h_type .eq. 2) then
c           MTS
c     Things that must be greater than zero:
c           tau_hat, g_o, b, p, q, boltz,
c           eps_do_o, mu_o, t_o, theta_o, tau_o
c
c     Handle errors locally
      if (c_array(cnum)%theta_o .le. 0.0) then
            write(out,9001) cnum, 'theta_o'
            go to 1001
      end if

      if (c_array(cnum)%k_o .lt. 0.0) then
            write (out,9003) cnum, 'k_o'
            go to 1001
      end if

      if (c_array(cnum)%tau_a .lt. 0.0) then
            write(out,9001) cnum, 'tau_a'
            go to 1001
      end if

      if (c_array(cnum)%tau_hat_y .le. 0.0) then
            write(out,9001) cnum, 'tau_hat_y'
            go to 1001
      end if

      if (c_array(cnum)%g_o_y .le. 0.0) then
            write(out,9001) cnum, 'g_o_y'
            go to 1001
      end if

      if (c_array(cnum)%tau_hat_v .le. 0.0) then
            write(out,9001) cnum, 'tau_hat_v'
            go to 1001
      end if

      if (c_array(cnum)%g_o_v .le. 0.0) then
            write(out,9001) cnum, 'g_o_v'
            go to 1001
      end if

      if (c_array(cnum)%b .le. 0.0) then
            write(out,9001) cnum, 'b'
            go to 1001
      end if

      if (c_array(cnum)%p_v .le. 0.0) then
            write(out,9001) cnum, 'p_v'
            go to 1001
      end if

      if (c_array(cnum)%q_v .le. 0.0) then
            write(out,9001) cnum, 'q_v'
            go to 1001
      end if

      if (c_array(cnum)%p_y .le. 0.0) then
            write(out,9001) cnum, 'p_y'
            go to 1001
      end if

      if (c_array(cnum)%q_y .le. 0.0) then
            write(out,9001) cnum, 'q_y'
            go to 1001
      end if

      if (c_array(cnum)%boltz .le. 0.0) then
            write(out,9001) cnum, 'boltz'
            go to 1001
      end if

      if (c_array(cnum)%eps_dot_o_v .le. 0.0) then
            write(out,9001) cnum, 'eps_dot_o_v'
            go to 1001
      end if

      if (c_array(cnum)%eps_dot_o_y .le. 0.0) then
            write(out,9001) cnum, 'eps_dot_o_y'
            go to 1001
      end if

      if (c_array(cnum)%mu_o .le. 0.0) then
            write(out,9001) cnum, 'mu_o'
            go to 1001
      end if

      if (c_array(cnum)%D_o .lt. 0.0) then
            write(out,9001) cnum, 'D_o'
            go to 1001
      end if

      if (c_array(cnum)%t_o .le. 0.0) then
            write(out,9001) cnum, 't_o'
            go to 1001
      end if

      else
c           User, just exit

      end if



      return

 1001 continue
      call die_gracefully

 9001 format(/1x,'>>>> Fatal error in crystal ', i3, '. Property ',
     &            a8, ' must be greater than zero.'/)
 9002 format(/1x,'>>>> Fatal error in crystal ', i3, '. Nu must be ',
     &            'less than 0.5.'/)
 9003 format(/1x,'>>>> Fatal error in crystal ', i3, '. Property ',
     &            a8, ' must not be less than zero.'/)

      end subroutine
