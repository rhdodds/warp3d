c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouhrks                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 8/1/2103 rhd               *
c     *                                                              *
c     *     this subroutine drives the computaion of the necessary   *
c     *     stress or strain hard copy output quantities for a block *
c     *     of non-conflicting, similar elements.                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks( span, blk, felem, type, order,
     &                   ngp, nnode,
     &                   geonl, long, nodpts, do_stresses,
     &                   mat_type, center_output, num_short_stress,
     &                   num_short_strain )
      use global_data ! old common.main
      implicit integer (a-z)
      logical do_stresses, geonl, long, nodpts, center_output
      double precision ::  nowtime

c
c                       get the element stresses or strains
c                       into a final form.
c                       for geonl: the cauchy stresses are computed
c                       from the unrotated cauchy stresses; the
c                       strains are accumulated spatial deformation
c                       increments.
c

      matnum = iprops(38,felem)  ! common.main
      kout = out                 !   "
      nowtime = total_model_time !   "
      nowstep = ltmstp           !   "
c
      call ougts1( span, blk, felem, do_stresses, ngp,
     &             geonl, mat_type, matnum, kout, nowtime, nowstep )
c
c                       for output at element center, make
c                       all gauss points have the average values.
c
      if ( center_output ) then
        call oumkcv( span, ngp, do_stresses, num_short_stress,
     &               num_short_strain )
      end if
c
c                       extrapolate gauss point results to element
c                       nodes. the 6 primary
c                       components, energy density, mises equiv.
c                       stress and material specific output values
c                       are extrapolated from gauss points to nodes.
c
      if( nodpts ) then
         call ounds1( span, type, order, nnode, ngp,
     &                do_stresses, num_short_stress, num_short_strain )
      end if
c
c                       for strain output, compute equivalent
c                       strain. if nodal, we are using extrapolated
c                       components. for center-output, we are using
c                       the averaged components. for stress output,
c                       we get the mises equivalent stress rather than
c                       equivalement strain.
c
      op_code = 1
      call ouext1( span, do_stresses, nnode, ngp, nodpts,
     &             num_short_stress, num_short_strain, op_code,
     &             mxvl )
c
c                       compute extra element quantities for long
c                       output (invariants, principal values, directions).
c                       for nodal output, the 6 primary strain/stress
c                       components are used. for nodal output, the
c                       extrapolated gp values are being used.
c
      if ( long ) then
        call ouext1( span, do_stresses, nnode, ngp, nodpts,
     &               num_short_stress, num_short_strain, 2, mxvl )
      end if
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine ouhrks_link                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/20/2017 rhd              *
c     *                                                              *
c     *               drive output for a single link element         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_link( elem, blk, felem, elem_type, num_enodes,
     &                       do_stress, mat_type, wide, eform, prec,
     &                       lnum, pgnum, lbltyp, strlbl, hedtyp,
     &                       noheader, out_packet_now, geo_non_flg  )
      use global_data ! old common.main
c
      use elblk_data, only : urcs_blk_n, elestr, ddtse
c
      implicit none
c
      integer :: elem, blk, felem, elem_type, num_enodes, mat_type,
     &           pgnum, lnum, lbltyp
      logical :: do_stress, wide, eform, prec, noheader,
     &           out_packet_now, geo_non_flg
      character(len=8) :: strlbl(30), hedtyp*30
c
c             locals
c
      integer :: nvals_to_output, bele
      logical :: newel
      logical, parameter :: local_debug = .false.
c
      if( local_debug ) write(out,9000) elem, blk, elem_type
c
      if( do_stress ) then
        nvals_to_output = 3
        elestr(1,1,1) = urcs_blk_n(1,1,1)
        elestr(1,2,1) = urcs_blk_n(1,2,1)
        elestr(1,3,1) = urcs_blk_n(1,3,1)
      else ! strain output
        nvals_to_output = 3
        elestr(1,1,1) = ddtse(1,1,1)
        elestr(1,2,1) = ddtse(1,2,1)
        elestr(1,3,1) = ddtse(1,3,1)
      end if
c
c             set labels for link element as needed.
c             then output values.
c
      call oulbst( do_stress, lbltyp, elem_type, elem, strlbl,
     &             .false., hedtyp, .false., 0 )
c
      bele  = 1
      newel = .true.
c
      call ouhel( bele, elem, hedtyp, 2, 1,
     &            .false., nvals_to_output, wide,
     &            eform, prec, strlbl, pgnum, lnum, newel,
     &            .false., noheader, out_packet_now )
c
      return

 9000 format(/,'.... entered ouhkrs_link ....',
     &  /,10x,'elem, blk, etype: ', i8, i7, i5 )
 9010 format(15x,'node1:2, x1,y1,z1:', 2i7,3f10.3,
     & /,15x,'x2,y2,z2, len0:',4f10.3,
     & /,15x,'u1,v1,w1: ',3f15.6,
     & /,15x,'u2,v2,w2: ',3f15.6 )
 9100 format(/,">>>> FATAL ERROR. oulbst_bar. val: ",i10,
     & /,      "                  job aborted...",//)
c
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine ouhrks_bar                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/4/2017 rhd               *
c     *                                                              *
c     *               drive output for a single bar element          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_bar( elem, blk, felem, elem_type, num_enodes,
     &                       do_stress, mat_type, wide, eform, prec,
     &                       lnum, pgnum, lbltyp, strlbl, hedtyp,
     &                       noheader, out_packet_now, geo_non_flg  )
      use global_data ! old common.main
c
      use elblk_data, only : urcs_blk_n, ddtse, elestr
      use main_data, only : incmap, incid, crdmap
c
      implicit none
c
      integer :: elem, blk, felem, elem_type, num_enodes, mat_type,
     &           pgnum, lnum, lbltyp
      logical :: do_stress, wide, eform, prec, noheader,
     &           out_packet_now, geo_non_flg
      character(len=8) :: strlbl(30), hedtyp*30
c
c             locals
c
      integer :: nvals_to_output, pos, node1, node2, bele
      double precision :: x1, x2, y1, y2, z1, z2, len0, u1, u2, v1, v2,
     &                    w1, w2, len, area
      double precision, parameter :: zero = 0.0d0
      logical :: newel
      logical, parameter :: local_debug = .false.
c
      if( local_debug ) write(out,9000) elem, blk, elem_type
c
      if( do_stress ) then  !  must get area @ n1 for geonl
        nvals_to_output = 2
        elestr(1,1,1) = urcs_blk_n(1,1,1)
        pos   = incmap(elem)
        node1 = incid(pos)
        node2 = incid(pos+1)
        x1 = c(crdmap(node1))
        y1 = c(crdmap(node1)+1)
        z1 = c(crdmap(node1)+2)
        x2 = c(crdmap(node2))
        y2 = c(crdmap(node2)+1)
        z2 = c(crdmap(node2)+2)
        len0 = sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
        u1 = u(node1)
        u2 = u(node2)
        v1 = u(node1+nonode)
        v2 = u(node2+nonode)
        w1 = u(node1+2*nonode)
        w2 = u(node2+2*nonode)
c        if( local_debug ) write(out,9010) node1, node2, x1, y1, z1,
c     &      x2, y2, z2, len0, u1, v1, w1, u2, v2, w2
        y1 = y1 + v1
        z1 = z1 + w1
        x2 = x2 + u2
        y2 = y2 + v2
        z2 = z2 + w2
        len = sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
        area = props(43,elem)
        if( geo_non_flg ) area = area * len0 / len
        elestr(1,2,1) = urcs_blk_n(1,1,1) * area
      else ! strain output
        nvals_to_output = 1
        elestr(1,1,1) = ddtse(1,1,1)
      end if
c
c             set labels for bar element as needed. then output values.
c
      call oulbst( do_stress, lbltyp, elem_type, elem, strlbl,
     &             .false., hedtyp, geo_non_flg, 0 )
c
      bele  = 1
      newel = .true.
c
      call ouhel( bele, elem, hedtyp, 2, 1,
     &            .false., nvals_to_output, wide,
     &            eform, prec, strlbl, pgnum, lnum, newel,
     &            .false., noheader, out_packet_now )
c
      return

 9000 format(/,'.... entered ouhkrs_bar ....',
     &  /,10x,'elem, blk, etype: ', i8, i7, i5 )
 9010 format(15x,'node1:2, x1,y1,z1:', 2i7,3f10.3,
     & /,15x,'x2,y2,z2, len0:',4f10.3,
     & /,15x,'u1,v1,w1: ',3f15.6,
     & /,15x,'u2,v2,w2: ',3f15.6 )
 9100 format(/,">>>> FATAL ERROR. oulbst_bar. val: ",i10,
     & /,      "                  job aborted...",//)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                subroutine ouhrks_cohesive                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/1/2013 rhd               *
c     *                                                              *
c     *     drive output for a single interface-cohesive element     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_cohesive( elem, blk, felem, elem_type,
     1    int_order, num_int_points, num_enodes, center_output,
     2    do_stress, mat_type, cohesive_type,wide,eform,
     3    prec, lnum, pgnum, lbltyp, strlbl, hedtyp, noheader,
     4    out_packet_now, geo_non_flg  )
      use global_data ! old common.main
c
      use elblk_data, only : elem_hist, urcs_blk_n, ddtse, elestr

      implicit integer (a-z)
      logical :: do_stress, center_output, wide, eform, newel,
     &        prec, noheader, out_packet_now, geo_non_flg
      character(len=8) :: strlbl(30), hedtyp*30
c
c             locals also visible by contains
c
      double precision :: zero,  avgs(mxoupr)
      logical :: local_debug, exp_model, ppr_model
      integer :: kout
      data zero / 0.0d00 /
c
c             get the element tractions or displacement jumps
c             with suplemental values in final form for output.
c
      local_debug = .false.
      kout = out
c
      if( local_debug ) write(kout,9000) elem, blk, elem_type,
     &       num_int_points, center_output, cohesive_type
c
      avgs(1:mxoupr) = zero
c
      select case( cohesive_type )
c
      case( 1, 2, 3, 5 )
c
        if( do_stress ) then
          call ouhrks_simple_tractions
          nvals_to_output = 6
        else
          call ouhrks_simple_displ_jumps
          nvals_to_output = 4
        end if
c
      case( 4 )  ! simple exponential model
c
        if( do_stress ) then
          call ouhrks_exp1_tractions
          nvals_to_output = 8
        else
          call ouhrks_exp1_displ_jumps
          nvals_to_output = 6
        end if
c
      case( 6 )  ! ppr
c
        if( do_stress ) then
          call ouhrks_ppr_tractions
          nvals_to_output = 8
        else
          call ouhrks_ppr_displ_jumps
          nvals_to_output = 6
        end if
c
      case( 7 )  ! cavit
c
        if( do_stress ) then
          call ouhrks_cavit_tractions
          nvals_to_output = 13
        else
          call ouhrks_simple_displ_jumps
          nvals_to_output = 4
        end if
c
      case default
c
        write(out,9100) cohesive_type
        call die_abort
c
       end select
c
c             set the traction-displacement jump labels for
c             the current element. then output values.
c
      call oulbst( do_stress, lbltyp, elem_type, elem, strlbl,
     &             long, hedtyp, geo_non_flg, cohesive_type )

      bele            = 1
      newel           = .true.
c
      call ouhel( bele, elem, hedtyp, num_enodes, num_int_points,
     &            .false., nvals_to_output, wide,
     &            eform, prec, strlbl, pgnum, lnum, newel,
     &            center_output, noheader, out_packet_now )
c
      return

 9000 format(/,'.... entered ouhkrs_cohesive ....',
     &  /,10x,'elem, blk, etype, npts, cntrout, c_type: ',
     &  i8,i4,i3,i3,l3,i4 )
 9100 format(/,">>>> FATAL ERROR. oulbst_cohesive. val: ",i10,
     & /,      "                  job aborted...",//)


      contains
c     ****************************************************************
c     *                                                              *
c     *              subroutine ouhrks_cavit_tractions               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/4/2015 rhd               *
c     *                                                              *
c     *     pull together traction-type output values for            *
c     *     the cavity cohesive material type                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_cavit_tractions
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   t1, t2, tn, ts, gamma, gamma_ur, n_bar, a_bar, b_bar,
     &   v_bar, criticality, max_Tn, max_Dn, omega, lambda_4,
     &   a_over_l_nr
c
      if( local_debug ) write(kout,9020)
c
      do gpn = 1, num_int_points
        t1               = urcs_blk_n(1,1,gpn)
        t2               = urcs_blk_n(1,2,gpn)
        tn               = urcs_blk_n(1,3,gpn)
        gamma            = urcs_blk_n(1,7,gpn)
        gamma_ur         = urcs_blk_n(1,8,gpn)
        ts               = sqrt( t1*t1 + t2*t2 )
        elestr(1,1,gpn)  = t1
        elestr(1,2,gpn)  = t2
        elestr(1,3,gpn)  = ts
        elestr(1,4,gpn)  = tn
        elestr(1,5,gpn)  = gamma
        elestr(1,6,gpn)  = gamma_ur
        n_bar            = elem_hist(1,1,gpn)
        a_bar            = elem_hist(1,2,gpn)
        b_bar            = elem_hist(1,3,gpn)
        lambda_4         = elem_hist(1,11,gpn)
        omega            = lambda_4 * a_bar / b_bar ! a/b
        v_bar            = elem_hist(1,4,gpn)
        criticality      = elem_hist(1,8,gpn)
        max_Tn           = elem_hist(1,9,gpn)
        max_Dn           = elem_hist(1,10,gpn)
        a_over_l_nr      = elem_hist(1,12,gpn)
        elestr(1,7,gpn)  = n_bar
        elestr(1,8,gpn)  = omega
        elestr(1,9,gpn)  = a_bar
        elestr(1,10,gpn) = criticality
        elestr(1,11,gpn) = max_Tn
        elestr(1,12,gpn) = max_Dn
        elestr(1,13,gpn) = a_over_l_nr
        avgs(1:13)     = avgs(1:13) + elestr(1,1:13,gpn)
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:6,gpn)
      end do
c
      if( center_output ) then
        avgs(1:13)      = avgs(1:13) / dble( num_int_points )
        t1              = avgs(1)
        t2              = avgs(2)
        ts              = sqrt( t1**2 + t2**2 )
        tn              = avgs(4)
        gamma           = avgs(5)
        gamma_ur        = avgs(6)
        n_bar           = avgs(7)
        omega           = avgs(8)  ! a / b average over element
        a_bar           = avgs(9)
        criticality     = avgs(10) ! need to think what this means
        max_Tn          = avgs(11)
        max_Dn          = avgs(12)
        a_over_l_nr     = avgs(13)
c
        elestr(1,1,1) = t1
        elestr(1,2,1) = t2
        elestr(1,3,1) = ts
        elestr(1,4,1) = tn
        elestr(1,5,1) = gamma
        elestr(1,6,1) = gamma_ur
        elestr(1,7,1) = n_bar
        elestr(1,8,1) = omega
        elestr(1,9,1) = a_bar
        elestr(1,10,1) = criticality
        elestr(1,11,1) = max_Tn
        elestr(1,12,1) = max_Dn
        elestr(1,13,1) = a_over_l_nr
c
        if( local_debug ) write(kout,9032) elestr(1,1:6,1)
      end if
c
      return
c
 9020 format(/,"... cavity cohesive model: tractions ...",/,
     & 10x,'int pt',5x,'t1',9x,'t2',10x,'ts',10x,'tn',
     &  5x,'E(tot)',6x,'E(plas)' )
 9030 format(10x,i2,6f12.4)
 9032 format(9x,'avg',6f12.4)
c
      end subroutine ouhrks_cavit_tractions
c     ****************************************************************
c     *                                                              *
c     *              subroutine ouhrks_simple_tractions              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
c     *                                                              *
c     *     pull together tractions output values for other          *
c     *     cohesive types                                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_simple_tractions
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   t1, t2, tn, ts, gamma, gamma_ur
c
      if( local_debug ) write(kout,9020)
c
      do gpn = 1, num_int_points
        t1       = urcs_blk_n(1,1,gpn)
        t2       = urcs_blk_n(1,2,gpn)
        tn       = urcs_blk_n(1,3,gpn)
        gamma    = urcs_blk_n(1,7,gpn)
        gamma_ur = urcs_blk_n(1,8,gpn)
        ts       = sqrt( t1*t1 + t2*t2 )
        elestr(1,1,gpn) = t1
        elestr(1,2,gpn) = t2
        elestr(1,3,gpn) = ts
        elestr(1,4,gpn) = tn
        elestr(1,5,gpn) = gamma
        elestr(1,6,gpn) = gamma_ur
        avgs(1:6)       = avgs(1:6) + elestr(1,1:6,gpn)
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:6,gpn)
      end do
c
      if( center_output ) then
        avgs(1:6)       = avgs(1:6) / dble( num_int_points )
        t1              = avgs(1)
        t2              = avgs(2)
        ts              = sqrt( t1**2 + t2**2 )
        tn              = avgs(4)
        gamma           = avgs(5)
        gamma_ur        = avgs(6)
c
        elestr(1,1,1) = t1
        elestr(1,2,1) = t2
        elestr(1,3,1) = ts
        elestr(1,4,1) = tn
        elestr(1,5,1) = gamma
        elestr(1,6,1) = gamma_ur
        if( local_debug ) write(kout,9032) elestr(1,1:6,1)
      end if
c
      return
c
 9020 format(/,"... other cohesive model: tractions ...",/,
     & 10x,'int pt',5x,'t1',9x,'t2',10x,'ts',10x,'tn',
     &  5x,'E(tot)',6x,'E(plas)' )
 9030 format(10x,i2,6f12.4)
 9032 format(9x,'avg',6f12.4)
c
      end subroutine ouhrks_simple_tractions
c

c     ****************************************************************
c     *                                                              *
c     *          subroutine ouhrks_simple_displ_jumps                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
c     *                                                              *
c     *     pull together displacement jumps output values for       *
c     *     other cohesive material type                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_simple_displ_jumps
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   d1, d2, dn, ds
c
      if( local_debug ) write(kout,9020)
c
      do gpn = 1, num_int_points
c
        d1              = ddtse(1,1,gpn)
        d2              = ddtse(1,2,gpn)
        dn              = ddtse(1,3,gpn)
        ds              = sqrt( d1*d1 + d2*d2 )
c
        elestr(1,1,gpn) = d1
        elestr(1,2,gpn) = d2
        elestr(1,3,gpn) = ds
        elestr(1,4,gpn) = dn
        avgs(1:4)       = avgs(1:4) + elestr(1,1:4,gpn)
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:4,gpn)
c
      end do
c
      if( center_output ) then
        avgs(1:4)       = avgs(1:4) / dble( num_int_points )
        d1              = avgs(1)
        d2              = avgs(2)
        ds              = sqrt( d1**2 + d2**2 )
        dn              = avgs(4)
c
        elestr(1,1,1) = d1
        elestr(1,2,1) = d2
        elestr(1,3,1) = ds
        elestr(1,4,1) = dn
        if( local_debug ) write(kout,9032) elestr(1,1:4,1)
      end if
c
      return
c
 9020 format(/,"... other cohesive model: displ jumps ...",/,
     & 10x,'int pt',5x,'d1',9x,'d2',10x,'ds',10x,'dn' )
 9030 format(10x,i2,4f12.4)
 9032 format(9x,'avg',4f12.4)
c
       end subroutine ouhrks_simple_displ_jumps

c     ****************************************************************
c     *                                                              *
c     *             subroutine ouhrks_ppr_displ_jumps                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
c     *                                                              *
c     *     pull together displacement jumps output values for       *
c     *     ppr cohesive                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_ppr_displ_jumps
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   d1, d2, dn, ds, dn_limit, ds_limit,
     &   ratio_normal, ratio_shear, dn_at_peak, ds_at_peak
c
      if( local_debug ) write(kout,9020)
c
      ratio_normal       = props(37,elem)
      ratio_shear        = props(39,elem)
c
      do gpn = 1, num_int_points
c
        dn_limit   = elem_hist(1,1,gpn)
        ds_limit   = elem_hist(1,2,gpn)
        ds_at_peak = ds_limit * ratio_shear
        dn_at_peak = dn_limit * ratio_normal
c
        d1              = ddtse(1,1,gpn)
        d2              = ddtse(1,2,gpn)
        dn              = ddtse(1,3,gpn)
        ds              = sqrt( d1*d1 + d2*d2 )
c
        elestr(1,1,gpn) = d1
        elestr(1,2,gpn) = d2
        elestr(1,3,gpn) = ds
        elestr(1,4,gpn) = dn
        elestr(1,5,gpn) = ds / ds_at_peak
        elestr(1,6,gpn) = dn / dn_at_peak
        avgs(1:6)       = avgs(1:6) + elestr(1,1:6,gpn)
        avgs(7)         = avgs(7) + ds_at_peak
        avgs(8)         = avgs(8) + dn_at_peak
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:6,gpn)
c
      end do
c
      if( center_output ) then
        avgs(1:8)       = avgs(1:8) / dble( num_int_points )
        d1              = avgs(1)
        d2              = avgs(2)
        ds              = sqrt( d1**2 + d2**2 )
        dn              = avgs(4)
        ds_at_peak      = avgs(7)
        dn_at_peak      = avgs(8)
c
        elestr(1,1,1) = d1
        elestr(1,2,1) = d2
        elestr(1,3,1) = ds
        elestr(1,4,1) = dn
        elestr(1,5,1) = ds / ds_at_peak
        elestr(1,6,1) = dn / dn_at_peak
        if( local_debug ) write(kout,9032) elestr(1,1:6,1)
      end if
c
      return
c
 9020 format(/,"... ppr model: displ jumps ...",/,
     & 10x,'int pt',5x,'d1',9x,'d2',10x,'ds',10x,'dn',
     &  9x,'ds/ds-peak',5x,'dn/dn-peak' )
 9030 format(10x,i2,6f12.4)
 9032 format(9x,'avg',6f12.4)
c
       end subroutine ouhrks_ppr_displ_jumps

c     ****************************************************************
c     *                                                              *
c     *             subroutine ouhrks_exp1_displ_jumps               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
c     *                                                              *
c     *     pull together displacement jumps output values for       *
c     *     exponential cohesive                                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_exp1_displ_jumps
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   deff_peak, beta, d1, d2, dn, ds, deff,
     &   normalized_deff
c
      if( local_debug ) write(kout,9020)
c
      deff_peak          = props(22,elem)
      beta               = props(23,elem)
c
      do gpn = 1, num_int_points
        d1              = ddtse(1,1,gpn)
        d2              = ddtse(1,2,gpn)
        dn              = ddtse(1,3,gpn)
        ds              = sqrt( d1*d1 + d2*d2 )
        deff            = sqrt( beta*beta*ds*ds + dn*dn )
        normalized_deff = deff / deff_peak
        elestr(1,1,gpn) = d1
        elestr(1,2,gpn) = d2
        elestr(1,3,gpn) = ds
        elestr(1,4,gpn) = dn
        elestr(1,5,gpn) = deff
        elestr(1,6,gpn) = normalized_deff
        avgs(1:6)       = avgs(1:6) + elestr(1,1:6,gpn)
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:6,gpn)
      end do
c
      if( center_output ) then
        avgs(1:6)       = avgs(1:6) / dble( num_int_points )
        d1              = avgs(1)
        d2              = avgs(2)
        ds              = sqrt( d1**2 + d2**2 )
        dn              = avgs(4)
        deff            = sqrt( beta*beta*ds*ds + dn*dn )
        normalized_deff = deff / deff_peak
c
        elestr(1,1,1) = d1
        elestr(1,2,1) = d2
        elestr(1,3,1) = ds
        elestr(1,4,1) = dn
        elestr(1,5,1) = deff
        elestr(1,6,1) = normalized_deff
        if( local_debug ) write(kout,9032) elestr(1,1:6,1)
      end if
c
      return
c
 9020 format(/,"... exponential model: displ jumps ...",/,
     & 10x,'int pt',5x,'d1',9x,'d2',10x,'ds',10x,'dn',
     &  9x,'deff',5x,'deff/dpeak' )
 9030 format(10x,i2,6f12.4)
 9032 format(9x,'avg',6f12.4)
c
       end subroutine ouhrks_exp1_displ_jumps

c     ****************************************************************
c     *                                                              *
c     *                subroutine ouhrks_exp1_tractions              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
c     *                                                              *
c     *     pull together tractions output values for exponential    *
c     *     cohesive                                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_exp1_tractions
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   t1, t2, tn, gamma, gamma_ur,
     &   teff, teff_peak, beta,ts, normalized_teff
c
      if( local_debug ) write(kout,9020)
c
      teff_peak = props(13,elem)
      beta      = props(23,elem)
c
      do gpn = 1, num_int_points
        t1       = urcs_blk_n(1,1,gpn)
        t2       = urcs_blk_n(1,2,gpn)
        tn       = urcs_blk_n(1,3,gpn)
        gamma    = urcs_blk_n(1,7,gpn)
        gamma_ur = urcs_blk_n(1,8,gpn)
        teff     = elem_hist(1,3,gpn)
        ts       = sqrt( t1*t1 + t2*t2 )
        normalized_teff = teff / teff_peak
        elestr(1,1,gpn) = t1
        elestr(1,2,gpn) = t2
        elestr(1,3,gpn) = ts
        elestr(1,4,gpn) = tn
        elestr(1,5,gpn) = teff
        elestr(1,6,gpn) = normalized_teff
        elestr(1,7,gpn) = gamma
        elestr(1,8,gpn) = gamma_ur
        avgs(1:8)     = avgs(1:8) + elestr(1,1:8,gpn)
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:8,gpn)
      end do
c
      if( center_output ) then
        avgs(1:8)       = avgs(1:8) / dble( num_int_points )
        t1              = avgs(1)
        t2              = avgs(2)
        ts              = sqrt( t1**2 + t2**2 )
        tn              = avgs(4)
        teff            = sqrt( ts*ts/beta/beta + tn**2 )
        normalized_teff = teff / teff_peak
        gamma           = avgs(7)
        gamma_ur        = avgs(8)
c
        elestr(1,1,1) = t1
        elestr(1,2,1) = t2
        elestr(1,3,1) = ts
        elestr(1,4,1) = tn
        elestr(1,5,1) = teff
        elestr(1,6,1) = normalized_teff
        elestr(1,7,1) = gamma
        elestr(1,8,1) = gamma_ur
        if( local_debug ) write(kout,9032) elestr(1,1:8,1)
      end if
c
      return
c
 9020 format(/,"... exponential model: tractions ...",/,
     & 10x,'int pt',5x,'t1',9x,'t2',10x,'ts',10x,'tn',
     &  9x,'teff',5x,'teff/tmax',5x,'E(tot)',6x,'E(plas)' )
 9030 format(10x,i2,8f12.4)
 9032 format(9x,'avg',8f12.4)
c
       end subroutine ouhrks_exp1_tractions
c     ****************************************************************
c     *                                                              *
c     *                subroutine ouhrks_ppr_tractions               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/13 rhd                *
c     *                                                              *
c     *     pull together traction output values for ppr cohesive    *
c     *     model                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhrks_ppr_tractions
      implicit integer (a-z)
c
c             locals (see also contains main)
c
      double precision
     &   tn_peak, ts_peak, t1, t2, tn, gamma,
     &   gamma_ur, ts
c
      tn_peak = props(13,elem)
      ts_peak = props(14,elem)
c
      if( local_debug ) write(kout,9020)

c
      do gpn = 1, num_int_points
        t1       = urcs_blk_n(1,1,gpn)
        t2       = urcs_blk_n(1,2,gpn)
        tn       = urcs_blk_n(1,3,gpn)
        gamma    = urcs_blk_n(1,7,gpn)
        gamma_ur = urcs_blk_n(1,8,gpn)
        ts       = sqrt( t1*t1 + t2*t2 )
        elestr(1,1,gpn) = t1
        elestr(1,2,gpn) = t2
        elestr(1,3,gpn) = ts
        elestr(1,4,gpn) = tn
        elestr(1,5,gpn) = ts / ts_peak
        elestr(1,6,gpn) = tn / tn_peak
        elestr(1,7,gpn) = gamma
        elestr(1,8,gpn) = gamma_ur
        avgs(1:8)     = avgs(1:8) + elestr(1,1:8,gpn)
        if( local_debug ) write(kout,9030) gpn, elestr(1,1:8,gpn)
      end do
c
      if( center_output ) then
        avgs(1:8)       = avgs(1:8) / dble( num_int_points )
        t1              = avgs(1)
        t2              = avgs(2)
        ts              = sqrt( t1**2 + t2**2 )
        tn              = avgs(4)
        gamma           = avgs(7)
        gamma_ur        = avgs(8)
c
        elestr(1,1,1) = t1
        elestr(1,2,1) = t2
        elestr(1,3,1) = ts
        elestr(1,4,1) = tn
        elestr(1,5,1) = ts / ts_peak
        elestr(1,6,1) = tn / tn_peak
        elestr(1,7,1) = gamma
        elestr(1,8,1) = gamma_ur
        if( local_debug ) write(kout,9032) elestr(1,1:8,1)
      end if
c
      return
c
 9020 format(/,"... ppr model: tractions ...",/,
     & 10x,'int pt',5x,'t1',9x,'t2',10x,'ts',10x,'tn',
     &  5x,'ts/tpeak',6x,'tn/tpeak',5x,'E(tot)', 6x,'E(plas)' )
 9030 format(10x,i2,8f12.4)
 9032 format(9x,'avg',8f12.4)
c
      end subroutine ouhrks_ppr_tractions
c
      end subroutine ouhrks_cohesive
