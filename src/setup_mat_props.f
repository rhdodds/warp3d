c     ****************************************************************
c     *                                                              *
c     *                 subroutine set_up_seg_type                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/29/2011                 *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_segmental_type( curve_set, segmental_type,
     &                               eps_curve_pts )
c
      use segmental_curves, only : seg_curve_table, seg_curves_type,
     &                             active_curve_set, num_seg_points,
     &                             curve_plseps_values
c
      implicit integer (a-z)
c
      double precision
     & eps_curve_pts(*)
c
c          set the type of segmental stress-strain curve: curve_set
c
c           =0 temperature and strain-rate independent
c           =1 temperature depdendent
c           =2 strain-rate dependent
c
c          set the active curve set in the module so it does not
c          have to be passed everywhere
c
      first_curve           = seg_curve_table(2,curve_set)
      segmental_type        = seg_curves_type(first_curve)
      num_points_on_curves  = num_seg_points(first_curve)
      active_curve_set      = curve_set
c
c          copy the plastic strain values for points on the
c          segmental curves into module vector so they don't have to
c          be copied everywhere.

      do i = 1, num_points_on_curves
       curve_plseps_values(i) = eps_curve_pts(i)
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine set_up_segmental                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/29/2011                 *
c     *                                                              *
c     *     set up the segmental stress-strain curve for each        *
c     *     element in this block for the gauss point being          *
c     *     processed. interpolate as necessary to produce a unique  *
c     *     curve for each element when response is temperature or   *
c     *     strain-rate dependent. also set temperature dependent    *
c     *     e, nu and isotropic thermal expansion values. Now        *
c     *     expanded to accommodate temperature varying props for    *
c     *     generalized_plasticity option of the cyclic model        *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_up_segmental( span, gp_temps, e_block, nu_block,
     & alpha_block, e_block_n, nu_block_n, alpha_block_n,
     & gp_sig_0_block, gp_h_u_block, gp_beta_u_block, gp_delta_u_block,
     & gp_sig_0_block_n, gp_h_u_block_n, gp_beta_u_block_n,
     & gp_delta_u_block_n,
     & gp_dtemps, gp_eps_rates, gpn, nrowd )
c
      use segmental_curves
c
      implicit integer (a-z)
c
c                      parameter declarations
c
      double precision
     & gp_temps(*), e_block(*), nu_block(*),
     & alpha_block(nrowd,*), gp_dtemps(*), gp_eps_rates(*),
     & e_block_n(*), nu_block_n(*), alpha_block_n(nrowd,*),
     & gp_sig_0_block(*), gp_h_u_block(*), gp_beta_u_block(*),
     & gp_delta_u_block(*), gp_sig_0_block_n(*), gp_h_u_block_n(*),
     & gp_beta_u_block_n(*), gp_delta_u_block_n(*)
c
c                      local declarations
c
      double precision
     & linear_interpolate,  stress_val, big_num, point_temp_np1,
     & point_temp_n, alpha, point_rate
c
      external linear_interpolate
      logical local_debug
      data local_debug, big_num / .false., 1.0e30 /
c
c
c          if stress-strain curve is temperature or strain rate
c          independent, just copy the first curve in the set
c          for use by all elements in the block. the previously
c          loaded e and nu values are constant for elements in block.
c
      curve_set             = active_curve_set
      first_curve           = seg_curve_table(2,curve_set)
      curve_set_type        = seg_curves_type(first_curve)
      num_points_on_curves  = num_seg_points(first_curve)
      num_curves_in_set     = seg_curve_table(1,curve_set)
c
      if ( gpn .eq. 1 ) then
         if ( allocated( sigma_curves ) ) deallocate( sigma_curves )
         if ( allocated( sigma_curve_min_values ) )
     &                       deallocate( sigma_curve_min_values )
         if ( allocated( curve_rates ) ) deallocate( curve_rates )
         if ( allocated( curve_temps ) ) deallocate( curve_temps )
         if ( allocated( curve_e_values ) ) deallocate( curve_e_values )
         if ( allocated( curve_nu_values ) )
     &                       deallocate( curve_nu_values )
         if ( allocated( curve_alpha_values ) )
     &                       deallocate( curve_alpha_values )
         if ( allocated( curve_gp_sig_0_values ) )
     &                       deallocate( curve_gp_sig_0_values )
         if ( allocated( curve_gp_h_u_values ) )
     &                       deallocate( curve_gp_h_u_values )
         if ( allocated( curve_gp_beta_u_values ) )
     &                       deallocate( curve_gp_beta_u_values )
         if ( allocated( curve_gp_delta_u_values ) )
     &                       deallocate( curve_gp_delta_u_values )
         if ( allocated( sigma_inter_table ) )
     &                       deallocate( sigma_inter_table  )
         ksize = num_curves_in_set
         allocate(
     &      sigma_curves(num_points_on_curves,span),
     &      sigma_curve_min_values(span), curve_temps(ksize),
     &      curve_e_values(ksize), curve_nu_values(ksize),
     &      curve_alpha_values(ksize),
     &      curve_gp_sig_0_values(ksize),
     &      curve_gp_h_u_values(ksize), curve_gp_beta_u_values(ksize),
     &      curve_gp_delta_u_values(ksize),
     &      curve_rates(ksize),
     &      sigma_inter_table(ksize,num_points_on_curves) )
      end if
c
      if ( local_debug ) then
         write(*,9000) first_curve, num_points_on_curves,
     &                 num_curves_in_set, curve_set_type
      end if
c
      select case ( curve_set_type )

      case( 0 )
c
c                temperature and rate independent curve. just use
c                the first curve in the curve set for all elements in
c                the block. the modulus, nu, and thermal expansion
c                coefficients have been already loaded into block
c                data structures.
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
       sigma_curves(1:num_points_on_curves,i) =
     &          seg_curves(1:num_points_on_curves,2,first_curve)
       sigma_curve_min_values(i) = seg_curves_min_stress(first_curve)
      end do
c
c                temperature dependent stress-strain curves.
c                interpolate to build a unique stress-strain curve for
c                this gauss point for each element in the block. first
c                build a transposed table so that a column has
c                the temperature dependent stress values for a plastic
c                strain point on the curve. we also build temperature
c                dependent elastic modulus, nu and isotropic thermal
c                expansion coefficients at end of step and at start
c                of step. note that temperature dependent alpha
c                values must be isotropic.
c
      case( 1 )
c
       if ( gpn .eq. 1 ) then
!DIR$ VECTOR ALIGNED
         do i = 1, num_curves_in_set
          curve_no                   = seg_curve_table(i+1,curve_set)
          curve_temps(i)             = seg_curves_value(curve_no)
          curve_e_values(i)          = seg_curves_ym(curve_no)
          curve_nu_values(i)         = seg_curves_nu(curve_no)
          curve_alpha_values(i)      = seg_curves_alpha(curve_no)
          curve_gp_sig_0_values(i)   = seg_curves_gp_sigma_0(curve_no)
          curve_gp_h_u_values(i)     = seg_curves_gp_h_u(curve_no)
          curve_gp_beta_u_values(i)  = seg_curves_gp_beta_u(curve_no)
          curve_gp_delta_u_values(i) = seg_curves_gp_delta_u(curve_no)
          do j = 1, num_points_on_curves
            sigma_inter_table(i,j) = seg_curves(j,2,curve_no)
          end do
         end do
       end if
c
       if ( local_debug ) write(*,9200)
c
       do ielem = 1, span
          point_temp_np1    = gp_temps(ielem)
          point_temp_n      = gp_temps(ielem) - gp_dtemps(ielem)
          sigma_curve_min_values(ielem) = big_num
c
          e_block(ielem) = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_e_values )
          e_block_n(ielem) = linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_e_values )
          nu_block(ielem) = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_nu_values )
          nu_block_n(ielem) = linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_nu_values )
          alpha =  linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_alpha_values )
          alpha_block(ielem,1:3) = alpha
          alpha =  linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_alpha_values )
          alpha_block_n(ielem,1:3) = alpha
c
          gp_sig_0_block(ielem) = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_sig_0_values )
          gp_sig_0_block_n(ielem) = linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_sig_0_values )
c
          gp_h_u_block(ielem) = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_h_u_values )
          gp_h_u_block_n(ielem) = linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_h_u_values )
c
          gp_beta_u_block(ielem) = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_beta_u_values )
          gp_beta_u_block_n(ielem) = linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_beta_u_values )
c
          gp_delta_u_block(ielem) = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_delta_u_values )
          gp_delta_u_block_n(ielem) = linear_interpolate( point_temp_n,
     &                        num_curves_in_set, curve_temps,
     &                        curve_gp_delta_u_values )
c
          do curve_pt = 1, num_points_on_curves
            stress_val = linear_interpolate( point_temp_np1,
     &                        num_curves_in_set, curve_temps,
     &                        sigma_inter_table(1,curve_pt) )
            sigma_curves(curve_pt,ielem) = stress_val
            sigma_curve_min_values(ielem) =
     &           min( sigma_curve_min_values(ielem), stress_val )
          end do
c
          if ( local_debug ) then
           write(*,9210) ielem, gp_temps(ielem),
     &                   sigma_curve_min_values(ielem)
           write(*,9215) e_block(ielem), nu_block(ielem)
           write(*,9220) (curve_pt, sigma_curves(curve_pt,ielem),
     &                    curve_pt=1,num_points_on_curves)
          end if
c
       end do
c
       return
c
c                strain-rate dependent stress-strain curves.
c                interpolate to build a unique stress-strain curve for
c                this gauss point for each element in the block. first
c                build a transposed table so that a column has
c                the rate dependent stress values for a plastic
c                strain point on the curve. the elastic modulus,
c                nu, thermal expansion are rate independent.

      case ( 2 )
c
       if ( gpn .eq. 1 ) then
         do i = 1, num_curves_in_set
           curve_no              = seg_curve_table(i+1,curve_set)
           curve_rates(i)        = seg_curves_value(curve_no)
           do j = 1, num_points_on_curves
             sigma_inter_table(i,j) = seg_curves(j,2,curve_no)
           end do
         end do
       end if
c
       if ( local_debug ) write(*,9240)
c
       do ielem = 1, span
          point_rate    = gp_eps_rates(ielem)
          sigma_curve_min_values(ielem) = big_num
          do curve_pt = 1, num_points_on_curves
            stress_val = linear_interpolate( point_rate,
     &                        num_curves_in_set, curve_rates,
     &                        sigma_inter_table(1,curve_pt) )
            sigma_curves(curve_pt,ielem) = stress_val
            sigma_curve_min_values(ielem) =
     &           min( sigma_curve_min_values(ielem), stress_val )
          end do
c
          if ( local_debug ) then
           write(*,9230) ielem,point_rate,sigma_curve_min_values(ielem)
           write(*,9220) (curve_pt, sigma_curves(curve_pt,ielem),
     &                   curve_pt=1,num_points_on_curves)
          end if
c
       end do
c
       if ( local_debug ) write(*,*) ' '
       return
c
      case default
        write(*,9100) 1
        stop
      end select
c
 9000 format('>> inside set_up_segmental...',
     & /,10x,'first_curve, num_points_on_curves: ',2i6,
     & /,10x,'num_curves_in_set, curve_set_type: ',2i6 )
 9100 format('>> FATAL ERROR:  routine set_up_segmental @ :',i3,
     & /,    '                 job aborted...' )
 9200 format(
     & '>> interpolated temperature dependent stress vs. plastic',
     & ' strain curve for elements in block....')
 9210 format(3x,'block element:',i7,2x,'gp temperature: ',f15.6,
     & ' min value on curve: ',f15.6)
 9215 format(10x,'ym, nu: ',f15.3,f10.3)
 9220 format(10x,i4,f15.6)
 9230 format(3x,'block element:',i7,2x,'gp strain rate: ',f15.6,
     & ' min value on curve: ',f15.6)
 9240 format(
     & '>> interpolated strain-rate dependent stress vs. plastic',
     & ' strain curve for elements in block....')
c
      end



c     ****************************************************************
c     *                                                              *
c     *                  function linear_interpolate                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 04/6/2018 rhd               *
c     *                                                              *
c     *    execute linear iterpolation on a tabular function         *
c     *    of a single variable where the x values are sorted in     *
c     *    increasing value but are not req'd to be uniformly spaced *
c     *                                                              *
c     ****************************************************************
c

      function linear_interpolate( xvalue, n, x, y ) result( ans )
      implicit none
      integer :: n, point
      double precision :: xvalue, x(n), y(n), x1, x2, y1, y2, ans
c
      if ( xvalue .le. x(1) ) then
        ans = y(1)
        return
      end if
c
      if ( xvalue .ge. x(n) ) then
        ans = y(n)
        return
      end if
c
      do point = 2, n
        if ( xvalue .gt. x(point) ) cycle
        x1 = x(point-1)
        x2 = x(point)
        y1 = y(point-1)
        y2 = y(point)
        ans = y1 + (xvalue-x1)*(y2-y1)/(x2-x1)
        return
      end do
c
      write(*,*) '>> FATAL ERROR: linear_interpolate'
      write(*,*) '                job aborted'
      stop
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine set_up_h_prime                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/3/00                    *
c     *                                                              *
c     *     set up h-prime for this gauss point for all elements     *
c     *     in the block. use the first 2 points on the segemetnal   *
c     *     stress-strain curve to define h-prime. this could be     *
c     *     temperature dependent. the interpolated curve at the     *
c     *     gauss point temperature is passed in                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_up_h_prime( span, h_block, sigyld_vec, felem )
      use segmental_curves, only : sigma_curves, curve_plseps_values
c
      implicit integer (a-z)
      include 'param_def'
c
c                      parameter declarations
c
      double precision
     & h_block(*), sigyld_vec(*)
c
c                      local declarations
c
      double precision
     &  dsigma, depspls
      logical local_debug
      data local_debug / .false. /
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
       dsigma        = sigma_curves(2,i) - sigma_curves(1,i)
       depspls       = curve_plseps_values(2) - curve_plseps_values(1)
       h_block(i)    = dsigma / depspls
       sigyld_vec(i) = sigma_curves(1,i)
      end do
c
      if ( local_debug ) then
       do i = 1, span
         write(*,9100) i+felem, dsigma, depspls, h_block(i),
     &                 sigyld_vec(i)
       end do
      end if
c
      return
c
 9100 format(i8, f10.2, f15.8, 2f10.2)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine set_e_nu_for_block                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/25/00                   *
c     *                                                              *
c     *     set up the temperature dependent young's moduli and      *
c     *     poisson's ratio for this gauss point for all elements    *
c     *     in the block                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_e_nu_for_block(
     &  span, nu_block, e_block, segmental, felem, etype, int_order,
     &  gpn, nnodel, temps_node_to_process, temps_node_blk )
      use global_data ! old common.main
c
      use segmental_curves, only :
     &      seg_curve_table, curve_temps, curve_e_values,
     &      curve_alpha_values, seg_curves_value, seg_curves_type,
     &      seg_curves_ym, seg_curves_nu, curve_nu_values
c
      implicit integer (a-z)
c
c                parameter declarations
c
      double precision
     &   e_block(*), nu_block(*), temps_node_blk(mxvl,*)
      logical segmental, temps_node_to_process
c
c                local declarations
c
      logical local_debug
      double precision
     &   sf(mxndel), xi, eta, zeta, weight, zero, gp_temps(mxvl),
     &   point_temp, linear_interpolate
      data local_debug, zero / .false., 0.0 /
c
c                we can only have temperature dependent elastic
c                properties if the material has segmental stress-
c                strain curves. if not, the current nu_block and
c                e_block values are constant for all elements in the
c                block.
c
      if ( .not. segmental ) return
c
c                the material for elements in the block is defined
c                by one or more segmental curves. if they are not
c                temperature dependent, the current nu_block and
c                e_block values are constant for all elements in the
c                block.
c
      curve_set         = iprops(21,felem)
      first_curve       = seg_curve_table(2,curve_set)
      num_curves_in_set = seg_curve_table(1,curve_set)
c
      if ( seg_curves_type(first_curve) .ne. 1 ) return
c
c                elastic properties are temperature dependent.
c                extract e, nu values and corresponding temperatures
c                from segmental curves (in order of increasing temp).
c                this is done only for the first gauss point of the
c                element block and saved for reuse during processing
c                of subsequent gauss points for element block.
c                all elements of the block must use the same set of
c                segmental stress-strain curves.
c
      if ( gpn .eq. 1 ) then
         if ( allocated( curve_temps ) )
     &     deallocate( curve_temps, curve_e_values,
     &                 curve_nu_values, curve_alpha_values )
         ksize = num_curves_in_set
         allocate(
     &      curve_temps(ksize), curve_e_values(ksize),
     &      curve_nu_values(ksize), curve_alpha_values(ksize) )
         do i = 1, num_curves_in_set
          curve_no           = seg_curve_table(i+1,curve_set)
          curve_temps(i)     = seg_curves_value(curve_no)
          curve_e_values(i)  = seg_curves_ym(curve_no)
          curve_nu_values(i) = seg_curves_nu(curve_no)
         end do
      end if
c
      if ( local_debug ) then
        write(*,9000) felem, curve_set, first_curve, num_curves_in_set
      end if
c
c                the material for elements in the block is described by
c                a set of temperature dependent, segmental stress-
c                strain curves. poisson's ratio and young's
c                modulus can be temperature dependent as well.
c
c                get the temperature at this gauss point for each
c                element in the block. interpolate from node
c                temperatures unless they are all zero.
c
!DIR$ VECTOR ALIGNED
      gp_temps(1:span) = zero
      if ( temps_node_to_process ) then
c
        call getgpts( etype, int_order, gpn, xi, eta, zeta, weight )
        call shapef( etype, xi, eta, zeta, sf(1) )
c
        do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
          do i = 1, span
            gp_temps(i) = gp_temps(i) +
     &                    sf(enode) * temps_node_blk(i,enode)
          end do
        end do
c
      end if
c
c                given young's modulus and poisson's ratio for
c                each curve in the set, use linear interpolation to
c                get the young's modulus and poisson ratio for
c                each element in the block at this gauss point.
c
      do ielem = 1, span
        point_temp     = gp_temps(ielem)
        e_block(ielem) = linear_interpolate( point_temp,
     &                        num_curves_in_set, curve_temps,
     &                        curve_e_values )
        nu_block(ielem) = linear_interpolate( point_temp,
     &                        num_curves_in_set, curve_temps,
     &                        curve_nu_values )
      end do
c
      if ( local_debug ) then
       write(*,9100) gpn
       do ielem = 1, span
         write(*,9200) ielem, gp_temps(ielem), e_block(ielem),
     &                 nu_block(ielem)
       end do
      end if
c
      return
c
 9000 format('>> inside set_e_nu_for_block. material is temperature',
     & ' dependent.',
     & /,4x,'felem, curve_set, first_curve, num_curves_in_set: ',i8,
     & 3i3 )
 9100 format('>> e and nu interpolated at elements for block, gpn:',i2,
     & /,'   rel. elem      temp            e          nu')
 9200 format(10x,i4, f12.4, f12.2, f10.4)
c
      end
c     ****************************************************************
c     *                                                              *
c     *           subroutine set_fgm_solid_props_for_block           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/15/01 M. Walters         *
c     *                                                              *
c     *     some elements in the block may require specific material *
c     *     properties to be interpolated from values specified at   *
c     *     model nodes via the fgm properties capability            *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_fgm_solid_props_for_block(
     &          span, felem, elem_type, gpn, nnode,
     &          blk_vec, shape, enode_mat_props, matcol, blk_fgm_flags )
      implicit integer (a-z)
      include 'param_def'
c
c                parameter declarations
c
      double precision
     &   blk_vec(*), shape(*), enode_mat_props(mxndel,mxvl,*),
     &   blk_fgm_flags(*)
c
c                local declarations
c
      logical local_debug, fgm
      double precision
     &   zero, fgm_mark, fgm_tol
      data local_debug, zero, fgm_mark, fgm_tol
     &     / .false., 0.0, 99.0, 1.0 /
c
      do i = 1, span
c
c               build vector of element properties at this gpn
c               interpolating from nodal values if necessary.
c               blk_fgm_flags(i) = -99 indicates that element
c               in block has an fgm material property requiring
c               interpolation.
c
        fgm = .false.
        if ( blk_fgm_flags(i) .lt. zero )
     &      fgm = abs(blk_fgm_flags(i)+fgm_mark) .le. fgm_tol
        if ( fgm ) then
            blk_vec(i) = zero
            do j = 1, nnode
              blk_vec(i) = blk_vec(i) +
     &                     shape(j)*enode_mat_props(j,i,matcol)
            end do
        end if
c
      end do
c
      if ( local_debug) then
         do node = 1, nnode
            write(*, 1000) node, shape(node)
         end do
         write(*,2000) gpn
         do ielem = 1, span
            write(*,3000) ielem, blk_vec(ielem)
         end do
      end if
c
      return
c
 1000 format('**shape function for node ',i4,2x,f15.8)
 2000 format('>> value for elements in block, gauss point # ',i2,
     & /,'          elem     block value')
 3000 format(10x,i7,2x,f12.4)
c
c
      end
