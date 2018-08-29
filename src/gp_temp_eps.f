c     ****************************************************************
c     *                                                              *
c     *                      subroutine gp_temps                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/1/2018 rhd               *
c     *                                                              *
c     *     compute increment of temperature at a gauss point for    *
c     *     all elements in the block.                               *
c     *                                                              *
c     *     compute current temperature at end of step for           *
c     *     all elements of the block.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gauss_pt_temps(
     &       dtemps_node_blk, gpn, etype, span, int_order,
     &       nnodel, gp_dtemps, temps_node_blk, gp_temps,
     &       temper_increment, temps_node_to_process,
     &       temperatures_ref, temps_ref_node_blk, gp_rtemps  )
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      integer :: gpn, etype, span, int_order, nnodel
      double precision ::
     &  dtemps_node_blk(mxvl,*), gp_dtemps(*), gp_temps(*),
     &  temps_node_blk(mxvl,*), temps_ref_node_blk(mxvl,*),
     &  gp_rtemps(*)
      logical :: temper_increment, temps_node_to_process,
     &           temperatures_ref
c
c                     locally defined arrays-variables
c
      integer :: i, k, enode
      double precision :: sf(mxndel), xi, eta, zeta, weight
      double precision, parameter :: zero = 0.0d0
      logical, parameter :: local_debug = .false.
c
c                     if there are incremental temperature changes
c                     imposed in this load step, interpolate temperature
c                     change at this gauss point for all elements
c                     in the block from nodal values.
c
c                     get the parametric coordinates for this
c                     integration point. then get the nodal
c                     shape functions evaluated at the point.
c
      if( temper_increment .or. temps_node_to_process .or.
     &     temperatures_ref ) then
         call getgpts( etype, int_order, gpn, xi, eta, zeta, weight )
         call shapef( etype, xi, eta, zeta, sf(1) )
      end if
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         gp_dtemps(i) = zero
         gp_temps(i)  = zero
         gp_rtemps(i) = zero
      end do
c
      if( temper_increment ) then
         do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
           do i = 1, span
             gp_dtemps(i) = gp_dtemps(i) +
     &                      sf(enode) * dtemps_node_blk(i,enode)
           end do
         end do
      end if
c
c                     interpolate the current temperature at the end
c                     of this load step for this gauss point for all elements
c                     in the block from nodal values.
c
      if( temps_node_to_process ) then
        do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
          do i = 1, span
            gp_temps(i) = gp_temps(i) +
     &                    sf(enode) * temps_node_blk(i,enode)
          end do
        end do
      end if
c
c                     interpolate the reference temperature
c                     for this gauss point for all elements
c                     in the block from nodal values.
c
      if( temperatures_ref ) then
        do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
          do i = 1, span
            gp_rtemps(i) = gp_rtemps(i) +
     &                     sf(enode) * temps_ref_node_blk(i,enode)
          end do
        end do
      end if
c
      if ( .not. local_debug ) return
         write(*,*) '>> in  gauss_pt_temps:'
         write(*,*) '>> temps_ref_node_block:'
         do i = 1, span
          write(*,*) '   element: ',i
          write(*,*) (k,temps_ref_node_blk(i,k),k=1,nnodel)
         end do
         write(*,*) '>> gp_temps...'
         write(*,9000) (i,gp_temps(i),i=1,span)
         write(*,*) '>> gp_dtemps...'
         write(*,9000) (i,gp_dtemps(i),i=1,span)
         write(*,*) '>> gp_rtemps...'
         write(*,9000) (i,gp_rtemps(i),i=1,span)
      return
c
 9000 format(1x,i3,f15.6 )
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine gp_temp_eps                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/1/2018 rhd               *
c     *                                                              *
c     *     compute incremental thermal strains at gauss point for   *
c     *     all elements of the block. subtract them from incr.      *
c     *     strains due to displacements.                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gp_temp_eps( span, deps, alpha_n1, gp_dtemps,
     &                        gp_temps, gp_rtemps, alpha_n, etype  )
      use main_data,  only : bar_types, link_types
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      integer :: span, etype
      double precision :: deps(mxvl,*), alpha_n1(mxvl,6),
     &                    alpha_n(mxvl,6), gp_dtemps(*),
     &                    gp_temps(*), gp_rtemps(*)
c
c                     locals
c
      integer :: i, j
      double precision :: theta_n, dtn, dtn1, deps_theta(6)
      double precision, parameter :: zero = 0.0d0
      logical :: is_bar_elem, is_link_elem
      logical, parameter :: local_debug = .false.
c
c            compute incremental thermal strain and subtract from
c            strain increment due to displacements.
c
      if( local_debug ) write(*,9100) (i,(deps(i,j),j=1,6),i=1,span)
      is_bar_elem  = bar_types(etype)
      is_link_elem = link_types(etype)
c
c           \Delta\veps_\theta = \alpha_{n+1) T_{n+1)  -
c                                \alpha_n T_n
c           adjust for reference temp
c
c           \alpha here can be anisotropic
!DIR$ VECTOR ALIGNED
      do i = 1, span
c
       theta_n = gp_temps(i) -  gp_dtemps(i)
       dtn     = theta_n     -  gp_rtemps(i)
       dtn1    = gp_temps(i) -  gp_rtemps(i)
c
       deps_theta(1) = alpha_n1(i,1)*dtn1 - alpha_n(i,1)*dtn
       deps_theta(2) = alpha_n1(i,2)*dtn1 - alpha_n(i,2)*dtn
       deps_theta(3) = alpha_n1(i,3)*dtn1 - alpha_n(i,3)*dtn
       deps_theta(4) = alpha_n1(i,4)*dtn1 - alpha_n(i,4)*dtn
       deps_theta(5) = alpha_n1(i,5)*dtn1 - alpha_n(i,5)*dtn
       deps_theta(6) = alpha_n1(i,6)*dtn1 - alpha_n(i,6)*dtn
       deps(i,1) = deps(i,1) -   deps_theta(1)
       deps(i,2) = deps(i,2) -   deps_theta(2)
       deps(i,3) = deps(i,3) -   deps_theta(3)
       deps(i,4) = deps(i,4) -   deps_theta(4)
       deps(i,5) = deps(i,5) -   deps_theta(5)
       deps(i,6) = deps(i,6) -   deps_theta(6)
c
      end do
c
      if( is_bar_elem .or. is_link_elem ) deps(1:span,2:6) = zero
c
      if ( .not. local_debug ) return
      write(*,*) '>> thermal strain increments...'
      do i = 1, span
       theta_n = gp_temps(i) -  gp_dtemps(i)
       dtn     = theta_n -  gp_rtemps(i)
       dtn1    = gp_temps(i) -  gp_rtemps(i)
       deps_theta(1) = alpha_n1(i,1)*dtn1 - alpha_n(i,1)*dtn
       deps_theta(2) = alpha_n1(i,2)*dtn1 - alpha_n(i,2)*dtn
       deps_theta(3) = alpha_n1(i,3)*dtn1 - alpha_n(i,3)*dtn
       deps_theta(4) = alpha_n1(i,4)*dtn1 - alpha_n(i,4)*dtn
       deps_theta(5) = alpha_n1(i,5)*dtn1 - alpha_n(i,5)*dtn
       deps_theta(6) = alpha_n1(i,6)*dtn1 - alpha_n(i,6)*dtn
       if( is_bar_elem ) deps_theta(2:6) = zero
       write(*,9200) i, deps_theta
      end do
      write(*,9000) (i,(deps(i,j),j=1,6),i=1,span)
      return
c
 9000 format(1x,'>> leaving  gp_temp_eps...',
     & /,'element strain increments after subtracting thermal',
     & ' increment...',/
     & (10x,i3,6e17.6) )
 9100 format(1x,'>> entering gp_temp_eps...',
     & /,'element strain increments before subtracting thermal',
     & ' increment...',/
     & (10x,i3,6e17.6) )
 9200 format(10x,i3,6e17.6)
c
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine gp_temp_eps_n1               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/27/2018 rhd              *
c     *                                                              *
c     *     compute total thermal strains at gauss point for         *
c     *     all elements of the block.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gp_temp_eps_n1( span, eps_theta_n1, alpha_n1, 
     &                           gp_temps, gp_rtemps, etype  )
      use main_data,  only : bar_types, link_types
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      integer :: span, etype
      double precision :: eps_theta_n1(span,6), alpha_n1(mxvl,6),
     &                    gp_temps(*), gp_rtemps(*)
c
c                     locals
c
      integer :: i
      double precision :: dtn1
      double precision, parameter :: zero = 0.0d0
      logical :: is_bar_elem, is_link_elem
      logical, parameter :: local_debug = .false.
c
      if( local_debug ) write(*,9100)
c
c            compute (total) thermal strain at n+1
c
      is_bar_elem  = bar_types(etype)
      is_link_elem = link_types(etype)
c
c           instantaneous thermal strain is always just
c           current temperature * current  CTEs
c           adjust for reference temp.
c           \alpha here can be anisotropic
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         dtn1 = gp_temps(i) -  gp_rtemps(i)
         eps_theta_n1(i,1) =  alpha_n1(i,1)*dtn1
         eps_theta_n1(i,2) =  alpha_n1(i,2)*dtn1
         eps_theta_n1(i,3) =  alpha_n1(i,3)*dtn1
         eps_theta_n1(i,4) =  alpha_n1(i,4)*dtn1
         eps_theta_n1(i,5) =  alpha_n1(i,5)*dtn1
         eps_theta_n1(i,6) =  alpha_n1(i,6)*dtn1
      end do
c
      if( is_bar_elem .or. is_link_elem ) 
     &     eps_theta_n1(1:span,2:6) = zero
c
      if ( .not. local_debug ) return
      write(*,*) '>> thermal strains @ n+1...'
      do i = 1, span
       write(*,9200) i, eps_theta_n1(i,1:6)
      end do
      return
c
 9100 format(1x,'>> entering gp_temp_eps_n1...' )
 9200 format(10x,i3,6e17.6)
c
      end






















