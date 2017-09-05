c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gpifv1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/29/2017 rhd             *          
c     *                                                              *          
c     *     computes the internal resisting force                    *          
c     *     vectors for a block of similar elements in uniform       *          
c           global coordinates for a integration point               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gpifv1( eleifv, nrow_ifv, gpn, weight, dj, b,                  
     &                   urcs_blk_n1, local_work, element_volumes )             
c                                                                               
      implicit none                                                   
      include 'param_def'                                                       
      include 'include_sig_up'                                                  
c                                                                               
c                       parameter declarations                                  
c             
      integer :: nrow_ifv ! same as span
      integer :: gpn                                                                  
      double precision :: eleifv(nrow_ifv,*), weight, dj(*), 
     &                    element_volumes(*), b(mxvl,mxedof,*), 
     &                    urcs_blk_n1(mxvl,*)                                   
c                                                                               
c                       local declarations                                      
c      
      integer :: i, j, span, felem, type, nnode, totdof, iter                                                                         
      logical :: geonl, bbar                                                       
      double precision :: qtn1(mxvl,nstr,nstr), cs_blk_n1(mxvl,nstr),                            
     &                    eps_bbar, w, scalar                                                    
c                                                                               
      span            = local_work%span                                         
      felem           = local_work%felem                                        
      type            = local_work%elem_type                                    
      nnode           = local_work%num_enodes                                   
      geonl           = local_work%geo_non_flg                                  
      bbar            = local_work%bbar_flg                                     
      eps_bbar        = local_work%eps_bbar                                     
      totdof          = local_work%totdof                                       
      iter            = local_work%iter                                         
c                                                                               
c                       if element is triangle, wedge, tet sacle                
c                       weights to correctly integrate B*sigma.                 
c                                                                               
      scalar = 1.0d00                                                           
      if ( local_work%adjust_const_elem  ) then                                 
        call adjust_scalar_weights( type, scalar )                              
      end if                                                                    
c                                                                               
      w = weight*scalar                                                         
c                                                                               
c                       compute the strain-displacement matrices                
c                       for the gauss point. for large                          
c                       displacement analysis we are computing [B]              
c                       linear using n+1 configuration and                      
c                       multiplying into Cauchy stresses.                       
c                       the cohesive materials have only 3 tractions            
c                       located in 1st 3 cols of urcs..                         
c                                                                               
      if ( local_work%is_cohes_elem ) then                                      
        call blcmp_cohes( span, b, local_work%cohes_rot_block,                  
     &                    local_work%shape(1,gpn), type, nnode )                
c                                                                               
          do j = 1, totdof                                                      
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
            do i = 1, span                                                      
               eleifv(i,j) = eleifv(i,j) +                                      
     &                       b(i,j,1) * urcs_blk_n1(i,1) * w * dj(i) +          
     &                       b(i,j,2) * urcs_blk_n1(i,2) * w * dj(i) +          
     &                       b(i,j,3) * urcs_blk_n1(i,3) * w * dj(i)            
            end do                                                              
          end do                                                                
          return                                                                
      end if                                                                    
c                                                                               
c                       solid elements                                          
c                                                                               
      call blcmp1( span, b,                                                     
     &               local_work%gama(1,1,1,gpn),                                
     &               local_work%nxi(1,gpn), local_work%neta(1,gpn),             
     &               local_work%nzeta(1,gpn), nnode )                           
c                                                                               
      if ( bbar .and. type .eq. 2 ) then                                        
            call bmod ( b, local_work%vol_block,                                
     &                  span, mxvl, eps_bbar, mxedof )                          
      end if                                                                    
c                                                                               
c                       compute gauss point contributions to                    
c                       the internal force vectors for the element              
c                       block. treat [B] as full even though it is              
c                       not. b-bar fills upper part. eliminating                
c                       zero multiplications makes code unreadable.             
c                       for geonl, convert urcs -> cauchy stress at             
c                       n+1 using [R,n+1] from polar decompositions             
c                       during strain update.                                   
c                                                                               
       if ( geonl ) then                                                        
          call getrm1( span, qtn1, local_work%rot_blk_n1(1,1,gpn), 2 )          
          call qmply1( span, mxvl, nstr, qtn1,                                  
     &                 local_work%urcs_blk_n1(1,1,gpn),                         
     &                 cs_blk_n1 )                                              
          do j = 1, totdof                                                      
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
            do i = 1, span                                                      
               eleifv(i,j) = eleifv(i,j) +                                      
     &                       b(i,j,1)*cs_blk_n1(i,1)*w*dj(i) +                  
     &                       b(i,j,2)*cs_blk_n1(i,2)*w*dj(i) +                  
     &                       b(i,j,3)*cs_blk_n1(i,3)*w*dj(i) +                  
     &                       b(i,j,4)*cs_blk_n1(i,4)*w*dj(i) +                  
     &                       b(i,j,5)*cs_blk_n1(i,5)*w*dj(i) +                  
     &                       b(i,j,6)*cs_blk_n1(i,6)*w*dj(i)                    
            end do                                                              
          end do                                                                
          go to 9000                                                            
      end if                                                                    
c                                                                               
c                       small displacement formulation. the urcs                
c                       are stresses. no transformation needed.                 
c                                                                               
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
        do i = 1, span                                                          
           eleifv(i,j) = eleifv(i,j) +                                          
     &                   b(i,j,1)*urcs_blk_n1(i,1)*w*dj(i) +                    
     &                   b(i,j,2)*urcs_blk_n1(i,2)*w*dj(i) +                    
     &                   b(i,j,3)*urcs_blk_n1(i,3)*w*dj(i) +                    
     &                   b(i,j,4)*urcs_blk_n1(i,4)*w*dj(i) +                    
     &                   b(i,j,5)*urcs_blk_n1(i,5)*w*dj(i) +                    
     &                   b(i,j,6)*urcs_blk_n1(i,6)*w*dj(i)                      
        end do                                                                  
      end do                                                                    
      go to 9000                                                                
c                                                                               
c                       update gauss point contribution to volume               
c                       of element.                                             
c                                                                               
 9000 continue                                                                  
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span                                                            
        element_volumes(i) = element_volumes(i)   +   w * dj(i)                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gpifv3                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/29/2017 rhd             *          
c     *                                                              *          
c     *     computes the internal resisting force                    *          
c     *     vectors for a block of bar elements in global            *
c     *     coordinates                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gpifv3( span, mxvl, felem, iout, etype, geonl, 
     &                   urcs_blk_n1, eleifv, ce, areas_0, areas_n1 )             
c                                                                               
      implicit none                                                    
c                                                                               
c                       parameter declarations                                  
c                
      integer :: span, mxvl, felem, iout, etype
      logical :: geonl  
      double precision ::  eleifv(span,*), urcs_blk_n1(mxvl,*),
     &                     ce(mxvl,*), areas_0(*), areas_n1(*)                                   
c                                                                               
c                       local declarations                                      
c 
      integer :: i
      logical, parameter :: local_debug = .false.
      double precision :: x1, x2, y1, y2, z1, z2, dx, dy, dz, len,
     &                    cos_l, cos_m, cos_n, area, force 
c
c                        compute force at 2 end of bar then global
c                        components at each end.
c 
c                        ce = ce_0 for small displacements
c                        ce = ce_n1 for large displacements
c
c                        urcs.. is engineering stress or Cauchy stress   
c 
c                        force ordering: Px1, Px2, Py1, Py2, Pz1, Pz2   
c      
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span   
        x1 = ce(i,1)
        x2 = ce(i,2)
        y1 = ce(i,3)
        y2 = ce(i,4)
        z1 = ce(i,5)
        z2 = ce(i,6)
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        len = sqrt( dx*dx + dy*dy + dz*dz )
        cos_l = dx / len
        cos_m = dy / len
        cos_n = dz / len
        area = areas_0(i)
        if( geonl ) area = areas_n1(i)
        force = urcs_blk_n1(i,1) * area
        eleifv(i,1) = -force * cos_l
        eleifv(i,2) =  force * cos_l
        eleifv(i,3) = -force * cos_m
        eleifv(i,4) =  force * cos_m
        eleifv(i,5) = -force * cos_n
        eleifv(i,6) =  force * cos_n
      end do  
c
      if( .not. local_debug ) return
      write(iout,9000) felem
      do i = 1, span
        x1 = ce(i,1)
        x2 = ce(i,2)
        y1 = ce(i,3)
        y2 = ce(i,4)
        z1 = ce(i,5)
        z2 = ce(i,6)
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        len = sqrt( dx*dx + dy*dy + dz*dz )
        cos_l = dx / len
        cos_m = dy / len
        cos_n = dz / len
        area = areas_0(i)
        if( geonl ) area = areas_n1(i)
        force = urcs_blk_n1(i,1) * area
        write(iout,9010) i + felem - 1, len, cos_l, cos_m, cos_n
        write(iout,9020) area, force
        write(iout,9030) eleifv(i,1:6)
      end do                
c   
      return
 9000 format(3x,'... internal bar forces for blk starting @ elem: ',i8)  
 9010 format(8x,'elem, len, l, m, n: ',i8,4f15.6)
 9020 format(8x,'area, force: ',2f15.6)
 9030 format(8x,'node 1, 2 ifv: ',3f15.5,/,
     &       8x,'               ',3f15.5 )    
      end                                                                            
                                                                                
                                                                                
