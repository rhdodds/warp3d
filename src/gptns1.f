c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gptns1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *               last modified : 9/19/2016 rhd                  *          
c     *                                                              *          
c     *     computes the contributon to the tangent                  *          
c     *     stiffnes matrices for a block of similar elements in     *          
c     *     uniform global coordinates for a single gauss point.     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gptns1( cp, icp, gpn, props, iprops, glb_ek_blk,               
     &                   local_work )                                           
      use main_data, only: asymmetric_assembly                                  
      use elem_block_data, only : global_cep_blocks => cep_blocks               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn                                                            
      real props(mxelpr,*)                                                      
      integer cp(*), icp(mxutsz,*), iprops(mxelpr,*)                            
c                                                                               
      double precision :: glb_ek_blk(*)                                         
c                                                                               
c                     locals                                                    
c                                                                               
      integer :: etype, span, felem, utsz, nnode, totdof, mat_type,             
     &           iter, int_order, local_iout, nrow_ek                           
      double precision ::                                                       
     & eps_bbar, weight, rad(mxvl), dummy, factors(mxvl), one                   
      logical include_qbar, geonl, bbar, first, qbar_flag,                      
     &        temps_to_process, iscp, symmetric_assembly                        
      data one  / 1.0d00 /                                                      
c                                                                               
c                       set local versions of the data structure                
c                       scalars. set logical to include/not include the         
c                       so-called qbar modification of the                      
c                       constitutive matrix. see addtional                      
c                       comments in routine ctran1                              
c                                                                               
      etype            = local_work%elem_type                                   
      span             = local_work%span                                        
      felem            = local_work%felem                                       
      utsz             = local_work%utsz                                        
      eps_bbar         = local_work%eps_bbar                                    
      nnode            = local_work%num_enodes                                  
      totdof           = local_work%totdof                                      
      geonl            = local_work%geo_non_flg                                 
      utsz             = local_work%utsz                                        
      bbar             = local_work%bbar_flg                                    
      mat_type         = local_work%mat_type                                    
      weight           = local_work%weights(gpn)                                
      first            = local_work%first                                       
      iter             = local_work%iter                                        
      qbar_flag        = local_work%qbar_flag                                   
      int_order        = local_work%int_order                                   
      temps_to_process = local_work%temps_node_to_process                       
      include_qbar     = qbar_flag                                              
      symmetric_assembly = .not. asymmetric_assembly                            
c                                                                               
c                      for axisymmetric elements, compute the radius to         
c                      the current gauss point for each element in the          
c                      block.                                                   
c                                                                               
      if( local_work%is_axisymm_elem ) then                                     
        call get_radius( rad, nnode, span, local_work%shape(1,gpn),             
     &                   local_work%ce, mxvl )                                  
      end if                                                                    
c                                                                               
c                     compute [B] matrices at this gauss                        
c                     point for all elements in block. for geonl,               
c                      a) [B] is the linear form evaluated using                
c                         the updated coordinates at n+1                        
c                      b) the unrotated cauchy stresses at n+1                  
c                         are rotated to define the cauchy stresses             
c                         at n+1 for initial stress stiffness                   
c                         and stress modification of [D].                       
c                         the [R,n+1] is retrieved from block storage           
c                         by getrm1 and used to form the 6x6                    
c                         transformation: {T} = [qn1] * {uc}.                   
c                         [qn1] is also used to rotated [D] matrices            
c                         from unrotated to spatial configuration.              
c                         [D*] = trans[ qn1 ] * [D] * [qn1 ]                    
c                     routine bmod applies the b-bar modifications              
c                     for both small and geonl formulations.                    
c                     note: [D] and [cep] mean the same here.                   
c                     the ctran1 routine performs the [D] rotation.             
c                                                                               
      if( local_work%is_cohes_elem ) then                                       
           call blcmp_cohes( span, local_work%b_block,                          
     &                       local_work%cohes_rot_block,                        
     &                       local_work%shape(1,gpn), etype, nnode )            
      else                                                                      
           call gptns1_a  ! set up solid elements                               
      end if                                                                    
c                                                                               
c                 branch on material type:                                      
c                      1 = simple mises model- linear hardening                 
c                          (isotropic or mixed kinematic), has geonl            
c                          option                                               
c                          rate effects on flow properties. temperature         
c                          dependent flow properties, modulus, nu,              
c                          alpha.                                               
c                      2 = nonlinear elastic model (deformation                 
c                          plasticity                                           
c                          with linear + power-law stress-strain                
c                          relation. no geonl option.                           
c                      3 = general mises and gurson model -                     
c                          rate dependent/independent, linear (iso)             
c                          hardening or power-law hardening. gurson             
c                          model with w/o nucleation (strain                    
c                          controlled).                                         
c                          temperature dependent flow props, modulus,           
c                          nu, alpha                                            
c                      4 = interface constitutive models                        
c                          supports several cohesive zone models.               
c                          no geometric stiffness matrix. set geonl             
c                          false to bypass those additional                     
c                          computations                                         
c                      5 = adv. cyclic plasticity model                         
c                      6 = creep                                                
c                      7 = mises + hydrogen                                     
c                      8 = Abaqus UMAT                                          
c                     10 = CP model                                             
c                                                                               
c                    Always see funtion: warp3d_matl_num                        
c                    for the latets updates.                                    
c                                                                               
c                    [Dt] computed in unrotated configuration.                  
c                                                                               
      local_iout = local_work%iout                                              
      select case ( mat_type )                                                  
      case ( 1 )                                                                
        call drive_01_cnst( gpn, local_iout, local_work )                       
      case ( 2 )                                                                
        call drive_02_cnst( gpn, local_iout, local_work )                       
      case ( 3 )                                                                
        call drive_03_cnst( gpn, local_iout, local_work )                       
      case ( 4 )                                                                
        call drive_04_cnst( gpn, local_iout, local_work )                       
        geonl = .false.                                                         
      case ( 5 )                                                                
        call drive_05_cnst( gpn, local_iout, local_work )                       
      case ( 6 )                                                                
        call drive_06_cnst( gpn, local_iout, local_work )                       
      case ( 7 )                                                                
        call drive_07_cnst( gpn, local_iout, local_work )                       
      case ( 8 )                                                                
        call drive_umat_cnst( gpn, local_iout, local_work )                     
      case (10 )                                                                
        call drive_10_cnst( gpn, local_iout, local_work )                       
      case default                                                              
          write(local_iout,9500)                                                
          call die_abort                                                        
      end select                                                                
c                                                                               
c                       Check to see if we're actually a damaged                
c                       interface model.  If we are, apply damage to            
c                       the tangent.                                            
c                                                                               
c                       commented out until Mark M wants to purse this          
c                       effort                                                  
c                                                                               
c      if( local_work%is_inter_dmg ) then                                       
c        call drive_11_cnst(gpn, local_iout, local_work)                        
c      end if                                                                   
c                                                                               
c                       convert [Dt] from unrotated cauchy to cauchy            
c                       at current deformed configuration for geometric         
c                       nonlinear analysis. no computations                     
c                       for cohesive or deformation plasticity.                 
c                       for UMAT with hyperelastic formulations which           
c                       use [F] to get strains, the                             
c                       [Dt] stored in WARP3D is really for Cauchy              
c                       stress - not unrotated Cauchy stress. The               
c                       code below skips the rotation but may include           
c                       the [Q] modification as requested in user               
c                       input.                                                  
c                                                                               
      if( geonl .and. local_work%is_deform_plas ) then                          
            write(local_iout,9000)                                              
            call die_abort                                                      
            stop                                                                
      end if                                                                    
      if( geonl .and. local_work%is_solid_matl )                                
     &  call ctran1( span, local_work%cep, local_work%qn1,                      
     &               local_work%cs_blk_n1,                                      
     &               include_qbar, local_work%det_jac_block(1,gpn),             
     &               weight, local_work%is_umat,                                
     &               local_work%umat_stress_type,                               
     &               local_work%is_crys_pls )                                   
c                                                                               
c                      include other required scalars with the material         
c                      matrix for axisymmetric (triangle or                     
c                      quadrillateral) or planar triangle elements.             
c                      A scalar to adjust the jacobian determinant              
c                      is also needed for tetrahedral elements.                 
c                                                                               
      if( local_work%adjust_const_elem ) then                                   
        call adjust_cnst( etype, nstr, mxvl, span, rad,                         
     &                    local_work%cep )                                      
      end if                                                                    
c                                                                               
c                     compute each part of the element tangent                  
c                     stiffness matrix and add it in to the                     
c                     total. for geonl, we include initial stress               
c                     stiffness. separate code (unrolled) is used for           
c                     the 8-node brick for speed.                               
c                                                                               
c                                                                               
      if( asymmetric_assembly ) then                                            
        call bdbt_asym( span, local_work%b_block,                               
     &                local_work%bd_block,                                      
     &                local_work%cep, local_work%ek_full, mxvl,                 
     &                mxedof, totdof*totdof, totdof )                           
      else   ! symmetric assembly                                               
         call bdbtgen( span, icp, local_work%b_block,                           
     &                 local_work%bd_block,                                     
     &                 local_work%cep, local_work%ek_symm, mxvl,                
     &                 mxedof, utsz, nstr, totdof, mxutsz )                     
      end if                                                                    
c                                                                               
      if( geonl ) then  ! add tans([G]) M [G] to [Ke]                           
        call kgstiff( span, cp, icp, local_work%gama_block(1,1,1,gpn),          
     &            local_work%nxi(1,gpn), local_work%neta(1,gpn),                
     &            local_work%nzeta(1,gpn), nnode,                               
     &            local_work%cs_blk_n1,                                         
     &            local_work%det_jac_block(1,gpn),                              
     &            weight, local_work%ek_full, local_work%ek_symm,               
     &            local_work%vol_block, bbar, totdof )                          
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('FATAL ERROR: the nonlinear elastic material model',               
     &   /,  '             is not valid for use with elements that',            
     &   /,  '             have geometric nonlinearity.',                       
     &   /,  '             job terminated.' )                                   
 9500 format(1x,'>> Fatal Error: gptns1. invalid material type..',              
     &    /, 1x,'                job terminated' )                              
c                                                                               
c                                                                               
      contains                                                                  
c     ========                                                                  
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gptns1_a                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/21/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gptns1_a                                                       
c                                                                               
      logical :: axisymm                                                        
c                                                                               
      axisymm = etype .eq. 10  .or.  etype .eq. 11                              
c                                                                               
      if( geonl ) then                                                          
          call getrm1( span, local_work%qn1,                                    
     &                 local_work%rot_blk_n1(1,1,gpn), 2 )                      
          call qmply1( span, mxvl, nstr, local_work%qn1,                        
     &                 local_work%urcs_blk_n1(1,1,gpn),                         
     &                 local_work%cs_blk_n1 )                                   
      end if                                                                    
c                                                                               
      if( axisymm )                                                             
     &   call blcmp1_axisymm( span, local_work%b_block,                         
     &                 local_work%gama_block(1,1,1,gpn),                        
     &                 local_work%nxi(1,gpn), local_work%neta(1,gpn),           
     &                 local_work%nzeta(1,gpn),                                 
     &                 local_work%shape(1,gpn),                                 
     &                 local_work%ce, rad, etype, nnode )                       
      if( .not. axisymm )                                                       
     &   call blcmp1(  span, local_work%b_block,                                
     &                 local_work%gama_block(1,1,1,gpn),                        
     &                 local_work%nxi(1,gpn), local_work%neta(1,gpn),           
     &                 local_work%nzeta(1,gpn),                                 
     &                 nnode )                                                  
c                                                                               
      if( bbar ) then                                                           
          call bmod ( local_work%b_block, local_work%vol_block,                 
     &                span, mxvl, eps_bbar, mxedof )                            
      end if                                                                    
c                                                                               
      return                                                                    
      end subroutine gptns1_a                                                   
      end subroutine gptns1                                                     
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_01_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 1/10/2016 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_01_cnst( gpn, iout, local_work )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, now_blk, ielem, k, felem                                 
      logical :: local_debug                                                    
c                                                                               
      span    = local_work%span                                                 
      now_blk = local_work%blk                                                  
      felem   = local_work%felem                                                
c                                                                               
      local_debug = .false. ! felem .eq. 1 .and. gpn .eq. 3                     
c                                                                               
      call extract_D_symmetric( gpn, local_work )                               
c                                                                               
      if( local_debug ) then                                                    
        write(iout,9000) now_blk, felem, gpn, span                              
        do ielem = 1, span                                                      
          write(iout,9010) felem + ielem - 1                                    
          do k = 1, 6                                                           
            write(iout,9020) local_work%cep(ielem,k,1:6)                        
          end do                                                                
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(1x,'.... debug cnst1. now_blk, felem, gpn, span: ', 4i8)           
 9010 format(10x,'[D] for element: ',i7)                                        
 9020 format(15x,6es14.6)                                                       
                                                                                
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_02_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  1/10/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_02_cnst( gpn, iout, local_work )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, now_blk, ielem, k, felem                                 
      logical :: local_debug                                                    
c                                                                               
      span    = local_work%span                                                 
      now_blk = local_work%blk                                                  
      felem   = local_work%felem                                                
c                                                                               
      local_debug = .false. ! felem .eq. 1 .and. gpn .eq. 3                     
c                                                                               
      call extract_D_symmetric( gpn, local_work )                               
c                                                                               
      if( local_debug ) then                                                    
        write(iout,9000) now_blk, felem, gpn, span                              
        do ielem = 1, span                                                      
          write(iout,9010) felem + ielem - 1                                    
          do k = 1, 6                                                           
            write(iout,9020) local_work%cep(ielem,k,1:6)                        
          end do                                                                
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(1x,'.... debug cnst2. now_blk, felem, gpn, span: ', 4i8)           
 9010 format(10x,'[D] for element: ',i7)                                        
 9020 format(15x,6es14.6)                                                       
                                                                                
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_03_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  1/10/2016 rhd             *          
c     *                                                              *          
c     *                      mises and gurson                        *          
c     ****************************************************************          
c                                                                               
      subroutine drive_03_cnst( gpn, iout, local_work )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
      call extract_D_symmetric( gpn, local_work )                               
      return                                                                    
c                                                                               
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine drive_04_cnst                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 10/26/2015 rhd             *          
c     *                                                              *          
c     *     drive [D] consistent update for cohesive model           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_04_cnst( gpn, iout, local_work )                         
c                                                                               
      use main_data, only : matprp, lmtprp                                      
      use elem_block_data, only : gbl_cep_blocks => cep_blocks                  
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      double precision :: weight, symm_part_cep(6), f                           
      logical :: ldebug                                                         
      integer :: span, felem, now_blk, ielem, k, i                              
c                                                                               
      ldebug  = .false.                                                         
      span    = local_work%span                                                 
      felem   = local_work%felem                                                
      weight  = local_work%weights(gpn)                                         
      now_blk = local_work%blk                                                  
c                                                                               
c              pull 6 terms of lower-triangle for this element                  
c              global cep block is 6 x span x num integration                   
c              points                                                           
c                                                                               
c              expand to 3x3 symmetric [D] scaled by gpn                        
c              weight factor                                                    
c                                                                               
!DIR$ IVDEP                                                                     
      do ielem = 1, span                                                        
       f = weight * local_work%det_jac_block(ielem,gpn)                         
       k = ( 6 * span * (gpn-1) ) + 6 * (ielem-1)                               
       symm_part_cep(1) = f * gbl_cep_blocks(now_blk)%vector(k+1)               
       symm_part_cep(2) = f * gbl_cep_blocks(now_blk)%vector(k+2)               
       symm_part_cep(3) = f * gbl_cep_blocks(now_blk)%vector(k+3)               
       symm_part_cep(4) = f * gbl_cep_blocks(now_blk)%vector(k+4)               
       symm_part_cep(5) = f * gbl_cep_blocks(now_blk)%vector(k+5)               
       symm_part_cep(6) = f * gbl_cep_blocks(now_blk)%vector(k+6)               
       local_work%cep(ielem,1,1) = symm_part_cep(1)                             
       local_work%cep(ielem,2,1) = symm_part_cep(2)                             
       local_work%cep(ielem,2,2) = symm_part_cep(3)                             
       local_work%cep(ielem,3,1) = symm_part_cep(4)                             
       local_work%cep(ielem,3,2) = symm_part_cep(5)                             
       local_work%cep(ielem,3,3) = symm_part_cep(6)                             
       local_work%cep(ielem,1,2) = symm_part_cep(2)                             
       local_work%cep(ielem,1,3) = symm_part_cep(4)                             
       local_work%cep(ielem,2,3) = symm_part_cep(5)                             
      end do                                                                    
c                                                                               
      if( ldebug ) then                                                         
         write(iout,9000)                                                       
         do i = 1, span                                                         
           f = 1.0d00 / (weight * local_work%det_jac_block(i,gpn) )             
           write(iout,9100) felem+i-1, gpn                                      
           write(iout,9110) local_work%cep(i,1:3,1:3) * f                       
         end do                                                                 
      end if                                                                    
                                                                                
      return                                                                    
c                                                                               
9000  format('... drive_04_cnst. returned [D]s')                                
9100  format(10x,'...element, gpn: ',i7,i3)                                     
9110  format(3(15x,3f10.3,/))                                                   
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_05_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  1/10/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_05_cnst( gpn, iout, local_work )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, now_blk, ielem, k, felem                                 
      logical :: local_debug                                                    
c                                                                               
      span    = local_work%span                                                 
      now_blk = local_work%blk                                                  
      felem   = local_work%felem                                                
c                                                                               
      local_debug = .false. ! felem .eq. 1 .and. gpn .eq. 3                     
c                                                                               
      call extract_D_symmetric( gpn, local_work )                               
c                                                                               
      if( local_debug ) then                                                    
        write(iout,9000) now_blk, felem, gpn, span                              
        do ielem = 1, span                                                      
          write(iout,9010) felem + ielem - 1                                    
          do k = 1, 6                                                           
            write(iout,9020) local_work%cep(ielem,k,1:6)                        
          end do                                                                
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(1x,'.... debug cnst5. now_blk, felem, gpn, span: ', 4i8)           
 9010 format(10x,'[D] for element: ',i7)                                        
 9020 format(15x,6es14.6)                                                       
                                                                                
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_06_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  1/10/2016 rhd             *          
c     *                                                              *          
c     *     drive [D] consistent update for creep                               
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_06_cnst( gpn, iout, local_work )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
      call extract_D_symmetric( gpn, local_work )                               
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_07_cnst                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  1/10/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_07_cnst( gpn, iout, local_work )                         
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, now_blk, ielem, k, felem                                 
      logical :: local_debug                                                    
c                                                                               
      span    = local_work%span                                                 
      now_blk = local_work%blk                                                  
      felem   = local_work%felem                                                
c                                                                               
      local_debug = .false. ! felem .eq. 1 .and. gpn .eq.                       
c                                                                               
      call extract_D_symmetric( gpn, local_work )                               
c                                                                               
      if( local_debug ) then                                                    
        write(iout,9000) now_blk, felem, gpn, span                              
        do ielem = 1, span                                                      
          write(iout,9010) felem + ielem - 1                                    
          do k = 1, 6                                                           
            write(iout,9020) local_work%cep(ielem,k,1:6)                        
          end do                                                                
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(1x,'.... debug cnst7. now_blk, felem, gpn, span: ', 4i8)           
 9010 format(10x,'[D] for element: ',i7)                                        
 9020 format(15x,6es14.6)                                                       
                                                                                
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine extract_D_symmetric                  *          
c     *           - service routine. should be inlined -             *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified :  1/10/2016 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine extract_D_symmetric( gpn, local_work )                         
c                                                                               
      use elem_block_data, only : gbl_cep_blocks => cep_blocks                  
      implicit none                                                             
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      double precision :: weight, symm_part_cep(mxvl,21), f                     
      integer :: span, now_blk, ielem, sloc, k, felem                           
c                                                                               
      span    = local_work%span                                                 
      weight  = local_work%weights(gpn)                                         
      now_blk = local_work%blk                                                  
      felem   = local_work%felem                                                
c                                                                               
c              pull 21 terms of lower-triangle for this element                 
c              global cep block is 21 x span x num integration                  
c              points                                                           
c                                                                               
c              expand to 6 x 6 symmetric [D] and scale by                       
c              integration weight factor                                        
c                                                                               
c      do ielem = 1, span                                                       
c        sloc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)                        
c        f = weight * local_work%det_jac_block(ielem,gpn)                       
c!DIR$ IVDEP                                                                    
c!DIR$ VECTOR ALIGNED                                                           
c        do k = 1, 21                                                           
c         symm_part_cep(ielem,k) = f *                                          
c     &        gbl_cep_blocks(now_blk)%vector(sloc+k)                           
c        end do                                                                 
c      end do                                                                   
                                                                                
      do k = 1, 21                                                              
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
        do ielem = 1, span                                                      
        sloc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)                         
        f = weight * local_work%det_jac_block(ielem,gpn)                        
         symm_part_cep(ielem,k) = f *                                           
     &        gbl_cep_blocks(now_blk)%vector(sloc+k)                            
        end do                                                                  
      end do                                                                    
                                                                                
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
       do ielem = 1, span                                                       
        local_work%cep(ielem,1,1) = symm_part_cep(ielem,1)                      
        local_work%cep(ielem,2,1) = symm_part_cep(ielem,2)                      
        local_work%cep(ielem,2,2) = symm_part_cep(ielem,3)                      
        local_work%cep(ielem,3,1) = symm_part_cep(ielem,4)                      
        local_work%cep(ielem,3,2) = symm_part_cep(ielem,5)                      
        local_work%cep(ielem,3,3) = symm_part_cep(ielem,6)                      
        local_work%cep(ielem,4,1) = symm_part_cep(ielem,7)                      
        local_work%cep(ielem,4,2) = symm_part_cep(ielem,8)                      
        local_work%cep(ielem,4,3) = symm_part_cep(ielem,9)                      
        local_work%cep(ielem,4,4) = symm_part_cep(ielem,10)                     
        local_work%cep(ielem,5,1) = symm_part_cep(ielem,11)                     
        local_work%cep(ielem,5,2) = symm_part_cep(ielem,12)                     
        local_work%cep(ielem,5,3) = symm_part_cep(ielem,13)                     
        local_work%cep(ielem,5,4) = symm_part_cep(ielem,14)                     
        local_work%cep(ielem,5,5) = symm_part_cep(ielem,15)                     
        local_work%cep(ielem,6,1) = symm_part_cep(ielem,16)                     
        local_work%cep(ielem,6,2) = symm_part_cep(ielem,17)                     
        local_work%cep(ielem,6,3) = symm_part_cep(ielem,18)                     
        local_work%cep(ielem,6,4) = symm_part_cep(ielem,19)                     
        local_work%cep(ielem,6,5) = symm_part_cep(ielem,20)                     
        local_work%cep(ielem,6,6) = symm_part_cep(ielem,21)                     
        local_work%cep(ielem,1,2) = symm_part_cep(ielem,2)                      
        local_work%cep(ielem,1,3) = symm_part_cep(ielem,4)                      
        local_work%cep(ielem,1,4) = symm_part_cep(ielem,7)                      
        local_work%cep(ielem,1,5) = symm_part_cep(ielem,11)                     
        local_work%cep(ielem,1,6) = symm_part_cep(ielem,16)                     
        local_work%cep(ielem,2,3) = symm_part_cep(ielem,5)                      
        local_work%cep(ielem,2,4) = symm_part_cep(ielem,8)                      
        local_work%cep(ielem,2,5) = symm_part_cep(ielem,12)                     
        local_work%cep(ielem,2,6) = symm_part_cep(ielem,17)                     
        local_work%cep(ielem,3,4) = symm_part_cep(ielem,9)                      
        local_work%cep(ielem,3,5) = symm_part_cep(ielem,13)                     
        local_work%cep(ielem,3,6) = symm_part_cep(ielem,18)                     
        local_work%cep(ielem,4,5) = symm_part_cep(ielem,14)                     
        local_work%cep(ielem,4,6) = symm_part_cep(ielem,19)                     
        local_work%cep(ielem,5,6) = symm_part_cep(ielem,20)                     
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine drive_umat_cnst                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 10/1/2015 rhd              *          
c     *                                                              *          
c     *          drive [D] consistent update for warp3d umat         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_umat_cnst( gpn, iout, local_work )                       
c                                                                               
      use elem_block_data, only : gbl_cep_blocks => cep_blocks                  
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      include 'include_tan_ek'  ! has local_work definition                     
c                                                                               
c                     local variables                                           
c                                                                               
      double precision                                                          
     & weight, symm_part_cep(21), factor                                        
      logical local_debug, debug_now                                            
c                                                                               
c           1. pull a few values from work space for block                      
c                                                                               
      span              = local_work%span                                       
      felem             = local_work%felem                                      
      weight            = local_work%weights(gpn)                               
      now_blk           = local_work%blk                                        
      local_debug       =  .false.                                              
c                                                                               
c           2. the tangent [D] matrices are stored in the                       
c              global_cep_blocks. only lower symmetric terms                    
c              stored in a vector. stored there in rstgp1.                      
c                                                                               
c              for this integration point (gpn), process all                    
c              elements in block.                                               
c                                                                               
      do ielem = 1, span                                                        
c                                                                               
        noel = felem + ielem - 1                                                
        debug_now = local_debug                                                 
        if( debug_now ) write(iout,9100) ielem, noel, gpn                       
c                                                                               
c                2 (a) pull 21 terms of lower-triangle for this element         
c                      global cep block is 21 x span x num integration          
c                      points                                                   
c                                                                               
        start_loc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)                    
!DIR$ IVDEP                                                                     
        do k = 1, 21                                                            
          symm_part_cep(k) = gbl_cep_blocks(now_blk)%vector(start_loc+k)        
        end do                                                                  
c                                                                               
c                2 (b) expand to 6 x 6 symmetric [D] and scale by               
c                      integration weight factor                                
c                                                                               
        factor = weight * local_work%det_jac_block(ielem,gpn)                   
        k = 1                                                                   
        do i = 1, 6                                                             
!DIR$ IVDEP                                                                     
         do j = 1, i                                                            
           local_work%cep(ielem,i,j) = symm_part_cep(k) * factor                
           local_work%cep(ielem,j,i) = symm_part_cep(k) * factor                
           k = k + 1                                                            
         end do                                                                 
        end do                                                                  
c                                                                               
        if( debug_now ) write(iout,9110) symm_part_cep(1:4)                     
      end do                                                                    
c                                                                               
      if( debug_now ) write(iout,9900)                                          
      return                                                                    
c                                                                               
 9000 format('>> Enter UMAT cnst driver...')                                    
 9001 format(5x,'span, felem, gpn, iter: ',i4,i10,i3,i4 )                       
 9002 format(5x,'integration weight: ',f10.5)                                   
 9006 format(5x,'num history terms: ',i4 )                                      
 9100 format(5x,'... processing i, elem, gpn: ',i4,i10, i3)                     
 9110 format(5x,'symd 1-4: ',4e14.6)                                            
 9900 format('>> Leave UMAT cnst driver...')                                    
c                                                                               
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_10_cnst                *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 12/8/2015 rhd              *          
c     *                                                              *          
c     *              drive [D] consistent update for CP model        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_10_cnst( gpn, iout, local_work )                         
c                                                                               
      use main_data, only : matprp, lmtprp                                      
      use elem_block_data, only : gbl_cep_blocks => cep_blocks                  
      use mm10_defs, only : indexes_common, index_crys_hist                     
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                     parameter declarations                                    
c                                                                               
      include 'include_tan_ek'                                                  
      integer :: gpn, iout                                                      
c                                                                               
c                     local variables                                           
c                                                                               
      integer :: span, felem, iter, now_blk,                                    
     &           start_loc, k, i, j, eh, sh                                     
      double precision ::                                                       
     & weight, f, cep(6,6), cep_vec(36), tol                                    
      logical :: local_debug                                                    
      equivalence ( cep, cep_vec )                                              
c                                                                               
      span             = local_work%span                                        
      felem            = local_work%felem                                       
      weight           = local_work%weights(gpn)                                
      iter             = local_work%iter                                        
      now_blk          = local_work%blk                                         
      sh  = indexes_common(1,1) ! first index of cep tangent                    
      eh  = indexes_common(1,2) ! last index of cep tangent                     
c                                                                               
      local_debug =  .false.                                                    
      if( local_debug ) then                                                    
        write(iout,9100) now_blk, span, felem,                                  
     &                   local_work%hist_size_for_blk                           
      end if                                                                    
c                                                                               
c                     the consistent [D] for each integration point is          
c                     stored in the 1st 36 positions of history                 
c                                                                               
      do i = 1, span                                                            
         f = weight * local_work%det_jac_block(i,gpn)                           
         cep_vec(1:36) = local_work%elem_hist1(i,sh:eh,gpn)                     
         local_work%cep(i,1:6,1:6) = cep(1:6,1:6) * f                           
      end do                                                                    
c                                                                               
c                     code to optionally check symmetry of the                  
c                     [D] linear                                                
c                                                                               
      if( local_debug ) then                                                    
      tol = 0.01d00                                                             
      do i = 1, span                                                            
         cep_vec(1:36) = local_work%elem_hist1(i,sh:eh,gpn)                     
         do j = 1, 6                                                            
            if( cep(j,j) .lt. tol ) then                                        
               write(iout,*) ' .. fatal @ 1'                                    
               call die_abort                                                   
            end if                                                              
        end do                                                                  
        do j = 1, 6                                                             
          do k = 1, 6                                                           
            if( abs( cep(j,k) -cep(k,j) )                                       
     &           .gt. 1.0d-8 ) then                                             
               write(iout,*) ' .. fatal @ 2'                                    
               call die_abort                                                   
            end if                                                              
          end do                                                                
       end do                                                                   
      end do                                                                    
      end if                                                                    
c                                                                               
      if( local_debug ) then                                                    
        write(iout,*) ".... linear elastic [D] for CP"                          
        do k = 1, span                                                          
           write(iout,*) '    ... element: ', felem+k-1                         
           do i = 1, 6                                                          
             write(iout,9000) cep(i,1:6)                                        
           end do                                                               
        end do                                                                  
      end if                                                                    
      return                                                                    
c                                                                               
 9000 format(10x,6e14.5)                                                        
 9100 format(2x,".. debug for routine drive_10_cnst",                           
     & /,10x,"now_blk, span, felem, history size:", 4i6)                        
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine drive_11_cnst                *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 03/11/14                   *          
c     *                                                              *          
c     *     drive damage update to the tangent matrix                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine drive_11_cnst( gpn, iout, local_work )                         
c                                                                               
      use main_data, only : matprp, lmtprp                                      
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      include 'include_tan_ek'                                                  
c                                                                               
c                     parameter declarations                                    
c                                                                               
c                     local variables                                           
c                                                                               
      double precision                                                          
     & weight                                                                   
      logical first                                                             
c                                                                               
      span             = local_work%span                                        
      felem            = local_work%felem                                       
      weight           = local_work%weights(gpn)                                
      first            = local_work%first                                       
      iter             = local_work%iter                                        
c                                                                               
      dmg_loc = 1+local_work%macro_sz                                           
c                                                                               
c                Note: commented out until Mark resumes work on this            
c                                                                               
c      call cnst11( gpn, span, local_work%cep(1:span,1:6,1:6),                  
c     &            local_work%elem_hist(1:span,dmg_loc,gpn),                    
c     &            local_work%sv, local_work%lv, local_work%tv, iout)           
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine bdbt_asym                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                last modified : 5/1/2017 rhd                  *          
c     *                                                              *          
c     *     trans(B) * D * B for the asymmetric case.  Handles all   *          
c     *     element types                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine bdbt_asym( span, b, bd, d, ek_full, mxvl, mxedof,              
     &                      ncol_ek, totdof )                                   
      implicit none                                                             
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, mxvl, mxedof, ncol_ek, nstr, totdof                      
      double precision ::                                                       
     &   b(mxvl,mxedof,6), ek_full(span,ncol_ek), d(mxvl,6,6),                  
     &   bd(mxvl,mxedof,6)                                                      
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i, j, k, row, col                                              
c                                                                               
c              perform multiplication of [D]*[B-mat].                           
c              assumes full [D] and [B]full matrix.                             
c              [D] for solids is 6x6. Also for cohesive,                        
c              i: 4->6 and j: 4->6 should be zeroed by                          
c              cnst.. routine.                                                  
c                                                                               
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
       do i = 1, span                                                           
           bd(i,j,1) = d(i,1,1) * b(i,j,1)                                      
     &               + d(i,2,1) * b(i,j,2)                                      
     &               + d(i,3,1) * b(i,j,3)                                      
     &               + d(i,4,1) * b(i,j,4)                                      
     &               + d(i,5,1) * b(i,j,5)                                      
     &               + d(i,6,1) * b(i,j,6)                                      
c                                                                               
           bd(i,j,2) = d(i,1,2) * b(i,j,1)                                      
     &               + d(i,2,2) * b(i,j,2)                                      
     &               + d(i,3,2) * b(i,j,3)                                      
     &               + d(i,4,2) * b(i,j,4)                                      
     &               + d(i,5,2) * b(i,j,5)                                      
     &               + d(i,6,2) * b(i,j,6)                                      
c                                                                               
           bd(i,j,3) = d(i,1,3) * b(i,j,1)                                      
     &               + d(i,2,3) * b(i,j,2)                                      
     &               + d(i,3,3) * b(i,j,3)                                      
     &               + d(i,4,3) * b(i,j,4)                                      
     &               + d(i,5,3) * b(i,j,5)                                      
     &               + d(i,6,3) * b(i,j,6)                                      
c                                                                               
           bd(i,j,4) = d(i,1,4) * b(i,j,1)                                      
     &               + d(i,2,4) * b(i,j,2)                                      
     &               + d(i,3,4) * b(i,j,3)                                      
     &               + d(i,4,4) * b(i,j,4)                                      
     &               + d(i,5,4) * b(i,j,5)                                      
     &               + d(i,6,4) * b(i,j,6)                                      
c                                                                               
           bd(i,j,5) = d(i,1,5) * b(i,j,1)                                      
     &               + d(i,2,5) * b(i,j,2)                                      
     &               + d(i,3,5) * b(i,j,3)                                      
     &               + d(i,4,5) * b(i,j,4)                                      
     &               + d(i,5,5) * b(i,j,5)                                      
     &               + d(i,6,5) * b(i,j,6)                                      
c                                                                               
           bd(i,j,6) = d(i,1,6) * b(i,j,1)                                      
     &               + d(i,2,6) * b(i,j,2)                                      
     &               + d(i,3,6) * b(i,j,3)                                      
     &               + d(i,4,6) * b(i,j,4)                                      
     &               + d(i,5,6) * b(i,j,5)                                      
     &               + d(i,6,6) * b(i,j,6)                                      
c                                                                               
       end do                                                                   
      end do                                                                    
c                                                                               
c              perform multiplication of tran(B)*([D][B]).                      
c              do for all entries since [D] is not symmetric.                   
c                                                                               
c              the full element stiffness is                                    
c              stored. but the row/column ordering corresponds to               
c              element dof: u1 u2 u3, .. v1 v2 v3, .. w1 w2 w3, ..              
c              not the typical u1 v1 w1, u2 v2 w2, ...                          
c                                                                               
      j = 0                                                                     
      do col = 1, totdof                                                        
       do row = 1, totdof                                                       
             j = j + 1                                                          
!DIR$ IVDEP                                                                     
        do i = 1, span                                                          
         ek_full(i,j) = ek_full(i,j)                                            
     &         +   b(i,col,1) * bd(i,row,1)                                     
     &         +   b(i,col,2) * bd(i,row,2)                                     
     &         +   b(i,col,3) * bd(i,row,3)                                     
     &         +   b(i,col,4) * bd(i,row,4)                                     
     &         +   b(i,col,5) * bd(i,row,5)                                     
     &         +   b(i,col,6) * bd(i,row,6)                                     
        end do                                                                  
       end do                                                                   
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine bdbtgen   (symmetric only)          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/1/2017 rhd               *          
c     *                                                              *          
c     *     this subroutine performs the multiplication of the       *          
c     *     transpose of the strain-displacement matrix by the       *          
c     *     constituitive matrix by the strain-displacement matrix   *          
c     *     that is a building block of a stiffness matrix. this     *          
c     *     routine handles any type of element.                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine bdbtgen( span, icp, b, bd, d, ek_symm,                         
     &                    mxvl, mxedof, utsz, nstr, totdof, mxutsz )            
      implicit none                                                             
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, mxvl, mxedof, utsz, nstr, totdof, mxutsz,                
     &           icp(mxutsz,*)                                                  
      double precision ::                                                       
     &   b(mxvl,mxedof,*), ek_symm(span,utsz), d(mxvl,nstr,*),                  
     &   bd(mxvl,mxedof,*)                                                      
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i, j, k, row, col                                              
c                                                                               
c              perform multiplication of [D]*[B-mat].                           
c              assumes full [D] and [B]full matrix.                             
c              [D] for solids is 6x6. Also for cohesive,                        
c              i: 4->6 and j: 4->6 should be zeroed by                          
c              cnst.. routine.                                                  
c                                                                               
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
       do i = 1, span                                                           
           bd(i,j,1) = d(i,1,1) * b(i,j,1)                                      
     &               + d(i,2,1) * b(i,j,2)                                      
     &               + d(i,3,1) * b(i,j,3)                                      
     &               + d(i,4,1) * b(i,j,4)                                      
     &               + d(i,5,1) * b(i,j,5)                                      
     &               + d(i,6,1) * b(i,j,6)                                      
c                                                                               
           bd(i,j,2) = d(i,1,2) * b(i,j,1)                                      
     &               + d(i,2,2) * b(i,j,2)                                      
     &               + d(i,3,2) * b(i,j,3)                                      
     &               + d(i,4,2) * b(i,j,4)                                      
     &               + d(i,5,2) * b(i,j,5)                                      
     &               + d(i,6,2) * b(i,j,6)                                      
                                                                                
c                                                                               
           bd(i,j,3) = d(i,1,3) * b(i,j,1)                                      
     &               + d(i,2,3) * b(i,j,2)                                      
     &               + d(i,3,3) * b(i,j,3)                                      
     &               + d(i,4,3) * b(i,j,4)                                      
     &               + d(i,5,3) * b(i,j,5)                                      
     &               + d(i,6,3) * b(i,j,6)                                      
c                                                                               
           bd(i,j,4) = d(i,1,4) * b(i,j,1)                                      
     &               + d(i,2,4) * b(i,j,2)                                      
     &               + d(i,3,4) * b(i,j,3)                                      
     &               + d(i,4,4) * b(i,j,4)                                      
     &               + d(i,5,4) * b(i,j,5)                                      
     &               + d(i,6,4) * b(i,j,6)                                      
c                                                                               
           bd(i,j,5) = d(i,1,5) * b(i,j,1)                                      
     &               + d(i,2,5) * b(i,j,2)                                      
     &               + d(i,3,5) * b(i,j,3)                                      
     &               + d(i,4,5) * b(i,j,4)                                      
     &               + d(i,5,5) * b(i,j,5)                                      
     &               + d(i,6,5) * b(i,j,6)                                      
c                                                                               
           bd(i,j,6) = d(i,1,6) * b(i,j,1)                                      
     &               + d(i,2,6) * b(i,j,2)                                      
     &               + d(i,3,6) * b(i,j,3)                                      
     &               + d(i,4,6) * b(i,j,4)                                      
     &               + d(i,5,6) * b(i,j,5)                                      
     &               + d(i,6,6) * b(i,j,6)                                      
c                                                                               
       end do                                                                   
      end do                                                                    
c                                                                               
c              perform multiplication of tran(B)*([D][B]).                      
c              do only for upper triangular entries.                            
c                                                                               
      do j = 1, utsz                                                            
        row = icp(j,1)                                                          
        col = icp(j,2)                                                          
!DIR$ IVDEP                                                                     
        do i = 1, span                                                          
         ek_symm(i,j) = ek_symm(i,j)                                            
     &         +   b(i,col,1) * bd(i,row,1)                                     
     &         +   b(i,col,2) * bd(i,row,2)                                     
     &         +   b(i,col,3) * bd(i,row,3)                                     
     &         +   b(i,col,4) * bd(i,row,4)                                     
     &         +   b(i,col,5) * bd(i,row,5)                                     
     &         +   b(i,col,6) * bd(i,row,6)                                     
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ctran1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 4/18/2016 rhd              *          
c     *                                                              *          
c     *     transform [Dt] from a form relating the unrotated stress *          
c     *     rate and the unrotated rate of deformation tensor to     *          
c     *     one relating the cauchy stress rate and the rate of      *          
c     *     (spatial) deformation strain for a block of elements.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ctran1( span, cep, qn1, cs, qbar, dj, w, is_umat,              
     &                   umat_stress_type, is_crys_pls  )                       
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision ::                                                       
     &     cep(mxvl,nstr,*), qn1(mxvl,nstr,*), tc(mxvl,nstr,nstr),              
     &     cs(mxvl,*), half, two, dj(*), w, wf, halfw                           
      logical :: qbar, is_umat, is_crys_pls, do_transform                       
      data half, two / 0.5d00, 2.0d00 /                                         
c                                                                               
c             [cep] (mxvl x 6 x 6) relates increments                           
c             of unrotated cauchy stress to increments                          
c             of the unrotated deformation. transform [cep]                     
c             so it relates increments of cauchy stress to                      
c             increments of the deformation, both on the                        
c             spatial coordinates.                                              
c                                                                               
c             [cep*] = [qn1] * [cep] * trans([qn1])                             
c                                                                               
c             [qn1] is a rotation matrix constructed from the                   
c             [R] obtained by polar decomposition of the deformation            
c             gradient, [F] =[R][U].                                            
c                                                                               
c             for UMATs the UMAT computed [cep] may already refer to            
c             the Cauchy stress. No rotation to be done.                        
c                                                                               
c             For crystal plasticity model,                                     
c                                                                               
      do_transform = .true.                                                     
      if( is_crys_pls ) do_transform = .true.                                   
      if( is_umat ) then                                                        
        if( umat_stress_type .eq. 1 ) do_transform = .false.                    
      end if                                                                    
      if( do_transform ) then                                                   
c                                                                               
c             perform multiplication of [tc] = [qn1] * [cep]                    
c             use 6 as number of stress components to expose                    
c             value to compiler                                                 
c                                                                               
      do j = 1, 6                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
c                                                                               
            tc(i,j,1)= (qn1(i,j,1)*cep(i,1,1)+                                  
     &                  qn1(i,j,2)*cep(i,2,1)+                                  
     &                  qn1(i,j,3)*cep(i,3,1)+                                  
     &                  qn1(i,j,4)*cep(i,4,1)+                                  
     &                  qn1(i,j,5)*cep(i,5,1)+                                  
     &                  qn1(i,j,6)*cep(i,6,1))                                  
c                                                                               
            tc(i,j,2)= (qn1(i,j,1)*cep(i,1,2)+                                  
     &                  qn1(i,j,2)*cep(i,2,2)+                                  
     &                  qn1(i,j,3)*cep(i,3,2)+                                  
     &                  qn1(i,j,4)*cep(i,4,2)+                                  
     &                  qn1(i,j,5)*cep(i,5,2)+                                  
     &                  qn1(i,j,6)*cep(i,6,2))                                  
c                                                                               
            tc(i,j,3)= (qn1(i,j,1)*cep(i,1,3)+                                  
     &                  qn1(i,j,2)*cep(i,2,3)+                                  
     &                  qn1(i,j,3)*cep(i,3,3)+                                  
     &                  qn1(i,j,4)*cep(i,4,3)+                                  
     &                  qn1(i,j,5)*cep(i,5,3)+                                  
     &                  qn1(i,j,6)*cep(i,6,3))                                  
c                                                                               
            tc(i,j,4)= (qn1(i,j,1)*cep(i,1,4)+                                  
     &                  qn1(i,j,2)*cep(i,2,4)+                                  
     &                  qn1(i,j,3)*cep(i,3,4)+                                  
     &                  qn1(i,j,4)*cep(i,4,4)+                                  
     &                  qn1(i,j,5)*cep(i,5,4)+                                  
     &                  qn1(i,j,6)*cep(i,6,4))                                  
c                                                                               
            tc(i,j,5)= (qn1(i,j,1)*cep(i,1,5)+                                  
     &                  qn1(i,j,2)*cep(i,2,5)+                                  
     &                  qn1(i,j,3)*cep(i,3,5)+                                  
     &                  qn1(i,j,4)*cep(i,4,5)+                                  
     &                  qn1(i,j,5)*cep(i,5,5)+                                  
     &                  qn1(i,j,6)*cep(i,6,5))                                  
c                                                                               
            tc(i,j,6)= (qn1(i,j,1)*cep(i,1,6)+                                  
     &                  qn1(i,j,2)*cep(i,2,6)+                                  
     &                  qn1(i,j,3)*cep(i,3,6)+                                  
     &                  qn1(i,j,4)*cep(i,4,6)+                                  
     &                  qn1(i,j,5)*cep(i,5,6)+                                  
     &                  qn1(i,j,6)*cep(i,6,6))                                  
c                                                                               
         end do                                                                 
      end do                                                                    
c                                                                               
c                       perform multiplication of                               
c                       [cep*] =  [tc] * transpose([qn1])                       
c                                                                               
      do j = 1, 6                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
c                                                                               
            cep(i,j,1)= tc(i,j,1)*qn1(i,1,1)+                                   
     &                  tc(i,j,2)*qn1(i,1,2)+                                   
     &                  tc(i,j,3)*qn1(i,1,3)+                                   
     &                  tc(i,j,4)*qn1(i,1,4)+                                   
     &                  tc(i,j,5)*qn1(i,1,5)+                                   
     &                  tc(i,j,6)*qn1(i,1,6)                                    
c                                                                               
            cep(i,j,2)= tc(i,j,1)*qn1(i,2,1)+                                   
     &                  tc(i,j,2)*qn1(i,2,2)+                                   
     &                  tc(i,j,3)*qn1(i,2,3)+                                   
     &                  tc(i,j,4)*qn1(i,2,4)+                                   
     &                  tc(i,j,5)*qn1(i,2,5)+                                   
     &                  tc(i,j,6)*qn1(i,2,6)                                    
c                                                                               
            cep(i,j,3)= tc(i,j,1)*qn1(i,3,1)+                                   
     &                  tc(i,j,2)*qn1(i,3,2)+                                   
     &                  tc(i,j,3)*qn1(i,3,3)+                                   
     &                  tc(i,j,4)*qn1(i,3,4)+                                   
     &                  tc(i,j,5)*qn1(i,3,5)+                                   
     &                  tc(i,j,6)*qn1(i,3,6)                                    
c                                                                               
            cep(i,j,4)= tc(i,j,1)*qn1(i,4,1)+                                   
     &                  tc(i,j,2)*qn1(i,4,2)+                                   
     &                  tc(i,j,3)*qn1(i,4,3)+                                   
     &                  tc(i,j,4)*qn1(i,4,4)+                                   
     &                  tc(i,j,5)*qn1(i,4,5)+                                   
     &                  tc(i,j,6)*qn1(i,4,6)                                    
c                                                                               
            cep(i,j,5)= tc(i,j,1)*qn1(i,5,1)+                                   
     &                  tc(i,j,2)*qn1(i,5,2)+                                   
     &                  tc(i,j,3)*qn1(i,5,3)+                                   
     &                  tc(i,j,4)*qn1(i,5,4)+                                   
     &                  tc(i,j,5)*qn1(i,5,5)+                                   
     &                  tc(i,j,6)*qn1(i,5,6)                                    
c                                                                               
            cep(i,j,6)= tc(i,j,1)*qn1(i,6,1)+                                   
     &                  tc(i,j,2)*qn1(i,6,2)+                                   
     &                  tc(i,j,3)*qn1(i,6,3)+                                   
     &                  tc(i,j,4)*qn1(i,6,4)+                                   
     &                  tc(i,j,5)*qn1(i,6,5)+                                   
     &                  tc(i,j,6)*qn1(i,6,6)                                    
c                                                                               
         end do                                                                 
      end do                                                                    
      end if ! on do_transform                                                  
c                                                                               
c            subtract the [Q-bar] matrix from the transformed                   
c            [cep]. this is the "initial stress" at the material                
c            point level. this remains an option indicated by qbar.             
c            note: we must multiply in the gauss weight factor and              
c            gauss point det[J] for the subtracted terms. the [cep]             
c            passed in had these factors included by the cnst...                
c            routines. the [Q-bar] 6x6 comes from the tensor                    
c            expression -2 (de.De):s, where, s is the stress tensor,            
c            de is the rate of deformation tensor and De is the virtual         
c            rate of deformation tensor. this expression in matrix form         
c            is: - trans([B]) * [Q-bar] * [B]. this modification of             
c            [cep] is essential for convergence of nearly homogeneous           
c            deformation problems.                                              
c                                                                               
      if( qbar ) then                                                           
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
        do i = 1, span                                                          
         wf    = dj(i) * w                                                      
         halfw = half * wf                                                      
         cep(i,1,1) = cep(i,1,1) - two * cs(i,1) * wf                           
         cep(i,2,2) = cep(i,2,2) - two * cs(i,2) * wf                           
         cep(i,3,3) = cep(i,3,3) - two * cs(i,3) * wf                           
         cep(i,4,1) = cep(i,4,1) - cs(i,4) * wf                                 
         cep(i,6,1) = cep(i,6,1) - cs(i,6) * wf                                 
         cep(i,4,2) = cep(i,4,2) - cs(i,4) * wf                                 
         cep(i,5,2) = cep(i,5,2) - cs(i,5) * wf                                 
         cep(i,5,3) = cep(i,5,3) - cs(i,5) * wf                                 
         cep(i,6,3) = cep(i,6,3) - cs(i,6) * wf                                 
         cep(i,4,4) = cep(i,4,4) - halfw * ( cs(i,1)+cs(i,2) )                  
         cep(i,5,5) = cep(i,5,5) - halfw * ( cs(i,2)+cs(i,3) )                  
         cep(i,6,6) = cep(i,6,6) - halfw * ( cs(i,1)+cs(i,3) )                  
         cep(i,5,4) = cep(i,5,4) - halfw * cs(i,6)                              
         cep(i,6,4) = cep(i,6,4) - halfw * cs(i,5)                              
         cep(i,6,5) = cep(i,6,5) - halfw * cs(i,4)                              
         cep(i,1,4) = cep(i,4,1)                                                
         cep(i,1,6) = cep(i,6,1)                                                
         cep(i,2,4) = cep(i,4,2)                                                
         cep(i,2,5) = cep(i,5,2)                                                
         cep(i,3,5) = cep(i,5,3)                                                
         cep(i,3,6) = cep(i,6,3)                                                
         cep(i,4,5) = cep(i,5,4)                                                
         cep(i,4,6) = cep(i,6,4)                                                
         cep(i,5,6) = cep(i,6,5)                                                
c                                                                               
c                      experiment with symmetrized version of the               
c                      nonsymmetric term. see Crisfield vol. 2, pg55.           
c                                                                               
c      z(1:6,1:6) = 0.0d0                                                       
C      z(1,1:3) = cs(i,1)                                                       
C      z(2,1:3) = cs(i,2)                                                       
C      z(3,1:3) = cs(i,3)                                                       
C      z(4,1:3) = cs(i,4)                                                       
C      z(5,1:3) = cs(i,5)                                                       
C      z(6,1:3) = cs(i,6)                                                       
c      zt = transpose( z)                                                       
c      do k = 1, 6                                                              
c      do l = 1, 6                                                              
c        z(k,l) = 0.5d0*(z(k,l)+zt(l,k))                                        
c        cep(i,k,l) = cep(i,k,l) + z(k,l)*wf                                    
c      end do                                                                   
c      end do                                                                   
c                                                                               
        end do                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *           subroutine bdbt_asym_mcm : not used. saved for     *          
c     *                                      future reference        *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 12/9/13                    *          
c     *  ===>  deprecated. retained for possible future reference    *          
c     *                                                              *          
c     *     B * D * B.T for the asymmetric case.  Handles all        *          
c     *     elements, at least for now.                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine bdbt_asym_mcm( span, icp, b, bt, bd, d, ek,                    
     &                    mxvl, mxedof, utsz, nstr, totdof, mxutsz )            
      implicit integer (a-z)                                                    
c                                                                               
c                       parameter declarations                                  
c                                                                               
      double precision                                                          
     &   b(mxvl,mxedof,*), ek(totdof*totdof,*), d(mxvl,nstr,*),                 
     &   bd(mxvl,mxedof,*), bt(mxvl,nstr,*)                                     
      integer icp(mxutsz,*)                                                     
c                                                                               
      do i=1,span                                                               
        ek(:,i) = ek(:,i) + reshape(matmul(b(i,1:totdof,1:6),                   
     &      matmul(d(i,1:6,1:6), transpose(b(i,1:totdof,1:6)))),                
     &      (/totdof*totdof/))                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end subroutine                                                            
c     ****************************************************************          
c     *                                                              *          
c     *           subroutine bdbt_asym2 : not used. saved for        *          
c     *                                   future reference                      
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/25/2015 rhd              *          
c     *                                                              *          
c     *     trans(B) * D * B for the asymmetric case.  Handles all   *          
c     *     elements, at least for now.                              *          
c     *     B and D are treated as full. D is non-symmetric.         *          
c     *     fastest of many variants of the operations.              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c      subroutine bdbt_asym2( span, b, d, ek,                                   
c     &                      mxvl, mxedof, nstr, totdof )                       
c      implicit none                                                            
c                                                                               
c                       parameter declarations                                  
c                                                                               
c      integer :: span, mxvl, mxedof, nstr, totdof                              
c      double precision ::                                                      
c     &   b(mxvl,mxedof,*), ek(totdof*totdof,*), d(mxvl,6,*)                    
c                                                                               
c      integer :: i, j, k                                                       
c                                                                               
c      double precision,                                                        
c     &  allocatable, dimension (:,:) ::                                        
c     &   local_b, local_bt, local_db, local_btdb                               
c      double precision :: local_d(6,6)                                         
c                                                                               
c                                                                               
c      allocate( local_b(6,totdof), local_bt(totdof,6),                         
c     &          local_db(6,totdof), local_btdb(totdof,totdof) )                
c                                                                               
c      do i = 1, span                                                           
c                                                                               
c        do k = 1, 6                                                            
c@!DIR$ IVDEP                                                                   
c          do j = 1, totdof                                                     
c             local_bt(j,k) = b(i,j,k)                                          
c             local_b(k,j)  = b(i,j,k)                                          
c          end do                                                               
c        end do                                                                 
c                                                                               
c        do j = 1, 6                                                            
c@!DIR$ IVDEP                                                                   
c          do k = 1, 6                                                          
c            local_d(k,j) = d(i,k,j)                                            
c          end do                                                               
c        end do                                                                 
c                                                                               
c                       matmul is faster than dgemm. reshape                    
c                       is faster than hand-coded loops                         
c                       allocates of local_.. really speeds up.                 
c                       reshape same speed as codes do loops                    
c                                                                               
c        local_db = matmul( local_d, local_b )                                  
c        local_btdb = matmul( local_bt, local_db )                              
c@!DIR$ IVDEP                                                                   
c        ek(:,i) = ek(:,i) + reshape(local_btdb,(/totdof*totdof/))              
c                                                                               
c      end do ! on span                                                         
c                                                                               
c      deallocate( local_b, local_bt, local_db, local_btdb )                    
c                                                                               
c      return                                                                   
c      end                                                                      
                                                                                
