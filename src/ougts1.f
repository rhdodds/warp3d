c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ougts1                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 1/27/2017 rhd              *          
c     *                                                              *          
c     *     transfer stresses or strains to simple data structure    *          
c     *     for direct output. Handle transformations for geometric  *          
c     *     nonlinear analyses. for stress output, material model    *          
c     *     has opportunity to insert additional values for output   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ougts1( span, blk, felem, do_stresses, ngp,                    
     &                   geonl, mat_type, matnum, iout, nowtime,                
     &                   nowstep )                                              
      use elblk_data, only : ddtse, elestr, rot_blk_n1, urcs_blk_n,             
     &                       elem_hist                                          
      use main_data, only : matprp, lmtprp, dmatprp                             
c                                                                               
      implicit none                                                             
      include 'param_def'                                                       
                                                                                
      integer :: span, blk, felem, ngp, mat_type, matnum, iout,                 
     &           nowstep                                                        
      logical :: do_stresses, geonl                                             
      double precision ::  nowtime                                              
c                                                                               
c                       local declarations                                      
c                                                                               
      integer :: iword(2)                                                       
      double precision :: dword                                                 
      equivalence ( iword, dword )                                              
c                                                                               
      integer :: info(10), gpn, i, ii, nprops, kinc, kout, nstatv,              
     &           ielem, noel, npt                                               
      logical :: do_strains                                                     
      logical ,parameter :: local_debug = .false.                               
      character(len=8) :: cmname                                                
      double precision :: qn(mxvl,nstr,nstr),  umat_props(50),                  
     &                    time, umat_statev(500), umat_stress(6),               
     &                    mat_vals(3)                                           
      double precision, parameter :: zero = 0.d0,                               
     &                               small_number = 1.0d-10,                    
     &                               root3 = dsqrt(3.0d0)                       
c                                                                               
      do_strains = .not. do_stresses                                            
c                                                                               
      if( do_strains ) then                                                     
c                                                                               
c                       do the element strains. for small and large             
c                       strains, the results are available in ddtse.            
c                       for geonl, these are the accumulated increments of      
c                       unrotated deformation tensors (\delta d).               
c                       rotate to current spatial axes to define strain values  
c                       for output. These are the equivalent of incremental     
c                       D tensors with accumulated values at n rotation to      
c                       n+1 before addition. Same approach used in Abaqus       
c                       explicit. Components are always in global XYZ system.   
c                                                                               
c                       Example: impose large axial deformation on an element in
c                       global X direction. Key strain is eps-xx. Now impose    
c                       rigid rotation of element +90-degrees around Z. Unrotate
c                       strains do not change. Rotated strains for output have  
c                       eps-yy now as key value. Output strains are thus "rotati
c                       corrected" logarithmic values (assuming the load        
c                       descritization into steps approximates adequately the   
c                       logarithmic integral. i.e., int dx/x.                   
c                                                                               
c                       The output strains are thus (approximately) work        
c                       conjugate  with the output Cauchy stresses.             
c                       This the same strain output as Abaqus - but             
c                       must continually update the rotated strain              
c                       fixed global coords at each time step.                  
c                                                                               
                                                                                
        if( .not. geonl ) then                                                  
c                                                                               
c                       small strains                                           
c                                                                               
           do gpn = 1, ngp                                                      
             do i = 1, span                                                     
               elestr(i,1,gpn) = ddtse(i,1,gpn)                                 
               elestr(i,2,gpn) = ddtse(i,2,gpn)                                 
               elestr(i,3,gpn) = ddtse(i,3,gpn)                                 
               elestr(i,4,gpn) = ddtse(i,4,gpn)                                 
               elestr(i,5,gpn) = ddtse(i,5,gpn)                                 
               elestr(i,6,gpn) = ddtse(i,6,gpn)                                 
            end do                                                              
           end do                                                               
           return                                                               
        end if                                                                  
c                                                                               
c                       large strains/rotations                                 
c                                                                               
        do gpn = 1, ngp                                                         
         call ou_get_spatial_from_material( span, elestr(1,1,gpn),              
     &               rot_blk_n1(1,1,gpn), ddtse(1,1,gpn) )                      
        end do                                                                  
        return                                                                  
c                                                                               
      end if     !   on do_strains                                              
c                                                                               
c                       build block structure for stresses for                  
c                       output. we insert the 6 stress                          
c                       components and the energy density. then                 
c                       zero a location for the mises equiv. stress.            
c                       Insert 3 material model specific values thru            
c                       a material model routine.                               
c                                                                               
c                       for geometric nonlinear, convert unrotated              
c                       cauchy stresses                                         
c                       to cauchy stresses at n+1 using [R,n+1].                
c                       {T} = [qn] * {urcs}.                                    
      do gpn = 1, ngp                                                           
        if ( geonl ) then                                                       
          call getrm1( span, qn, rot_blk_n1(1,1,gpn), 2 )                       
          call qmply1( span, mxvl, nstr, qn, urcs_blk_n(1,1,gpn),               
     &                 elestr(1,1,gpn) )                                        
          do ii = 1, span                                                       
             elestr(ii,7,gpn) = urcs_blk_n(ii,7,gpn)                            
          end do                                                                
        else    !  load in small-strain theory stresses                         
          do i = 1, span                                                        
             elestr(i,1,gpn)  = urcs_blk_n(i,1,gpn)                             
             elestr(i,2,gpn)  = urcs_blk_n(i,2,gpn)                             
             elestr(i,3,gpn)  = urcs_blk_n(i,3,gpn)                             
             elestr(i,4,gpn)  = urcs_blk_n(i,4,gpn)                             
             elestr(i,5,gpn)  = urcs_blk_n(i,5,gpn)                             
             elestr(i,6,gpn)  = urcs_blk_n(i,6,gpn)                             
             elestr(i,7,gpn)  = urcs_blk_n(i,7,gpn)                             
          end do                                                                
        end if                                                                  
c                                                                               
        do i = 1, span                                                          
           elestr(i,8,gpn)  = zero                                              
           elestr(i,9,gpn)  = zero                                              
           elestr(i,10,gpn) = zero                                              
           elestr(i,11,gpn) = zero                                              
        end do                                                                  
c                                                                               
      end do   ! on gpn                                                         
c                                                                               
c                       invoke the material model dependent                     
c                       routine to insert the 3 "mat_val"s.                     
c                                                                               
      select case( mat_type )                                                   
                                                                                
      case ( 1 )                                                                
c                                                                               
c                       bilinear mises model                                    
c                                                                               
           do gpn = 1, ngp                                                      
              do i = 1, span                                                    
                elestr(i,10,gpn) = elem_hist(i,2,gpn)*root3                     
                elestr(i,9,gpn)  = urcs_blk_n(i,9,gpn)                          
                dword            = elem_hist(i,4,gpn)                           
                elestr(i,11,gpn) = iword(1)                                     
               end do                                                           
           end do                                                               
c                                                                               
      case ( 2 )                                                                
c                                                                               
c                       power-law deformation plasticity model                  
c                                                                               
           do gpn = 1, ngp                                                      
              do i = 1, span                                                    
                elestr(i,9,gpn) =  urcs_blk_n(i,9,gpn)                          
               end do                                                           
           end do                                                               
c                                                                               
      case ( 3 )                                                                
c                                                                               
c                       general mises and gurson-tvergaard model                
c                                                                               
           do gpn = 1, ngp                                                      
              do i = 1, span                                                    
                elestr(i,9,gpn)  = urcs_blk_n(i,9,gpn)                          
                elestr(i,10,gpn) = elem_hist(i,2,gpn)                           
                elestr(i,11,gpn) = elem_hist(i,5,gpn)                           
              end do                                                            
           end do                                                               
c                                                                               
      case ( 5 )                                                                
c                                                                               
c                       adv. cyclic plasticity model                            
c                                                                               
           do gpn = 1, ngp                                                      
             call oumm05( gpn, mxvl, span, iout, elestr(1,1,gpn),               
     &                    urcs_blk_n(1,1,gpn), elem_hist(1,1,gpn) )             
           end do                                                               
c                                                                               
      case ( 6 )                                                                
c                                                                               
c                       creep model                                             
c                                                                               
           do gpn = 1, ngp                                                      
             call oumm06( gpn, mxvl, span, iout, elestr(1,1,gpn),               
     &                    urcs_blk_n(1,1,gpn), elem_hist(1,1,gpn) )             
           end do                                                               
c                                                                               
      case ( 7 )                                                                
c                                                                               
c                       adv. mises + hydrogen                                   
c                                                                               
           do gpn = 1, ngp                                                      
             call oumm07( gpn, mxvl, span, iout, elestr(1,1,gpn),               
     &                    urcs_blk_n(1,1,gpn), elem_hist(1,1,gpn) )             
           end do                                                               
c                                                                               
      case ( 8 )                                                                
c                                                                               
c                       Abaqus compatible umat                                  
c                                                                               
        cmname(1:) = "UMAT-WRP"                                                 
        umat_props(1:50) =  dmatprp(151:200,matnum)                             
        nprops = 50                                                             
        time = nowtime                                                          
        kinc = nowstep                                                          
        kout = iout                                                             
        call umat_set_features( info )                                          
        nstatv = info(1)                                                        
        do gpn = 1, ngp                                                         
         do ielem = 1, span                                                     
           umat_statev(1:nstatv) = elem_hist(ielem,1:nstatv,gpn)                
           noel = felem + ielem - 1                                             
           npt = gpn                                                            
           umat_stress(1:4) =  elestr(ielem,1:4,gpn)                            
           umat_stress(5) =  elestr(ielem,6,gpn)                                
           umat_stress(6) =  elestr(ielem,5,gpn)                                
           mat_vals(1:3) = zero                                                 
           call umat_output( umat_stress, mat_vals, umat_statev, nstatv,        
     &                       time, cmname, umat_props, nprops, noel,            
     &                       npt, kinc, kout )                                  
           elestr(ielem,9:11,gpn) = mat_vals(1:3)                               
        end do                                                                  
       end do                                                                   
c                                                                               
c           Crystal plasticity: c1 -> elastic work density                      
c                               c2 -> plastic work density                      
c                               c3 -> accumulated plastic strain                
c                   urcs_blk_n(i,7,gpn) enters with total work density.         
c                   that value is printed by WARP3D under "energy"              
c                   label of stresses                                           
c                                                                               
      case (10 )                                                                
c                                                                               
          do gpn = 1, ngp                                                       
              do i = 1, span                                                    
                elestr(i,9,gpn)  = urcs_blk_n(i,7,gpn) -                        
     &                             urcs_blk_n(i,8,gpn)                          
                elestr(i,10,gpn) = urcs_blk_n(i,8,gpn)                          
                elestr(i,11,gpn) = urcs_blk_n(i,9,gpn)                          
              end do                                                            
           end do                                                               
                                                                                
      case default                                                              
c                                                                               
c                 This is a kluge to avoid a patran bug -- in patran 2.5,       
c                 if a element value that you want to print is                  
c                 identically zero, then patran will not print the element.     
c                 This is a problem for displaying the porosity.  Thus          
c                 here we set the porosity to be equal to a small number        
c                 so that patran will print the non-gurson elements.            
c                                                                               
           do gpn = 1, ngp                                                      
              do i = 1, span                                                    
                elestr(i,11,gpn) = small_number                                 
              end do                                                            
           end do                                                               
c                                                                               
      end select                                                                
c                                                                               
      return                                                                    
c                                                                               
 9100 format(i4,6f15.6,/,4x,6f15.6)                                             
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine ou_get_spatial_from_material          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/30/2014 rhd              *          
c     *                                                              *          
c     *     transforms deformation rate d -> spatial rate D.         *          
c     *     can be used for increments \Delta d -> \Delta D          *          
c     *     and the acummulated \Delta d over all loading to get the *          
c     *     equivalent sum \Delta D where incremental rotations      *          
c     *     are performed during each time step (what Abaqus does)   *          
c     *     convert from vector form of d to tensor, rotate,         *          
c     *     extract vector form of spatial strain                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ou_get_spatial_from_material( span, spatial_strain,            
     &                                         R, material_strain )             
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c             R        -- 3x3 rotations at int pts from polar                   
c                         decompositions F = RU                                 
c      material_strain -- the unrotated strain ("d") in vector                  
c                         form (input)                                          
c      spatial_strain  -- the spatial strain "D" in vector                      
c                         form (output)                                         
c      span            -- number of elements in block                           
c                                                                               
c            spatial_strain = R * material_strain * trans( R )                  
c                                                                               
c            convert vector strain to strain tensor form, rotate,               
c            convert back to vector                                             
c                                                                               
      integer :: span                                                           
      double precision :: R(mxvl,3,3), spatial_strain(mxvl,6),                  
     &                    material_strain(mxvl,6)                               
c                                                                               
c            local work arrays                                                  
c                                                                               
      integer :: i                                                              
      double precision :: rt(mxvl,3,3), t(mxvl,3,3), d(mxvl,3,3),               
     &                    Deps(mxvl,3,3)                                        
      double precision, parameter ::  half = 0.5d0, two = 2.0d0                 
c                                                                               
      do i = 1, span                                                            
c                                                                               
c             pull trans (R ) for convenience                                   
c                                                                               
        rt(i,1,1) = R(i,1,1)                                                    
        rt(i,1,2) = R(i,2,1)                                                    
        rt(i,1,3) = R(i,3,1)                                                    
        rt(i,2,1) = R(i,1,2)                                                    
        rt(i,2,2) = R(i,2,2)                                                    
        rt(i,2,3) = R(i,3,2)                                                    
        rt(i,3,1) = R(i,1,3)                                                    
        rt(i,3,2) = R(i,2,3)                                                    
        rt(i,3,3) = R(i,3,3)                                                    
c                                                                               
c             form unrotated strain tensor from vector                          
c             material (unrotated) strains                                      
c                                                                               
        d(i,1,1) = material_strain(i,1)                                         
        d(i,1,2) = material_strain(i,4) * half                                  
        d(i,1,3) = material_strain(i,6) * half                                  
        d(i,2,1) = material_strain(i,4) * half                                  
        d(i,2,2) = material_strain(i,2)                                         
        d(i,2,3) = material_strain(i,5) * half                                  
        d(i,3,1) = material_strain(i,6) * half                                  
        d(i,3,2) = material_strain(i,5) * half                                  
        d(i,3,3) = material_strain(i,3)                                         
c                                                                               
c             temp = unrotated strain tensor d * trans(R)                       
c                                                                               
        t(i,1,1) =  d(i,1,1)*rt(i,1,1) + d(i,1,2)*rt(i,2,1) +                   
     &                 d(i,1,3)*rt(i,3,1)                                       
        t(i,1,2) =  d(i,1,1)*rt(i,1,2) + d(i,1,2)*rt(i,2,2) +                   
     &                 d(i,1,3)*rt(i,3,2)                                       
        t(i,1,3) =  d(i,1,1)*rt(i,1,3) + d(i,1,2)*rt(i,2,3) +                   
     &                 d(i,1,3)*rt(i,3,3)                                       
c                                                                               
        t(i,2,1) =  d(i,2,1)*rt(i,1,1) + d(i,2,2)*rt(i,2,1) +                   
     &                 d(i,2,3)*rt(i,3,1)                                       
        t(i,2,2) =  d(i,2,1)*rt(i,1,2) + d(i,2,2)*rt(i,2,2) +                   
     &                 d(i,2,3)*rt(i,3,2)                                       
        t(i,2,3) =  d(i,2,1)*rt(i,1,3) + d(i,2,2)*rt(i,2,3) +                   
     &                 d(i,2,3)*rt(i,3,3)                                       
c                                                                               
        t(i,3,1) =  d(i,3,1)*rt(i,1,1) + d(i,3,2)*rt(i,2,1) +                   
     &                 d(i,3,3)*rt(i,3,1)                                       
        t(i,3,2) =  d(i,3,1)*rt(i,1,2) + d(i,3,2)*rt(i,2,2) +                   
     &                 d(i,3,3)*rt(i,3,2)                                       
        t(i,3,3) =  d(i,3,1)*rt(i,1,3) + d(i,3,2)*rt(i,2,3) +                   
     &                 d(i,3,3)*rt(i,3,3)                                       
c                                                                               
c             spatial strain tensor D = R * temp                                
c                                                                               
        Deps(i,1,1) =  r(i,1,1)*t(i,1,1) + r(i,1,2)*t(i,2,1) +                  
     &                 r(i,1,3)*t(i,3,1)                                        
        Deps(i,1,2) =  r(i,1,1)*t(i,1,2) + r(i,1,2)*t(i,2,2) +                  
     &                 r(i,1,3)*t(i,3,2)                                        
        Deps(i,1,3) =  r(i,1,1)*t(i,1,3) + r(i,1,2)*t(i,2,3) +                  
     &                 r(i,1,3)*t(i,3,3)                                        
c                                                                               
        Deps(i,2,1) =  r(i,2,1)*t(i,1,1) + r(i,2,2)*t(i,2,1) +                  
     &                 r(i,2,3)*t(i,3,1)                                        
        Deps(i,2,2) =  r(i,2,1)*t(i,1,2) + r(i,2,2)*t(i,2,2) +                  
     &                 r(i,2,3)*t(i,3,2)                                        
        Deps(i,2,3) =  r(i,2,1)*t(i,1,3) + r(i,2,2)*t(i,2,3) +                  
     &                 r(i,2,3)*t(i,3,3)                                        
c                                                                               
        Deps(i,3,1) =  r(i,3,1)*t(i,1,1) + r(i,3,2)*t(i,2,1) +                  
     &                 r(i,3,3)*t(i,3,1)                                        
        Deps(i,3,2) =  r(i,3,1)*t(i,1,2) + r(i,3,2)*t(i,2,2) +                  
     &                 r(i,3,3)*t(i,3,2)                                        
        Deps(i,3,3) =  r(i,3,1)*t(i,1,3) + r(i,3,2)*t(i,2,3) +                  
     &                 r(i,3,3)*t(i,3,3)                                        
c                                                                               
c             return vector form of spatial strain tensor D                     
c                                                                               
        spatial_strain(i,1) = Deps(i,1,1)                                       
        spatial_strain(i,2) = Deps(i,2,2)                                       
        spatial_strain(i,3) = Deps(i,3,3)                                       
        spatial_strain(i,4) = two * Deps(i,2,1)                                 
        spatial_strain(i,5) = two * Deps(i,2,3)                                 
        spatial_strain(i,6) = two * Deps(i,1,3)                                 
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
