c                                                                               
c           user_routines_umat.f   Distribution version                         
c                                                                               
c           Updated: 8/20/2016 rhd                                              
c                                                                               
c                                                                               
c           The first routine here, umat_set_features, is called by WARP3D      
c           at various times to obtain information about the umat.              
c                                                                               
c           This file is compiled and included in the normal executable         
c           build of WARP3D.                                                    
c                                                                               
c           WARP3D calls umat routines in:                                      
c                                                                               
c              gplns1.f    -- linear stiffness computation                      
c              rstgp1.f    -- stress update and new [D] computation             
c                                                                               
c           These above two code fies have lots of comments about setting       
c           up data arrays, values for the umat from WARP3D data                
c           structures.                                                         
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine umat_set_features                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 12/14/14 rhd                *          
c     *                                                              *          
c     *  called by warp3d for the umat to obtain statev size and     *          
c     *  other characteristic information about the umat             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine umat_set_features( info_vector )                               
      implicit none                                                             
      integer info_vector(*)                                                    
c                                                                               
c        set info_data                                                          
c                                                                               
c         1        number of history values per integration                     
c                  point. Abaqus calls these "statev". Values                   
c                  double or single precsion based on hardware.                 
c                                                                               
c         2        number of values in the symmetric part of the                
c                  [D] for each integration point. for solid                    
c                  elements this is 21, for cohesive elements this 6.           
c                                                                               
c         3        = 0, the material model returns "unrotated"                  
c                       Cauchy stresses at n+1                                  
c                  = 1, the material model returns the standard                 
c                       Cauchy stresses at n+1                                  
c                                                                               
c         4        number of state variables per point to be output             
c                  when user requests this type of results                      
c                                                                               
c                                                                               
c                                                                               
c        example umat included here is for mises plasticity with                
c        bilinear kinematic hardening                                           
c                                                                               
      info_vector(1) = 20   ! includes back stresses                            
      info_vector(2) = 21                                                       
      info_vector(3) = 0                                                        
      info_vector(4) = 14                                                       
c                                                                              
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *   subroutine umat -- bilinear kinematic hardening            *          
c     *                                                              *          
c     *  Abaqus example umat for elastic-plastic response with       *          
c     *  linear kinematic hardening                                  *          
c     *                                                              *          
c     ****************************************************************          
      subroutine umat( stress, statev, ddsdde, sse, spd, scd, rpl,              
     1 ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,                
     2 dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv,                  
     3 props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1,             
     4 noel, npt, layer, kspt, kstep, kinc, kiter, kout,                        
     5 kthread, knumthreads, nonlocal_shared, nshared )                         
c                                                                               
      implicit none                                                             
                                                                                
      integer :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer,            
     1           kspt, kstep, kinc, kiter, kout, kthread,                       
     2           knumthreads, nshared                                           
      integer, parameter :: nprecd=2                                            
c                                                                               
      double precision :: stress(*), statev(nstatv),                            
     1 ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens), stran(ntens),         
     2 dstran(ntens), predef(1), dpred(1), props(nprops), coords(3),            
     3 drot(3,3), dfgrd0(3,3), dfgrd1(3,3), nonlocal_shared(nshared)
c
      logical :: yielding        
c                                                                               
      double precision :: time, dtime, temp, dtemp, sse, spd, scd, rpl,         
     1                    drpldt, celent, pnewdt                                
                                                                                
      character(len=8) :: cmname                                                
c                                                                               
c                                                                               
c           locals                                                              
c                                                                               
c              eelas  - elastic strains                                         
c              eplas  - plastic strains                                         
c              alpha  - shift tensor                                            
c              flow   - plastic flow directions                                 
c              olds   - stress at start of increment                            
c              oldpl  - plastic strains at start of increment                   
c                                                                               
      integer :: k1, k2                                                         
      double precision :: eelas(6), eplas(6), alpha(6), flow(6),                
     1                    olds(6), oldpl(6), avg_stress(6), 
     2                    deps_plas(6), dsig(6)                                     
      logical ::  debug                                    
      double precision :: cte, emod, enu, ebulk3, eg2, eg, eg3, elam,           
     1                    old_eqpl, deqpl, smises, syield, hard, shydro,        
     2                    effg, effg2, effg3, efflam, effhrd        
      double precision :: state_np1, dspd, dsst                                 
      double precision, parameter :: zero=0.d0, one=1.d0, two=2.d0,             
     1                               three=3.d0, six=6.d0,                      
     2                               enumax=.4999d0, toler=1.0d-6,
     3                               half = 0.5d0               
c                                                                               
                                                                                
c                                                                               
c ----------------------------------------------------------------              
c           umat for isotropic elasticity and mises plasticity                  
c           with kinematic hardening                                            
c                                                                               
c                 props(1) - e                                                  
c                 props(2) - nu                                                 
c                 props(3) - syield                                             
c                 props(4) - plastic hardening modulus (H')                     
c                 props(5) - coefficient of thermal expansion                   
c                                                                               
c                 statev requries 20 entries per material point                 
c                 for current model                                             
c                                                                               
c ----------------------------------------------------------------              
c                                                                               
c      
      debug =    .false. ! noel .eq. 2648 .and. npt .eq. 8                                                                      
      if( debug ) write(kout,*) ' .... umat '    
      if( debug ) write(kout,*) '... spd in: ', spd                                    
      cte = props(5)                                                            
      dstran(1) = dstran(1) - cte * dtemp                                       
      dstran(2) = dstran(2) - cte * dtemp                                       
      dstran(3) = dstran(3) - cte * dtemp                                       
c      if( kiter .eq. 0 ) return   ! for new WARP3D solution model              
c                                                                               
c           elastic material constants                                          
c                                                                               
      emod = props(1)                                                           
      enu = min(props(2), enumax)                                               
      ebulk3 = emod/(one-two*enu)                                               
      eg2 = emod/(one+enu)                                                      
      eg = eg2/two                                                              
      eg3 = three*eg                                                            
      elam = (ebulk3-eg2)/three                                                 
c                                                                               
c           linear-elastic material stiffness                                   
c                                                                               
      do k1 = 1, 3                                                              
        do k2 = 1, 3                                                            
          ddsdde(k2,k1) = elam                                                  
        end do                                                                  
        ddsdde(k1,k1) = eg2 + elam                                              
      end do                                                                    
      ddsdde(4,4) = eg                                                          
      ddsdde(5,5) = eg                                                          
      ddsdde(6,6) = eg   
      yielding = .false.                                                       
c                                                                               
c           recover elastic strain, plastic strain and shift tensor             
c           and rotate. note:                                                   
c           use code 1 for (tensor) stress, code 2 for (engineering) strain     
c           recover plastic strain at n. set current yielding state.            
c           0=never yielded, 1=active plasticity, 2=prior plasticity by         
c           currently linear elastic                                            
c                                                                               
      call rotsig( statev(1),  drot, eelas, 2, 3, 3)                            
      call rotsig( statev(7),  drot, eplas, 2, 3, 3)                            
      call rotsig( statev(13), drot, alpha, 1, 3, 3)                            
      old_eqpl   = statev(19)                                                   
      deqpl      = zero                                                         
      dspd       = zero                                                         
      dsst       = zero                                                         
      state_np1  = zero                                                         
      if( old_eqpl .gt. zero ) state_np1 = two                                  
c                                                                               
c           save stress and plastic strains and calculate predictor             
c           (trial elastic) stress and elastic strain. the original             
c           umat coding from Abaqus PPT slides for next few loops               
c           had an error.                                                       
c                                                                               
      do k1 = 1, 6                                                              
         olds(k1)  = stress(k1)                                                 
         oldpl(k1) = eplas(k1)                                                  
      end do                                                                    
c                                                                               
      do k1 = 1, 6                                                              
         eelas(k1) = eelas(k1) + dstran(k1)                                     
         do k2 = 1, 6                                                           
           stress(k2) = stress(k2) + ddsdde(k2,k1)*dstran(k1)                   
         end do                                                                 
      end do                                                                    
c                                                                               
c           calculate equivalent von mises stress                               
c                                                                               
      smises = (stress(1)-alpha(1)-stress(2)+alpha(2))**2                       
     1      + (stress(2)-alpha(2)-stress(3)+alpha(3))**2                        
     2      + (stress(3)-alpha(3)-stress(1)+alpha(1))**2                        
      smises = smises+six*(stress(4)-alpha(4))**2                               
      smises = smises+six*(stress(5)-alpha(5))**2                               
      smises = smises+six*(stress(6)-alpha(6))**2                               
      smises = sqrt(smises/two)                                                 
c                                                                               
c           get yield stress and plastic hardening modulus                      
c                                                                               
      syield = props(3)                                                         
      hard   = props(4)                                                         
c                                                                               
c           determine if actively yielding                                      
c                                                                               
      if( smises .gt. (one+toler)*syield ) then                                 
c                                                                               
c              actively yielding                                                
c                separate the hydrostatic from the deviatoric stress            
c                calculate the flow direction                                   
c                                                                               
         shydro  = (stress(1)+stress(2)+stress(3))/three                        
         flow(1) = (stress(1)-alpha(1)-shydro)/smises                           
         flow(2) = (stress(2)-alpha(2)-shydro)/smises                           
         flow(3) = (stress(3)-alpha(3)-shydro)/smises                           
         flow(4) = (stress(4)-alpha(4))/smises                                  
         flow(5) = (stress(5)-alpha(5))/smises                                  
         flow(6) = (stress(6)-alpha(6))/smises                                  
c                                                                               
c                 solve for equivalent plastic strain increment                 
c                                                                               
         deqpl = (smises-syield)/(eg3+hard)                                     
         state_np1 = one
         yielding = .true.                                                        
c                                                                               
c                 update shift tensor, elastic and plastic strains              
c                 and stress                                                    
c                                                                               
         do k1 = 1, 3                                                           
           alpha(k1)  = alpha(k1)+hard*flow(k1)*deqpl                           
           eplas(k1)  = eplas(k1)+three/two*flow(k1)*deqpl                      
           eelas(k1)  = eelas(k1)-three/two*flow(k1)*deqpl                      
           stress(k1) = alpha(k1)+flow(k1)*syield+shydro                        
         end do                                                                 
         do k1 =  4, 6                                                          
           alpha(k1)  = alpha(k1)+hard*flow(k1)*deqpl                           
           eplas(k1)  = eplas(k1)+three*flow(k1)*deqpl                          
           eelas(k1)  = eelas(k1)-three*flow(k1)*deqpl                          
           stress(k1) = alpha(k1)+flow(k1)*syield                               
         end do                                                                 
c                                                                               
c                 formulate the jacobian (material tangent)                     
c                 first calculate effective moduli                              
c                                                                               
         effg  = eg*(syield+hard*deqpl)/smises                                  
         effg2  = two*effg                                                      
         effg3  = three*effg                                                    
         efflam = (ebulk3-effg2)/three                                          
         effhrd = eg3*hard/(eg3+hard)-effg3                                     
         do k1 = 1, 3                                                           
           do k2 = 1, 3                                                         
             ddsdde(k2,k1) = efflam                                             
           end do                                                               
           ddsdde(k1,k1) = effg2 + efflam                                       
         end do                                                                 
         ddsdde(4,4) = effg                                                     
         ddsdde(5,5) = effg                                                     
         ddsdde(6,6) = effg                                                     
         do k1 = 1, 6                                                           
           do k2 = 1, 6                                                         
             ddsdde(k2,k1) = ddsdde(k2,k1) + effhrd*flow(k2)*flow(k1)           
           end do                                                               
         end do                                                                 
      end if                                                                    
c                                                                               
c          store elastic strains, plastic strains, shift tensor                 
c          and updated (scalar) plastic strain in state variable array.         
c                                                                               
      do k1 = 1, 6                                                              
         statev(k1)    = eelas(k1)                                              
         statev(k1+6)  = eplas(k1)                                              
         statev(k1+12) = alpha(k1)                                              
         avg_stress(k1) = ( stress(k1) + olds(k1) ) / two 
         dsig(k1) = stress(k1) - olds(k1)                        
         dsst = dsst + avg_stress(k1) * dstran(k1) 
      end do   
      if( yielding ) then
         deps_plas(1) = dstran(1) -
     &                  (dsig(1) - enu*(dsig(2)+dsig(3)))/emod
         deps_plas(2) = dstran(2) -
     &                  (dsig(2) - enu*(dsig(1)+dsig(3)))/emod
         deps_plas(3) = dstran(3) -
     &                  (dsig(3) - enu*(dsig(1)+dsig(2)))/emod
         deps_plas(4) = dstran(4) - dsig(4) / eg
         deps_plas(5) = dstran(5) - dsig(5) / eg
         deps_plas(6) = dstran(6) - dsig(6) / eg
         dspd = half * dot_product( deps_plas, avg_stress )
      end if
                                                                 
c                                                                               
c          store updated scalar plastic strain. material state                  
c          = 0 (never yielded), = 1 actively yielding, = 2                      
c          (previously yielded but currently linear-elastic)                    
c          these values can be output by WARP3D. see                            
c          optional ou_umat routine.                                            
c                                                                               
      statev(19) = old_eqpl + deqpl                                             
      statev(20) = state_np1                                                    
c                                                                               
c          update specific elastic energy and dissipatation.                    
c                                                                               
      sse = sse + (dsst - dspd)                                                 
      spd = spd + dspd                                                          
      if( debug ) write(kout,*) '... spd out: ', spd                                    
c                                                                               
c          fill nonlocal shared values for integration point with               
c          the element number for testing purposes.                             
c                                                                               
      nonlocal_shared(1:nshared) = dble( noel )                                 
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine umat_states_labels                    *          
c     *                                                              *          
c     *   call by WARP3D to get the number of states variables for   *          
c     *   output, an 8 character id for each variable and an         *          
c     *   descriptor string for each variable                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine umat_states_labels( size_state,                                
     &      num_states, state_labels, state_descriptors, out,                   
     &      comment_lines, max_comment_lines, num_comment_lines )               
      implicit none                                                             
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: size_state, num_states, out, max_comment_lines,                
     &           num_comment_lines                                              
      character(len=8)  :: state_labels(size_state)                             
      character(len=60) :: state_descriptors(size_state)                        
      character(len=80) :: comment_lines(max_comment_lines)                     
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i                                                              
c                                                                               
      num_states = 14                                                           
      num_comment_lines = 5                                                     
c                                                                               
      state_labels(1) = "epspls"                                                
      state_labels(2) = "state"                                                 
      state_labels(3) = "epls_xx"                                               
      state_labels(4) = "epls_yy"                                               
      state_labels(5) = "epls_zz"                                               
      state_labels(6) = "epls_xy"                                               
      state_labels(7) = "epls_xz"                                               
      state_labels(8) = "epls_yz"                                               
      state_labels(9) = "alpha_xx"                                              
      state_labels(10) = "alpha_yy"                                             
      state_labels(11) = "alpha_zz"                                             
      state_labels(12) = "alpha_xy"                                             
      state_labels(13) = "alpha_xz"   ! abaqus ordering                         
      state_labels(14) = "alpha_yz"   ! abaqus ordering                         
c                                                                               
      state_descriptors(1) = "Equivalent plastic strain"                        
      state_descriptors(2) = "=1, 2, 3 see states_header"                       
      state_descriptors(2:7) = "Plastic strain component"                       
      state_descriptors(8:14) = "Backstress component"                          
c                                                                               
      comment_lines(1) = "Notes on 'state' quantity"                            
      comment_lines(2) = "  = 1 never yielded"                                  
      comment_lines(3) = "  = 2 active yielding"                                
      comment_lines(4) = "  = 3 previous plasticity. "                          
     &                       // " now linear-elastic"                           
      comment_lines(5) = "  value is average of int points"                     
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *             subroutine umat_output_states                       *           
c *                                                                 *           
c *     routine called just before output of states. return         *           
c *     state values for output at this integration point           *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine umat_output_states( statev, output_statev, nvals )             
      implicit none                                                             
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: nvals                                                          
      double precision ::                                                       
     & statev(*), output_statev(*)                                              
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: k, k1                                                          
c                                                                               
      nvals = 14                                                                
c                                                                               
      output_statev(1) = statev(19)                                             
      output_statev(2) = statev(20)                                             
      output_statev(3:14) = statev(7:18)                                        
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *  optional UMAT routine to set additional stress output values   *           
c *                                                                 *           
c *   set up to 3 material model dependent output values            *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine umat_output( stress, mat_vals, statev, nstatv, time,           
     &                        cmname, props, nprops, noel,                      
     &                        npt, kinc, kout )                                 
      implicit double precision (a-h,o-z)                                       
c                                                                               
c                   parameter declarations                                      
c                   ----------------------                                      
c                                                                               
      character(len=8) :: cmname                                                
      double precision                                                          
     & stress(6), mat_vals(3), statev(nstatv), props(nprops)                    
c                                                                               
c               description of parameters                                       
c               -------------------------                                       
c                                                                               
c     stress            : 6 x 1 stresses for point. Cauchy stresses for         
c                         GEONL solution. Abaqus ordering (see below)           
c (*) mat_vals          : 3x1. these are the c1, c2, c3 values that may be      
c                         set by this routine for including in the stress       
c                         output below. Values are floating point. WARP3D       
c                         does nothing but print values.                        
c     statev            : current state vector for point                        
c     nstatv            : umat decalerd number of state variables/point.        
c                         see routine umat_set_features                         
c     time              : current simultation time (sum of all load step        
c                         time increments)                                      
c     cmname            : always "UMAT-WRP"                                     
c     props             : values of umat properties specified by user           
c                         in model definition (um_1, um_2, um_3, ...)           
c     nprops            : umat spedified number of properties                   
c                         (see umat_set_features)                               
c     noel              : current element number in model                       
c     npt               : current material point number                         
c     kinc              : current analysis (WARP3D) load step number            
c     kout              : Fortran device number for debugging messages          
c                                                                               
c    =>   only mat_vals should be modified by this routine                      
c                                                                               
c   stress ordering                                                             
c     (1) sig-xx                                                                
c     (2) sig-yy                                                                
c     (3) sig-zz                                                                
c     (4) tau-xy                                                                
c     (5) tau-xz                                                                
c     (6) tau-yz                                                                
c                                                                               
c            this example implementation supports the umat for                  
c            bilinear plasticity with kinematic hardening.                      
c            see umat code for details.                                         
c                                                                               
c            **** modify code below for your umat -- or delete the   ****       
c            **** assignment statements to mat_vals                  ****       
c                                                                               
c            statev(1:6)   = elastic strains                                    
c            statev(7:12)  = plastic strains                                    
c            statev(13:18) = backstress vector                                  
c            statev(19)    = scalar plastic strain (eqpl)                       
c            statev(20)    = curretn loading state at point                     
c                            = 0.0 never yielded                                
c                            = 1.0 actively yielded this step                   
c                            = 2.0 previosuly yielding, but not this step       
c                                                                               
c                                                                               
c                   local declarations                                          
c                   ------------------                                          
c                                                                               
c                                                                               
      logical debug                                                             
c                                                                               
      debug = .false.                                                           
c      debug =  noel .eq. 10 .and. npt .eq. 1                                   
c                                                                               
      if( debug ) then                                                          
          call getnumcpus( numranks )                                           
          call getrank( myrank )                                                
          write(kout,9000) noel, npt, kinc, time, numranks, myrank              
          write(kout,9010) stress(1:6)                                          
          write(kout,9020) statev(1:nstatv)                                     
          write(kout,9030) props(1:5)                                           
          write(kout,9040)                                                      
      end if                                                                    
c                                                                               
      mat_vals(1) = statev(19)                                                  
      mat_vals(2) = statev(20)                                                  
      mat_vals(3) = 0.0                                                         
c                                                                               
      return                                                                    
 9000 format(/,'.... debugging from umat_output ....',                          
     & /,10x,'noel:                ',i7,                                        
     & /,10x,'npt:                 ',i7,                                        
     & /,10x,'kinc (load step):    ',i7,' simulation time: ',e14.6,             
     & /,10x,'numranks:            ',i7,                                        
     & /,10x,'myrank:              ',i7 )                                       
 9010 format(/,'   stresses: ',6e14.6 )                                         
 9020 format(/,'   state variables: ', 30(/10x,6 e14.6) )                       
 9030 format(/,'   1st 5 umat properties: ',                                    
     & /,10x,5e14.6 )                                                           
 9040 format(//)                                                                
c                                                                               
       end                                                                      
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine uhard  (called by UMATs)          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 3/22/12                     *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
      subroutine uhard( syield, hard, eqplas, eqplasrt, time, dtime,            
     &     temp, dtemp, noel, npt, layer, kspt, kstep, kinc,                    
     &     cmname, nstatv, statev, numfieldv,                                   
     &     predef, dpred, nvalue, table, kout )                                 
      implicit double precision (a-h,o-z)                                       
      parameter (nprecd=2)                                                      
c                                                                               
      character(len=80) :: cmname                                               
      dimension hard(3), statev(nstatv), time(*),                               
     1          predef(numfieldv), dpred(*)                                     
c                                                                               
      dimension table(2,nvalue)                                                 
c                                                                               
      parameter(zero = 0.d0)                                                    
c                                                                               
c            set yie      ld stress to last value of table,                     
c            hardening to zero                                                  
c                                                                               
      syield  = table(1,nvalue)                                                 
      hard(1) = zero                                                            
c                                                                               
c            if more than one entry, search table                               
c                                                                               
      if( nvalue .gt. 1 ) then                                                  
         do k1 = 1, nvalue-1                                                    
           eqpl1 = table(2,k1+1)                                                
           if( eqplas .lt. eqpl1 ) then                                         
              eqpl0 = table(2,k1)                                               
              if( eqpl1 .le .eqpl0 ) then                                       
                  write(kout,100)                                               
                  call xit                                                      
              endif                                                             
c                                                                               
c            current yield stress and hardening                                 
c                                                                               
              deqpl = eqpl1 - eqpl0                                             
              syiel0 = table(1,k1)                                              
              syiel1 = table(1,k1+1)                                            
              dsyiel = syiel1 - syiel0                                          
              hard(1) = dsyiel / deqpl                                          
              syield = syiel0 + ( eqplas-eqpl0 ) * hard(1)                      
              goto 10                                                           
            endif                                                               
           end do                                                               
10         continue                                                             
      endif                                                                     
      return                                                                    
c                                                                               
 100  format(//,  30x,  '***error - plastic strain must be',                    
     1         ' entered in ascending order' )                                  
c                                                                               
      end                                                                       
