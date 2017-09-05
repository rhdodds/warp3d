c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine rknmas                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/11/2017 rhd             *          
c     *                                                              *          
c     *     compute effective diagonal mass matrix for block of      *          
c     *     elements. consistent mass scaled to diagonal             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine rknmas( span, felem, type, order, ngp, nnode, totdof,          
     &                   props, mel, beta_fact, totvol, fgm_props,              
     &                   cohesive_elem, linear_displ_elem,                      
     &                   axisymm_elem, implemented_elem, ce_block,              
     &                   rho_block, iout )                                            
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c              parameter declarations                                           
c                                                                               
      integer :: span, felem, type, order, ngp, nnode, totdof, iout                  
      double precision :: mel(totdof,*), beta_fact, totvol, 
     &                    ce_block(mxvl,*), rho_block(mxndel,*)                                                      
      real :: props(mxelpr,*)                                                   
      logical :: fgm_props, cohesive_elem, linear_displ_elem,                   
     &           axisymm_elem, implemented_elem                                 
c                                                                               
c              locals                                                
c                                                                               
      integer :: i, j, gpn                                               
      double precision, allocatable, dimension(:) :: rho, volume,
     &         emass, rho_fgm_flags, area
      double precision, parameter :: zero=0.0d0, one=1.0d0                                                     
c                                                                               
c              initialize block of mass matrices                                
c              initialize element volumes; element masses.                      
c              set default mass densities.                                      
c                                                                               
c              local_work%fgm_flags(i,j) = -99 indicates that                   
c              element in block has an fgm material property                    
c              requiring interpolation.                                         
c   
      allocate(rho(mxvl), volume(mxvl), emass(mxvl),                         
     &                    rho_fgm_flags(mxvl), area(mxvl) )
c                                                                                  
!DIR$ VECTOR ALIGNED                                                            
      mel(1:totdof,1:span)= zero                                                
!DIR$ VECTOR ALIGNED
      do i = 1, span                                                            
        volume(i)        = zero                                                 
        rho(i)           = props(10,i)                                          
        rho_fgm_flags(i) = props(10,i)                                          
        emass(i)         = zero 
        area(i)          = props(43,i)  
      end do  
c                                                                               
c              for linear displacement elements, nodal densities are            
c              averaged to be consistent with the averaging                     
c              of nodal temperature, e, nu, and alpha values.                   
c                                                                               
      if( linear_displ_elem .and. fgm_props ) then                              
         call average_fgm_properties( rho_block, mxndel, mxvl, 1,               
     &                                nnode, span )                             
      end if                                                                    
c                                                                               
c              an interface element has zero mass and zero volume               
c                                                                               
      if( cohesive_elem ) return                                                
c                                                                               
c              compute lumped mass terms from diagonals of the                  
c              consistent mass for each element of the block.                   
c              each dof at an element node has the same lumped                  
c              mass. this routine just computes the mass value                  
c              one intergration point at a time.                                
c       
      do gpn = 1, ngp                                                           
        call gpmas1( span, felem, type, order, gpn, nnode, volume,              
     &               mel, totdof, fgm_props, area, rho, emass, iout,                  
     &               rho_fgm_flags, axisymm_elem, implemented_elem,             
     &               ce_block, rho_block )                                      
      end do                                                                    
c                                                                               
c              set mass value at other dof v, w for each element node           
c                                                                               
      call elmas1( span, nnode, emass, volume, mel, totdof )                    
c                                                                               
c                             multiply by thickness ratio if not equal to 1     
c                             this ratio is really only for 2d work             
c                                                                               
      if( beta_fact .ne. one )  mel(1:totdof,1:span) =
     &                          beta_fact * mel(1:totdof,1:span)
c                                                                               
c              add volumes of elements in this block to total                   
c              volume of the model                                              
c                                                                                     
      totvol = totvol + sum( volume(1:span) )
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
