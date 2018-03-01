c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gpmas1                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/11/2017 rhd             *          
c     *                                                              *          
c     *     this subroutine computes the lumped mass matrix terms    *          
c     *     for elements in the block. uses undeformed coordinates.  *          
c     *     compute contribution just for the specified gauss point. *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gpmas1( span, felem, type, order, gpn, nnode,                  
     &                   volume, mel, totdof, fgm_props, area, rho,            
     &                   emass, iout, rho_fgm_flags, axisymm_elem,                     
     &                   implemented_elem, ce_block, rho_block )
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, felem, type, order, gpn, nnode, totdof, iout                                             
      double precision ::                                                       
     &  volume(*), mel(totdof,*), rho(*), emass(*), rho_fgm_flags(*),           
     &  ce_block(mxvl,*), rho_block(mxndel,*), area(*)
      logical :: fgm_props, axisymm_elem, implemented_elem                      
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i, j                                                           
      double precision ::                                                       
     &   shape(mxndel), radius(mxvl), dummy(mxvl,3,3),                          
     &   jac(mxvl,ndim,ndim),  gama(mxvl,ndim,ndim),                            
     &   xi, eta ,zeta, wfactor, dj(mxvl),                                      
     &   nxi(mxndel), neta(mxndel), nzeta(mxndel)
      double precision, parameter :: one = 1.0d0                           
c                                                                               
c                       message for unsupported elements                        
c  
      if( .not. implemented_elem ) then                                         
        write(iout,9500) type,span,felem,order,gpn,nnode,totdof                 
        do j = 1, nnode                                                         
          do i = 1, span                                                        
            mel(j,i) = one                                                      
          end do                                                                
        end do                                                                  
        do i = 1, span                                                          
          volume(i) = one                                                       
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c                       compute the shape functions and the                     
c                       shape function derivatives.                             
c                                                                               
      call getgpts( type, order, gpn, xi, eta, zeta, wfactor )                  
      call derivs(  type, xi, eta, zeta, nxi, neta, nzeta )                     
      call shapef(  type, xi, eta, zeta, shape )                                
c                                                                               
c                       for axisymmetric elements, compute the radius           
c                       to the current gauss point for each element             
c                       in the block.                                           
c                                                                               
      if( axisymm_elem ) then                                                   
       call get_radius( radius, nnode, span, shape, ce_block, mxvl )            
      end if                                                                    
c                                                                               
c                       compute the determinate of the coordinate               
c                       jacobian (use initial coordinates).                     
c                                                                               
      call jacob1( type, span, felem, gpn, jac, dj, gama, dummy,                
     &             nxi, neta, nzeta, ce_block, nnode )                          
c                                                                               
c                       some elements may have the mass density                 
c                       defined using the fgm capabilities, i.e,                
c                       thru nodal values. interpolate density                  
c                       at this gauss point for all elements in                 
c                       the block.                                              
c                                                                               
      if( fgm_props )                                                          
     &    call set_fgm_solid_props_for_block(                                   
     &                span, felem, type, gpn, nnode, rho, shape,                
     &                rho_block, 1, rho_fgm_flags )                             
c                                                                               
c                       compute the gauss point contribution to the             
c                       the nodal mass term.                                    
c                                                                               
      call mass_sum( type, totdof, nnode, mxvl, span, radius, shape,                  
     &               dj, wfactor, ce_block, mel, volume, area, rho,
     &               emass )              
c                                                                               
c                                                                               
      return                                                                    
c                                                                               
9500    format(1x,//,'>>>>>  in gpmas1, Mass matrix not yet defined',/,         
     &               '       for element type = ',i3,/,                         
     &               '   Skipping mass matrix computation, set',/,              
     &               '   unit values for now.',//,                              
     &               '     span =',i7,/,                                        
     &               '    felem =',i7,/,                                        
     &               '    order =',i7,/,                                        
     &               '      gpn =',i7,/,                                        
     &               '    nnode =',i7,/,                                        
     &               '   totdof =',i6)                                          
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine mass_sum                      *          
c     *                                                              *          
c     *                       written by: rau                        *          
c     *                    last modified: 8/21/2017 rhd              *          
c     *                                                              *          
c     *     this routine computes the current gauss point            *          
c     *     contribution to the terms of the consistent              *          
c     *     diagonal mass matrix for one (any) dof. the              *          
c     *     matrix will be scaled later to include mass              *          
c     *     density.                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         for triangular elements, the correct area is given by                 
c         0.5*|J|. for tetrahedral elements, the correct volume is              
c         given by (1/6)*|J|.                                                   
c                                                                               
c         For axisymmetric elements, include the 2*pi*radius scalar term        
c         in the element stiffness summation.                                   
c                                                                               
c                                                                               
c         Variables:                                                            
c                                                                               
c         elem_type = integer flag for the element type                         
c         totdof    = total number of degrees of freedom,                       
c                     dimension for mel()                                       
c         mxvl      = maximum number of elements in a block                     
c         span      = number of elements in the block                           
c         radius()  = radius to the current Gauss point for each element        
c                     in the current block.                                     
c         shape()   = shape function array evaluated at current Gauss point     
c         dj()      = Jacobian determinant for each element in the block        
c         w         = Gauss integration weight                                  
c         mel()     = mass contribution at a node for each element              
c                     in the block at the current Gauss integration point       
c         volume()  = volume contribution at a node for each element in the     
c                     block at the current Gauss integration point              
c                                                                               
c                                                                               
      subroutine mass_sum( elem_type, totdof, nnode, mxvl, span, radius,        
     &                     shape, dj, w, ce_0, mel, volume, area, 
     &                     rho, emass )              
      implicit none                                                             
c                                                                               
c               parameter declarations                                          
c                                                                               
      integer :: elem_type, totdof, nnode, mxvl, span                           
      double precision :: radius(*), shape(*), dj(*), w, mel(totdof,*),                                                   
     &         volume(*), rho(*), emass(*), ce_0(mxvl,*), area(*)
c                                                                               
c               local variables                                                 
c                                                                               
      integer :: i,j                                                            
      double precision :: scalar, weight, x1, x2, y1, y2, z1, z2, len
      double precision, parameter :: two = 2.0d0, half = 0.5d0,
     &                               pi =  3.14159265358979323846d0                                                          
c                                                                               
c               the element type determines the calculations                    
c               performed. compute the current mass and volume                  
c               contributions for the different element types.                  
c                                                                               
      select case( elem_type )                           
c                                                                               
c               element numbers 1-5, 7, 9, 12: hex elements                     
c                                                                               
       case( 1,2,3,4,5, 7, 9, 12 )                                                                  
c                                                                               
        do j = 1, nnode                                                           
!DIR$ VECTOR ALIGNED
          do i = 1, span                                                         
             mel(j,i) = mel(j,i) + shape(j)*shape(j)*dj(i)*w*rho(i)               
         end do                                                                 
        end do                                                                    
c                                                                               
!DIR$ VECTOR ALIGNED
        do i = 1, span                                                            
          volume(i) = volume(i) + dj(i)*w                                        
          emass(i) = emass(i) + dj(i)*w*rho(i)                                   
        end do                                                                    
c                                                                               
c               element number 6, 13: tet10, tet4                               
c                                                                               
       case( 6, 13 )                                                                 
c                                                                               
         scalar = 1.0D0 / 6.0D0                                                     
         weight = scalar * w                                                        
         do j = 1, nnode                                                          
!DIR$ VECTOR ALIGNED
          do i = 1, span                                                        
            mel(j,i) = mel(j,i) + shape(j)*shape(j)*dj(i)*weight*rho(i)         
          end do                                                                
         end do                                                                   
c                                                                               
!DIR$ VECTOR ALIGNED
         do i = 1, span                                                           
          volume(i) = volume(i) + dj(i)*weight                                  
          emass(i) = emass(i) + dj(i)*weight*rho(i)                             
         end do                                                                   
c
c               element number 8: tri6                                          
c                                                                               
        case( 8 )                                                                
c                                                                               
         scalar = half                                                            
         weight = scalar*w                                                        
         do j = 1, nnode                                                          
!DIR$ VECTOR ALIGNED
          do i = 1, span                                                        
            mel(j,i) = mel(j,i) + shape(j)*shape(j)*dj(i)*weight*rho(i)         
          end do                                                                
         end do                                                                   
c                                                                               
!DIR$ VECTOR ALIGNED
         do i = 1, span                                                           
          volume(i) = volume(i) + dj(i)*weight                                  
          emass(i)  = emass(i) + dj(i)*weight*rho(i)                             
         end do                                                                   
c                                                                               
c               element numbers 10: axisymmetric quad                           
c                                                                               
        case( 10 ) 
c                                                                               
c               axisymmetric elements, multiply by 2*pi*radius                  
c                                                                               
         scalar = two * pi                                                         
         do j = 1, nnode                                                          
          do i = 1, span                                                        
            mel(j,i) = mel(j,i) +                                               
     &              scalar*radius(i)*shape(j)*shape(j)*dj(i)*w*rho(i)        
          end do                                                                
         end do                                                                   
         do i = 1, span                                                           
          volume(i) = volume(i) + scalar*radius(i)*dj(i)*w                      
          emass(i) = emass(i) + dj(i)*w*rho(i)                                  
         end do                                                                   
c                                                                               
c               element number 11: axisymmetric tri                             
c               note that for the axisymmetric triangle 2*pi*0.5 = pi.          
c                                                                               
        case( 11 )
c                                                                               
         scalar = pi                                                              
         do j = 1, nnode                                                          
          do i = 1, span                                                        
            mel(j,i) = mel(j,i) +                                               
     &               scalar*radius(i)*shape(j)*shape(j)*dj(i)*w*rho(i)        
          end do                                                                
         end do                                                                   
         do i = 1, span                                                           
          volume(i) = volume(i) + scalar*radius(i)*dj(i)*w                      
          emass(i)  = emass(i) + dj(i)*w*rho(i)                                 
         end do  
c                                                                               
c               element number 18: 2-node element
c                                                                               
       case( 18 )
c              
         nnode = 2                                                                 
         do j = 1, nnode                                                 
!DIR$ VECTOR ALIGNED
          do i = 1, span   
            x1 = ce_0(i,1); y1 = ce_0(i,nnode+1)
            z1 = ce_0(i,2*nnode+1)      
            x2 = ce_0(i,2); y2 = ce_0(i,nnode+2)
            z2 = ce_0(i,2*nnode+2)      
            len = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )                                               
            mel(j,i) = mel(j,i) + rho(i) * area(i) * len * half    
          end do                                                                
         end do                                                                   
!DIR$ VECTOR ALIGNED
         do i = 1, span                                                           
          volume(i) = area(i) * len                   
          emass(i)  = emass(i) + area(i) * len * rho(i)                                 
         end do    
c                                                                               
c               element number 19: 2-node link element. has no 
c               volume. user set rho is total element mass
c                                                                               
       case( 19 )
c              
         nnode = 2                                                                 
         do j = 1, nnode                                                 
          do i = 1, span   
            mel(j,i) = mel(j,i) + rho(i) * half    
          end do                                                                
         end do
         do i = 1, span                                                           
          volume(i) = 0.0d0                   
          emass(i)  = emass(i) + rho(i)
         end do    

      case default 
         write(*,9000) elem_type
         call die_abort
      end select
c                                                                           
      return 
 9000 format('>> FATAL ERROR: mass_sum. no element type:',i8,
     & /,    '                job terminated' )                                                                         
       end                                                                      
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine get_radius                    *          
c     *                                                              *          
c     *                       written by: gvt                        *          
c     *                   last modified : 09/28/2015 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c         for axisymmetric elements, compute the radius to the current          
c         gauss point for each element in the block.  The radius is             
c         used for the hoop strain = u/r, and in the summation of the           
c         element stiffness matrix. Use the first nnode values in ce()          
c         for x=radial coordinates.                                             
c                                                                               
c         Variables:                                                            
c                                                                               
c          nnode = number of nodes                                              
c          span  = number of elements in the block                              
c          mxvl  = maximum number of elements in a block                        
c          rad() = radius to the current Gauss point for each element           
c                  in the current block. Using the shape functions              
c                  and node coordinates, the radius is given by:                
c                  r = Sum[n(i)*x(i)], i=1,nnode where n(i) =                   
c                  shape functions, x(i) = x=radial coordinate                  
c          n()   = element shape functions, already evaluated at the curent     
c                  gauss integration point                                      
c          ce()  = matrix of node coordinates, the first nnode columns          
c                  contain the x=radial coordinates                             
c                                                                               
      subroutine get_radius( rad, nnode, span, n, ce, mxvl )                    
      implicit none                                                             
      integer :: nnode, span, mxvl                                              
      double precision ::                                                       
     &         rad(*), n(*), ce(mxvl,*)                                         
c                                                                               
      integer i,j                                                               
      double precision                                                          
     &         zero                                                             
      data zero / 0.0d00 /                                                      
c                                                                               
      do i = 1, span                                                            
         rad(i) = zero                                                          
      end do                                                                    
c                                                                               
      do j = 1, nnode                                                           
          do i = 1, span                                                        
            rad(i) = rad(i) + n(j)*ce(i,j)                                      
          end do                                                                
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                     Subroutine adjust_cnst                   *          
c     *                                                              *          
c     *                       written by: gvt                        *          
c     *                   last modified : 08/25/98                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c         Include other scalars with the linear material matrix required for    
c         axisymmetric or planar elements.                                      
c                                                                               
c         For triangular elements, the correct area is given by                 
c         0.5*|J|.                                                              
c                                                                               
c         For axisymmetric elements, include the 2*pi*radius scalar term        
c         in the element stiffness summation.                                   
c                                                                               
c         Determine the correct scalar value and then multiply all the          
c         terms in the material matrix [D] = cep_block.                         
c                                                                               
c         Variables:                                                            
c                                                                               
c         elem_type = integer flag for the element type, should get one of      
c                     8=tri6, 10=axiquad8, 11=axitri6                           
c         nstr  = number of strains, should be 6                                
c         mxvl  = maximum number of elements in a block                         
c         span  = number of elements in the block                               
c         rad() = radius to the current Gauss point for each element in the     
c                 current block.                                                
c         cep() = material matrix, should already have been updated, extra      
c                 scalar values in the Gauss integration                        
c                 loop are included here into the material matrix.              
c                                                                               
      subroutine adjust_cnst( elem_type, nstr, mxvl, span, rad, cep )           
      implicit none                                                             
c                                                                               
c                  parameters                                                   
c                                                                               
      integer elem_type, nstr, mxvl, span                                       
      double precision                                                          
     &         rad(*),                                                          
     &         cep(mxvl,nstr,*)                                                 
c                                                                               
c                   local variables                                             
c                                                                               
      integer i, row                                                            
c                                                                               
      double precision                                                          
     &         scalar                                                           
c                                                                               
c                    the element type determines the scalar                     
c                    multiple to get the correct adjustment to |J|.             
c                                                                               
      call adjust_scalar_weights( elem_type, scalar )                           
      do row = 1, nstr                                                          
          do i = 1, span                                                        
              cep(i,row,1) = scalar*cep(i,row,1)                                
              cep(i,row,2) = scalar*cep(i,row,2)                                
              cep(i,row,3) = scalar*cep(i,row,3)                                
              cep(i,row,4) = scalar*cep(i,row,4)                                
              cep(i,row,5) = scalar*cep(i,row,5)                                
              cep(i,row,6) = scalar*cep(i,row,6)                                
          end do                                                                
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
