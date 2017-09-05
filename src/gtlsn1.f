c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gtlsn1                       *          
c     *                                                              *          
c     *             -- incremental strain for solid elements --      *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 4/30/2016 rhd              *          
c     *                                                              *          
c     *      processes a material (gauss) point for a block of       *          
c     *      identical type solid elements. compute the [B] matrix   *          
c     *      at the point for each element in the block. multiply    *          
c     *      [B] into incremental displacement vector for element    *          
c     *      nodes: (n+1) - n to define a strain increment for       *          
c     *      n->n+1. used for both small and finite strain           *          
c     *      theory. for finite strains, the [B] is evaluated using  *          
c     *      mid-step (n+1/2) element configuration. the strain      *          
c     *      increment is then most often noted D (or the rate D     *          
c     *      * dt).                                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gtlsn1( span, nnode,                                           
     &                   due, deps, gama, nxi, neta,                            
     &                   nzeta, vol_block, bbar, eps_bbar, b )                  
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision ::                                                       
     & due(mxvl,*), deps(mxvl,nstr), gama(mxvl,3,3),                            
     & nxi(*), neta(*), nzeta(*), vol_block(mxvl,8,*), eps_bbar,                
     & b(mxvl,mxedof,*)                                                         
      logical :: bbar                                                           
c                                                                               
c                      locals                                                   
c                                                                               
c                       compute linear strain-displacement                      
c                       [B] matrix for this material (gauss) point at           
c                       all elements in the block. modify                       
c                       for B-bar as needed.                                    
                                                                                
      call blcmp1( span, b, gama, nxi, neta, nzeta, nnode )                     
      if ( bbar ) call bmod( b, vol_block, span, mxvl, eps_bbar,                
     &                       mxedof )                                           
c                                                                               
c                       multiply [B] x displacement increment. this is          
c                       done for all elements in the block for this             
c                       material (gauss) point nuber. take advantage of         
c                       sparsity in [B] during multiply.                        
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      deps = 0.0d00                                                             
c                                                                               
      bpos1 = nnode                                                             
      bpos2 = 2*nnode                                                           
c                                                                               
      do j = 1, nnode                                                           
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
         do i = 1,span                                                          
           deps(i,1) = deps(i,1) +  b(i,j,1) * due(i,j) +                       
     &                              b(i,bpos1+j,1) * due(i,bpos1+j) +           
     &                              b(i,bpos2+j,1) * due(i,bpos2+j)             
           deps(i,2) = deps(i,2) +  b(i,j,2) * due (i,j) +                      
     &                              b(i,bpos1+j,2) * due(i,bpos1+j)+            
     &                              b(i,bpos2+j,2) * due(i,bpos2+j)             
           deps(i,3) = deps(i,3) +  b(i,j,3) * due (i,j) +                      
     &                              b(i,bpos1+j,3) * due(i,bpos1+j)+            
     &                              b(i,bpos2+j,3) * due(i,bpos2+j)             
           deps(i,4) = deps(i,4) +  b(i,j,4)*due(i,j)+                          
     &                              b(i,bpos1+j,4)*due(i,bpos1+j)               
           deps(i,5) = deps(i,5) +  b(i,bpos1+j,5)*due(i,bpos1+j)+              
     &                              b(i,bpos2+j,5)*due(i,bpos2+j)               
           deps(i,6) = deps(i,6) +  b(i,j,6)*due(i,j)+                          
     &                              b(i,bpos2+j,6)*due(i,bpos2+j)               
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gtlsn2                       *          
c     *                                                              *          
c     *         -- incremental "strain" for cohesive elements --     *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                                                              *          
c     *                   last modified : 06/12/99                   *          
c     *                                                              *          
c     *    this subroutine computes the incremental linear           *          
c     *    "strain" of a given gauss point for use in the stress     *          
c     *    recovery routine for an element in a block of similar,    *          
c     *    non-conflicting cohesive elements excluding geometric     *          
c     *    nonlinearity. the accumulated linear strains are updated  *          
c     *    "strains" are just displacement jumps across the          *          
c     *    interface for cohesive elements                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gtlsn2( span, nnode, due, dgstrn, dgstrs, rot, shape,          
     &                   etype, gpn, felem, iout )                              
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision ::                                                       
     & due(mxvl,*), dgstrn(mxvl,nstr), dgstrs(mxvl,*), rot(mxvl,3,3),           
     & shape(*), b(mxvl,mxedof,nstr)                                            
c                                                                               
c                      locals                                                   
c                                                                               
      double precision, parameter :: zero = 0.0d00                              
      logical :: local_debug                                                    
      data local_debug / .false. /                                              
c                                                                               
c         compute linear relative displacement jump                             
c         [b] matrix for this gauss point at all elements in                    
c         the block.                                                            
c                                                                               
      call blcmp_cohes( span, b, rot, shape, etype,  nnode )                    
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      dgstrn = zero                                                             
c                                                                               
      bpos1 = nnode                                                             
      bpos2 = 2*nnode                                                           
c                                                                               
      do j = 1, nnode                                                           
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
           dgstrn(i,1)= dgstrn(i,1)+b(i,j,1) * due(i,j) +                       
     &                              b(i,bpos1+j,1) * due(i,bpos1+j) +           
     &                              b(i,bpos2+j,1) * due(i,bpos2+j)             
           dgstrn(i,2)= dgstrn(i,2) + b(i,j,2) * due (i,j) +                    
     &                                b(i,bpos1+j,2) * due(i,bpos1+j)+          
     &                                b(i,bpos2+j,2) * due(i,bpos2+j)           
           dgstrn(i,3)= dgstrn(i,3)+ b(i,j,3) * due (i,j) +                     
     &                               b(i,bpos1+j,3) * due(i,bpos1+j)+           
     &                               b(i,bpos2+j,3) * due(i,bpos2+j)            
         end do                                                                 
      end do                                                                    
c                                                                               
c                       update the accumulated displacement jumps.              
c                       dgstrs is total "strain" at end of step.                
c                       dgstrn is total "strain" increment over step.           
c                                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span                                                            
         dgstrs(i,1) = dgstrs(i,1) + dgstrn(i,1)                                
         dgstrs(i,2) = dgstrs(i,2) + dgstrn(i,2)                                
         dgstrs(i,3) = dgstrs(i,3) + dgstrn(i,3)                                
         dgstrs(i,4) = zero                                                     
         dgstrs(i,5) = zero                                                     
         dgstrs(i,6) = zero                                                     
      end do                                                                    
c                                                                               
      if ( local_debug .and. gpn .eq. 1 ) then                                  
        write(iout,9000)                                                        
        do i = 1, span                                                          
          write(iout,9100) felem + i - 1                                        
          do j = 1, nnode                                                       
            write(iout,9200) j, due(i,j), due(i,bpos1+j), due(i,bpos2+j)        
          end do                                                                
        end do                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>>>> Debug in gtlsn2 (cohesive element)...')                      
 9100 format('      Element: ',i10)                                             
 9200 format(15x,i3,3f20.10)                                                    
      end                                                                       
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gtlsn3                       *          
c     *                                                              *          
c     *         -- incremental "strain" for bar2 element --          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/29/2017 rhd             *          
c     *                                                              *          
c     *    computes the incremental strain for the bar2 element.     *          
c     *    for small displacements, ce has original coordinates.     *
c     *    for geonl, ce has coords at t = n+1/2                     *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gtlsn3( span, due, deps, strain_np1, ce,
     &                   mxvl, nstr, etype, felem, iout )                              
      implicit none                                                    
c                                                                               
c                      parameter declarations                                   
c                    
      integer :: span, etype, gpn, felem, iout, nstr, mxvl                                                         
      double precision :: due(mxvl,*), deps(mxvl,nstr), 
     &                    strain_np1(mxvl,nstr), ce(mxvl,*)
c                                                                               
c                      locals                                                   
c     
      integer :: i
      double precision :: x1, x2, y1, y2, z1, z2, dx, dy, dz, du_1,
     &                    du_2, len, cos_l, cos_m, cos_n                                                                     
      double precision, parameter :: zero = 0.0d00                              
      logical, parameter :: local_debug = .false.                                              
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      deps = zero                                                             
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
        du_2  = due(i,2) * cos_l +  due(i,4) * cos_m +
     &          due(i,6) * cos_n
        du_1  = due(i,1) * cos_l +  due(i,3) * cos_m +
     &          due(i,5) * cos_n
        deps(i,1) = ( du_2 - du_1 ) / len   
        strain_np1(i,1) = strain_np1(i,1) + deps(i,1)                                                    
      end do
c                                                                               
      if ( local_debug  ) then                                  
        write(iout,9000)                                                        
        do i = 1, span                                                          
          write(iout,9100) felem + i - 1  
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
          du_2  = due(i,2) * cos_l +  due(i,4) * cos_m +
     &            due(i,6) * cos_n
          du_1  = due(i,1) * cos_l +  due(i,3) * cos_m +
     &            due(i,5) * cos_n
          deps(i,1) = ( du_2 - du_1 ) / len   
          write(iout,9200) dx, dy, dz, len
          write(iout,9210) du_2, du_1, ( du_2 - du_1 ) / len         
        end do  
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>>>> Debug in gtlsn3 (bar2 element)...')                      
 9100 format('      Element: ',i10)                                             
 9200 format(15x,'dx, dy, dz, len: ',4f15.6)
 9210 format(15x,'du_2, du_1, deps: ',3f15.6)                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gtlsn3_vols                  *          
c     *                                                              *          
c     *               -- volumes, areas bar2 element --              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/4/2017 rhd              *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gtlsn3_vols( span, mxvl, felem, iout, bar_areas_0,
     &                        bar_areas_nx, ce_0, ce_nx, bar_volumes )
      implicit none
c
      integer :: span, mxvl, iout, felem
      double precision :: bar_areas_0(mxvl), bar_areas_nx(mxvl),
     &                    ce_0(mxvl,*), ce_nx(mxvl,*), bar_volumes(mxvl)  
c
      integer :: i
      logical, parameter :: local_debug = .false.
      double precision :: x1, x2, y1, y2, z1, z2, dx, dy, dz, len0, len
c
c              volume of each bar in undeformed coordinates
c              areas in a deformed configuration t = nx (n+1/2 or n+1)
c              nx be t = 0 for small displacement solutions.
c              assumes incompressibility.
c
      if( local_debug ) write(iout,*) '.. compute bar volumes, areas'         
c
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span
        x1 = ce_0(i,1)
        x2 = ce_0(i,2)
        y1 = ce_0(i,3)
        y2 = ce_0(i,4)
        z1 = ce_0(i,5)
        z2 = ce_0(i,6)
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        len0 = sqrt( dx*dx + dy*dy + dz*dz )
        bar_volumes(i) = len0 * bar_areas_0(i)
        x1 = ce_nx(i,1)
        x2 = ce_nx(i,2)
        y1 = ce_nx(i,3)
        y2 = ce_nx(i,4)
        z1 = ce_nx(i,5)
        z2 = ce_nx(i,6)
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        len = sqrt( dx*dx + dy*dy + dz*dz )
        bar_areas_nx(i) = bar_volumes(i) / len
c        if( local_debug ) 
c     &      write(iout,9000) i+felem-1, len0, len, bar_volumes(i),
c     &                       bar_areas_n1(i)
      end do 
c
      return
 9000 format(10x,'element, len0, len,vol,area(n1):', i8,4f15.6)      
      end
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gtlsn4                       *          
c     *                                                              *          
c     *   -- updated strain and stress for link2 element --          *
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/20/2017 rhd             *          
c     *                                                              *          
c     *    computes updated strain (relative displacement) and       *
c     *    stress (link forces) for small displacements              *
c     *    link element is always linear-elastic                     *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gtlsn4( span, mxelpr, props, due, deps, strain_np1, 
     &                   stress_n, stress_np1, mxvl, nstr, nstrs,
     &                   etype, felem, iout )                               
      implicit none                                                    
c                                                                               
c                      parameter declarations                                   
c                    
      integer :: span, mxelpr, etype, gpn, felem, iout, nstr, nstrs,
     &           mxvl  
      real :: props(mxelpr,span)                                                       
      double precision :: due(mxvl,*), deps(mxvl,nstr), 
     &                    strain_np1(mxvl,nstr), stress_n(mxvl,nstrs),
     &                    stress_np1(mxvl,nstrs)
c                                                                               
c                      locals                                                   
c     
      integer :: i
      double precision ::  du_x, du_y, du_z, kx, ky, kz
      double precision, parameter :: zero = 0.0d0, half = 0.5d0                              
      logical, save :: local_debug = .false.    
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span
        du_x  = due(i,2) - due(i,1)
        du_y  = due(i,4) - due(i,3)
        du_z  = due(i,6) - due(i,5)
        deps(i,1) = du_x
        deps(i,2) = du_y
        deps(i,3) = du_z
        deps(i,4) = zero
        deps(i,5) = zero
        deps(i,6) = zero
        strain_np1(i,1) = strain_np1(i,1) + deps(i,1)                                                    
        strain_np1(i,2) = strain_np1(i,2) + deps(i,2)                                                    
        strain_np1(i,3) = strain_np1(i,3) + deps(i,3)    
        strain_np1(i,4) = zero
        strain_np1(i,5) = zero
        strain_np1(i,6) = zero
        kx = props(7,i)                                                
        ky = props(8,i)                                                
        kz = props(9,i)  
        stress_np1(i,1) = stress_n(i,1) + du_x * kx
        stress_np1(i,2) = stress_n(i,2) + du_y * ky
        stress_np1(i,3) = stress_n(i,3) + du_z * kz
        stress_np1(i,4) = zero
        stress_np1(i,5) = zero
        stress_np1(i,6) = zero
        stress_np1(i,7) = half * ( strain_np1(i,1)*stress_np1(i,1) +
     &                             strain_np1(i,2)*stress_np1(i,2) + 
     &                             strain_np1(i,3)*stress_np1(i,3) )
      end do
c                                                                               
      if ( local_debug  ) then                                  
        write(iout,9000)                                                        
        do i = 1, span                                                          
          du_x  = due(i,2) - due(i,1)
          du_y  = due(i,4) - due(i,3)
          du_z  = due(i,6) - due(i,5)
          kx    = props(7,i)                                                
          ky    = props(8,i)                                                
          kz    = props(9,i)  
          write(iout,9100) felem + i - 1  
          write(iout,9200) du_x, du_y, du_z
          write(iout,9210) kx, ky, kz         
          write(iout,9220) du_x*kx, du_y*ky, du_z*kz
          write(iout,9230) stress_np1(i,1), stress_np1(i,2),
     &                     stress_np1(i,3)
          write(iout,9235) strain_np1(i,1), strain_np1(i,2),
     &                     strain_np1(i,3)
          write(iout,9240) stress_np1(i,7)
        end do  
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>>>> Debug in gtlsn4 (link2 element)...')                      
 9100 format('      Element: ',i10)                                             
 9200 format(15x,'du_x, du_y, du_z: ',3f15.6)
 9210 format(15x,'kx, ky, kz:       ',3f15.1)  
 9220 format(15x,'force increments: ',3f15.1)                                                  
 9230 format(15x,'forces(n+1):      ',3f15.1)      
 9235 format(15x,'rel-dis(n+1):     ',3d15.6)      
 9240 format(15x,'work(n+1):        ',f15.6)      
c
      end                                             
