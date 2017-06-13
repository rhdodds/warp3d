c                                                                               
c           user_routines_support.f   Distribution version                      
c                                                                               
c           Updated:  6/13/2017 rhd                                            
c                                                                               
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *         integer function warp3d_get_device_number()          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/26/12                    *          
c     *                                                              *          
c     *     find a non-connected Fortran device number               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      integer function warp3d_get_device_number()                               
      implicit none                                                             
      integer :: iunit                                                             
      logical :: connected                                                         
c                                                                               
c        1.  find an available device number, open the neutral                  
c            file                                                               
c                                                                               
      warp3d_get_device_number = -1                                             
      do iunit = 11, 5000                                                       
        inquire( unit=iunit, opened=connected )                                 
        if ( connected ) cycle                                                  
        warp3d_get_device_number = iunit                                        
        return                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c *********************************************************************         
c **                                                                 **         
c **              Abaqus Compatible Support Routines                 **         
c **                                                                 **         
c *********************************************************************         
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine xit  (called by UMATs)            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 3/22/12                     *          
c     *                                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine xit                                                            
      implicit integer (a-z)                                                    
c                                                                               
      write(*,*) '>>> UMAT called to abort execution'                           
      call die_abort                                                            
      stop                                                                      
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine rotsig                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 03/29/12                       *          
c     *                                                              *          
c     *     tensor rotation routine for Abaqus compatible UMAT       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine rotsig( in_vec, drot, out_vec, type, nrow, ncol )              
      implicit none                                                             
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision :: in_vec(6), drot(3,3), out_vec(6)                                        
      integer :: type, nrow, ncol                                                  
c                                                                               
c                     locally defined arrays-variables                          
c                                                                               
      double precision :: factor, one, half, two, a(3,3), t(3,3), 
     &                    c(3,3)                          
      logical :: local_debug                                                       
      data one, half, two, local_debug / 1.0d00, 0.5d00,                        
     &                                       2.0d00, .true. /                   
c                                                                               
c                                                                               
c                     out = drot * in * trans(drot)                             
c                                                                               
c                     put input tensor into 3x3 matrix form.                    
c                     type = 1 input tensor has engr strain                     
c                     type = 2 input tensor has stress                          
c                     put into 3 x 3 form                                       
c                     Abaqus umat ordering: x,y,z,xy,xz,yz                      
c                                                                               
      factor = one                                                              
      if( type .eq. 1 ) factor = half                                           
c                                                                               
      a(1,1) = in_vec(1)                                                        
      a(2,1) = in_vec(4) * factor                                               
      a(3,1) = in_vec(5) * factor                                               
      a(1,2) = a(2,1)                                                           
      a(2,2) = in_vec(2)                                                        
      a(3,2) = in_vec(6) * factor                                               
      a(1,3) = a(3,1)                                                           
      a(2,3) = a(3,2)                                                           
      a(3,3) = in_vec(3)                                                        
c                                                                               
c                     t = a * trans(drot)                                       
c                                                                               
      t(1,1) = a(1,1)*drot(1,1) + a(1,2)*drot(1,2)  + a(1,3)*drot(1,3)          
      t(2,1) = a(2,1)*drot(1,1) + a(2,2)*drot(1,2)  + a(2,3)*drot(1,3)          
      t(3,1) = a(3,1)*drot(1,1) + a(3,2)*drot(1,2)  + a(3,3)*drot(1,3)          
c                                                                               
      t(1,2) = a(1,1)*drot(2,1) + a(1,2)*drot(2,2)  + a(1,3)*drot(2,3)          
      t(2,2) = a(2,1)*drot(2,1) + a(2,2)*drot(2,2)  + a(2,3)*drot(2,3)          
      t(3,2) = a(3,1)*drot(2,1) + a(3,2)*drot(2,2)  + a(3,3)*drot(2,3)          
c                                                                               
      t(1,3) = a(1,1)*drot(3,1) + a(1,2)*drot(3,2)  + a(1,3)*drot(3,3)          
      t(2,3) = a(2,1)*drot(3,1) + a(2,2)*drot(3,2)  + a(2,3)*drot(3,3)          
      t(3,3) = a(3,1)*drot(3,1) + a(3,2)*drot(3,2)  + a(3,3)*drot(3,3)          
c                                                                               
c                     c = drot * t                                              
c                                                                               
      c(1,1) = drot(1,1)*t(1,1) + drot(1,2)*t(2,1) + drot(1,3)*t(3,1)           
      c(2,1) = drot(2,1)*t(1,1) + drot(2,2)*t(2,1) + drot(2,3)*t(3,1)           
      c(3,1) = drot(3,1)*t(1,1) + drot(3,2)*t(2,1) + drot(3,3)*t(3,1)           
c                                                                               
      c(1,2) = drot(1,1)*t(1,2) + drot(1,2)*t(2,2) + drot(1,3)*t(3,2)           
      c(2,2) = drot(2,1)*t(1,2) + drot(2,2)*t(2,2) + drot(2,3)*t(3,2)           
      c(3,2) = drot(3,1)*t(1,2) + drot(3,2)*t(2,2) + drot(3,3)*t(3,2)           
c                                                                               
      c(1,3) = drot(1,1)*t(1,3) + drot(1,2)*t(2,3) + drot(1,3)*t(3,3)           
      c(2,3) = drot(2,1)*t(1,3) + drot(2,2)*t(2,3) + drot(2,3)*t(3,3)           
      c(3,3) = drot(3,1)*t(1,3) + drot(3,2)*t(2,3) + drot(3,3)*t(3,3)           
c                                                                               
      factor = one                                                              
      if( type .eq. 1 ) factor = two                                            
      out_vec(1) = c(1,1)                                                       
      out_vec(2) = c(2,2)                                                       
      out_vec(3) = c(3,3)                                                       
      out_vec(4) = c(2,1) * two                                                 
      out_vec(5) = c(3,1) * two                                                 
      out_vec(6) = c(3,2) * two                                                 
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine sinv                         *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 05/14/12                       *          
c     *                                                              *          
c     *              invariants 1 & 2 of 3D stress tensor            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine sinv( sig, sinv1, sinv2, ndi, nshr )                           
      implicit none                                                             
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision :: sig(*), sinv1, sinv2                                                    
      integer :: ndi, nshr                                                         
c                                                                               
c                     locally defined arrays-variables                          
c                                                                               
      double precision :: sig_dev(6), t1, t2, one_third, two, oneptfive                           
      data one_third, two, oneptfive                                            
     &    /  0.3333333333333333d00, 2.0d00, 1.5d00 /                            
c                                                                               
      sinv1 = one_third * ( sig(1) + sig(2) + sig(3) )                          
c                                                                               
      sig_dev(1) = sig(1) - sinv1                                               
      sig_dev(2) = sig(2) - sinv1                                               
      sig_dev(3) = sig(3) - sinv1                                               
      sig_dev(4) = sig(4)                                                       
      sig_dev(5) = sig(5)                                                       
      sig_dev(6) = sig(6)                                                       
c                                                                               
      t1 = sig_dev(1)**2 + sig_dev(2)**2 + sig_dev(3)**2                        
      t2 = sig_dev(4)**2 + sig_dev(5)**2 + sig_dev(6)**2                        
      sinv2 = sqrt( oneptfive * (t1 + two * t2 ) )                              
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine sprinc                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/14/12                    *          
c     *                                                              *          
c     *     compute principal strain or stresses compatible with     *          
c     *     UMAT specifications                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine sprinc( s, ps, lstr, ndi, nshr )                               
      implicit none  
      integer :: lstr, ndi, nshr, ier                                                 
      double precision :: s(*), ps(*)                                                          
c                                                                               
c                    locally allocated                                          
c                                                                               
      double precision :: temp(6), wk(3), evec(3,3), half, one, factor                            
c                                                                               
      data  half, one                                                           
     &   / 0.5d00, 1.0d00 /                                                     
c                                                                               
c        calculate the principal strains or stresses for UMAT support.          
c        Abaqus stress/strain ordering: xx, yy, zz, xy, xz, yz                  
c        lstr = 1 (stresses), lstr = 2 (strains)                                
c                                                                               
      factor = one                                                              
      if( lstr .eq. 2 ) factor = half                                           
c                                                                               
      temp(1) = s(1)                                                            
      temp(2) = s(4) * factor                                                   
      temp(3) = s(2)                                                            
      temp(4) = s(5) * factor                                                   
      temp(5) = s(6) * factor                                                   
      temp(6) = s(3)                                                            
      call ou3dpr( temp, 3, 0, ps, evec, 3, wk, ier )   ! warp3d routine        
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine sprind                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/14/12                    *          
c     *                                                              *          
c     *     compute principal strain or stresses compatible with     *          
c     *     UMAT specifications                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine sprind( s, ps, an, lstr, ndi, nshr )                           
      implicit none
      integer :: lstr, ndi, nshr, ier                                                    
      double precision :: s(*), ps(*), an(3,3)                                                 
c                                                                               
c                    locally allocated                                          
c                                                                               
      double precision :: temp(6), wk(3), half, one, factor                                       
c                                                                               
      data  half, one  / 0.5d00, 1.0d00 /                                                     
c                                                                               
c        calculate the principal strains or stresses and eigenvectors           
c        for UMAT support.                                                      
c        Abaqus stress/strain ordering: xx, yy, zz, xy, xz, yz                  
c        lstr = 1 (stresses), lstr = 2 (strains)                                
c                                                                               
      factor = one                                                              
      if( lstr .eq. 2 ) factor = half                                           
c                                                                               
      temp(1) = s(1)                                                            
      temp(2) = s(4) * factor                                                   
      temp(3) = s(2)                                                            
      temp(4) = s(5) * factor                                                   
      temp(5) = s(6) * factor                                                   
      temp(6) = s(3)                                                            
      call ou3dpr( temp, 3, 1, ps, an, 3, wk, ier )   ! warp3d routine          
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine getnumcpus                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/13/2017                  *          
c     *                                                              *          
c     *     return the number of ranks during MPI execution          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getnumcpus( nranks )                                           
      use global_data ! old common.main
      implicit none
      integer :: nranks                                                   
c                                                                               
      nranks = 1                                                                
      if( use_mpi ) nranks = numprocs                                           
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine getrank                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/13/2017                  *          
c     *                                                              *          
c     *     return the MPI rank number for this process              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getrank( thisrank )                                            
      use global_data ! old common.main
      implicit none
      integer :: thisrank                                                   
c                                                                               
      thisrank = 0                                                              
      if( use_mpi ) thisrank = myid                                             
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine getmodelsizes                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/13/2017                  *          
c     *                                                              *          
c     *     return number of nodes and elements                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine getmodelsizes( num_model_nodes, num_model_elements )           
      use global_data ! old common.main
      implicit none
      integer :: num_model_nodes, num_model_elements                                                    
c                                                                               
      num_model_nodes    = nonode                                               
      num_model_elements = noelem                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
