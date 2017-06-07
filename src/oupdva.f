c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupdva                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *               last modified : 4/1/2017 rhd                   *          
c     *                                                              *          
c     *     write nodal displacements, velocities, accelerations,    *          
c     *     reactions, or temperatures                               *          
c     *     to (1) patran in either binary or formatted forms;       *          
c     *     or (2) a flat file in text or stream (binary)            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine oupdva( dva, oubin, ouasc, flat_file, stream_file,             
     &                   text_file, compressed  )                               
      use global_data ! old common.main
c                                                                               
      use main_data, only : trn, trnmat, inverse_incidences                     
      implicit none                                                             
c                                                                               
      integer :: dva                                                            
      logical :: oubin, ouasc, flat_file, stream_file,                          
     &           text_file, compressed                                          
c                                                                               
c                       local declarations                                      
c                                                                               
      integer :: nodmax, nod, elem, ndof, dof, i, nwidth                        
      double precision, allocatable :: edva(:,:), trnmte(:,:,:)                 
      double precision :: defmax                                                
      double precision, parameter :: zero = 0.0d0                               
      logical, allocatable :: trne(:,:)                                         
c                                                                               
c                       if it is displacements to be output, then a             
c                       pass must be made over all nodes to determine           
c                       the maximum absolute displacement and the               
c                       node at which it occurs.                                
c                                                                               
      allocate( edva(mxvl,mxndof), trnmte(mxvl,mxedof,mxndof) )                 
      allocate( trne(mxvl,mxndel) )                                             
c                                                                               
      nodmax = 0                                                                
      defmax = zero                                                             
c                                                                               
      if( dva .eq. 1 ) then                                                     
         defmax = zero                                                          
         nodmax = 0                                                             
         do nod = 1, nonode                                                     
            elem = inverse_incidences(nod)%element_list(1)                      
            ndof = iprops(4,elem)                                               
            do dof = 1, ndof                                                    
               edva(1,dof) = u(dstmap(nod)+dof-1)                               
            end do                                                              
c                                                                               
c                       transform displacement vector at the current            
c                       node to uniform global coordinates.                     
c                                                                               
            trne(1,1) = trn(nod)                                                
            if( trne(1,1) ) then                                                
               trnmte(1,1:3,1:3) = trnmat(nod)%mat(1:3,1:3)                     
               call trnvec( edva, trnmte, trne, ndof, 1, 1, 2 )                 
            end if                                                              
c                                                                               
c                       loop over the displacement vector to check for          
c                       the maximum absolute displacement.                      
c                                                                               
            do i = 1, ndof                                                      
              if( abs(edva(1,i)) .gt. defmax ) then                             
                  defmax = abs(edva(1,i))                                       
                  nodmax = nod                                                  
               end if                                                           
            end do                                                              
         end do  !  on nodes                                                    
      end if  !  on dva == 1                                                    
c                                                                               
c                       set the width of nodal output records for               
c                       dis/vel/acc/reactions or temperatures.                  
c                                                                               
      nwidth = mxndof                                                           
      if( dva .eq. 5 ) nwidth = 2  ! temperatures                               
c                                                                               
c                       output the nodal results for dis/vel/acc/               
c                       reactions/temperatures to                               
c                       (1) patran results file or (2) flat                     
c                       file                                                    
c                                                                               
      call ouddpa( dva, oubin, ouasc, nodmax, defmax, nwidth,                   
     &             flat_file, stream_file, text_file, compressed )              
c                                                                               
      return                                                                    
      end                                                                       
