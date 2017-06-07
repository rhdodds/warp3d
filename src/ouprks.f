c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouprks                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 3/3/2016 rhd                *          
c     *                                                              *          
c     *   drives stress or strain output quantities for a block of   *          
c     *   similar elements in preparation for patran or flat         *          
c     *   results file                                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ouprks( span, blk, felem, type, order, ngp,                    
     &                   nnode, geonl,                                          
     &                   do_stresses, mat_type, center_output,                  
     &                   num_short_stress, num_short_strain,                    
     &                   element_output )                                       
      use global_data ! old common.main
      use elblk_data, only : elestr                                             
c                                                                               
      implicit none                                                             
c                                                                               
      integer :: span, blk, felem, type, order, ngp, nnode, mat_type,           
     &           num_short_stress, num_short_strain                             
      logical :: do_stresses, geonl, center_output, element_output              
c                                                                               
      integer :: matnum, kout, nowstep                                          
      double precision :: nowtime                                               
c                                                                               
c                                                                               
c                       technically, the output configuration flag              
c                       is not a criteria for similarity. but for               
c                       patran output to make any sense, all ele-               
c                       ments must have the same output config-                 
c                       uration. thus, the flag is being treated                
c                       like a similarity criteria.                             
c                                                                               
c                                                                               
c                       get the element or strains                              
c                       into a final form.                                      
c                       for geonl: the cauchy stresses are computed             
c                       from the unrotated cauchy stresses; the                 
c                       strains are accumulated spatial deformation             
c                       increments.                                             
c                                                                               
      matnum = iprops(38,felem)  ! common.main                                  
      kout = out                 !   "                                          
      nowtime = total_model_time !   "                                          
      nowstep = ltmstp           !   "                                          
c                                                                               
      call ougts1( span, blk, felem, do_stresses, ngp,                          
     &             geonl, mat_type, matnum, kout, nowtime, nowstep )            
c                                                                               
c                                                                               
c                       for output at element center, make                      
c                       all gauss points have the average values.               
c                       Note: mises stress and equivalent strain                
c                       will be replaced with values computed from              
c                       the averaged computavalues later.                       
c                                                                               
      if( element_output ) call oumkcv( span, ngp, do_stresses,                 
     &                     num_short_stress + 1, num_short_strain + 1 )         
c                                                                               
c                       if we are generating nodal results,                     
c                       extrapolate gauss point results to element              
c                       nodes. the 6 primary components, energy density,        
c                       mises equiv. stress and material specific               
c                       output values are extrapolated from gauss               
c                       points to nodes. for strains, we only                   
c                       use the 6 extrapolated components and the equivalent    
c                       strain. Note: the mises stress and equivalent strain    
c                       extrapolated to nodes will be subsequently overwritten  
c                       with values computed from the extrapolated and          
c                       averaged components of the stress or strain tensor.     
c                       this will keep extrapolated mises values and            
c                       equivalent strains from ever becoming negative.         
c                                                                               
      if( .not. element_output ) call ounds1( span, type, order, nnode,         
     &   ngp, do_stresses, num_short_stress + 1, num_short_strain + 1 )         
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
