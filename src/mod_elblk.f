c     ****************************************************************          
c     *                                                              *          
c     *                     f-90 data modules:                       *          
c     *                        elblk_data                            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 05/20/2017 rhd             *          
c     *                                                              *          
c     *     define the data structures for a block of elements       *          
c     *     used in various places.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      module elblk_data                                                         
c     -----------------                                                         
      implicit none                                                             
c                                                                               
      include 'param_def'                                                       
c                                                                               
c                stresses, strains, histories, etc. for a block                 
c                                                                               
c                alignment needed only for non-allocatables                     
c                                                                               
!dir$ attributes align: 64 :: ddtse, urcs_blk_n, urcs_blk_n1, rot_blk_n1        
      double precision, save :: ddtse(mxvl,nstr,mxgp),                          
     &   urcs_blk_n(mxvl,nstrs,mxgp), urcs_blk_n1(mxvl,nstrs,mxgp),             
     &   rot_blk_n1(mxvl,9,mxgp)                                                
c                                                                               
                                                                                
      double precision, dimension(:,:,:), save, allocatable ::                  
     &        elem_hist, elem_hist1                                             
c                                                                               
c                internal forces, transformations for block of elements         
c                                                                               
!dir$ attributes align: 64 :: elestr                                            
      double precision, save :: elestr(mxvl,mxoupr,mxoupt)                      
c                                                                               
      integer, save :: blk_size_hist, blk_size_gp                               
c                                                                               
c      intrinsic allocated                                                      
c                                                                               
      end module                                                                
