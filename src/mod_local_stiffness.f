c                                                                               
c     ***********************************************************************   
c     *                                                                     *   
c     *     Module local_stiffness_mod                                      *   
c     *                                                                     *   
c     *     created by: mcm 2/11                                            *   
c     *     last modified: mcm 5/11                                         *   
c     *                                                                     *   
c     *     Contains a data structure for holding the local part of         *   
c     *     a row-block distributed matrix-vector system and functions      *   
c     *     to take a fully  assembled stiffness + vector system            *   
c     *     on 1 processor and distribute it (by contiguous row blocks)     *   
c     *     to all the processors on a communicator                         *   
c     *                                                                     *   
c     *     Notes: Will fail if a diagonal entry is not represented in the  *   
c     *     sparsity structure.  Of course this should never happen         *   
c     *     in WARP (b/c we store the full diagonal) but it is something    *   
c     *     to keep in mind if you want to reuse the module.                *   
c     *                                                                     *   
c     ***********************************************************************   
c                                                                               
      module local_stiffness_mod                                                
            implicit none                                                       
c                                                                               
      end module local_stiffness_mod                                            
