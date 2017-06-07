c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine wmpi_handle_slaves           *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : mcm 05/11                  *          
c     *                                                              *          
c     *     The MPI implementation of warp3d follows a master-slave  *          
c     *     approach for most of the code.  The root processor       *          
c     *     (processor 0) runs warp3d in full, notifying the slave   *          
c     *     processors when they have work to do.                    *          
c     *                                                              *          
c     *     This routine serves as the waiting point for all of the  *          
c     *     slave processors, where they wait for a message from     *          
c     *     the root telling them what to do next.                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine wmpi_handle_slaves                                             
c                                                                               
      end                                                                       
                                                                                
