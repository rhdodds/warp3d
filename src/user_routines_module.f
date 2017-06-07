c     ****************************************************************          
c     *                                                              *          
c     *  f-90 module user_nodal_load_data -- distribution version    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *              last modified : 11/7/2015 rhd                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      module mod_user_routines                                                  
c                                                                               
c234567890123456                                                                
c              available as needed by the Abaqus compatible                     
c              user_routines_....                                               
c                                                                               
c              WARP3D is unaware of the contents here                           
c                                                                               
c              The WARP3D makefile system insures that this module.f            
c              is compiled before the user_routines_...f.                       
c                                                                               
      integer, save :: xxx                                                      
c                                                                               
      end module mod_user_routines                                              
