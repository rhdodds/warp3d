c     ****************************************************************          
c     *                                                              *          
c     *                      f-90 module crack_growth_data           *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 10/10/97                   *          
c     *                                                              *          
c     *     define the data structures for crack growth              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c           module elem_extinct_data contains allocated arrays for              
c           the element extinction algorithms                                   
c                                                                               
c                                                                               
      module elem_extinct_data                                                  
c                                                                               
c                                                                               
      integer, allocatable, dimension(:)  :: dam_node_elecnt                    
      integer, allocatable, dimension(:)  :: dam_state                          
      allocatable                         :: dam_ifv(:,:)                       
      allocatable                         :: dam_dbar_elems(:,:)                
      logical, allocatable, dimension(:)  :: dam_blk_killed                     
      integer, allocatable, dimension(:)  :: dam_print_list                     
      integer, allocatable, dimension(:)  :: kill_order_list                    
      integer, allocatable                :: dam_face_nodes(:,:)                
      allocatable                         :: old_porosity(:)                    
      allocatable                         :: old_deff(:)                        
      allocatable                         :: old_plast_strain(:)                
      allocatable                         :: old_mises(:)                       
      allocatable                         :: old_mean(:)                        
c                                                                               
      double precision                                                          
     &  dam_ifv, dam_dbar_elems, old_porosity, old_plast_strain,                
     &  old_mean, old_mises, old_deff                                           
c                                                                               
c                                                                               
      end module                                                                
c                                                                               
c                                                                               
c           module node_release_data contains allocated arrays for              
c           the node release algorithms                                         
c                                                                               
c                                                                               
      module node_release_data                                                  
c                                                                               
c                                                                               
      integer, allocatable, dimension(:)  :: crack_plane_nodes                  
      integer, allocatable, dimension(:)  :: inv_crkpln_nodes                   
      integer, allocatable, dimension(:)  :: num_neighbors                      
      integer, allocatable                :: neighbor_nodes(:,:)                
      integer, allocatable                :: crack_front_nodes(:,:)             
      allocatable                         :: crkpln_nodes_react(:)              
      integer, allocatable, dimension(:)  :: crkpln_nodes_state                 
      allocatable                         :: node_release_frac(:)               
      allocatable                         :: old_angles_at_front(:,:)           
      integer, allocatable, dimension(:)  :: master_nodes                       
      integer, allocatable, dimension(:,:):: crack_front_list                   
      integer, allocatable, dimension(:,:):: master_lines                       
c                                                                               
      double precision                                                          
     &  crkpln_nodes_react, node_release_frac, old_angles_at_front              
c                                                                               
c                                                                               
      end module                                                                
