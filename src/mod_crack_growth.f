c     ****************************************************************          
c     *                                                              *          
c     *                      f-90 module crack_growth_data           *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 4/28/2019 rhd              *          
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
      double precision, allocatable       :: dam_ifv(:,:)                       
      double precision, allocatable       :: dam_dbar_elems(:,:)                
      logical, allocatable, dimension(:)  :: dam_blk_killed                     
      integer, allocatable, dimension(:)  :: dam_print_list                     
      integer, allocatable, dimension(:)  :: kill_order_list                    
      integer, allocatable                :: dam_face_nodes(:,:)                
      double precision, allocatable       :: old_porosity(:)                    
      double precision, allocatable       :: old_deff(:)                   
      double precision, allocatable       :: old_plast_strain(:)                
      double precision, allocatable       :: old_mises(:)                       
      double precision, allocatable       :: old_mean(:)    
      double precision, allocatable       :: smcs_weighted_T(:),
     &                                       smcs_old_epsplas(:)                   
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
