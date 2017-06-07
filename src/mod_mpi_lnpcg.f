c     ****************************************************************          
c     *                                                              *          
c     *                      f-90 module mpi_lnpcg                   *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 02/04/98                   *          
c     *                                                              *          
c     *     define the data structures for domain decomposed         *          
c     *     version of lnpcg.                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      module mpi_lnpcg                                                          
c                                                                               
c                                                                               
c               -----------------------------------------------                 
c                data structure for storing node ownership info                 
c               -----------------------------------------------                 
c                                                                               
c                                                                               
        type :: list_type                                                       
c                                                                               
           integer, dimension(:), pointer   :: ptr                              
c                                                                               
        end type                                                                
c                                                                               
c                                                                               
        type :: node_owner_type                                                 
c                                                                               
c                storage for global info                                        
c                                                                               
           integer                          :: num_local_nodes                  
           integer, dimension(:), pointer   :: local2global                     
           integer, dimension(:), pointer   :: global2local                     
c                                                                               
c                storage for nodes owned locally and not shared                 
c                                                                               
           integer                          :: num_private                      
           integer, dimension(:), pointer   :: private                          
c                                                                               
c                storage for nodes owned locally but shared                     
c                                                                               
           integer                          :: num_own_shared                   
           integer, dimension(:), pointer   :: own_shared                       
           integer, dimension(:), pointer   :: sharing_count                    
           type(list_type), dimension(:), pointer :: sharing_procs              
           integer, dimension(:), pointer   :: MPI_sharing_type                 
c                                                                               
c                storage for nodes only shared (not owned locally)              
c                                                                               
           integer                          :: num_shared                       
           integer, dimension(:), pointer   :: shared                           
           integer, dimension(:), pointer   :: shared_count                     
           type(list_type), dimension(:), pointer :: shared_owner               
           integer, dimension(:), pointer   :: MPI_shared_type                  
c                                                                               
c               storage for ordering info for ebe preonditioner                 
c                                                                               
         integer                          :: num_int_blks                       
           integer, dimension(:), pointer   :: internal_blks                    
c                                                                               
        end type                                                                
c                                                                               
        type (node_owner_type), save        :: local_nodes                      
c                                                                               
        type (node_owner_type), save, dimension(:),                             
     &       allocatable  :: proc_nodes(:)                                      
c                                                                               
        type :: map_type                                                        
c                                                                               
           integer :: num_dof                                                   
           integer, dimension(:), pointer   :: dof                              
c                                                                               
        end type                                                                
c                                                                               
      type (map_type), dimension(:), allocatable,                               
     &            save :: procdof2glob                                          
c                                                                               
c                                                                               
c               -------------------------------------------                     
c                local version of edest indexes for blocks                      
c               -------------------------------------------                     
c                                                                               
c                                                                               
        type :: int_blocks_ptr_type                                             
              integer, dimension(:,:), pointer :: ptr                           
        end type                                                                
c                                                                               
        type (int_blocks_ptr_type), save, dimension(:),                         
     &         allocatable :: ledest_blocks                                     
        integer, dimension (:), allocatable, save :: ledest_blk_list            
c                                                                               
c                                                                               
c               ------------------------------------------------                
c                data structure for graph nodes for parallelism                 
c                in ebe preconditioner                                          
c               ------------------------------------------------                
c                                                                               
c                                                                               
        type :: link_type                                                       
c                                                                               
           integer                          :: id                               
           integer                          :: owner                            
           integer                          :: order                            
           integer                          :: num_dofs                         
           integer, dimension(:), pointer   :: dofs                             
           integer                          :: MPI_type                         
c                                                                               
        end type                                                                
c                                                                               
        type :: graph_type                                                      
c                                                                               
           integer                          :: id                               
           integer                          :: owner                            
           integer                          :: order                            
           integer                          :: num_blks                         
           integer, dimension(:), pointer   :: blks                             
           integer                          :: num_elem                         
           integer                          :: num_proc_share                   
           integer, dimension(:), pointer   :: proc_share                       
           integer                          :: num_int_links                    
           integer, dimension(:), pointer   :: int_links                        
           integer                          :: num_ext_links                    
           integer, dimension(:), pointer   :: ext_links                        
           logical                          :: allo_link                        
           type(link_type), dimension(:), pointer   :: link                     
c                                                                               
        end type                                                                
c                                                                               
c                                                                               
        type (graph_type), save, dimension(:),                                  
     &       allocatable  :: local_graph(:)                                     
c                                                                               
        integer, save    :: num_groups                                          
        integer, save    :: num_graph                                           
c                                                                               
c                                                                               
c               ------------------------------------------------                
c                data structure for ordering execution of parts                 
c                for ebe preconditioner                                         
c               ------------------------------------------------                
c                                                                               
c                                                                               
c                                                                               
        integer, save, allocatable          ::  order_list(:,:)                 
        integer, save                       ::  order_ptr                       
c                                                                               
      end module                                                                
c                                                                               
