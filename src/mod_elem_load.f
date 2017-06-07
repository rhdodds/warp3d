c     ****************************************************************          
c     *                                                              *          
c     *                      f-90 module elem_load_data              *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 8/1/13 rhd                 *          
c     *                                   chg pointer to allocatable *          
c     *                                                              *          
c     *     define the data structures for element face and body     *          
c     *     loadings                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      module elem_load_data                                                     
c                                                                               
c        ** Note **                                                             
c        Input data describing load floating point values are                   
c        stored as - single precision -                                         
c                                                                               
         type :: elem_load_ptr_type                                             
            real, dimension(:), allocatable      :: vals                        
            integer, dimension(:), allocatable   :: piston_tabnum               
            integer, dimension(:,:), allocatable :: data                        
            integer, dimension(:), allocatable   :: thread_number               
            integer                          :: size                            
         end type                                                               
c                                                                               
         type (elem_load_ptr_type), save, dimension(:),                         
     &        allocatable :: elem_loads                                         
c                                                                               
         real, allocatable, save, dimension(:,:)     :: temp_body               
         real, allocatable, save, dimension(:,:,:)   :: temp_face               
         real, allocatable, save, dimension(:,:)     :: temp_press              
         real, allocatable, save, dimension(:)       :: temp_temper             
         integer, allocatable, save, dimension(:,:)  :: temp_piston             
c                                                                               
         logical, save :: temp_allocated                                        
c                                                                               
         integer, parameter :: numfaces = 6                                     
c                                                                               
         intrinsic allocated                                                    
c                                                                               
      end module                                                                
c                                                                               
