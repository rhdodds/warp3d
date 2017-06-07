c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine innum                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 6/5/2017 rhd               *          
c     *                                                              *          
c     *             translate input of model sizes                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine innum( sbflg1, sbflg2 )                                        
      use global_data ! old common.main
      use main_data, only : elstor, incmap, crdmap                              
      use erflgs                                                                
c                                                                               
      implicit none                                                             
c                                                                               
      logical ::  sbflg1, sbflg2                                                
c                                                                               
      integer :: i, dum, param                                                  
      double precision :: dumd                                                  
      real :: dumr                                                              
      character :: dums*1                                                       
      logical :: matchs, integr, endcrd, true                                   
c                                                                               
c                       branch on whether nodes or elements are                 
c                       to be input.                                            
c                                                                               
 310  if( matchs('of',2)       ) call splunj                                    
      if( matchs('nodes',4)    ) go to 320                                      
      if( matchs('elements',4) ) go to 330                                      
c                                                                               
c                       there is no match. check for end of card.               
c                       if not, print error message.                            
c                                                                               
      if( endcrd(dum) ) then                                                    
         go to 9999                                                             
      else                                                                      
         call errmsg(6,dum,dums,dumr,dumd)                                      
         if( true(dum) ) go to 310                                              
      end if                                                                    
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     nodes input.                                   *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 320  continue                                                                  
c                                                                               
c                       check to make sure the number of nodes                  
c                       has not already been input.                             
c                                                                               
      if( numnod ) then                                                         
         call errmsg(2,dum,'nodes',dumr,dumd)                                   
         go to 310                                                              
      end if                                                                    
c                                                                               
c                       input number of nodes. check to make sure               
c                       the maximum number of nodes allowed is not              
c                       exceeded and that the number of nodes is                
c                       greater than zero.                                      
c                                                                               
      numnod = .true.                                                           
      if( .not. integr(nonode) ) then                                           
         call errmsg(7,dum,dums,dumr,dumd)                                      
         numnod = .false.                                                       
      else                                                                      
c                                                                               
         if( nonode .le. 0 ) then                                               
            param = 1                                                           
            call errmsg(111,param,dums,dumr,dumd)                               
            fatal  = .true.                                                     
            numnod = .false.                                                    
            go to 9999                                                          
         end if                                                                 
c                                                                               
c                       got valid number of nodes                               
c                       allocate and initialize  arrays dependent               
c                       only upon the number of nodes   
c
c                       0. displacements, velocities, ... coordinates
c
         call mem_allocate( 4 )                        
c                                                                               
c                       1. accumulated and step incremental temperatures.       
c                                                                               
         call mem_allocate( 1 )                                                 
c                                                                               
c                       2. vectors which record the presence of non-            
c                          global constraint transformations.                   
c                                                                               
         call mem_allocate( 3 )                                                 
c                                                                               
c                       3. array of integers used to store nodal load           
c                          descriptors for each loading conditions.             
c                                                                               
         call mem_allocate( 5 )                                                 
c                                                                               
c                       4. diagonal mass and                                    
c                          effective dynamic load, crdmap.  also                
c                          allocate diagonal stiffness if using serial          
c                          version.                                             
c                                                                               
         if ( .not. use_mpi ) call mem_allocate(11)                             
         call mem_allocate( 12 )                                                
         call mem_allocate( 13 )                                                
         call mem_allocate( 14 )                                                
c                                                                               
c                       5. vectors for constraint definitions                   
c                                                                               
         call mem_allocate( 15 )                                                
c                                                                               
c                       6. vectors for residual loads and step                  
c                          incremental forces                                   
c                                                                               
         call mem_allocate( 16 )                                                
c                                                                               
c                       initialize the global coordinate vectors.               
c                       the initialization of the coordinate                    
c                       pointer vector crdmap is done so that                   
c                       error checking is facilitated                           
c                                                                               
         do i = 1, nonode                                                       
            crdmap(i) = 0                                                       
            dstmap(i) = 0                                                       
         end do                                                                 
      end if                                                                    
      go to 310                                                                 
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     elements input.                                *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 330  continue                                                                  
c                                                                               
c                       check to make sure the number of elements               
c                       has not already been input.                             
c                                                                               
      if( numel ) then                                                          
         call errmsg(2,dum,'elems',dumr,dumd)                                   
         go to 310                                                              
      end if                                                                    
c                                                                               
c                       input number of elements. check to make sure            
c                       the maximum number of elements allowed is not           
c                       exceeded and that the number of elements is             
c                       greater than zero.                                      
c                                                                               
                                                                                
c                                                                               
      numel = .true.                                                            
      if( .not. integr(noelem) ) then                                           
         call errmsg(9,dum,dums,dumr,dumd)                                      
         numel = .false.                                                        
      else                                                                      
c                                                                               
         if( noelem .le. 0 ) then                                               
            param = 2                                                           
            call errmsg(111,param,dums,dumr,dumd)                               
            fatal = .true.                                                      
            numel = .false.                                                     
            go to 9999                                                          
         end if                                                                 
c                                                                               
         if( noelem .gt. mxel ) then                                            
            call errmsg(10,dum,dums,dumr,dumd)                                  
            fatal = .true.                                                      
            numel = .false.                                                     
            go to 9999                                                          
         end if                                                                 
      end if                                                                    
c                                                                               
c                       got valid number of elements.                           
c                       allocate and initialize  arrays dependent               
c                       only upon the number of elements                        
c                                                                               
c                       1. accumulated and step incremental temperatures.       
c                                                                               
        call mem_allocate( 2 )                                                  
c                                                                               
c                       2. temporary element definitions for input              
c                          processing                                           
c                                                                               
        call mem_allocate( 6 )                                                  
c                                                                               
c                       3. temporary element incmap and incidences.             
c                          incid is resized once final data know,               
c                                                                               
        call mem_allocate( 8 )                                                  
c                                                                               
c                       initialize the element property array                   
c                       and the element incidence mapping vector                
c                       to facilitate error checking. also initialize           
c                       the plasticity paramters map and the stress/            
c                       strain map.                                             
c                                                                               
      do i = 1, noelem                                                          
         elstor(1,i) = 0                                                        
         incmap(i)   = 0                                                        
      end do                                                                    
c                                                                               
c                       4. integer flags for damage state in elements,          
c                          initialized to zero                                  
c                                                                               
          call allocate_damage( 12 )                                            
c                                                                               
      go to 310                                                                 
c                                                                               
c                                                                               
c                                                                               
c **********************************************************************        
c **********************************************************************        
c                                                                               
c                                                                               
 9999 sbflg1 = .false.                                                          
      sbflg2 = .false.                                                          
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
