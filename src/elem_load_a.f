c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine driv_eload                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 11/15/11 jcs               *          
c     *                                                              *          
c     *     this subroutine drives the calculation of the            *          
c     *     equivalent nodal loads for element loadings.             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c           Here is how the face/body/pressure loading storage                  
c           strategy works:                                                     
c                                                                               
c                 There is one 3 by 'eload_size' integer array and one          
c                 'eload_size' long real array that store the loading           
c                 information:                                                  
c                                                                               
c                   eload_data(i,1) -- the element with the loading             
c                                       description for entry i                 
c                   eload_data(i,2) -- the type of the loading for entry i:     
c                               1 to  6 : face number for traction loading      
c                                  0    : body force loading                    
c                              -1 to -6 : pressure loading                      
c                              -7 to -12: -face # for piston loading -6         
c                   eload_data(i,3)  -- the dof for the loading for entry i     
c                               (zero for pressure loading)                     
c                   eload_val(i)     -- the value of the force for entry i      
c                   eload_pist(i)   -- piston table number for entry i          
c                   thread_number(i) -- thread number for entry i               
c                                                                               
c                                                                               
c                                                                               
      subroutine driv_eload( loadnum, mult_fact, rload )                        
      use global_data ! old common.main
      use elem_load_data, only : elem_loads                                     
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
c                                                                               
c                global variables                                               
c                                                                               
      double precision                                                          
     &     mult_fact, rload(*)                                                  
c                                                                               
c                                                                               
c                   originally, this routine was to govern element              
c                   dependencies for load calculations. the current             
c                   code structure deals with different element types           
c                   in later routines. simply call elem_load to start           
c                   calculations for all elements.                              
c                                                                               
c                                                                               
      if ( elem_loads(loadnum)%size .gt. 0 )                                    
     &      call elem_load ( elem_loads(loadnum)%data(1,1),                     
     &      elem_loads(loadnum)%vals(1), elem_loads(loadnum)%size,              
     &      elem_loads(loadnum)%piston_tabnum(1),                               
     &      elem_loads(loadnum)%thread_number(1), mult_fact, rload )            
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
