c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                     module stiffness_data                    *          
c     *                                                              *          
c     *                       written by : bjb                       *          
c     *                                                              *          
c     *                    last modified : 01/19/04 bjb              *          
c     *                                                              *          
c     *  Module for the stiffness matrix data structures; use a      *          
c     *  module to permit use of global values in the solvers.       *          
c     *                                                              *          
c     *  NOTE: this module is currently only used by the mpc         *          
c     *        implementation routines                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      module  stiffness_data                                                    
c                                                                               
c           Variables in the 'stiffness_data' module                            
c                                                                               
c           ncoeff = # nonzero coefficients (not counting diagonal              
c               terms in stiffness matrix                                       
c           big_ncoeff = max # ncoeff during mpc implementation--this           
c               number is kept to allow for quick re-dimensioning during        
c               subsequent calls to the mpc solver routines                     
c           newcount = # new terms added to stiffness matrix by mpc             
c               implementation routines--length of new_locations and            
c               new_indexes vectors                                             
c           dep_len = length of dep_locations vector                            
c           ind_len = length of ind_locations vector                            
c           dia_len = length of diag_locations vector                           
c                                                                               
      integer, save ::  ncoeff, big_ncoeff, newcount, dep_len,                  
     &                  ind_len, dia_len, temp_len, new_len,                    
     &                  ncoeff_from_assembled_profile                           
c                                                                               
c           Arrays in the 'stiffness_data' module                               
c                                                                               
c           k_indexes = the column numbers on each row containing a             
c               nonzero term in the stiffness matrix                            
c           new_locations = location in expanded k_indexes vector where         
c               the mpc solver routines added terms                             
c           new_indexes = the column numbers inserted in the locations          
c               pointed to by 'new_locations'                                   
c           new_ptrs = number of terms on each row of stiffness matrix          
c               after mpc implementation (replaces k_ptrs)                      
c           dep_locations = locations in expanded 'k_indexes' of the terms      
c               in the dependent equations (set to 0.0 after modification)      
c           ind_locations = locations in expanded 'k_indexes' of the terms      
c               in the independent equations                                    
c           diag_locations = locations of terms that will update the ind        
c               equation diagonal terms                                         
c           dep_loc = temporary storage (deallocated after use)                 
c           ind_loc = temporary storage (deallocated after use)                 
c           dia_loc = temporary storage (deallocated after use)                 
c           k_coeffs = actual terms of stiffness matrix                         
c           cof_temp = temporary storage (deallocated after use)                
c           ind_temp = temporary storage (deallocated after use)                
c                                                                               
      integer, save, allocatable, dimension (:) :: k_indexes,                   
     &                                       ind_temp,                          
     &                                       new_locations,                     
     &                                       new_indexes,                       
     &                                       new_ptrs,                          
     &                                       dep_locations,                     
     &                                       ind_locations,                     
     &                                       diag_locations,                    
     &                                       dep_loc, ind_loc, dia_loc,         
     &                                       new_loc, new_ind                   
c                                                                               
      double precision,                                                         
     &          allocatable, save, dimension (:) ::                             
     &          k_coeffs, cof_temp                                              
c                                                                               
c           supporting vectors to store nodal forces that impose                
c           the multipoint and tied contact facilities (just MPCs for           
c           short). These may be                                                
c           considered the same as Lagrange multiplier forces acting            
c           on model nodes.                                                     
c                                                                               
c           We do not actually use the Lagrange multiplier method which         
c           increases the number of equations to solve. We use a                
c           method that modifies the assembled K and rhs vectors                
c           to impose the constraints using the concept of dependent            
c           (dep) and independent (ind) structural dofs. The dep dofs           
c           are effectively condensed out using a very efficient                
c           procedure with row-by-row operations on K.                          
c                                                                               
c           The number of assembled equilibrium equations is                    
c           always = 3 * nnode - number of dof w/ absolute                      
c           constraints.                                                        
c                                                                               
c           We enforce K and pvec changes to reflect the MPCs                   
c           and do not remove the dep dofs. Rather we zero their                
c           rows/columns, put 1.0 on the diagonal and make the rhs              
c           term = 0. The dep dofs thus = 0 from equation solver                
c           with the ind dofs computed. values of the dep                       
c           dof are found from the MPC definition along with the                
c           nodal forces acting on both dep and ind dof                         
c           (we call those Lagrange forces even though are strictly             
c           not them in our implementation).                                    
c                                                                               
c           our implementation is motivated by the formulation                  
c           described in: IJNME 1979, pp. 464-467 by J. Abel                    
c           and M. Shepard                                                      
c                                                                               
c           all 3 vectors below always contain set of self-                     
c           equilibrating nodal forces                                          
c                                                                               
c           total_lagrange_forces:                                              
c            accumulated values over solution history. always contains          
c            forces at time step n. Must be zeroed at start of solution         
c            and save/restored across restarts                                  
c                                                                               
c           d_lagrange_forces:                                                  
c            zeroed at start of solution for n -> n+1                           
c            accumulated change in Lagrange forces over n -> n+1.               
c            At completion of Newton iterations for a step, update              
c            total_lagrange_forces <- total_lagrange_forces +                   
c                                     d_lagrange_forces                         
c           i_lagrange_forces:                                                  
c            change in Lagrange forces caused only by the just                  
c            completed Newton iteration. update                                 
c            d_lagrange_forces <- d_lagrange_forces +                           
c                                 i_lagrange_forces                             
c                                                                               
c            These 3 vectors have direct analogies with those for the           
c            nodal displacements: u, du, idu.                                   
c            They will all have size nodof.                                     
c                                                                               
      double precision,                                                         
     &          allocatable, save, dimension (:) ::                             
     &          total_lagrange_forces,                                          
     &          d_lagrange_forces,                                              
     &          i_lagrange_forces                                               
c                                                                               
      end module  stiffness_data                                                
