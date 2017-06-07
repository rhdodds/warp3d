c     ****************************************************************          
c     *                                                              *          
c     *              f-90 module adaptive_steps                      *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                    last modified : 11/11/2015 rhd            *          
c     *                                                              *          
c     *     define the variables and data structures to support      *          
c     *     adaptive load sub-sizing of solution during a user       *          
c     *     defined load step to enhance global newton convergence   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      module adaptive_steps                                                     
      implicit none                                                             
c                                                                               
      integer adapt_rows, adapt_cols                                            
      parameter ( adapt_rows=5, adapt_cols=110 )                                
c                                                                               
      integer, save ::  adapt_level, adapt_result, adapt_divisions              
c                                                                               
      double precision, save ::                                                 
     &   adaptive_stack(adapt_rows,adapt_cols),                                 
     &   adapt_disp_fact,                                                       
     &   adapt_load_fact,                                                       
     &   predict_disp_fact,                                                     
     &   adapt_temper_fact,                                                     
     &   adapt_min_fact                                                         
c                                                                               
      end module                                                                
                                                                                
