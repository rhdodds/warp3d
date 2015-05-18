c     ****************************************************************
c     *                                                              *
c     *              f-90 module adaptive_steps                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                    last modified : 5/2/04                    *
c     *                                                              *
c     *     define the variables and data structures to support      *
c     *     adaptive load sub-sizing of solution during a user       *
c     *     defined load step to enhance global newton convergence   *
c     *                                                              *
c     ****************************************************************
c
c
      module adaptive_steps
c
      integer adapt_rows, adapt_cols
      parameter ( adapt_rows=5, adapt_cols=30 )
c
      integer, save ::
     &  adapt_level,
     &  adapt_result,
     &  adapt_divisions
c
#dbl      double precision, save :: 
#sgl      real, save ::
     &   adaptive_stack(adapt_rows,adapt_cols),
     &   adapt_disp_fact,
     &   adapt_load_fact, 
     &   predict_disp_fact,
     &   adapt_temper_fact
c
      end module

