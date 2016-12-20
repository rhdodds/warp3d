c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouext2                       *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 5/15/01 rhd                *
c     *                                                              *
c     *     this subroutine computes derived stress/strain values    *
c     *     after the primary values have been averaged              * 
c     *                                                              *
c     ****************************************************************
c
      subroutine ouext2 ( results, nrowd, noval, stress )
      implicit integer (a-z)
      double precision
     &     results(nrowd,*) 
      logical stress
c
      if ( stress ) then
         call princ_inv_stress ( results, nrowd, noval )
         call princ_stress ( results, nrowd, noval )
         call yield_function( results, nrowd, noval )
      else
         call princ_inv_strain ( results, nrowd, noval )
         call princ_strain ( results, nrowd, noval )
         call equiv_strain( results, nrowd, noval )
      end if
c
      return
      end
        
