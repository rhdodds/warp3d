c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine smcs_cut_step                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/07/96                   *          
c     *                                                              *          
c     *         This routine checks if the load step size is too     *          
c     *         large. If the plastic strain in any killable mises   *          
c     *         element has grown more than max_plast_strain_change, *          
c     *         then permanently cut the load step size by half.     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine smcs_cut_step ( debug )                                        
      use global_data ! old common.main
      use elem_extinct_data, only : old_plast_strain                            
      use damage_data                                                           
      implicit integer (a-z)                                                    
      logical debug                                                             
c                                                                               
c                                                                               
c           local declarations                                                  
c                                                                               
      logical not_cut, duml                                                     
      double precision                                                          
     &     new_plast_strain, two, dum1, dumd2, dumd3, dumd4,                    
     &     dumd5, dumd6, dumd7                                                  
      character(len=1) :: dums                                                  
      real dumr                                                                 
      data two / 2.0 /                                                          
c                                                                               
      if ( debug ) write ( out, * ) '>>>>>> in smcs_cut_step'                   
c                                                                               
c           loop over all killable mises elements                               
c                                                                               
      not_cut = .true.                                                          
      do elem = 1, noelem                                                       
         elem_ptr = dam_ptr( elem )                                             
         if ( elem_ptr .eq. 0 ) cycle                                           
c                                                                               
c              calculate new plast_strain in element                            
c                                                                               
         call dam_param( elem, duml, debug, dum1, new_plast_strain,             
     &        dumd2, dumd3, dumd4, dumd5, dumd6, dumd7 )                        
c                                                                               
c              compare old plast_strain with new plast_strain -- if             
c              change is larger than the acceptible max, cut load               
c              step size.                                                       
c                                                                               
         if ( new_plast_strain - old_plast_strain(elem_ptr) .gt.                
     &        max_plast_strain_change ) then                                    
          if ( not_cut ) then                                                   
               perm_load_fact = perm_load_fact / two                            
               call errmsg ( 273, dum, dums, dumr, perm_load_fact )             
             not_cut = .false.                                                  
            end if                                                              
         end if                                                                 
       old_plast_strain(elem_ptr) = new_plast_strain                            
      end do                                                                    
c                                                                               
      if ( debug ) write(out,*) '<<<<< leaving smcs_cut_step'                   
c                                                                               
      return                                                                    
      end                                                                       
