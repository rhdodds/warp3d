c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine kill_element                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 05/9/2019 rhd              *          
c     *                                                              *          
c     *        This routine kills an element in the model. the       *          
c     *        material properties and history data for the element  *          
c     *        are modified so that subsequent element [k]s will be  *          
c     *        zero. the existing stresses are zeroed. the props and *          
c     *        history changes are such that stresses will remain    *          
c     *        zero for all subsequent steps.                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine kill_element( ielem, debug )                                   
      use global_data ! old common.main
      use main_data                                                             
      use elem_block_data, only : history_blocks, urcs_n_blocks,                
     &                            urcs_n1_blocks, history_blk_list              
      implicit none   
c   
      integer :: ielem
      logical :: debug
c 
c              locals
c
      integer :: elem_type, material_type, ngp, blk, rel_elem,
     &           hist_size, hist_offset, gp, hisloc, sig_start, 
     &           sig_end, i                                              
      double precision :: porosity
      double precision, parameter :: zero = 0.d0
      double precision, dimension(:), pointer :: history, urcs_n1,
     &                                           urcs_n 
      logical :: gurson, regular_material, cohesive_elem
c                                                                               
      if( debug ) write (*,*) '>>>> in kill_element'                           
c                                                                               
c              for regular elements (hex, tet, et.c),                           
c              set young's modulus and poisson's ratio to zero for              
c              element so that any future request for a linear stiffness        
c              computes zero. for cohesive elements, we set a row of            
c              element props table = 1 so the element and cohesive model        
c              routines know element is killed.                                 
c                                                                               
      elem_type     = iprops(1,ielem)   
      material_type = iprops(25,ielem)      
      gurson        = material_type .eq. 3  ! also could just be mises                                      
      regular_material = .not. gurson
      cohesive_elem = cohesive_ele_types(elem_type)    
c                         
      if( .not. cohesive_elem ) then                                           
         props(7,ielem) = 0.0    !   E and nu                                                
         props(8,ielem) = 0.0                                                  
      else  ! cohesive elements                                                                    
         props(7:9,ielem) = 0.0
         props(20,ielem)  = 0.0
         props(21,ielem)  = 0.0
         iprops(32,ielem) = 1                                                   
      end if                                                                    
c                                                                               
c              for each gauss point of the element:                             
c               (a)  reset history data to a null state characteristic          
c                    of an unstressed element                                   
c               (b)  set stresses and n and n+1 to zero.                        
c                                                                               
      ngp           = iprops(6,ielem)                                           
      blk           = elems_to_blocks(ielem,1)                                  
      rel_elem      = elems_to_blocks(ielem,2)                                  
      hist_size     = history_blk_list(blk)                                     
      hist_offset   = (rel_elem-1)*hist_size*ngp                                
      history       => history_blocks(blk)%ptr                                  
      urcs_n1       => urcs_n1_blocks(blk)%ptr                                  
      urcs_n        => urcs_n_blocks(blk)%ptr                                   
c    
      if( gurson ) then ! don't zero porosity
         do gp = 1, ngp                                                            
           hisloc = hist_offset + (gp-1)*hist_size 
           porosity = history(hisloc+5)                          
           history(hisloc+1:hisloc+hist_size) = zero 
           history(hisloc+5) = porosity  
         end do                                         
      end if
c
      if( cohesive_elem ) then
         do gp = 1, ngp                                                            
            hisloc = hist_offset + (gp-1)*hist_size                                
            history(hisloc+1) = zero                                            
            history(hisloc+2) = zero                                            
            history(hisloc+3) = zero                                            
         end do                                                                    
      end if
c
      if( regular_material ) then
         do gp = 1, ngp                                                            
           hisloc = hist_offset + (gp-1)*hist_size 
           history(hisloc+1:hisloc+hist_size) = zero 
         end do
      end if
c                                                                               
      sig_start = (rel_elem-1)*nstrs*ngp + 1                                    
      sig_end   = sig_start + ngp*nstrs - 1                                     
      do i = sig_start, sig_end                                                
        urcs_n1(i) = zero                                                      
        urcs_n(i)  = zero                                                      
      end do                                                                    
c                                                                               
      if( debug ) write (*,*) '<<<< leaving kill_element'                      
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine update_killed_energy         *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 05/9/2019 rhd              *          
c     *                                                              *          
c     *        An element is being killed. We need to save the       *          
c     *        internal and plastic work right now so it can be      *          
c     *        included in future work totals. this is because all   *          
c     *        such data is about to be zeroed.                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine update_killed_energy( ielem )                                  
      use global_data ! old common.main
      use main_data, only : elems_to_blocks, cohesive_ele_types                 
      use elem_block_data, only : urcs_n_blocks, element_vol_blocks             
      implicit none
c
      integer :: ielem
c
      integer :: elem_type, ngp, blk, rel_elem, sig_start, gp
      double precision :: sum_internal_energy, sum_plastic_energy,                       
     &                    element_volume, rgp                                                  
      double precision, pointer :: urcs_n(:) 
      double precision, parameter :: zero = 0.d0                                       
      logical :: cohesive_elem   
      logical, parameter :: debug = .false.                                           
c                                                                               
      if( debug ) write (*,*) '>>>> in update_killed_energy'                   
c                                                                               
c              for regular elements (hex, tet, et.c),                           
c              compute the average internal energy and plastic energy           
c              at the gauss points. these values are multipled by element       
c              volume and added to accumluated totals for killed                
c              elements.                                                        
c                                                                               
      elem_type     = iprops(1,ielem)                                           
      cohesive_elem = cohesive_ele_types(elem_type)                             
      if( cohesive_elem ) return                                               
      ngp           = iprops(6,ielem)                                           
      blk           = elems_to_blocks(ielem,1)                                  
      rel_elem      = elems_to_blocks(ielem,2)                                  
      urcs_n        => urcs_n_blocks(blk)%ptr                                   
      sig_start     = (rel_elem-1) * nstrs * ngp                                
c                                                                               
      sum_internal_energy = zero                                                
      sum_plastic_energy  = zero                                                
c                                                                               
      do gp = 1, ngp                                                            
        sum_internal_energy = sum_internal_energy + urcs_n(sig_start+7)         
        sum_plastic_energy  = sum_plastic_energy  + urcs_n(sig_start+8)         
        sig_start           = sig_start + nstrs                                 
      end do                                                                    
c                                                                               
      rgp                 = dble(ngp)                                           
      sum_internal_energy = sum_internal_energy / rgp                           
      sum_plastic_energy  = sum_plastic_energy  / rgp                           
c                                                                               
      element_volume      = element_vol_blocks(blk)%ptr(rel_elem)               
      killed_ele_int_work = killed_ele_int_work +                               
     &                      element_volume * sum_internal_energy                
      killed_ele_pls_work = killed_ele_pls_work +                               
     &                      element_volume * sum_plastic_energy                 
c                                                                               
      if ( debug ) write (*,*) '<<<< leaving update_killed_energy'              
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                 
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
