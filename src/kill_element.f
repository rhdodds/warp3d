c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine kill_element                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/03/94                   *          
c     *                   last modified : 08/13/94 rhd               *          
c     *                   last modified : 04/03/94 asg -- don't zero *          
c     *                                                 the porosity *          
c     *                   last modified : 03/10/04 (rhd)             *          
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
      implicit integer (a-z)                                                    
      double precision                                                          
     &     zero,word                                                            
      double precision,                                                         
     &   dimension(:), pointer :: history, urcs_n1, urcs_n                      
c                                                                               
      integer iword(2)                                                          
      logical debug, cohesive_elem                                              
      equivalence ( iword , word )                                              
      data zero / 0.0 /                                                         
      data iword(1), iword(2) /-1,0/                                            
c                                                                               
c                                                                               
      if ( debug ) write (*,*) '>>>> in kill_element'                           
c                                                                               
c              for regular elements (hex, tet, et.c),                           
c              set young's modulus and poisson's ratio to zero for              
c              element so that any future request for a linear stiffness        
c              computes zero. for cohesive elements, we set a row of            
c              element props table = 1 so the element and cohesive model        
c              routines know element is killed.                                 
c                                                                               
      elem_type = iprops(1,ielem)                                               
      cohesive_elem = cohesive_ele_types(elem_type)                             
      if ( .not. cohesive_elem ) then                                           
         props(7,ielem) = zero                                                  
         props(8,ielem) = zero                                                  
      else                                                                      
         props(7,ielem) = zero                                                  
         props(8,ielem) = zero                                                  
         props(9,ielem) = zero                                                  
         props(20,ielem) = zero                                                 
         props(21,ielem) = zero                                                 
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
      do gp = 1, ngp                                                            
         hisloc = hist_offset + (gp-1)*hist_size                                
         if ( .not. cohesive_elem ) then                                        
c                                                                               
c             do not zero the porosity -- keep it's killed value                
c             after it is killed                                                
c                                                                               
            history(hisloc+1) = zero                                            
            history(hisloc+4) = zero                                            
            history(hisloc+6)  = word                                           
            history(hisloc+7)  = zero                                           
            history(hisloc+8)  = zero                                           
            history(hisloc+9)  = zero                                           
            history(hisloc+10) = zero                                           
         else                                                                   
            history(hisloc+1) = zero                                            
            history(hisloc+2) = zero                                            
            history(hisloc+3) = zero                                            
         end if                                                                 
      end do                                                                    
c                                                                               
      sig_start = (rel_elem-1)*nstrs*ngp + 1                                    
      sig_end   = sig_start + ngp*nstrs - 1                                     
      do ii = sig_start, sig_end                                                
        urcs_n1(ii) = zero                                                      
        urcs_n(ii)  = zero                                                      
      end do                                                                    
c                                                                               
      if ( debug ) write (*,*) '<<<< leaving kill_element'                      
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine update_killed_energy         *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/26/99                   *          
c     *                                 : 01/02/01 sushovan          *          
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
      implicit integer (a-z)                                                    
      double precision                                                          
     &     zero, sum_internal_energy, sum_plastic_energy,                       
     &     element_volume, rgp                                                  
      double precision,                                                         
     &   dimension(:), pointer :: urcs_n                                        
c                                                                               
      logical debug, cohesive_elem                                              
      data zero, debug / 0.0, .false. /                                         
c                                                                               
c                                                                               
      if ( debug ) write (*,*) '>>>> in update_killed_energy'                   
c                                                                               
c              for regular elements (hex, tet, et.c),                           
c              compute the average internal energy and plastic energy           
c              at the gauss points. these values are multipled by element       
c              volume and added to accumluated totals for killed                
c              elements.                                                        
c                                                                               
      elem_type     = iprops(1,ielem)                                           
      cohesive_elem = cohesive_ele_types(elem_type)                             
      if ( cohesive_elem ) return                                               
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
      rgp                 = real(ngp)                                           
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
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
