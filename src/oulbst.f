c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulbst                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 8/21/2017 rhd              *          
c     *                                                              *          
c     *     sets the label type and corresponding                    *          
c     *     dof labels for an element for stress/strain output       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oulbst( do_stress, lbltyp, elem_type, elem,                    
     &                   strlbl, long, hedtyp, geonl, cohesive_type )           
      use main_data, only: bar_types, link_types                                      
      implicit none
c
      integer :: lbltyp, elem_type, elem, cohesive_type
      logical :: long, do_stress, geonl     
      character(len=8) :: strlbl(*)                                             
      character(len=*) :: hedtyp
c                                                      
      logical :: bar_type, link_type  
c                                                                               
c                       solids and interface-cohesive are separate              
c                                                                               
      bar_type  = bar_types(elem_type)
      link_type = link_types(elem_type)
c
      if( bar_type ) then
         call oulbst_bar
         return
      end if  
c             
      if( link_type ) then
         call oulbst_link
         return
      end if         
c
      if( cohesive_type. ne. 0 ) then                                           
         call oulbst_cohesive                                                   
      else                                                                      
         call oulbst_solid                                                      
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
      contains                                                                  
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine oulbst_cohesive                *          
c     *                                                              *          
c     *                       written by : rhs                       *          
c     *                                                              *          
c     *                   last modified : 3/4/2015 rhd               *          
c     *                                                              *          
c     *  sets output header labels for interface-cohesive element    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine oulbst_cohesive                                                
      implicit none                                                    
c                                                                               
      select case ( cohesive_type )                                             
c                                                                               
      case( 1, 2, 3, 5 )                                                        
c                                                                               
       if( do_stress ) then                                                     
          hedtyp(1:) = 'tractions'                                              
          if ( geonl ) hedtyp(1:) = 'tractions (large displacement)'            
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
          strlbl(5)  = 'gamma'                                                  
          strlbl(6)  = 'gamma-ur'                                               
       else                                                                     
          hedtyp(1:) = 'displacement jumps'                                     
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
       end if                                                                   
c                                                                               
      case ( 4 )  ! simple exponential                                          
                                                                                
       if( do_stress ) then                                                     
          hedtyp(1:) = 'tractions'                                              
          if ( geonl ) hedtyp(1:) = 'tractions (large displacement)'            
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
          strlbl(5)  = 'eff'                                                    
          strlbl(6)  = 'eff/peak'                                               
          strlbl(7)  = 'gamma'                                                  
          strlbl(8)  = 'gamma-ur'                                               
       else                                                                     
          hedtyp(1:) = 'displacement jumps'                                     
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
          strlbl(5)  = 'eff'                                                    
          strlbl(6)  = 'eff/peak'                                               
       end if                                                                   
c                                                                               
      case ( 6 )   ! ppr                                                        
                                                                                
       if( do_stress ) then                                                     
          hedtyp(1:) = 'tractions'                                              
          if ( geonl ) hedtyp(1:) = 'tractions (large displacement)'            
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
          strlbl(5)  = 'shr/peak'                                               
          strlbl(6)  = 'nml/peak'                                               
          strlbl(7)  = 'gamma'                                                  
          strlbl(8)  = 'gamma-ur'                                               
       else                                                                     
          hedtyp(1:) = 'displacement jumps'                                     
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
          strlbl(5)  = 'shr/peak'                                               
          strlbl(6)  = 'nml/peak'                                               
       end if                                                                   
c                                                                               
c                                                                               
      case ( 7 )   ! cavit                                                      
                                                                                
       if( do_stress ) then                                                     
          hedtyp(1:) = 'tractions'                                              
          if ( geonl ) hedtyp(1:) = 'tractions (large displacement)'            
c                       12345678                                                
          strlbl(1)  = ' shear-1'                                               
          strlbl(2)  = ' shear-2'                                               
          strlbl(3)  = '   shear'                                               
          strlbl(4)  = '  normal'                                               
          strlbl(5)  = '   gamma'                                               
          strlbl(6)  = 'gamma-ur'                                               
          strlbl(7)  = '   N/N_I'                                               
          strlbl(8)  = '     a/b'                                               
          strlbl(9)  = '   a/a_0'                                               
          strlbl(10) = 'crit val'                                               
          strlbl(11) = '   max_T'                                               
          strlbl(12) = 'mxdeltac'                                               
          strlbl(13) = '   a/Lnr'                                               
       else                                                                     
          hedtyp(1:) = 'displacement jumps'                                     
          strlbl(1)  = 'shear-1'                                                
          strlbl(2)  = 'shear-2'                                                
          strlbl(3)  = 'shear'                                                  
          strlbl(4)  = 'normal'                                                 
       end if                                                                   
c                                                                               
      case default                                                              
c                                                                               
        write(*,9000) cohesive_type                                             
        call die_abort                                                          
c                                                                               
      end select                                                                
c                                                                               
      return                                                                    
9000  format(/,">>>> FATAL ERROR. oulbst_cohesive. val: ",i10,                  
     & /,      "                  job aborted...",//)                           
c                                                                               
      end subroutine oulbst_cohesive                                            
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine oulbst_solid                   *          
c     *                                                              *          
c     *                       written by : rhs                       *          
c     *                                                              *          
c     *                   last modified : 03/11/13  rhd              *          
c     *                                                              *          
c     *       sets output header labels for solid elements           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
                                                                                
      subroutine   oulbst_solid                                                 
      implicit none                                                    
c                                                                               
c                       short output                                            
c                                                                               
      if( do_stress ) then                                                      
         hedtyp(1:) = 'stresses'                                                
         if ( geonl ) hedtyp(1:) = 'cauchy stresses'                            
         strlbl(1)  = 'sigma_xx'                                                
         strlbl(2)  = 'sigma_yy'                                                
         strlbl(3)  = 'sigma_zz'                                                
         strlbl(4)  = 'sigma_xy'                                                
         strlbl(5)  = 'sigma_yz'                                                
         strlbl(6)  = 'sigma_xz'                                                
         strlbl(7)  = 'eng_dens'                                                
         strlbl(8)  = 'mises   '                                                
         strlbl(9)  = 'c1      '                                                
         strlbl(10) = 'c2      '                                                
         strlbl(11) = 'c3      '                                                
      else                                                                      
         hedtyp(1:)= 'strains '                                                 
         if ( geonl ) hedtyp(1:) = 'logarithmic strains'                        
         strlbl(1) = '  eps_xx'                                                 
         strlbl(2) = '  eps_yy'                                                 
         strlbl(3) = '  eps_zz'                                                 
         strlbl(4) = '  gam_xy'                                                 
         strlbl(5) = '  gam_yz'                                                 
         strlbl(6) = '  gam_xz'                                                 
         strlbl(7) = ' eff_eps'                                                 
      end if                                                                    
c                                                                               
      if( long ) then                                                           
         if( do_stress ) then                                                   
            strlbl(12) = '   inv_1'                                             
            strlbl(13) = '   inv_2'                                             
            strlbl(14) = '   inv_3'                                             
            strlbl(15) = '   sig_1'                                             
            strlbl(16) = '   sig_2'                                             
            strlbl(17) = '   sig_3'                                             
            strlbl(18) = 'angle1_x'                                             
            strlbl(19) = 'angle1_y'                                             
            strlbl(20) = 'angle1_z'                                             
            strlbl(21) = 'angle2_x'                                             
            strlbl(22) = 'angle2_y'                                             
            strlbl(23) = 'angle2_z'                                             
            strlbl(24) = 'angle3_x'                                             
            strlbl(25) = 'angle3_y'                                             
            strlbl(26) = 'angle3_z'                                             
         else                                                                   
            strlbl(8)= '   inv_1'                                               
            strlbl(9)= '   inv_2'                                               
            strlbl(10)= '   inv_3'                                              
            strlbl(11)= '   eps_1'                                              
            strlbl(12)= '   eps_2'                                              
            strlbl(13)= '   eps_3'                                              
            strlbl(14)= 'angle1_x'                                              
            strlbl(15)= 'angle1_y'                                              
            strlbl(16)= 'angle1_z'                                              
            strlbl(17)= 'angle2_x'                                              
            strlbl(18)= 'angle2_y'                                              
            strlbl(19)= 'angle2_z'                                              
            strlbl(20)= 'angle3_x'                                              
            strlbl(21)= 'angle3_y'                                              
            strlbl(22)= 'angle3_z'                                              
         end if                                                                 
      end if                                                                    
      return                                                                    
      end subroutine oulbst_solid                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine oulbst_bar                     *          
c     *                                                              *          
c     *                       written by : rhs                       *          
c     *                                                              *          
c     *                   last modified : 03/3/2017 rhd              *          
c     *                                                              *          
c     *       sets output header labels for bar elements             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
                                                                                
      subroutine   oulbst_bar                                                 
      implicit none                                                    
c       
      if( do_stress ) then                                                      
         hedtyp(1:) = 'stresses'                                                
         if ( geonl ) hedtyp(1:) = 'cauchy stresses'                            
         strlbl(1)  = 'sigma'                                                
         strlbl(2)  = 'force'                                                
      else                                                                      
         hedtyp(1:)= 'strains '                                                 
         if ( geonl ) hedtyp(1:) = 'logarithmic strains'                        
         strlbl(1) = '  eps'                                                 
      end if                                                                    
c
      return                                                                    
      end subroutine oulbst_bar                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine oulbst_link                    *          
c     *                                                              *          
c     *                       written by : rhs                       *          
c     *                                                              *          
c     *                   last modified : 08/20/2017 rhd             *          
c     *                                                              *          
c     *       sets output header labels for link elements            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
                                                                                
      subroutine   oulbst_link                                                 
      implicit none                                                    
c                                                                               
      if( do_stress ) then                                                      
         hedtyp(1:) = 'forces'                                                
         strlbl(1)  = 'force_x'                                                
         strlbl(2)  = 'force_y'
         strlbl(3)  = 'force_z'
      else                                                                      
         hedtyp(1:)= 'rel-dis'                                                 
         strlbl(1) = '   du'                                                 
         strlbl(2) = '   dv'                                                 
         strlbl(3) = '   dw'                                                 
      end if                                                                    
c
      return                                                                    
      end subroutine oulbst_link                                               
c                                                                               
      end subroutine oulbst                                                     
