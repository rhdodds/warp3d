c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine infgm                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/26/01 mcw               *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *     material properties at model nodes to define fgms        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine infgm( sbflg1, sbflg2 )                                        
      use global_data ! old common.main
      use main_data, only: fgm_node_values, fgm_node_values_defined,            
     &                     fgm_node_values_cols                                 
      implicit integer (a-z)                                                    
      double precision                                                          
     & dumd                                                                     
      real dumr, young_mod, poisson_ratio, alpha, vol_fract_ductile,            
     &     rho, tan_e, n_power, yld_pt, zero                                    
      character(len=1) :: dums                                                  
      logical sbflg1, sbflg2                                                    
      logical matchs, matchs_exact, endcrd, true, numr,                         
     &        prop_flags(fgm_node_values_cols)                                  
      dimension intlst(mxlsz)                                                   
      data zero / 0.d0 /                                                        
c                                                                               
      if ( .not. fgm_node_values_defined )  call mem_allocate( 20 )             
c                                                                               
 505  continue                                                                  
      young_mod         = zero                                                  
      poisson_ratio     = zero                                                  
      alpha             = zero                                                  
      vol_fract_ductile = zero                                                  
      rho               = zero                                                  
      tan_e             = zero                                                  
      yld_pt            = zero                                                  
      n_power           = zero                                                  
      prop_flags(1:fgm_node_values_cols) = .false.                              
c                                                                               
      call readsc                                                               
      if ( matchs('dump',4) ) then                                              
       call infgm_dump                                                          
       go to 505                                                                
      end if                                                                    
                                                                                
c                                                                               
c                       translate the list of nodes input on                    
c                       this line.                                              
c                                                                               
      if ( matchs('node',4) ) then                                              
           call splunj                                                          
      else                                                                      
           call reset                                                           
      end if                                                                    
      call scan                                                                 
      call trlist( intlst, mxlsz, nonode, lenlst, errnum )                      
c                                                                               
c                       branch on the return code from trlist. a                
c                       value of 1 indicates no error. a value of               
c                       2 indicates that the parse rules failed in              
c                       the list. a value of 3 indicates that the               
c                       list overflowed its maximum length of mxlsz.            
c                       in these last two cases, the rest of the card           
c                       will be ignored and a new card will be sought.          
c                       a value of 4 indicates that no list was found.          
c                       in this case, either thrm. expansion input has          
c                       ceased.                                                 
c                                                                               
      if( errnum .eq. 2 ) then                                                  
         param = 1                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         go to 505                                                              
      else if( errnum .eq. 3 ) then                                             
         param = 2                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         go to 505                                                              
      else if( errnum .eq. 4) then                                              
         go to 9999                                                             
      else                                                                      
         if( errnum .eq. 1 ) then                                               
            call backsp(1)                                                      
            prop_flags(1:fgm_node_values_cols) = .false.                        
            go to 511                                                           
         end if                                                                 
         param = 3                                                              
         call errmsg(24,param,dums,dumr,dumd)                                   
         go to 505                                                              
      end if                                                                    
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *           a node list exists. parse the material property values   *        
c *           for these nodes. read values until end of line           *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 511  continue                                                                  
      if ( true(dumr) ) call splunj                                             
      if ( matchs_exact('e') ) then                                             
        if ( .not. numr(young_mod) ) then                                       
         call errmsg( 5,dum,'e   ',dumr,dumd )                                  
        else                                                                    
         prop_flags(1) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if ( matchs_exact('nu') ) then                                            
        if ( .not. numr(poisson_ratio) ) then                                   
         call errmsg( 5,dum,'nu  ',dumr,dumd )                                  
        else                                                                    
         prop_flags(2) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if ( matchs_exact('tan_e') ) then                                         
        if ( .not. numr(tan_e) ) then                                           
         call errmsg( 5,dum,'tan_e',dumr,dumd )                                 
        else                                                                    
         prop_flags(6) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if ( matchs_exact('yld_pt') ) then                                        
        if ( .not. numr(yld_pt) ) then                                          
         call errmsg( 5,dum,'yld_pt',dumr,dumd )                                
        else                                                                    
         prop_flags(7) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if ( matchs_exact('n_power') ) then                                       
        if ( .not. numr(n_power) ) then                                         
         call errmsg( 5,dum,'n_power',dumr,dumd )                               
        else                                                                    
         prop_flags(8) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if( matchs_exact('alpha') ) then                                          
        if ( .not. numr(alpha) ) then                                           
         call errmsg( 5,dum,'alpha',dumr,dumd )                                 
        else                                                                    
         prop_flags(3) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if( matchs('vol_fract_ductile',8) ) then                                  
        if ( .not. numr(vol_fract_ductile) ) then                               
         call errmsg( 5,dum,'vol_fract_ductile',dumr,dumd )                     
        else                                                                    
         prop_flags(4) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
c                                                                               
      if( matchs_exact('rho') ) then                                            
        if ( .not. numr(rho) ) then                                             
         call errmsg( 5,dum,'rho',dumr,dumd )                                   
        else                                                                    
         prop_flags(5) = .true.                                                 
         go to 511                                                              
        end if                                                                  
      end if                                                                    
      if( matchs(',',1) ) go to 525                                             
c                                                                               
c                       there is no match for a material property.              
c                       check for end of card. if not, print error              
c                       message.                                                
c                                                                               
      if( endcrd(dum) ) then                                                    
         go to 590                                                              
      else                                                                      
         call errmsg2(29,dum,dums,dumr,dumd)                                    
         if( true(dum) ) go to 511                                              
      end if                                                                    
c                                                                               
c                       comma separator. if at the end of a line,               
c                       denotes continuation. otherwise, ignore.                
c                                                                               
 525  continue                                                                  
      if( endcrd(dum) ) then                                                    
         call readsc                                                            
         go to 511                                                              
      else                                                                      
         go to 511                                                              
      end if                                                                    
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *                     a valid list of material properties at nodes   *        
c *                     have been input. store values for nodes        *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
 590  continue                                                                  
      call instore_fgm( intlst, lenlst, young_mod, poisson_ratio,               
     &                  alpha, vol_fract_ductile, rho, tan_e,                   
     &                  yld_pt, n_power, prop_flags )                           
c                                                                               
c                       all processed. examine another card for fgm             
c                       node data                                               
c                                                                               
      go to 505                                                                 
c                                                                               
c                                                                               
 9999 continue                                                                  
      sbflg1 = .true.                                                           
      sbflg2 = .true.                                                           
      call reset                                                                
      if ( true(dum) ) return                                                   
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine instore_fgm                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/26/01 mcw               *          
c     *                                                              *          
c     *     this subroutine stores the fgm material constants for    *          
c     *     a list of nodes                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine instore_fgm( intlst, lenlst, young_mod, poisson_ratio,         
     &                        alpha, vol_fract_ductile, rho, tan_e,             
     &                        yld_pt, n_power, prop_flags)                      
      use global_data ! old common.main
      use main_data, only: fgm_node_values, fgm_node_values_defined,            
     &                     fgm_node_values_cols                                 
      implicit integer (a-z)                                                    
      real dumr, young_mod, poisson_ratio, alpha, vol_fract_ductile,            
     &     rho, tan_e, yld_pt, n_power, zero                                    
      double precision                                                          
     &     dumd                                                                 
      logical prop_flags(*)                                                     
      character :: dums                                                         
      dimension intlst(*)                                                       
c                                                                               
c                       for each node in the list, set the                      
c                       node temporary storage array.                           
c                                                                               
      icn    = 0                                                                
      iplist = 1                                                                
 20   call trxlst( intlst, lenlst, iplist, icn, node )                          
c                                                                               
c                       check that the list element is not negative.            
c                                                                               
      if( node .lt. 0 ) then                                                    
         call errmsg( 58, node, dums, dumr, dumd )                              
         go to 30                                                               
      end if                                                                    
c                                                                               
c                       check that the list node does not exceed                
c                       the number of nodes in the structure.                   
c                                                                               
      if( node .gt. nonode ) then                                               
         call errmsg( 16, node, dums, dumr, dumd )                              
         go to 30                                                               
      end if                                                                    
c                                                                               
c                       store the node values of fgm material proprties         
c                                                                               
      if( prop_flags(1)) fgm_node_values(node,1) = young_mod                    
      if( prop_flags(2)) fgm_node_values(node,2) = poisson_ratio                
      if( prop_flags(3)) fgm_node_values(node,3) = alpha                        
      if( prop_flags(4)) fgm_node_values(node,4) = vol_fract_ductile            
      if( prop_flags(5)) fgm_node_values(node,5) = rho                          
      if( prop_flags(6)) fgm_node_values(node,6) = tan_e                        
      if( prop_flags(7)) fgm_node_values(node,7) = yld_pt                       
      if( prop_flags(8)) fgm_node_values(node,8) = n_power                      
c                                                                               
 30   if( iplist .ne. 0 ) go to 20                                              
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine infgm_dump                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/26/00                   *          
c     *                                                              *          
c     *     this subroutine prints the fgm material properties       *          
c     *     at the model nodes                                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine infgm_dump                                                     
      use global_data ! old common.main
      use main_data, only: fgm_node_values, fgm_node_values_defined             
      implicit integer (a-z)                                                    
      real young_mod, poisson_ratio, alpha, vol_fract_ductile,                  
     &     rho, tan_e, yld_pt, n_power, zero                                    
c                                                                               
c                       for each node, print the fgm material values            
c                                                                               
      if ( .not. fgm_node_values_defined ) then                                 
        write(out,9200)                                                         
        return                                                                  
      end if                                                                    
c                                                                               
      write(out,9100)                                                           
      do node = 1, nonode                                                       
        young_mod         = fgm_node_values(node,1)                             
        poisson_ratio     = fgm_node_values(node,2)                             
        alpha             = fgm_node_values(node,3)                             
        vol_fract_ductile = fgm_node_values(node,4)                             
        rho               = fgm_node_values(node,5)                             
        tan_e             = fgm_node_values(node,6)                             
        yld_pt            = fgm_node_values(node,7)                             
        n_power           = fgm_node_values(node,8)                             
        write(out,9000) node, young_mod, poisson_ratio, alpha,                  
     &                  vol_fract_ductile, rho, tan_e, yld_pt,                  
     &                  n_power                                                 
      end do                                                                    
c                                                                               
 9000 format(2x,i7,8e14.4)                                                      
 9100 format(//,'Dump of FGM properties defined at model nodes...',             
     & //,4x,'node',8x,'e',12x,'nu',12x,'alpha',                                
     &     4x,'vol. fract. ductile',3x, 'rho',10x,'tan_e',7x,                   
     &      'yld_pt',8x,'n_power' )                                             
 9200 format(//,'>>> no FGM properties defined at the model nodes....')         
c                                                                               
      return                                                                    
      end                                                                       
