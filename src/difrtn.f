c *******************************************************************           
c *                                                                 *           
c *   support routine: compute tangent to crack front               *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine difrtn( nfnode, fnodes, caseno, scoord, coord_map,             
     &                   tvec, fnctyp, crack_front_tangent,                     
     &                   tangent_vector_defined, iout, debug, status  )         
      implicit integer (a-z)                                                    
c                                                                               
c               parameter declarations                                          
c                                                                               
      double precision                                                          
     & scoord(*), tvec(*), crack_front_tangent(*)                               
      dimension coord_map(*), fnodes(*)                                         
      logical tangent_vector_defined                                            
c                                                                               
c               local declarations                                              
c                                                                               
      double precision                                                          
     &  dsf(4), sf(4), tvec1(3), tvec2(3),                                      
     &  dx, dy, dz, xsi, zero, one, half                                        
      logical debug                                                             
      data zero, one, half / 0.0d00, 1.0d00, 0.5d00 /                           
c                                                                               
c             compute or load unit tangent vector to the crack                  
c             front according to the situation indicated by the                 
c             case number or the user-defined flag.                             
c             for cases 2,5,6 we average two tangents                           
c             computed for adjacent line segments. for cases 1,3,4 we           
c             compute a single tangent at start or end of line segment.         
c             fnodes is list of structure nodes defining the segment(s).        
c                                                                               
c             tangent_vector_defined = .T. --> the user specified the           
c                                              tangent vector to use.           
c                                              components are in                
c                                              crack_front_tangent              
c                                              skip computations.               
c                                                                               
c                                                                               
c             case no.           situation                                      
c             --------           ---------                                      
c                1               linear; 2-nodes on front                       
c                2               linear; 3-nodes on front                       
c                3               quadratic; 3-nodes on front                    
c                4               cubic; 4-nodes on front                        
c                5               quadratic; 5-nodes on front                    
c                6               cubic; 7-nodes on front                        
c                                                                               
c             fnctyp: 1-4, used for cases 3,4 to determine if                   
c                     tangent is needed for first or last node                  
c                     on front. = 1 (first node), =3 (last node).               
c                                                                               
c             status: returned value =0 (ok) =1 (bad data)                      
c                                                                               
c                                                                               
c                                                                               
      if( tangent_vector_defined ) then                                         
        dx = crack_front_tangent(1)                                             
        dy = crack_front_tangent(2)                                             
        dz = crack_front_tangent(3)                                             
        call difrts( dx, dy, dz, tvec, status )                                 
        go to 1000                                                              
      end if                                                                    
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600 ), caseno                            
c                                                                               
c             case 1,3,4: 2, 3 or four nodes on front but only a single         
c                         element (linear, quadratic or cubic). compute         
c                         tangent at first or last node (2,3 or 4). use         
c                         utility routine to get derivatives of the 1-D         
c                         isoparametric shape functions.                        
c                                                                               
 100  continue                                                                  
 300  continue                                                                  
 400  continue                                                                  
      if ( fnctyp .eq. 1 ) then                                                 
         xsi = -one                                                             
       elseif ( fnctyp .eq. 3 ) then                                            
         xsi = one                                                              
       else                                                                     
         write(iout,9100) fnctyp                                                
         call die_abort                                                         
         stop                                                                   
      end if                                                                    
      call di1dsf( xsi, dsf, sf, nfnode )                                       
      dx = zero                                                                 
      dy = zero                                                                 
      dz = zero                                                                 
      do ii = 1, nfnode                                                         
        coord_loc  = coord_map(fnodes(ii)) - 1                                  
        dx = dx + dsf(ii)*scoord(coord_loc+1)                                   
        dy = dy + dsf(ii)*scoord(coord_loc+2)                                   
        dz = dz + dsf(ii)*scoord(coord_loc+3)                                   
      end do                                                                    
      call difrts( dx, dy, dz, tvec, status )                                   
      go to 1000                                                                
c                                                                               
c             case 2,5,6: two linear, quadratic or cubic segments               
c                         connecting 3, 5 or 6 front nodes.                     
c                         compute average tangent at center (common)            
c                         node. compute two unit tangents, average              
c                         components, then restore unit length.                 
c                                                                               
 200  continue                                                                  
 500  continue                                                                  
 600  continue                                                                  
      nfn = 2                                                                   
      if ( caseno .eq. 5) nfn = 3                                               
      if ( caseno .eq. 6) nfn = 4                                               
      xsi = one                                                                 
      call di1dsf( xsi, dsf, sf, nfn )                                          
      dx = zero                                                                 
      dy = zero                                                                 
      dz = zero                                                                 
      do ii = 1, nfn                                                            
        coord_loc = coord_map(fnodes(ii)) - 1                                   
        dx = dx + dsf(ii)*scoord(coord_loc+1)                                   
        dy = dy + dsf(ii)*scoord(coord_loc+2)                                   
        dz = dz + dsf(ii)*scoord(coord_loc+3)                                   
      end do                                                                    
      call difrts( dx, dy, dz, tvec1, status )                                  
      if ( status .eq. 1 ) return                                               
      xsi = -one                                                                
      call di1dsf( xsi, dsf, sf, nfn )                                          
      dx = zero                                                                 
      dy = zero                                                                 
      dz = zero                                                                 
      jj = nfn - 1                                                              
      do ii = 1, nfn                                                            
        coord_loc = coord_map(fnodes(jj+ii)) - 1                                
        dx = dx + dsf(ii)*scoord(coord_loc+1)                                   
        dy = dy + dsf(ii)*scoord(coord_loc+2)                                   
        dz = dz + dsf(ii)*scoord(coord_loc+3)                                   
      end do                                                                    
      call difrts( dx, dy, dz, tvec2, status )                                  
      if ( status .eq. 1 ) return                                               
      dx = (tvec1(1) + tvec2(1)) * half                                         
      dy = (tvec1(2) + tvec2(2)) * half                                         
      dz = (tvec1(3) + tvec2(3)) * half                                         
      call difrts( dx, dy, dz, tvec, status )                                   
      go to 1000                                                                
c                                                                               
c             all done.                                                         
c                                                                               
 1000 continue                                                                  
      if ( .not. debug ) return                                                 
      write(iout,9200)  tangent_vector_defined,  caseno, fnctyp,                
     &                  status, tvec(1), tvec(2), tvec(3)                       
      return                                                                    
c                                                                               
 9100 format('>>>>> SYSTEM ERROR:in difrtn. fnctyp has illegal',                
     &  /,   '                   value: ',i10,                                  
     &  /,   '                   job terminated.')                              
 9200 format('      > leaving difrtn.',                                         
     & /,20x,'tangent_vector_defined: ',l1,                                     
     & /,20x,'caseno:       ',i10,                                              
     & /,20x,'domain_type:  ',i10,                                              
     & /,20x,'status:       ',i10,                                              
     & /,20x,'tangent vec:  ',3f10.6)                                           
c                                                                               
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *   support routine: called by difrtn                             *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine difrts( dx, dy, dz, tvec, status )                             
      double precision                                                          
     & dx, dy, dz, tvec(*), length                                              
      integer status                                                            
c                                                                               
      length = sqrt( dx*dx + dy*dy + dz*dz )                                    
      if ( length .lt. 1.0e-07 ) then                                           
          status = 1                                                            
          return                                                                
      end if                                                                    
      tvec(1) = dx / length                                                     
      tvec(2) = dy / length                                                     
      tvec(3) = dz / length                                                     
      status = 0                                                                
c                                                                               
      return                                                                    
      end                                                                       
