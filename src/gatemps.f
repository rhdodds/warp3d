c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gartemps                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 10/25/00                   *          
c     *                                                              *          
c     *     this subroutine creates a table of nodal temperature     *          
c     *     values at the element nodes in the block for the user    *          
c     *     specified reference (initial) temperatures               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gartemps( rtemp_nodes, belinc, nnode,                          
     &                     span, felem, rtemps_node_blk, mxvl  )                
      implicit integer (a-z)                                                    
      double precision                                                          
     &      rtemp_nodes(*), rtemps_node_blk(mxvl,*)                             
      integer belinc(nnode,*)                                                   
      logical local_debug                                                       
      data local_debug / .false./                                               
c                                                                               
c           for each element in block gather the reference temperature          
c           at element nodes from structure vector.                             
c                                                                               
c                                                                               
      do j = 1, nnode                                                           
        do i = 1, span                                                          
         rtemps_node_blk(i,j) = rtemp_nodes(belinc(j,i))                        
        end do                                                                  
      end do                                                                    
c                                                                               
      if ( local_debug ) then                                                   
          write(*,*) ' '                                                        
          write(*,*) '>> gartemps...'                                           
          do i = 1, span                                                        
           write(*,9000) i+felem-1, (rtemps_node_blk(i,j),j=1,nnode)            
          end do                                                                
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(i5,10(4f10.2,/,5x))                                                
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gadtemps                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/9/00                    *          
c     *                                                              *          
c     *     this subroutine creates a table of nodal temperature     *          
c     *     increments for elements in this block. we add the        *          
c     *     constant element temperature change to the nodal changes *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gadtemps( dtemp_nodes, dtemp_elems, belinc, nnode,             
     &                     span, felem, dtemps_node_blk, mxvl  )                
      implicit integer (a-z)                                                    
      double precision                                                          
     &      dtemp_nodes(*), dtemp_elems(*), dtemps_node_blk(mxvl,*)             
      integer belinc(nnode,*)                                                   
      logical local_debug                                                       
      data local_debug / .false./                                               
c                                                                               
c           for each element in block:                                          
c             a) get element constant temperature change                        
c             b) gather the temperature change at element nodes                 
c                from structure vector.                                         
c             c) define total temperature change at each element                
c                node.                                                          
c                                                                               
c                                                                               
      do j = 1, nnode                                                           
        do i = 1, span                                                          
         dtemps_node_blk(i,j) = dtemp_nodes(belinc(j,i)) +                      
     &                          dtemp_elems(i)                                  
        end do                                                                  
      end do                                                                    
c                                                                               
      if ( local_debug ) then                                                   
          write(*,*) ' '                                                        
          write(*,*) '>> gatemps...'                                            
          do i = 1, span                                                        
           write(*,9000) i+felem-1, (dtemps_node_blk(i,j),j=1,nnode)            
          end do                                                                
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(i5,10(4f10.2,/,5x))                                                
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine gatemps                      *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/9/00                    *          
c     *                                                              *          
c     *     this subroutine creates a table of nodal temperatures    *          
c     *     at end of the step for elements in this block. includes  *          
c     *     initial temperatures, accumulated nodal and element      *          
c     *     temperatures (reference temps were loaded into nodal     *          
c     *     temps at input time)                                     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine gatemps( temp_nodes, temp_elems, belinc, nnode,                
     &                    span, felem, temps_node_blk, mxvl,                    
     &                    dtemps_node_blk, temps_node_to_process )              
      implicit integer (a-z)                                                    
c                                                                               
c                    parameter declarations                                     
c                                                                               
      double precision                                                          
     &      temp_nodes(*), temp_elems(*), temps_node_blk(mxvl,*),               
     &      dtemps_node_blk(mxvl,*)                                             
      integer belinc(nnode,*)                                                   
      logical temps_node_to_process                                             
c                                                                               
c                    local declarations                                         
c                                                                               
      double precision                                                          
     &  zero                                                                    
      logical local_debug                                                       
      data local_debug, zero / .false., 0.0 /                                   
c                                                                               
c           for each element in block:                                          
c             a) get element constant temperature                               
c             b) gather the temperature at element nodes                        
c                from structure vector.                                         
c             c) define total temperature at end of step for each element       
c                node.                                                          
c                                                                               
c                                                                               
      temps_node_to_process = .false.                                           
      do j = 1, nnode                                                           
        do i = 1, span                                                          
         temps_node_blk(i,j) = temp_nodes(belinc(j,i)) +                        
     &                         temp_elems(i) + dtemps_node_blk(i,j)             
         if ( abs( temps_node_blk(i,j) ) .gt. zero )                            
     &       temps_node_to_process = .true.                                     
        end do                                                                  
      end do                                                                    
c                                                                               
      if ( local_debug ) then                                                   
          write(*,*) ' '                                                        
          write(*,*) '>> gatemps...'                                            
          do i = 1, span                                                        
           write(*,9000) i+felem-1, (temps_node_blk(i,j),j=1,nnode)             
          end do                                                                
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(i5,10(4f10.2,/,5x))                                                
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
