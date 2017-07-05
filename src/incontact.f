c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine incontact                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/3/2017 rhd              *          
c     *                                                              *          
c     *        This routine processes the user input which describes *          
c     *        the contact surfaces.                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine incontact(sbflg1,sbflg2)                                       
      use global_data ! old common.main
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      logical numd, numi, sbflg1, sbflg2, rate, point_set, norm_set,            
     &    matchs, outside, inside                                                               
      character :: name*80, string*80, dums*1                                   
      double precision                                                          
     &     zero, dot, stiff, fric, point(3,3), mag1, mag2,                      
     &     tol_val, rate_val(3), depth, dumd, radius, length,                   
     &     ten_billion, norm(3)                                                 
      real dumr                                                                 
      data zero, tol_val, ten_billion /0.0, 0.0000001, 1.0e10 /                 
c                                                                               
c           read in new line                                                    
c                                                                               
 10   continue                                                                  
      call readsc                                                               
 20   continue                                                                  
c                                                                               
      if ( matchs('shape',4)) goto 100                                          
      if ( matchs('surface',4)) goto 100                                        
      if ( matchs('clear',5)) goto 1000                                         
      if ( matchs('dump',4)) goto 2000                                          
c                                                                               
      sbflg1         = .true.                                                   
      sbflg2         = .true.                                                   
c                                                                               
      goto 9999                                                                 
c                                                                               
c          input for contact shape definition                                   
c                                                                               
 100  continue                                                                  
c                                                                               
c             read shape number                                                 
c                                                                               
      if ( .not. numi(shape)) then                                              
         call errmsg (294, dumi, dums, dumr, dumd)                              
         goto 10                                                                
      endif                                                                     
c                                                                               
      if ( shape .lt. 1 .or. shape .gt. maxcontact) then                        
         call errmsg (295, maxcontact, dums, dumr, dumd)                        
         goto 10                                                                
      endif                                                                     
c                                                                               
c            initialize temporary shape variables                               
c                                                                               
      point_num = 0                                                             
      stiff = zero                                                              
      fric = zero                                                               
      depth = zero                                                              
      radius = zero                                                             
      length = zero                                                             
      point_set = .false.                                                       
      norm_set = .false.    
      outside = .true.  
      inside = .false.                                                  
      do i=1, 3                                                                 
         rate_val(i) = zero                                                     
         do j=1, 3                                                              
            point(i,j) = zero                                                   
         enddo                                                                  
      enddo                                                                     
c                                                                               
c            read shape type                                                    
c                                                                               
      if ( matchs ('plane',4)) goto 200                                         
      if ( matchs ('cylinder',3)) goto 300                                      
      if ( matchs ('sphere',3)) goto 400                                        
c                                                                               
      call errmsg (296, dumi, dums, dumr, dumd)                                 
      goto 10                                                                   
c                                                                               
c                                                                               
c          Case 1: input for contact plane definition                           
c          -------                                                              
c                                                                               
 200  continue                                                                  
c                                                                               
      call readsc                                                               
c                                                                               
      if ( matchs('point',5) ) goto 210                                         
      if ( matchs('contact',4) ) goto 200                                       
      if ( matchs('stiffness',5) ) goto 220                                     
      if ( matchs('friction',4) ) goto 230                                      
      if ( matchs('rate',4) ) goto 240                                          
      if ( matchs('velocity',3) ) goto 240                                      
      if ( matchs('depth',5) ) goto 250                                         
c                                                                               
c            if we don't match any keywords, we assume that the                 
c            plane definition is complete.                                      
c                                                                               
      goto 290                                                                  
c                                                                               
c        Input point for plane definition                                       
c                                                                               
 210  continue                                                                  
      point_num = point_num + 1                                                 
      if ( point_num .gt. 3) then                                               
         call errmsg (297, dumi, dums, dumr, dumd)                              
         goto 200                                                               
      endif                                                                     
c                                                                               
      do i = 1, 3                                                               
         if ( .not. numd (point(point_num,i))) then                             
            call errmsg (298, dumi, dums, dumr, dumd)                           
            point_num = point_num - 1                                           
            goto 200                                                            
         endif                                                                  
      enddo                                                                     
      goto 200                                                                  
c                                                                               
c        Input plane stiffness                                                  
c                                                                               
 220  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(stiff))                                                   
     &     call errmsg (299, dumi, dums, dumr, dumd)                            
      goto 200                                                                  
c                                                                               
c        Input plane friction coeffficient                                      
c                                                                               
 230  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(fric))                                                    
     &     call errmsg (300, dumi, dums, dumr, dumd)                            
      goto 200                                                                  
c                                                                               
c        Input rate of plane motion                                             
c                                                                               
 240  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      do i = 1, 3                                                               
         if ( .not. numd (rate_val(i))) then                                    
            call errmsg (301, dumi, dums, dumr, dumd)                           
            rate_val(1:3) = zero                                                
            goto 200                                                            
         endif                                                                  
      enddo                                                                     
      goto 200                                                                  
c                                                                               
c        Input depth of plane                                                   
c                                                                               
 250  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(depth))                                                   
     &     call errmsg (302, dumi, dums, dumr, dumd)                            
      goto 200                                                                  
c                                                                               
c          save this data as part of contact plane                              
c                                                                               
 290  continue                                                                  
c                                                                               
c             if we have too few points, then error and skip this plane.        
c             also, if stiffness is negative or zero, then skip this            
c             plane.                                                            
c                                                                               
      if ( point_num .lt. 3 ) then                                              
         call errmsg (303, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      else if ( stiff .le. zero ) then                                          
         call errmsg (304, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      endif                                                                     
c                                                                               
c             compute v1 and v2, edge vectors of rectangle.  Check to make      
c             sure the dot product is zero -- if not, produce error.            
c                                                                               
      dot = zero                                                                
      do i = 1, 3                                                               
         cshape_pnt(i,shape) = point(1,i)                                       
         cplane_vec(i,1,shape) = point(2,i) - point(1,i)                        
         cplane_vec(i,2,shape) = point(3,i) - point(1,i)                        
         dot = dot + cplane_vec(i,1,shape) * cplane_vec(i,2,shape)              
      enddo                                                                     
      if ( abs(dot) .gt. tol_val) then                                          
         call errmsg (305, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      endif                                                                     
c                                                                               
c             compute v3 -- normal vector                                       
c                                                                               
      call cross_prod ( cplane_vec(1,1,shape),cplane_vec(1,2,shape),            
     &     cshape_norm(1,shape) )                                               
      call normalize (cshape_norm(1,shape), dumd)                               
c                                                                               
c             store plane constats -- stiffness, friction, rate                 
c                                                                               
      contact_shape(shape) = 1                                                  
      contact_stiff(shape) = stiff                                              
      contact_fric(shape) = fric 
      contact_outside(shape) = .true.                                               
      if ( depth .ne. zero ) contact_depth(shape) = depth                       
      do i=1, 3                                                                 
         cshape_rate(i,shape) = rate_val(i)                                     
      enddo                                                                     
      use_contact = .true.                                                      
c                                                                               
      goto 20                                                                   
c                                                                               
c          Case 2: input for contact cylinder definition                        
c          -------                                                              
c                                                                               
 300  continue                                                                  
c                                                                               
      call readsc                                                               
c                                                                               
      if ( matchs('point',5) ) goto 310                                         
      if ( matchs('contact',4) ) goto 300                                       
      if ( matchs('center',4) ) goto 300                                        
      if ( matchs('stiffness',5) ) goto 320                                     
      if ( matchs('friction',4) ) goto 330                                      
      if ( matchs('rate',4) ) goto 340                                          
      if ( matchs('velocity',3) ) goto 340                                      
      if ( matchs('radius',3) ) goto 350                                        
      if ( matchs('length',3) ) goto 360                                        
      if ( matchs('direction',3) ) goto 370      
      if ( matchs('outside',3) ) go to  375
      if ( matchs('inside',3) ) go to 380                              
c                                                                               
c            if we don't match any keywords, we assume that the                 
c            cylinder definition is complete.                                      
c                                                                               
      goto 390                                                                  
c                                                                               
c        Input point on line along center of cylinder                           
c                                                                               
 310  continue                                                                  
c                                                                               
      point_set = .true.                                                        
      do i = 1, 3                                                               
         if ( .not. numd (point(i,1))) then                                     
            call errmsg (298, dumi, dums, dumr, dumd)                           
            point_set = .false.                                                 
            goto 300                                                            
         endif                                                                  
      enddo                                                                     
      goto 300                                                                  
c                                                                               
c        Input cylinder stiffness                                               
c                                                                               
 320  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(stiff))                                                   
     &     call errmsg (299, dumi, dums, dumr, dumd)                            
      goto 300                                                                  
c                                                                               
c        Input cylinder friction coeffficient                                   
c                                                                               
 330  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(fric))                                                    
     &     call errmsg (300, dumi, dums, dumr, dumd)                            
      goto 300                                                                  
c                                                                               
c        Input rate of motion of cylinder -- only translation is allowed        
c                                                                               
 340  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      do i = 1, 3                                                               
         if ( .not. numd (rate_val(i))) then                                    
            call errmsg (301, dumi, dums, dumr, dumd)                           
            rate_val(1:3) = zero                                                
            goto 300                                                            
         endif                                                                  
      enddo                                                                     
      goto 300                                                                  
c                                                                               
c        Input radius of cylinder                                               
c                                                                               
 350  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(radius))                                                  
     &     call errmsg (306, dumi, dums, dumr, dumd)                            
      goto 300                                                                  
c                                                                               
c        Input length of cylinder                                               
c                                                                               
 360  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(length))                                                  
     &     call errmsg (307, dumi, dums, dumr, dumd)                            
      goto 300                                                                  
c                                                                               
c        Input direction of center line of cylinder                             
c                                                                               
 370  continue                                                                  
      norm_set = .true.                                                         
      do i = 1, 3                                                               
         if ( .not. numd (norm(i))) then                                        
            call errmsg (308, dumi, dums, dumr, dumd)                           
            norm_set = .false.                                                  
            goto 300                                                            
         endif                                                                  
      enddo                                                                     
      goto 300    
c
c        enforce contact on outside (default) of cylinder or
c        inside                             
c                                                                               
 375  continue                                                                  
      outside = .true.
      inside = .false.
      go to 300
 380  continue
      inside = .true.
      outside = .false.
      go to 300                          
c                                                                               
c          save this data as part of cylinder                              
c                                                                               
 390  continue                                                                  
c                                                                               
c             if we don't have a valid point or normal, then skip this          
c             shape.                                                            
c             also, if stiffness is negative or zero, then skip this            
c             plane.                                                            
c                                                                               
      if (.not. point_set .or. .not. norm_set ) then                            
         call errmsg (309, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      else if ( stiff .le. zero ) then                                          
         call errmsg (304, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      else if ( radius .le. zero ) then                                         
         call errmsg (306, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      else if ( length .le. zero ) then                                         
         call errmsg (307, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      endif                                                                     
c                                                                               
c             store cylinder constants                                          
c                                                                               
      cshape_pnt (1:3,shape) = point(1:3,1)                                     
      cshape_norm (1:3,shape) = norm(1:3)                                       
      call normalize (cshape_norm (1:3,shape), dumd)                            
      cshape_param (1,shape) = radius                                           
      cshape_param (2,shape) = length                                           
      contact_shape(shape) = 2                                                  
      contact_stiff(shape) = stiff                                              
      contact_fric(shape) = fric  
      contact_outside(shape) = outside                                              
      if( depth .ne. zero ) contact_depth(shape) = depth                       
      cshape_rate(1:3,shape) = rate_val(1:3)                                    
      use_contact = .true.                                                      
c                                                                               
      goto 20                                                                   
c                                                                               
c          Case 3: input for contact sphere definition                          
c          -------                                                              
c                                                                               
 400  continue                                                                  
c                                                                               
      call readsc                                                               
c                                                                               
      if ( matchs('point',5) ) goto 410                                         
      if ( matchs('contact',4) ) goto 400                                       
      if ( matchs('center',4) ) goto 410                                        
      if ( matchs('stiffness',5) ) goto 420                                     
      if ( matchs('friction',4) ) goto 430                                      
      if ( matchs('rate',4) ) goto 440                                          
      if ( matchs('velocity',3) ) goto 440                                      
      if ( matchs('radius',3) ) goto 450                                        
      if ( matchs('outside',3) ) go to  460
      if ( matchs('inside',3) ) go to 465                             
c                                                                               
c            if we don't match any keywords, we assume that the                 
c            sphere definition is complete.                                     
c                                                                               
      goto 490                                                                  
c                                                                               
c        Input point on line along center of sphere                             
c                                                                               
 410  continue                                                                  
c                                                                               
      point_set = .true.                                                        
      do i = 1, 3                                                               
         if ( .not. numd (point(i,1))) then                                     
            call errmsg (297, dumi, dums, dumr, dumd)                           
            point_set = .false.                                                 
            goto 400                                                            
         endif                                                                  
      enddo                                                                     
      goto 400                                                                  
c                                                                               
c        Input sphere stiffness                                                 
c                                                                               
 420  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(stiff))                                                   
     &     call errmsg (299, dumi, dums, dumr, dumd)                            
      goto 400                                                                  
c                                                                               
c        Input sphere friction coeffficient                                     
c                                                                               
 430  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(fric))                                                    
     &     call errmsg (300, dumi, dums, dumr, dumd)                            
      goto 400                                                                  
c                                                                               
c        Input rate of motion of sphere -- only translation is allowed          
c                                                                               
 440  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      do i = 1, 3                                                               
         if ( .not. numd (rate_val(i))) then                                    
            call errmsg (301, dumi, dums, dumr, dumd)                           
            rate_val(1:3) = zero                                                
            goto 400                                                            
         endif                                                                  
      enddo                                                                     
      goto 400                                                                  
c                                                                               
c        Input radius of sphere                                                 
c                                                                               
 450  continue                                                                  
      if (matchs('=',1)) call splunj                                            
      if ( .not. numd(radius))                                                  
     &     call errmsg (306, dumi, dums, dumr, dumd)                            
      goto 400                                                                  
c
c        enforce contact on outside (default) of sphere or
c        inside                             
c                                                                               
 460  continue                                                                  
      outside = .true.
      inside = .false.
      go to 400
 465  continue
      inside = .true.
      outside = .false.
      go to 400                          
c                                                                               
c          save this data as part of contact sphere                             
c                                                                               
 490  continue                                                                  
c                                                                               
c             if we don't have a valid point or normal, then skip this          
c             shape.                                                            
c             also, if stiffness is negative or zero, then skip this            
c             sphere.                                                           
c                                                                               
      if (.not. point_set )then                                                 
         call errmsg (309, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      else if ( stiff .le. zero ) then                                          
         call errmsg (304, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      else if ( radius .le. zero ) then                                         
         call errmsg (306, dumi, dums, dumr, dumd)                              
         goto 20                                                                
      endif                                                                     
c                                                                               
c             store sphere constants                                            
c                                                                               
      cshape_pnt (1:3,shape) = point(1:3,1)                       
      cshape_param (1,shape) = radius                                           
      contact_shape(shape) = 3                                                  
      contact_stiff(shape) = stiff                                              
      contact_fric(shape) = fric                                                
      cshape_rate(1:3,shape) = rate_val(1:3)   
      contact_outside(shape) = outside                                 
      use_contact = .true.                                                      
c                                                                               
      goto 20                                                                   
c                                                                               
c          keyword "clear" was found.  zero out the contact plane               
c          definitions.                                                         
c                                                                               
 1000 continue                                                                  
      write (out,*) '>>> Clearing contact surfaces'                             
c                                                                               
      call contact_remove(.true.)                                               
c                                                                               
      do i=1, maxcontact                                                        
         cplane_vec(1:3,1:2,i) = zero                                           
         cshape_norm(1:3,i) = zero                                              
         cshape_pnt(1:3,i) = zero                                               
         cshape_rate(1:3,i) = zero                                              
         cshape_param(1:maxcntprm,i) = zero                                     
         contact_stiff(i) = zero                                                
         contact_fric(i) = zero                                                 
         contact_force(i) = zero                                                
         contact_shape(i) = 0                                                   
         contact_depth(i) = ten_billion   
         contact_outside(i) = .true.                                      
      end do                                                                     
c                                                                               
      use_contact = .false.                                                     
c                                                                               
      goto 10                                                                   
c                                                                               
c                                                                               
c          keyword "dump" was found.  zero out the contact plane                
c          definitions.                                                         
c                                                                               
 2000 continue                                                                  
      write (out,*) '>>> Dumping contact surfaces...'                           
      do shape = 1, maxcontact                                                  
c                                                                               
         if ( contact_shape(shape) .eq. 1 ) then                                
c                                                                               
            write (out,*) ' --> Shape ',shape,' is a sub-plane'                 
            write (out,*) '    corner point:'                                   
            write (out,'(7x,3e14.6)')(cshape_pnt(i,shape),i=1,3)                
            write (out,*) '    edge vectors:'                                   
            write (out,'(7x,3e14.6)')((cplane_vec(i,j,shape),i=1,3),            
     &              j=1,2)                                                      
            write (out,*) '    normal vector:'                                  
            write (out,'(7x,3e14.6)')(cshape_norm(i,shape),i=1,3)               
            write (out,*) '    velocity:'                                       
            write (out,'(7x,3e14.6)')(cshape_rate(i,shape),i=1,3)               
            write (out,*) '    stiffness:', contact_stiff(shape)                
            write (out,*) '    friction:', contact_fric(shape)                  
            write (out,*) '    depth:', contact_depth(shape)                    
c                                                                               
         else if ( contact_shape(shape) .eq. 2 ) then                           
c                                                                               
            write (out,*) ' --> Shape ',shape,' is a cylinder'                  
            write (out,*) '    center point:'                                   
            write (out,'(7x,3e14.6)')(cshape_pnt(i,shape),i=1,3)                
            write (out,*) '    center line direction:'                          
            write (out,'(7x,3e14.6)')(cshape_norm(i,shape),i=1,3)               
            write (out,*) '    velocity:'                                       
            write (out,'(7x,3e14.6)')(cshape_rate(i,shape),i=1,3)               
            write (out,*) '    radius:', cshape_param(1,shape)                  
            write (out,*) '    length:', cshape_param(2,shape)                  
            write (out,*) '    stiffness:', contact_stiff(shape)                
            write (out,*) '    friction:', contact_fric(shape)
            write (out,*) '    outside: ', contact_outside(shape)                  
c                                                                               
         else if ( contact_shape(shape) .eq. 3 ) then                           
c                                                                               
            write (out,*) ' --> Shape ',shape,' is a sphere'                    
            write (out,*) '    center point:'                                   
            write (out,'(7x,3e14.6)')(cshape_pnt(i,shape),i=1,3)                
            write (out,*) '    velocity:'                                       
            write (out,'(7x,3e14.6)')(cshape_rate(i,shape),i=1,3)               
            write (out,*) '    radius:', cshape_param(1,shape)                  
            write (out,*) '    stiffness:', contact_stiff(shape)                
            write (out,*) '    friction:', contact_fric(shape)                  
            write (out,*) '    outside: ', contact_outside(shape)                  
c                                                                               
         endif                                                                  
c                                                                               
      enddo                                                                     
c                                                                               
      goto 10                                                                   
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
c                                                                               
                                                                                
