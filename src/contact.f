c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine contact_find                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/3/2017 rhd              *          
c     *                                                              *          
c     *        This routine checks all nodes against all contact     *          
c     *        shapes to see if contact has occurred.  If contact is *          
c     *        found, then the contact force for the contacting      *          
c     *        node is calculated using the penetration and the      *          
c     *        shape stiffness. The contact force on all other       *          
c     *        dofs is zero.                                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine contact_find                                                   
      use global_data ! old common.main
c                                                                               
      use main_data, only : trn, trnmat                                         
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision ::                                                     
     &     pen_dist(maxcontact), zero, force, normal(3,maxcontact),             
     &     transmat(3,3), dumd, new_pen_sign   
      double precision :: pen_sign(maxcontact)                                         
      real dumr                                                                 
      character(len=1) :: dums                                                  
      logical debug, penetrated, allow_trn, referenced, owned                   
      data debug, zero / .false., 0.0 /                                         
c                                                                               
c         tell slaves we are about to assess contact                            
c                                                                               
      call wmpi_alert_slaves ( 29 )                                             
c                                                                               
c          zero out all the contact forces.                                     
c                                                                               
      do dof = 1, nodof                                                         
         contact_force (dof) = zero                                             
      enddo                                                                     
c                                                                               
c          skip this routine if we dont have any contact.                       
c                                                                               
      if ( .not. use_contact ) then                                             
         if (debug) write (*,*) '>>>>> Skipping contact check.'                 
         return                                                                 
      endif                                                                     
c                                                                               
c          if we have supplied any transformation matrices for nodes            
c          undergoing contact, un-transform those nodes so that                 
c          all references are in global coordinates.                            
c                                                                               
      call contact_remove(.false.)                                              
c                                                                               
c          loop over all nodes and check for contact                            
c                                                                               
      if (debug) then                                                           
         write (*,*) '>>>> checking for Contact!'                               
         write (*,*) '>>     dt is:', dt                                        
         write (*,*) '>>     I think coords of shape 1 is:'                     
         write (*,*) '    ',(cshape_pnt(i,1)+cshape_rate(i,1)*dt,i=1,3)         
         write (*,*) '>>     penetrating nodes, pen dist, force'                
      endif                                                                     
c                                                                               
      do node = 1, nonode !  big outside loop                                                      
c                                                                               
c             skip node if not referenced by this processor                     
c                                                                               
c                                                                               
         call wmpi_chknode ( node, referenced, owned )                          
c                                                                               
         if ( .not. referenced ) cycle                                          
c                                                                               
c             check if node contacts with one or more of the contact shapes.    
c             if it does, list the contacting surfaces ranked by decreasing     
c             penetration in contact_cause(xxx,node). Return the penetration    
c             and direction of contact for all shapes contacting this node.     
c                                                                               
         call contact_find_node ( node, pen_dist, pen_sign, normal )                       
c                                                                               
c             if no contact was found, cycle to next node.                      
c                                                                               
         if ( contact_cause(1,node).eq. 0 ) cycle                               
c                                                                               
c             contact has been found. We now must deal with the                 
c             transformation of the node.  In order to improve the              
c             conditioning of the stiffness matrix, we want to rotate           
c             the coordinate system of the node to be orthogonal to             
c             the normal of the contact shape.  There are only a few cases      
c             where this is possible/useful, however:                           
c                - if there is already a nodal coordinate transformation        
c                   applied to the node through the constraints input,          
c                   then we can't apply another transformation. Also,           
c                   computing the contact force becomes tricky. Just            
c                   print an error and stop execution of the program.           
c                - if the contact normal is already orthogonal to the           
c                   global coordinate axes, then there is no point              
c                   to applying a transformation; skip it and continue on.      
c                - if the normal is not orthognal to the global axes,           
c                   then we must consider the constraints.                      
c             If a transformation matrix is calculated, then apply the          
c             matrix to all data which stores dof information -- u, v, a,       
c             loads, etc.                                                       
c                                                                               
c                  assess whether or not we can apply the nodal coordinate      
c                  transformation.                                              
c                                                                               
         if ( trn(node) ) then                                                  
            call errmsg (306, node, dums, dumr, dumd)                           
            call die_abort                                                      
            stop                                                                
         endif                                                                  
c                                                                               
c                  compute a transformation matrix for the node,                
c                  if possible.                                                 
c                                                                               
         call contact_trnmat ( node, normal(1,contact_cause(1,node)),           
     &        transmat, allow_trn )                                             
c                                                                               
c                  if transformation was found, then                            
c                  put transformation matrix for contacting                     
c                  node into the general data structure for                     
c                  all transformation matricies (trnmat).                       
c                                                                               
         if ( allow_trn ) then                                                  
c                                                                               
            trn(node) = .true.                                                  
            call allo_trnmat(node,1, dum)                                       
            do i=1, 3                                                           
               do j = 1,3                                                       
                  trnmat(node)%mat(i,j) = transmat(i,j)                         
               enddo                                                            
            enddo                                                               
c                  Modify all vectors that are                                  
c                  effected by a change in the transformation                   
c                  matricies; u, du, idu, v, a, forces, etc.                    
c                                                                               
            call trn2all (node, 1)                                              
c                                                                               
         endif                                                                  
c                                                                               
c             Now we need to compute the contact force contribution             
c             from each valid contact shape. If we have more than one           
c             contact shape currently acting upon the node, then sum            
c             the force from each shape into the contact force for the          
c             node.                                                             
c                                                                               
c                 If this processor does not own the node, then                 
c                 don't calculate the contact force.                            
c                                                                               
         if ( .not. owned ) then                                                
            if (debug) write (*,*) myid,': ---->',node,                         
     &           ' but I dont own it.'                                          
            cycle                                                               
         endif                                                                  
c                                                                               
c                 Loop over the contact shapes.                                 
c                                                                               
         do loop = 1, num_contact                                               
c                                                                               
            cause = contact_cause(loop,node)                                    
            if ( cause .eq. 0 ) exit                                            
c                                                                               
c                    calculate the contact force normal to contact shape.       
c                                                                               
            force = pen_dist(cause) * contact_stiff(cause) *
     &              pen_sign(cause)                      
c                                                                               
c                   if transformation is allowed, then rotate                   
c                   force contribution to transformed coordinates.              
c                                                                               
            if ( allow_trn ) then                                               
c                                                                               
               do i = 1, 3                                                      
                  do j = 1, 3                                                   
                     contact_force( dstmap(node)+i-1 ) =                        
     &                    contact_force( dstmap(node)+i-1 ) +                   
     &                    force * normal(j,cause) * transmat(i,j)               
                  enddo                                                         
               enddo                                                            
c                                                                               
c                   if transformation is not allowed, then contact force        
c                   is applied in global coordinates.                           
c                                                                               
            else                                                                
c                                                                               
               do i=1, 3                                                        
                  contact_force( dstmap(node)+i-1 ) =                           
     &                 contact_force( dstmap(node)+i-1 ) +                      
     &                 force * normal(i,cause)                                  
               enddo                                                            
c                                                                               
            endif                                                               
c                                                                               
         enddo                                                                  
c                                                                               
c               if debug, output calculated contact force.                      
c                                                                               
         if (debug) write (*,'(i2,": ---->",i7,4e13.6)')myid,node,              
     &        pen_dist ( contact_cause (1,node)),                               
     &        (contact_force (dstmap(node)+i-1),i=1,3)                          
c                                                                               
      end do   ! over nonode                                                                
c                                                                               
c               now reduce the contact_force vector and the nodal               
c               transformation matrices caused by contact back to the           
c               root processor                                                  
c                                                                               
      call wmpi_contact_gthr                                                    
c                                                                               
      if ( debug ) then                                                         
         if ( root_processor) then                                              
            write (*,*) '>> NONZERO CONTACT FORCE TERMS'                        
            do i=1, nodof                                                       
               if ( contact_force(i) .ne. zero ) then                           
                  write (*,'(6x,i7,3x,e13.6)')i,contact_force(i)                
               endif                                                            
            enddo                                                               
         endif                                                                  
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine contact_remove               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/17/98                   *          
c     *                                                              *          
c     *        This routine un-transforms all nodes which were       *          
c     *        rotated into contact-compatible coordinates.  All     *          
c     *        vectors bearing dof data are rotated back to global   *          
c     *        coordinates.                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine contact_remove(call_others)                                    
      use global_data ! old common.main
c                                                                               
      use main_data, only : trn, trnmat, cnstrn                                 
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     dumvec(3)                                                            
      logical debug, ldum, call_others                                          
      data debug /.false./                                                      
c                                                                               
c          if this subroutine was called from a place other than                
c          contact_find, then alert the slaves to remove effects of             
c          nodal transformation matricies due to contact.                       
c                                                                               
      if ( call_others ) call wmpi_alert_slaves ( 11 )                          
c                                                                               
c          skip this routine if we dont have any contact.                       
c                                                                               
      if ( .not. use_contact ) return                                           
c                                                                               
c          loop over each node in the structure; cycle if it has no contact.    
c                                                                               
      do node = 1, nonode                                                       
c                                                                               
         if ( contact_cause(1,node) .ne. 0) then                                
c                                                                               
c                 if node was transformed into contact-compatible               
c                 coordinates, rotate it back to global coordinates             
c                 and deallocate the transformation matrix.                     
c                                                                               
            if (trn(node)) then                                                 
               if ( debug ) write (*,*) myid,':-> undoing trn on node ',        
     &              node                                                        
               call trn2all (node, 2)                                           
               call allo_trnmat(node,2,dum)                                     
               trn(node) = .false.                                              
            endif                                                               
c                                                                               
c                 set cause of contact for all contacting nodes to be           
c                 zero.                                                         
c                                                                               
            contact_cause (1:maxcontact,node) = 0                               
         endif                                                                  
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine contact_find_node            *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/17/98                   *          
c     *                                                              *          
c     *        This routine checks a node against all contact        *          
c     *        shapes to see if contact has occurred.  If more than  *          
c     *        one contacting shape is found, then this routine      *          
c     *        sorts them to identify which are important, then      *          
c     *        stores the list of contacting shapes in decreasing    *          
c     *        penetration order in contact_cause.                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine contact_find_node ( node, pen_dist, pen_sign, 
     &                               normal )                   
      use global_data ! old common.main
c                                                                               
      use main_data, only : crdmap, trn, trnmat                                 
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision :: pen_sign(*), new_pen_sign,                                                         
     &     pen_dist(*), zero, curr_coord(3),                                    
     &     normal(3,*), tmp_coord(3),                                           
     &     dumvec(3), new_dist, curr_normal(3), curr_dist,                      
     &     tmp_coord2(3)                                                        
      logical debug, penetrated, found_entry, rebuild_norm, first,              
     &     new_coords                                                           
      dimension pen_list(maxcontact)                                            
      integer i, node, j                                                        
      data debug, zero /.false., 0.0/                                           
c                                                                               
c             get coords of node, initialize vars for this node                 
c                                                                               
      call get_coords ( node, curr_coord )                                      
c                                                                               
      if ( debug) then                                                          
         write (*,*) '               >> check node ',node                       
         write (*,'(40x," init, u, du, coords:")')                              
         do i=1,3                                                               
            write (*,'(40x,4e13.6)') c(crdmap(node)+i-1),                       
     &           u(dstmap(node)+i-1),du(dstmap(node)+i-1),                      
     &           curr_coord(i)                                                  
         enddo                                                                  
      endif                                                                     
c                                                                               
      num_pen = 0                                                               
c                                                                               
      do i=1, maxcontact                                                        
         pen_list(i) = 0                                                        
      enddo                                                                     
c                                                                               
c             loop over all contact planes                                      
c                                                                               
      do shape = 1, num_contact                                                 
c                                                                               
c                find penetration for this contact shape.                       
c                if no penetration is found, then go to next shape.             
c                                                                               
         call find_shape_contact( curr_coord, shape, pen_dist(shape),           
     &        pen_sign(shape), penetrated, normal(1,shape), node )                               
         if (.not. penetrated) cycle                                            
c                                                                               
c                penetration has been found. Insert this shape into list of     
c                penetrated shapes.  Sort the list so that smallest             
c                penetration distance is first.                                 
c                                                                               
         num_pen = num_pen + 1                                                  
         if ( num_pen .eq. 1) then                                              
            pen_list(num_pen) = shape                                           
         else                                                                   
            found_entry = .false.                                               
            do i=1, num_pen - 1                                                 
               if ( pen_dist(pen_list(i)) .gt. pen_dist(shape)) then            
                  do j = num_pen, i+1 , -1                                      
                     pen_list(j) = pen_list(j-1)                                
                  enddo                                                         
                  pen_list(i) = shape                                           
                  found_entry = .true.                                          
                  exit                                                          
               endif                                                            
            enddo                                                               
            if ( .not. found_entry ) pen_list(num_pen) = shape                  
         endif                                                                  
c                                                                               
      enddo                                                                     
c                                                                               
c             If no contact was found, return                                   
c                                                                               
      if ( num_pen .eq. 0) then                                                 
         goto 9999                                                              
      endif                                                                     
c                                                                               
c             print order of penetration shapes for debugging                   
c                                                                               
      if (debug) then                                                           
         write (*,*) '   >>> penetration shapes for node ',node                 
         do i=1, num_pen                                                        
            write (*,*) '     shape:',pen_list(i),' pen:',                      
     &           pen_dist(pen_list(i))                                          
         enddo                                                                  
      endif                                                                     
c                                                                               
c             Now run through the penetration list to determine which           
c             of the penetrated shapes superceed others.                        
c             For instance, if a problem contains two                           
c             contact planes at an external 90 degree angle, then               
c             the node should be moved to the plane with the smaller            
c             contact distance.  However, if the planes are parallel,           
c             then the node should be moved to the furthest contact             
c             shape.                                                            
c                                                                               
c             To assess this, we compare each pair of nodes using the           
c             following scheme:                                                 
c                                                                               
c              - make the one with the larger penetration "deep", and           
c                   the other "shallow"                                         
c              - assume that "shallow" is the primary contact shape.            
c                   calculate the new location of the node on the surface       
c                   of the "shallow" shape.                                     
c              - now test if we still violate the "deep" contact shape.         
c              - if "deep" is no longer violated, then "shallow" is             
c                   the important contact shape, and deep can be ignored.       
c                   go to evaluation of next shape.                             
c              - if "deep" is still violated, then "deep" must still be         
c                   considered.  However, we still need to check the            
c                   "shallow" shape.                                            
c              - find new location for node on surface of "deep" that           
c                   no longer produces penetration.                             
c              - test if we still violate the "shallow" condition using the     
c                   new location.                                               
c              - if we still violate "shallow", then both shapes are            
c                   important.  Otherwise, remove "shallow" from pen_list.      
c                                                                               
c             If more than one contact shape is still valid after               
c             a comparison, then find the location which prevents               
c             penetration of any of the valid shapes.  This new location        
c             is then treated as the "shallow" location in future               
c             comparisions.                                                     
c                                                                               
c             See the users manual for additional description.                  
c                                                                               
c                Start by setting the current normal and penetration equal      
c                to the contact shape with the smallest penetration.            
c                                                                               
      do i=1, 3                                                                 
         curr_normal(i) = normal(i,pen_list(1))                                 
      enddo                                                                     
      curr_dist = pen_dist (pen_list(1))                                        
      new_coords = .true.                                                       
c                                                                               
c                Now loop over the remaining shapes.  Evaluate which contact    
c                shapes superceed others.  As we loop through the pen_list      
c                array, we set the corresponding pen_list entry to zero if      
c                the shape is superceeded by others.                            
c                                                                               
      do loop = 2, num_pen                                                      
c                                                                               
c                    Find where node would end up after application of          
c                    all of the currently approved contact surfaces. Note       
c                    that this location is calculated later in this routine,    
c                    using the subroutine correct_normal.                       
c                                                                               
         if ( new_coords) then                                                  
            do i = 1, 3                                                         
               tmp_coord(i) = curr_coord(i) + curr_dist * curr_normal(i)        
            enddo                                                               
            new_coords = .false.                                                
         endif                                                                  
c                                                                               
c                    find if new location still violates current                
c                    contact shape.                                             
c                                                                               
         call find_shape_contact( tmp_coord, pen_list(loop), new_dist,          
     &        new_pen_sign, penetrated, dumvec, node )                                        
c                                                                               
c                    if no penetration is found, then current contact           
c                    shape is no longer a valid contact surface. Set            
c                    corresponding entry in pen_list to zero and                
c                    move to next penetrated shape.                             
c                                                                               
         if ( .not. penetrated ) then                                           
            pen_list(loop) = 0                                                  
            cycle                                                               
         endif                                                                  
c                                                                               
c                    penetration is found, so the current contact shape         
c                    is still valid.  Loop through all the previous             
c                    contact shapes to see if the current shape superceeds      
c                    them or just combines with them.                           
c                                                                               
c                       compute location if we move to current shape.           
c                                                                               
         do i = 1, 3                                                            
            tmp_coord2(i) = curr_coord(i) +                                     
     &           pen_dist(pen_list(loop)) * normal(i,pen_list(loop))            
         enddo                                                                  
         new_coords = .true.                                                    
c                                                                               
c                       loop over previous contact shapes.                      
c                                                                               
         rebuild_norm = .false.                                                 
         do j = 1, loop - 1                                                     
c                                                                               
            if ( pen_list(j).eq.0) cycle                                        
c                                                                               
c                          check if new contact location still violates         
c                          this contact shape.                                  
c                                                                               
            call find_shape_contact( tmp_coord2, pen_list(j),                   
     &           new_pen_sign, new_dist, penetrated, dumvec, node )                           
c                                                                               
c                          if no penetration is found, then shape we are        
c                          checking is no longer valid.  Otherwise, do          
c                          nothing and the shape will still be considered.      
c                          Note that if we invalidate the contact shape, then   
c                          we will have to recalculate the nodal displacement   
c                          which satisfies all contact planes.                  
c                                                                               
            if ( .not. penetrated ) then                                        
               pen_list(j) = 0                                                  
               rebuild_norm = .true.                                            
               cycle                                                            
            endif                                                               
c                                                                               
         enddo                                                                  
c                                                                               
c                    if we invalidated any previously valid contact planes,     
c                    then we need to rebuild the vector which would move the    
c                    node to a location which satisifes all contact planes.     
c                    Loop over all the currently valid contact shapes and       
c                    rebuild the vector.                                        
c                                                                               
         if ( rebuild_norm ) then                                               
            first = .true.                                                      
            do j = 1, loop -1                                                   
               if (pen_list(j) .eq. 0 ) cycle                                   
               if ( first) then                                                 
                  do i=1,3                                                      
                     curr_normal(i) = normal(i,pen_list(j))                     
                  enddo                                                         
                  first = .false.                                               
               else                                                             
                  call correct_normal (curr_normal, curr_dist,                  
     &                 normal(1,pen_list(j)),pen_dist(pen_list(j)))             
               endif                                                            
            enddo                                                               
         endif                                                                  
c                                                                               
c                    compute vector which moves the node to location            
c                    which satisifes all the previous contact shapes            
c                    as well as the new shape.  curr_normal and curr_dist       
c                    hold the vector which satisifes all of the previous        
c                    contact shapes, while normal(1,pen_list(loop))             
c                    and pen_dist(pen_list(loop)) hold the vector which         
c                    satisifes the new contact shape.                           
c                                                                               
         call correct_normal ( curr_normal, curr_dist,                          
     &        normal(1,pen_list(loop)), pen_dist(pen_list(loop)))               
c                                                                               
      enddo                                                                     
c                                                                               
c             print order of penetration shapes for debugging                   
c                                                                               
      if (debug) then                                                           
         write (*,*) '   >>> penetration shapes for node ',node                 
         do i=1, num_pen                                                        
            write (*,*) '     shape:',pen_list(i),' pen:',                      
     &           pen_dist(pen_list(i))                                          
         enddo                                                                  
      endif                                                                     
c                                                                               
c             construct contact_cause for this node, listing the                
c             contacting shapes in order of decreasing penetration.             
c                                                                               
      entry = 0                                                                 
      do i = maxcontact, 1, -1                                                  
         if ( pen_list(i) .ne. 0 ) then                                         
            entry = entry + 1                                                   
            contact_cause(entry,node) = pen_list(i)                             
         endif                                                                  
      enddo                                                                     
c                                                                               
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine correct_normal               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/16/98                   *          
c     *                                                              *          
c     *     This routine takes the current vector (curr_norm and     *          
c     *     curr_dist) which moves the node to a location which      *          
c     *     satisifes one or more contact shapes and combines it     *          
c     *     with a new vector which satisfies a new contact shape.   *          
c     *     Essentially, this routine takes two vectors normal to    *          
c     *     two contact plane and finds the vector which points to   *          
c     *     the closest location where the planes intersect.         *          
c     *                                                              *          
c     *     See manual for additional description.                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine correct_normal ( curr_norm, curr_dist, new_norm,               
     &     new_dist )                                                           
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     curr_norm(3), new_norm(3), curr_dist, new_dist,                      
     &     alpha, dot, zero, corr_dist,                                         
     &     tolval, val, plane_norm(3), corr_norm(3), dumd                       
      logical debug                                                             
      data debug, zero, tolval  /.false., 0.0, .00001/                          
c                                                                               
c            compute angle between normals.                                     
c                                                                               
      call dot_prod (curr_norm, new_norm, dot )                                 
      alpha = acos(dot)                                                         
c                                                                               
c            if angle is almost zero, then make the further shape               
c            the primary shape.                                                 
c                                                                               
      if ( abs(alpha) .lt. tolval ) then                                        
         curr_dist = new_dist                                                   
         return                                                                 
      endif                                                                     
c                                                                               
c            now compute correction distance                                    
c                                                                               
      corr_dist = ( curr_dist - new_dist * cos(alpha) ) / sin(alpha)            
c                                                                               
c            compute direction normal to plane of new_norm                      
c            and curr_norm.                                                     
c                                                                               
      call cross_prod (curr_norm, new_norm, plane_norm)                         
      call normalize (plane_norm, dumd)                                         
c                                                                               
c            compute correction normal by finding a vector perpendicular to     
c            new_norm in direction of curr_norm.                                
c                                                                               
      call cross_prod (new_norm, plane_norm, corr_norm)                         
c                                                                               
c            we have correction distance and direction, now compute             
c            new current normal and the distance.                               
c                                                                               
      do i = 1, 3                                                               
         curr_norm(i) = new_dist*new_norm(i) + corr_dist*corr_norm(i)           
      enddo                                                                     
c                                                                               
      call normalize (curr_norm, curr_dist)                                     
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_shape_contact           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/3/2017 rhd              *          
c     *                                                              *          
c     *        This routine serves as a branching point for finding  *          
c     *        contact on all contact shapes.                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_shape_contact( curr_coord, shape, pen_dist,               
     &     pen_sign, penetrated, normal, node )                                           
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision :: pen_dist, pen_sign, curr_coord(3), normal(3)                                   
      logical :: debug, penetrated                                                 
      data debug /.false./                                                      
c                                                                               
c      if ( debug)  write (*,*) '                 >> looking for',              
c     &     ' contact: node, shape:',node, shape                                
c                                                                               
      if (contact_shape(shape) .eq. 1) then                                     
         call find_plane_contact( curr_coord, shape, pen_dist,                  
     &        penetrated, normal, node, debug )  
         pen_sign = 1.0d0                               
      else if (contact_shape(shape) .eq. 2) then                                
         call find_cyl_contact( curr_coord, shape, pen_dist,                    
     &        pen_sign, penetrated, normal, node, debug )                                 
      else if (contact_shape(shape) .eq. 3) then                                
         call find_sph_contact( curr_coord, shape, pen_dist,                    
     &        pen_sign, penetrated, normal, node, debug )                                 
      else                                                                      
         penetrated = .false.                                                   
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_plane_contact           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/17/98                   *          
c     *                                                              *          
c     *         This routine determines if contact has occurred      *          
c     *         on a subplane region.                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_plane_contact( curr_coord, shape, pen_dist,               
     &     penetrated, normal, node, debug )                                    
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     pen_dist, zero, curr_coord(3), vprime(3), dot1, dot2,                
     &     normal(3), mag1, mag2                                                
      logical debug, penetrated                                                 
      data zero /0.0/                                                           
c                                                                               
c                check if node is penetrating contact plane. Do this by the     
c                following method:                                              
c                                                                               
c                 - compute:                                                    
c                      v' = p' - p1' (v' is the vector between the current      
c                                     point and the contact plane corner.       
c                                     p' is coordinates of current point,       
c                                     p1' is current location of node 1         
c                                     of contact plane).                        
c                      -dot(v',v3)   (this is the penetration, where            
c                                     v3 is the normal to the plane )           
c                 - if penetration is positive, then we need to check the       
c                    the depth.  If the penetration is deeper than the depth,   
c                    we have no contact.                                        
c                 - if penetration is positive and inside depth, then we        
c                    need to check the sides. compute dot(v1,v') and            
c                    dot(v2,v').  If both are > 0 and < 1. then we have         
c                    contact.                                                   
c                                                                               
      penetrated = .false.                                                      
c                                                                               
c                      find v', dot(v',v3)                                      
c                                                                               
      pen_dist = zero                                                           
      do i = 1, 3                                                               
         vprime(i) = curr_coord(i) - ( cshape_pnt(i,shape) +                    
     &        cshape_rate(i,shape)*dt )                                         
         pen_dist = pen_dist - cshape_norm(i,shape) * vprime(i)                 
      enddo                                                                     
c                                                                               
c                      if we haven't penetrated, exit                           
c                                                                               
      if (pen_dist .le. zero) goto 9999                                         
c                                                                               
c                      check depth                                              
c                                                                               
      if (pen_dist .gt. contact_depth(shape)) goto 9999                         
c                                                                               
c                      we have penetrated plane. find dot1 and dot2 to          
c                      see if we are in subplane.                               
c                                                                               
      dot1 = zero                                                               
      dot2 = zero                                                               
      mag1 = zero                                                               
      mag2 = zero                                                               
      do i = 1, 3                                                               
         dot1 = dot1 + cplane_vec(i,1,shape) * vprime(i)                        
         dot2 = dot2 + cplane_vec(i,2,shape) * vprime(i)                        
         mag1 = mag1 + cplane_vec(i,1,shape) ** 2                               
         mag2 = mag2 + cplane_vec(i,2,shape) ** 2                               
      enddo                                                                     
      if ( dot1 .ge. zero .and. dot2 .ge. zero .and.                            
     &     dot1 .le. mag1  .and. dot2 .le. mag2 ) then                          
         penetrated = .true.                                                    
         do i=1, 3                                                              
            normal(i) = cshape_norm(i,shape)                                    
         enddo                                                                  
      endif                                                                     
c                                                                               
 9999 continue                                                                  
c                                                                               
      if ( penetrated .and. debug) then                                         
         write (*,*) '           -> penetrating node: ',node                    
         write (*,*) '                dist,dot ratios:',                        
     &        pen_dist,dot1/mag1,dot2/mag2                                      
         write (*,*) '                coords:',                                 
     &        (curr_coord(j),j=1,3)                                             
         write (*,*) '                point:',                                  
     &        (cshape_pnt(j,shape),j=1,3)                                       
         write (*,*) '                normal:',                                 
     &        (normal(j),j=1,3)                                                 
      endif                                                                     
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_cyl_contact             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/2/2017 rhd              *          
c     *                                                              *          
c     *         This routine determines if contact has occurred      *          
c     *         on a cylindrical region.                             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_cyl_contact( curr_coord, shape, pen_dist,               
     &      pen_sign, penetrated, normal, node, debug )                                    
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision :: pen_sign,                                                          
     &     pen_dist, zero, curr_coord(3), vprime(3), dot, proj_dist,            
     &     normal(3), mag, angle, point(3), dumd, dist, radius,                 
     &     length                                                               
      logical debug, penetrated, outside, inside                                                
      data zero /0.0/                                                           
c                                                                               
c      if( debug)  write (*,*) '                 >> looking for',              
c     &     ' contact: node,pln:',node, shape                                   
c                                                                               
c                check if node is penetrating cylinder. Do this by the          
c                following method:                                              
c                                                                               
c                  - find the vector perpendicular to center line which         
c                         runs thru given nodal location.                       
c                  - calc distance between node and center line along normal.   
c                         this gives penetration.                               
c                  - check if projection of node onto center line is within     
c                         the length of the center line.                        
c        
      outside = contact_outside(shape)
      inside = .not. outside                                                                       
      radius = cshape_param(1,shape)                                            
      length = cshape_param(2,shape)                                            
      penetrated = .false.                                                      
c                                                                               
c                calc vector between node and start point of cylinder           
c                                                                               
      do i = 1, 3                                                               
         vprime(i) = curr_coord(i) - ( cshape_pnt(i,shape) +                    
     &        cshape_rate(i,shape)*dt )                                         
      end do                                                                     
c                                                                               
c                now calc angle, and then the perpendiular distance             
c                                                                               
      call normalize (vprime,mag)                                               
      call dot_prod (vprime, cshape_norm(1,shape), dot)                         
      angle = acos (dot) 
      proj_dist = mag * sin(angle)                                              
      pen_dist = radius - proj_dist     
c
      if( outside ) then
        pen_sign = 1.0d0
        if( pen_dist <= zero) go to 9999   ! node still outside cylinder  
      end if
      if( inside ) then
        pen_sign = -1.0d0
        pen_dist = proj_dist - radius
        if( pen_dist <= zero ) go to 9999 ! node still inside cycliner                                       
      end if                                                                                        
c                      if we haven't penetrated, exit                           
c                                                                               
c                      we have penetrated cylinder. Find if we are              
c                      within the length of the cylinder.  First find           
c                      normal.                                                  
c                                                                               
      call find_shape_normal ( node, curr_coord, shape, normal )                
      dist = mag * cos(angle)                                                   
      if ( dist .gt. zero .and. dist .lt. length) penetrated = .true.           
c                                                                               
 9999 continue                                                                  
c                                                                               
      if ( penetrated .and. debug ) then                                        
         write (*,*) myid,':           -> penetration of cyl by node: ',        
     &        node                                                              
         write (*,*)  myid,':                dist ,length:',                    
     &        pen_dist, dist                                                    
         write (*,*)  myid,':                coords:',                          
     &        (curr_coord(j),j=1,3)                                             
         write (*,*)  myid,':                proj_point:',                      
     &        (point(j),j=1,3)                                                  
         write (*,*)  myid,':                normal:',                          
     &        (normal(j),j=1,3)                                                 
         write (*,*) myid,':          vprime, mag:',                            
     &        (vprime(i),i=1,3),mag                                             
         write (*,*) myid,':          dot:',dot,' angle:',angle,                
     &        ' proj_dist:', proj_dist                                          
         write (*,*) myid,':          radius:',radius,' tot length:',           
     &        length                                                            
         write (*,*) myid,':          dt:',dt                                   
      endif                                                                     
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_sph_contact             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 7/3/2017 rhd               *          
c     *                                                              *          
c     *         This routine determines if contact has occurred      *          
c     *         on a spherical region.                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_sph_contact ( curr_coord, shape, pen_dist,                
     &      pen_sign, penetrated, normal, node, debug )                                    
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision :: pen_sign,                                                         
     &     pen_dist, zero, curr_coord(3), vprime(3),                            
     &     normal(3), mag, angle, point(3), dumd, dist, radius,                 
     &     length                                                               
      logical debug, penetrated, outside, inside                                                 
      data zero /0.0/                                                           
c                                                                               
c      if ( debug)  write (*,*) '                 >> looking for',              
c     &     ' contact: node,pln:',node, shape                                   
c                                                                               
c                check if node is penetrating sphere. Do this by the            
c                following method:                                              
c                                                                               
c                  - find the vector between node and center point of           
c                         sphere                                                
c                  - if the length of the vector is less than the radius        
c                         we have penetration.                                  
c                  - normal is vector after normalization                       
c                                                                               
      outside = contact_outside(shape)
      inside = .not. outside                                                                       
      radius = cshape_param(1,shape)                                            
      penetrated = .false.     
      pen_sign = 1.0d0                                                 
c                                                                               
      do i = 1, 3                                                               
         vprime(i) = curr_coord(i) - ( cshape_pnt(i,shape) +                    
     &        cshape_rate(i,shape)*dt )                                         
      end do                                                                     
      call normalize (vprime,dist)                                              
c                                                                               
c                                                                               
      if( outside ) then
        if( dist .lt. radius ) then                                               
          penetrated = .true.                                                     
          pen_dist = radius - dist                                                
          normal(1:3) = vprime(1:3)  
          pen_sign = 1.0d0                                             
        end if 
      else ! chk penetration from inside sphere
        if( dist .gt. radius ) then
          penetrated = .true.                                                     
          pen_dist = abs( radius - dist )                                           
          normal(1:3) = vprime(1:3)  
          pen_sign = -1.0d0                                             
        end if 
      end if  
c                                                                               
c                                                                               
 9999 continue                                                                  
      if ( penetrated .and. debug) then                                         
         write (*,*) '           -> penetration of sph by node: ',node          
         write (*,*) '                dist:', pen_dist                          
         write (*,*) '                coords:',                                 
     &        (curr_coord(j),j=1,3)                                             
         write (*,*) '                normal:',                                 
     &        (normal(j),j=1,3)                                                 
         write (*,*) '          radius:',radius                                 
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_shape_normal            *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/18/98                   *          
c     *                                                              *          
c     *        This routine, given the coordinates of a point and    *          
c     *        a contact shape, returns the normal to the contact    *          
c     *        shape at that point.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_shape_normal ( node, curr_coord, shape, normal )          
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     normal(3), curr_coord(3), zero                                       
      data zero /0.0/                                                           
c                                                                               
c             branch on contact shape                                           
c                                                                               
      if ( contact_shape(shape) .eq. 1 ) then                                   
         call find_plane_normal ( shape, normal )                               
      else if ( contact_shape(shape) .eq. 2 ) then                              
         call find_cyl_normal ( curr_coord, shape, normal )                     
      else if ( contact_shape(shape) .eq. 3 ) then                              
         call find_sph_normal ( curr_coord, shape, normal )                     
      else                                                                      
         do i = 1, 3                                                            
            normal(i) = zero                                                    
         enddo                                                                  
      endif                                                                     
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_plane_normal            *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/18/98                   *          
c     *                                                              *          
c     *        This routine, given a contact plane, returns the      *          
c     *        normal to the plane.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_plane_normal( shape, normal)                              
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     normal(3)                                                            
c                                                                               
      do j = 1, 3                                                               
         normal(j) = cshape_norm(j,shape)                                       
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_cyl_normal              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/18/98                   *          
c     *                                                              *          
c     *        This routine, given the coordinates of a point and    *          
c     *        a contact cylinder, returns the normal to the         *          
c     *        cylinder at that point.                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_cyl_normal ( curr_coord, shape, normal )                  
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     curr_coord(3), vprime(3), normal(3), vec(3), dumd                    
      logical debug                                                             
c                                                                               
c                calc vector between node and start point of cylinder           
c                                                                               
      do i = 1, 3                                                               
         vprime(i) = curr_coord(i) - ( cshape_pnt(i,shape) +                    
     &        cshape_rate(i,shape)*dt )                                         
      enddo                                                                     
c                                                                               
c                now calc normal                                                
c                                                                               
      call cross_prod (cshape_norm(1,shape), vprime, vec)                       
      call cross_prod (vec, cshape_norm(1,shape), normal)                       
      call normalize (normal,dumd)                                              
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_sph_normal              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/18/98                   *          
c     *                                                              *          
c     *        This routine, given the coordinates of a point and    *          
c     *        a contact sphere, returns the normal to the contact   *          
c     *        sphere at that point.                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine find_sph_normal ( curr_coord, shape, normal )                  
      use global_data ! old common.main
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     curr_coord(3), normal(3), dumd                                       
c                                                                               
c                calc vector between node and start point of sphere             
c                                                                               
      do i = 1, 3                                                               
         normal(i) = curr_coord(i) - ( cshape_pnt(i,shape) +                    
     &        cshape_rate(i,shape)*dt )                                         
      enddo                                                                     
c                                                                               
c                now calc normal                                                
c                                                                               
      call normalize (normal, dumd)                                             
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine contact_trnmat               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/28/98                   *          
c     *                                                              *          
c     *        This routine computes the transformation matrix for   *          
c     *        a contacting node.  If a transformation matrix cannot *          
c     *        be computed, then the subroutine returns with         *          
c     *        success = .false.  Otherwise it returns with the      *          
c     *        transformation matrix.                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine contact_trnmat (node, curr_norm, transmat, success)            
      use global_data ! old common.main
c                                                                               
      use main_data, only : cnstrn, cnstrn_in                                   
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     zero, one, transmat(3,3),                                            
     &     d32460, dumd, workvec(3), vec1(3),                                   
     &     vec2(3), dot, tolval, curr_norm(3)                                   
      logical debug, done, undo, success, calc_trnmat, orthog,                  
     &     anticonst                                                            
      data debug, zero, one, d32460, tolval                                     
     &     /.false., 0.0, 1.0, 32460.0, 0.0001 /                                
c                                                                               
c           We want to compute a transformation matrix for the contacting       
c           node.                                                               
c                                                                               
c           We have two cases;                                                  
c           1) if the contact normal is perpendicular to the                    
c              coordinate axes, then we do not need to rotate                   
c              the constraints, nor do we need a transformation                 
c              matrix.                                                          
c           2) if the contact normal is not perpendicular to                    
c              the coordinate axes, then:                                       
c              - if we have 0 constraints, we are free to generate              
c                the transformation matrix                                      
c              - if we have 1 constraint, we can                                
c                make transformation matrix if the constraint is                
c                perpendicular to the contact normal.  If not,                  
c                we can't make a transformation matrix.                         
c                                                                               
c           First figure if contact normal is perpendicular to the              
c           coordinate axes and count the number of constraints.                
c                                                                               
         orthog = .false.                                                       
         num_const = 0                                                          
         do i=1, 3                                                              
            dof = dstmap(node) + i - 1                                          
            if ( cstmap(dof) .ne. 0 ) num_const = num_const + 1                 
            if (abs(abs(curr_norm(i))-one) .lt. tolval) orthog = .true.         
         enddo                                                                  
c                                                                               
c                if orthog, then contact plane is orthog to                     
c                the coordinate axes.  Don't form a                             
c                transformation matrix because the contact                      
c                "springs" will already be orthogonal                           
c                to one of the axes.                                            
c                                                                               
         if ( orthog ) then                                                     
c                                                                               
            success = .false.                                                   
c                                                                               
         else                                                                   
c                                                                               
c                if not orthog, then we have the following                      
c                cases:                                                         
c                - no constraints: form an appropriate                          
c                    transformation matrix                                      
c                - one constraint: If the constraint is                         
c                    perpendicular to the normal of the contact                 
c                    shape, then we can form a transformation                   
c                    matrix in which one of the transformed                     
c                    coordinate directions is parallel to the                   
c                    constraint.  If we can't do it, then                       
c                    don't form a transformation matrix.                        
c                - two or three constraints: cannot form                        
c                    a transformation matrix.                                   
c                                                                               
c                                                                               
            success = .true.                                                    
            if (debug) write (*,*) '     >>> NON-ORTHOGONAL COORD SYS'          
c                                                                               
c                     -------                                                   
c                     Case 1: No constraints                                    
c                     -------                                                   
c                                                                               
            if ( num_const .eq. 0 ) then                                        
c                                                                               
               if ( debug ) write (*,*) '      0 constraints'                   
c                                                                               
c                             form a vector not in plane with curr_norm         
c                                                                               
               if ( abs(abs(curr_norm(3)-one)) .lt. tolval ) then               
                  workvec(1) = zero                                             
                  workvec(2) = one                                              
                  workvec(3) = zero                                             
               else                                                             
                  workvec(1) = zero                                             
                  workvec(2) = zero                                             
                  workvec(3) = one                                              
               endif                                                            
c                                                                               
c                              cross curr_norm with other vec to get first      
c                              transmat vector                                  
c                                                                               
               call cross_prod (workvec, curr_norm, vec1)                       
               call normalize (vec1, dumd)                                      
c                                                                               
c                               now compute final vec                           
c                                                                               
               call cross_prod (curr_norm, vec1, vec2)                          
c                                                                               
c                               now put vecs into transmat.  note that          
c                               the normal to the contact shape is made         
c                               to be the third direction in transmat           
c                               (dof 3)                                         
c                                                                               
               transmat(1,1:3) = vec1(1:3)                                      
               transmat(2,1:3) = vec2(1:3)                                      
               transmat(3,1:3) = curr_norm(1:3)                                 
c                                                                               
            else if ( num_const .eq. 1 ) then                                   
c                                                                               
c                     -------                                                   
c                     Case 2: One constraint                                    
c                     -------                                                   
c                                                                               
               if (debug) write (*,*) '      1 constraint'                      
c                                                                               
c                               In this case, build a vector for the            
c                               constraint.                                     
c                                                                               
               entry = 0                                                        
               do i=1, 3                                                        
                  vec1(i) = zero                                                
                  dof = dstmap(node) + i - 1                                    
                  if ( cstmap(dof) .ne. 0 ) then                                
                     vec1(i) = one                                              
                     entry = i                                                  
                  endif                                                         
               enddo                                                            
c                                                                               
c                               find dot product of constraint vector           
c                               and contact normal.  If not zero,               
c                               then they are not perpendicular and             
c                               we can't make a transmat.                       
c                                                                               
               call dot_prod ( curr_norm, vec1, dot )                           
c                                                                               
               if ( abs(dot) .gt. tolval) then                                  
c                                                                               
                  if ( debug) write (*,*) '    >>> cant move it.'               
                  success = .false.                                             
                  return                                                        
c                                                                               
               endif                                                            
c                                                                               
c                               compute final vec                               
c                                                                               
               call cross_prod (vec1, curr_norm, vec2)                          
c                                                                               
c                               now put vecs into transmat.  Make sure that     
c                               the constraint vector corresponds to the        
c                               same degree of freedom so that we wont have     
c                               to modify the constraints table.                
c                                                                               
               ent2 = entry + 1                                                 
               if (ent2 .gt. 3) ent2 = 1                                        
               ent3 = entry - 1                                                 
               if (ent3 .lt. 1) ent3 = 3                                        
               transmat(entry,1:3) = vec1(1:3)                                  
               transmat(ent3,1:3) = vec2(1:3)                                   
               transmat(ent2,1:3) = curr_norm(1:3)                              
c                                                                               
            else                                                                
c                                                                               
c                     -------                                                   
c                     Case 2: More than one constraint                          
c                     -------                                                   
c                                                                               
c                               we can't compute a transmat for this case.      
c                                                                               
               if (debug) write (*,*) '      > 1 constraint'                    
               success = .false.                                                
               return                                                           
c                                                                               
            endif                                                               
c                                                                               
         endif                                                                  
c                                                                               
c                                                                               
         if ( debug) then                                                       
            if (success) then                                                   
               write (*,*) '     >> We have success! Here is transmat:'         
               write (*,'(10x,3e14.6)') ((transmat(k,j),j=1,3),k=1,3)           
               write (*,*) '     << end transmat.'                              
            else                                                                
               write (*,*) '     >> Failure! No transmat for you,',             
     &              ' my pretty...'                                             
            endif                                                               
         endif                                                                  
c                                                                               
         return                                                                 
         end                                                                    
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine contact_stfadd               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/28/98                   *          
c     *                                                              *          
c     *        This routine computes the appropriate stiffness       *          
c     *        contributions due to contact to the element           *          
c     *        stiffnesses for a block of elements.                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine contact_stfadd (span, felem, totdof, edest,                    
     &     ek, nrowek, nnode, belinc )                                          
      use global_data ! old common.main
c                                                                               
      use main_data, only : invdst, trn, trnmat, crdmap,                        
     &                      inverse_incidences,                                 
     &                      asymmetric_assembly                                 
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     zero, stiff, ek (nrowek, span) , normal(3), transmat(3,3),           
     &     mat(3,3), mat2(3,3), nodstf(3,3), elstf(24,24),                      
     &     curr_coord(3), factor                                                
      dimension edest(totdof,*), belinc(nnode,*)                                
      logical debug                                                             
      data debug, zero, factor /.false., 0.0, 1.001/                            
c                                                                               
c     For now, trigger an error if we try to do contact with asymmetric         
c     assembly                                                                  
c                                                                               
      if (use_contact .and. asymmetric_assembly) then                           
        write(*,*) "Contact not implemented for asym. assembly"                 
        call die_gracefully                                                     
      end if                                                                    
                                                                                
c                                                                               
c         skip if we aren't doing contact                                       
c                                                                               
      if ( .not. use_contact ) return                                           
c                                                                               
c         loop over elements in block.  For any dof of element that has         
c         contact caused by a contact shape (in contact_cause), we need         
c         to add the stiffness of the contact shape into the stiffness          
c         for the element.  Note that if the dof is shared by other             
c         elements, we must split the contact stiffness into equal              
c         parts to assign to each element stiffness matrix.  This               
c         means that when the element stiffnesses are added to get the          
c         structural stiffness matrix, then the dof will have the correct       
c         contact stiffness.  This approach also provides the correct           
c         stiffness for the LNPCG solver.                                       
c                                                                               
      do elem = 1, span                                                         
c                                                                               
         gelem = elem + felem - 1                                               
c                                                                               
         do node_loop = 1, nnode                                                
c                                                                               
            node = belinc(node_loop,elem)                                       
c                                                                               
c                skip if no contact                                             
c                                                                               
            if ( contact_cause(1,node) .eq. 0 ) cycle                           
c                                                                               
            if ( debug ) then                                                   
               write (*,*) '--> do node', node,'w/cause:',                      
     &              (contact_cause(i,node), i=1,num_contact)                    
               if ( trn(node) ) then                                            
                  write (*,*) '    > trnmat:'                                   
                  write (*,'(10x,3e13.6)')                                      
     &                 ((trnmat(node)%mat(k,j),j=1,3),k=1,3)                    
               else                                                             
                  write (*,*) '     > NO trnmat.'                               
               endif                                                            
            endif                                                               
c                                                                               
c                   Form the 3x3 nodal stiffness matrix to simplify             
c                   stiffness modifications.  Also get current                  
c                   coordinates of node.                                        
c                                                                               
            do j = 1, 3                                                         
               do i = 1, j                                                      
                  nodstf(j,i) = ek(dcp((j-1)*nnode + node_loop)                 
     &                 - nnode*(j-i),elem)                                      
               enddo                                                            
               do i = 1, j-1                                                    
                  nodstf(i,j) = nodstf(j,i)                                     
               enddo                                                            
            enddo                                                               
c                                                                               
            if ( debug) then                                                    
               write (*,*) '   > nodal stiffness for node ',node,               
     &              ' elem:', gelem                                             
               write (*,'(10x,3e13.6)') ((nodstf(i,j),j=1,3),i=1,3)             
            endif                                                               
c                                                                               
            call get_coords ( node, curr_coord )                                
c                                                                               
c                   loop over valid contact shapes                              
c                                                                               
            do cause = 1, num_contact                                           
c                                                                               
               if ( contact_cause(cause,node) .eq. 0 ) exit                     
c                                                                               
c                      Calculate the stiffness due to penetration of the        
c                      current contact shape.  Multiply by a factor larger      
c                      larger than one to make the contact shape appear         
c                      stiffer in the solution than it actually is.  This       
c                      prevents oscilations in which in alternating             
c                      iterations the contact force pushes the node             
c                      completely out of the contact shape, then the next       
c                      iteration pushes the node back into the contact shape.   
c                                                                               
               stiff = contact_stiff(contact_cause(cause,node)) /               
     &              dble(inverse_incidences(node)%element_count) *              
     &              factor                                                      
               if ( debug) write (*,*) '    > stiff for shape ',                
     &              contact_cause(cause,node),' is ',stiff                      
c                                                                               
c                      Compute normal vector for contact shape for node         
c                                                                               
               call find_shape_normal ( node, curr_coord,                       
     &              contact_cause(cause,node), normal )                         
c                                                                               
c                      Form stiff contribution in global coords                 
c                                                                               
               do i = 1, 3                                                      
                  do j = 1, 3                                                   
                     mat(i,j) = normal(i) * normal(j) * stiff                   
                  enddo                                                         
               enddo                                                            
c                                                                               
c                      if coordinate transformation is allowed, then            
c                      rotate the contribution into transformed coordinates.    
c                      if not, then we just leave the contribution in           
c                      global coordinates.                                      
c                                                                               
c                      to rotate to transformed coordinates, we need            
c                      to calculate:                                            
c                                                                        T      
c                             [nodstf] = [nodstf] + [trn] * [mat] * [trn]       
c                                                                               
               if ( trn(node) ) then                                            
c                                                                               
c                         first do:                                             
c                                                                               
c                             [mat2] = [trn] * [mat]                            
c                                                                               
                  do i = 1, 3                                                   
                     do j = 1, 3                                                
                        mat2(i,j) = zero                                        
                        do k = 1, 3                                             
                           mat2(i,j) = mat2(i,j) +                              
     &                          trnmat(node)%mat(i,k) * mat(k,j)                
                        enddo                                                   
                     enddo                                                      
                  enddo                                                         
c                                                                               
c                         now do:                                               
c                                                                T              
c                            [nodstf] = [nodstf] + [mat2] * [trn]               
c                                                                               
                  do i = 1, 3                                                   
                     do j = 1, 3                                                
                        mat(i,j) = zero                                         
                        do k = 1, 3                                             
                           mat(i,j) = mat(i,j) +                                
     &                          mat2(i,k) * trnmat(node)%mat(j,k)               
c     &                             mat2(i,k) * trnmat(node)%mat(k,j)           
                        enddo                                                   
                     enddo                                                      
                  enddo                                                         
c                                                                               
               endif                                                            
c                                                                               
c                      Add stiffness contribution to nodal stiffness matrix.    
c                                                                               
               do i = 1, 3                                                      
                  do j = 1, 3                                                   
                     nodstf(i,j) = nodstf(i,j) + mat(i,j)                       
                  enddo                                                         
               enddo                                                            
c                                                                               
            enddo                                                               
c                                                                               
c                                                                               
c                Now put nodal stiffness back into the element stiffness        
c                matrix.                                                        
c                                                                               
            do j = 1, 3                                                         
               do i = 1, j                                                      
                  ek(dcp((j-1)*nnode + node_loop) -                             
     &                 nnode*(j-i),elem) = nodstf(j,i)                          
               enddo                                                            
            enddo                                                               
c                                                                               
c                                                                               
            if ( debug) then                                                    
               write (*,*) '    > new nodal stiffness for node ',node,          
     &              ' elem:', gelem                                             
               write (*,'(10x,3e13.6)') ((nodstf(i,j),j=1,3),i=1,3)             
            endif                                                               
c                                                                               
         enddo                                                                  
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine trn2all                      *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/02/98                   *          
c     *                                                              *          
c     *        This routine applies a nodal transformation matrix    *          
c     *        to all the data structures that may be effected by    *          
c     *        the change of coordinates.  This allows the           *          
c     *        transformation matrices to be changed in the          *          
c     *        middle of a load step, which is needed by the         *          
c     *        contact algorithms.                                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine trn2all ( node, status )                                       
      use global_data ! old common.main
c                                                                               
      use main_data, only : trnmat, pbar                                        
      use contact                                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     trans(3,3)                                                           
c                                                                               
c            get the nodal transformation matrix.  If status is 1,              
c            then use the normal matix; if 2, then use the transpose.           
c                                                                               
c            note that we assume that each node has only 3 dofs.                
c                                                                               
      if ( status .eq. 1 ) then                                                 
         trans(1,1) = trnmat(node)%mat(1,1)                                     
         trans(1,2) = trnmat(node)%mat(1,2)                                     
         trans(1,3) = trnmat(node)%mat(1,3)                                     
         trans(2,1) = trnmat(node)%mat(2,1)                                     
         trans(2,2) = trnmat(node)%mat(2,2)                                     
         trans(2,3) = trnmat(node)%mat(2,3)                                     
         trans(3,1) = trnmat(node)%mat(3,1)                                     
         trans(3,2) = trnmat(node)%mat(3,2)                                     
         trans(3,3) = trnmat(node)%mat(3,3)                                     
      else                                                                      
         trans(1,1) = trnmat(node)%mat(1,1)                                     
         trans(1,2) = trnmat(node)%mat(2,1)                                     
         trans(1,3) = trnmat(node)%mat(3,1)                                     
         trans(2,1) = trnmat(node)%mat(1,2)                                     
         trans(2,2) = trnmat(node)%mat(2,2)                                     
         trans(2,3) = trnmat(node)%mat(3,2)                                     
         trans(3,1) = trnmat(node)%mat(1,3)                                     
         trans(3,2) = trnmat(node)%mat(2,3)                                     
         trans(3,3) = trnmat(node)%mat(3,3)                                     
      endif                                                                     
c                                                                               
c           now transform all of the dof vectors that need transformation       
c                                                                               
c                                                                               
      call trnnvec (trans, u, node, 3)                                          
      call trnnvec (trans, du, node, 3)                                         
      call trnnvec (trans, idu, node, 3)                                        
      call trnnvec (trans, ifv, node, 3)                                        
c                                                                               
      if (slave_processor) return                                               
c                                                                               
      call trnnvec (trans, pbar, node, 3)                                       
      call trnnvec (trans, v, node, 3)                                          
      call trnnvec (trans, a, node, 3)                                          
      call trnnvec (trans, load, node, 3)                                       
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine trnnvec                      *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/02/98                   *          
c     *                                                              *          
c     *        This routine applies a nodal transformation matrix    *          
c     *        to a vector which ranges over the dofs in the model.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine trnnvec (trans, vec, node, nod_dof)                            
      use global_data ! old common.main
c                                                                               
      use main_data, only : trn, trnmat                                         
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     trans(mxndof,mxndof), vec(*), newvec(mxndof), zero                   
      data zero /0.0/                                                           
c                                                                               
      do i = 1, nod_dof                                                         
         newvec(i) = zero                                                       
         do j = 1, nod_dof                                                      
            newvec(i)= newvec(i) + vec(dstmap(node)+j-1) * trans(i,j)           
         enddo                                                                  
      enddo                                                                     
c                                                                               
      do i = 1, nod_dof                                                         
         vec(dstmap(node)+i-1) = newvec(i)                                      
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dot_prod                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/02/98                   *          
c     *                                                              *          
c     *        This routine comp[utes the dot product of the two     *          
c     *        input vectors.                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dot_prod (vec1, vec2, dot)                                     
c                                                                               
c                                                                               
      double precision                                                          
     &     vec1(3), vec2(3), dot, zero                                          
      integer i, j                                                              
      data zero /0.0/                                                           
c                                                                               
      dot = zero                                                                
      do i = 1, 3                                                               
         dot = dot + vec1(i) * vec2(i)                                          
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine cross_prod                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/02/98                   *          
c     *                                                              *          
c     *        This routine computes the cross product of the two    *          
c     *        input vectors.                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine cross_prod (vec1, vec2, vec_out)                               
c                                                                               
c                                                                               
      double precision                                                          
     &     vec1(3), vec2(3), vec_out(3)                                         
      data zero /0.0/                                                           
c                                                                               
      vec_out(1) = vec1(2)*vec2(3) -                                            
     &     vec2(2)*vec1(3)                                                      
      vec_out(2) = vec2(1)*vec1(3) -                                            
     &     vec1(1)*vec2(3)                                                      
      vec_out(3) = vec1(1)*vec2(2) -                                            
     &     vec2(1)*vec1(2)                                                      
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine normalize                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/02/98                   *          
c     *                                                              *          
c     *        This routine normalizes an input vector.              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine normalize (vec, mag)                                           
c                                                                               
c                                                                               
      double precision                                                          
     &     vec(3), zero, mag                                                    
      integer i, j                                                              
      data zero /0.0/                                                           
c                                                                               
      mag = sqrt ( vec(1)**2 + vec(2)**2 + vec(3)**2 )                          
c                                                                               
      do i=1, 3                                                                 
         vec(i) = vec(i) / mag                                                  
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_coords                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/22/98                   *          
c     *                                                              *          
c     *        This routine computes the current location of the     *          
c     *        specified node in global coordinates.                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_coords ( node, curr_coord )                                
      use global_data ! old common.main
c                                                                               
      use main_data, only : crdmap, trn, trnmat                                 
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      double precision                                                          
     &     curr_coord(3), u_new(3), du_new(3), zero                             
      data zero /0.0/                                                           
c                                                                               
c               if the node is currently using nodal transformation             
c               coordinates, then rotate the current displacements and          
c               increase in displacements back to global coordinates.           
c                                                                               
      if ( trn(node) ) then                                                     
c                                                                               
         do i=1, 3                                                              
            u_new(i) = zero                                                     
            du_new(i) = zero                                                    
            do j=1,3                                                            
               u_new(i) = u_new(i) + trnmat(node)%mat(j,i) *                    
     &              u(dstmap(node)+j-1)                                         
               du_new(i) = du_new(i) + trnmat(node)%mat(j,i) *                  
     &              du(dstmap(node)+j-1)                                        
            enddo                                                               
         enddo                                                                  
c                                                                               
      else                                                                      
c                                                                               
         do i=1,3                                                               
            u_new(i) = u(dstmap(node)+i-1)                                      
            du_new(i) = du(dstmap(node)+i-1)                                    
         enddo                                                                  
c                                                                               
      endif                                                                     
c                                                                               
c               calculate current coordinates                                   
c                                                                               
      do i = 1, 3                                                               
         curr_coord (i) = c(crdmap(node)+i-1) +                                 
     &        u_new(i) + du_new(i)                                              
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine updt_contact                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/22/98                   *          
c     *                                                              *          
c     *        This routine updates the location of the center       *          
c     *        point for all contact surfaces after the completion   *          
c     *        of a load step.                                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine updt_contact                                                   
      use global_data ! old common.main
c                                                                               
      use contact                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
c                update goddamn center point.                                   
c                                                                               
      do shape = 1, num_contact                                                 
c                                                                               
         do i=1, 3                                                              
            cshape_pnt(i,shape) = cshape_pnt(i,shape) +                         
     &           cshape_rate(i,shape) * dt                                      
         enddo                                                                  
c                                                                               
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
