c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine elem_load                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 12/01/11 jcs               *          
c     *                                                              *          
c     *     this subroutine sets up the equivalent nodal loads       *          
c     *     calcs for all element loads on a 3-d isoparametrics      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c           Here is how the face/body/pressure/temperature loading              
c           storage strategy works:                                             
c                                                                               
c                 There is one 3 by 'eload_size' integer array and one          
c                 'eload_size' long real array that store the loading           
c                 information:                                                  
c                                                                               
c                   eload_data(i,1) -- the element with the loading             
c                                       description for entry i                 
c                   eload_data(i,2) -- the type of the loading for entry i:     
c                               1 to  6 : face number for traction loading      
c                                  0    : body force loading                    
c                              -1 to -6 : pressure loading                      
c                              -7 to -12: -face # for piston loading -6         
c                                -100   : element temperature                   
c                   eload_data(i,3)  -- the dof for the loading for entry i     
c                               (zero for pressure loading and                  
c                                temperature)                                   
c                   eload_val(i)     -- the value of the force for entry i      
c                                       value of element temperature            
c                   eload_pist(i)   -- piston table number for entry i          
c                   thread_number(i) -- thread number for entry i               
c                                                                               
c                                                                               
      subroutine elem_load( eload_data, eload_val, size, eload_pist,            
     &                      thread_number, mult_fact, eqloads )                 
      use global_data ! old common.main
c                                                                               
      use main_data                                                             
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                global variables                                               
c                                                                               
      dimension eload_data(size,3), thread_number(size),                        
     &          eload_pist(size)                                                
      real eload_val(size)                                                      
      double precision                                                          
     &     mult_fact, eqloads(*)                                                
c                                                                               
c                local variables                                                
c                                                                               
      logical debug                                                             
      double precision                                                          
     &      zero, nrmeq                                                         
c                                                                               
      double precision,                                                         
     & dimension(:,:), allocatable :: del_eqload                                
c                                                                               
      data debug, zero / .false., 0.0 /                                         
c                                                                               
      if ( debug ) then                                                         
         write (*,*) '>>>> entering elem_load.'                                 
         call dump_load( eload_data, eload_val, thread_number, size )           
      end if                                                                    
c                                                                               
c            allocate thread local arrays for eqload calculations.              
c            we're handling the reduction process aross                         
c            threads by allocating a (private) work vector                      
c            for each thread.                                                   
c                                                                               
      allocate( del_eqload(nodof,num_threads) )                                 
      del_eqload(:,:) = zero                                                    
c                                                                               
c            parallel section. each processor runs down the entire do           
c            loop, but only performs work on the elements assigned              
c            to it. otherwise the loop cycles.                                  
c            each change in the eqload vector is marked by processor.           
c                                                                               
c$OMP PARALLEL PRIVATE(erow,now_thread,element,exe_thread)                      
      now_thread = omp_get_thread_num()+1                                       
      do erow = 1,size                                                          
         exe_thread = thread_number(erow)                                       
         if ( now_thread .ne. exe_thread ) cycle                                
         element    = eload_data(erow,1)                                        
         if (debug) write(*,*) 'erow,elem,myid,exid',                           
     &        erow,element,now_thread,exe_thread                                
         call do_elem_load( element, eload_data(erow,2),                        
     &        eload_data(erow,3), eload_val(erow),                              
     &        eload_pist(erow), del_eqload(1,now_thread), mult_fact )           
      end do                                                                    
c$OMP END PARALLEL                                                              
c                                                                               
c            reduce individual components of the array eqloads                  
c                                                                               
      do j = 1, num_threads                                                     
         do i = 1, nodof                                                        
            eqloads(i) = eqloads(i) + del_eqload(i,j)                           
         end do                                                                 
      end do                                                                    
c                                                                               
      deallocate( del_eqload )                                                  
c                                                                               
      if ( debug ) then                                                         
         nrmeq = zero                                                           
         do i = 1,nodof                                                         
            nrmeq = nrmeq + eqloads(i)*eqloads(i)                               
         end do                                                                 
         nrmeq = sqrt(nrmeq)                                                    
         write(*,*) 'Norm of equivalent loading vector: ',nrmeq                 
      end if                                                                    
                                                                                
c                                                                               
 9999 continue                                                                  
      if (debug) write (*,*) '<<<< leaving elem_load.'                          
      return                                                                    
      end                                                                       
c  ****************************************************************             
c  *                                                              *             
c  *                  subroutine: do_elem_load                    *             
c  *                                                              *             
c  *                       written by : jcs                       *             
c  *                                                              *             
c  *                   last modified : 11/15/11 jcs               *             
c  *                                                              *             
c  *     this subroutine calculates the equivalent nodal loads    *             
c  *     for all element loads on a 3-d isoparametrics            *             
c  *                                                              *             
c  ****************************************************************             
c                                                                               
c                                                                               
      subroutine do_elem_load( element, eload_type, eload_dof,                  
     &     elval, eload_pist, del_eqload, mult_fact )                           
      use global_data ! old common.main
c                                                                               
      use main_data                                                             
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c     parameter declarations                                                    
c     ----------------------                                                    
c                                                                               
      integer element, eload_type, eload_dof, eload_pist                        
      real elval                                                                
      double precision                                                          
     &     del_eqload(*), mult_fact                                             
c                                                                               
c     declare local variables                                                   
c     -----------------------                                                   
c                                                                               
      double precision                                                          
     &      ecoord(3,mxndel), zero, body_intens, face_intens,                   
     &      equiv_loads_ele(mxndel), equiv_loads_face(mxndel,3),                
     &      face_intensities(mxndel,3), dummy,                                  
     &      p3, u3, m3, gam, fdirc(3)                                           
c                                                                               
      double precision,                                                         
     & dimension(:,:), pointer :: eforces                                       
c                                                                               
      dimension elem_nodes(mxndel), edest(mxedof)                               
      logical debug, geonl, body_flg, face_flg, nrml_flg, tet_elem,             
     &        failed                                                            
      data debug, zero / .false., 0.0 /                                         
c                                                                               
c                                                                               
      if (debug) write(*,*) '>> In do_elem_load'                                
c                                                                               
      if (debug) then                                                           
         write(*,9010) element                                                  
         write(*,9020) eload_type                                               
         write(*,9030) eload_dof                                                
         write(*,9040) elval                                                    
         write(*,9050) eload_pist                                               
         write(*,9070) mult_fact                                                
      end if                                                                    
c                                                                               
c                    initialize loading flags for current element row           
c                                                                               
      etype    = iprops(1,element)                                              
      nnode    = iprops(2,element)                                              
      iorder   = iprops(5,element)                                              
      ngp      = iprops(6,element)                                              
      totdof   = nnode * 3                                                      
      geonl    = lprops(18,element)                                             
      body_flg = .false.                                                        
      face_flg = .false.                                                        
      nrml_flg = .false.                                                        
      if (debug) write (*,*) '>>>>>> ELEMENT:' , element                        
      call get_single_edest_terms( edest, element )                             
c                                                                               
      tet_elem = etype .eq. 6 .or. etype .eq. 13                                
c                                                                               
c                    get the coordinates of nodes for the current element       
c                         elem_nodes stores the structure node numbers          
c                         for the element nodes; ecoord stores                  
c                         the coordinates for the element nodes.                
c                         if we have geometric nonlinearity, then the           
c                         current coordinates are the original coords           
c                         (stored in c) plus the displacements (stored in       
c                         u).                                                   
c                                                                               
      if ( debug ) write (*,*) '     getting new coords'                        
      if ( geonl ) then                                                         
         do node = 1, nnode                                                     
            elem_nodes(node) = incid(incmap(element)+(node-1))                  
            ecoord(1,node) =                                                    
     &           c(crdmap(elem_nodes(node)))                                    
     &           + u(edest(node))                                               
            ecoord(2,node) =                                                    
     &           c(crdmap(elem_nodes(node))+1)                                  
     &           + u(edest(nnode+node))                                         
            ecoord(3,node) =                                                    
     &           c(crdmap(elem_nodes(node))+2)                                  
     &           + u(edest(2*nnode+node))                                       
         end do                                                                 
      else                                                                      
         do node = 1, nnode                                                     
            elem_nodes(node) = incid(incmap(element)+(node-1))                  
            ecoord(1,node) =                                                    
     &           c(crdmap(elem_nodes(node)))                                    
            ecoord(2,node) =                                                    
     &           c(crdmap(elem_nodes(node))+1)                                  
            ecoord(3,node) =                                                    
     &           c(crdmap(elem_nodes(node))+2)                                  
         end do                                                                 
      end if                                                                    
      if ( debug ) write (*,*) '     done getting new coords'                   
c                                                                               
c     build loading data variables.                                             
c                                                                               
      if ( eload_type .eq. 0 ) then                                             
c                                                                               
c                        1) body force -- build applied load vector             
c                                                                               
         body_flg    = .true.                                                   
         body_dir    = eload_dof                                                
         body_intens = elval                                                    
         if ( debug ) then                                                      
            write(*,*) '    body force. element',element                        
            write(*,*) 'body_flg,body_dir,body_intens:'                         
            write(*,*)  body_flg,body_dir,body_intens                           
         end if                                                                 
c                                                                               
      else if ( eload_type .gt. 0 ) then                                        
c                                                                               
c                        2) face force constant direction                       
c                                                                               
         face        = eload_type                                               
         face_dir    = eload_dof                                                
         face_intens = elval                                                    
         face_flg    = .true.                                                   
         if ( debug ) then                                                      
            write(*,*) '    face load, const. dir'                              
            write(*,*) 'face,face_dir,face_intens:'                             
            write(*,*) face,face_dir,face_intens                                
         end if                                                                 
c                                                                               
      else if ( eload_type .eq. -100 ) then                                     
c                                                                               
c                        3) element temperature. save temperature               
c                           change in global vector                             
c                                                                               
         if ( debug ) write (*,*) '    element temper'                          
         dtemp_elems(element) = dtemp_elems(element)+                           
     &        elval*mult_fact                                                   
         temperatures = .true.                                                  
c                                                                               
      else if ( (eload_type .le. -7)                                            
     &        .and. (eload_type .ge. -12 ) ) then                               
c                                                                               
c                        4) piston theory pressure                              
c                                                                               
         face       = abs(eload_type) - 6                                       
         face_dir   = eload_dof                                                 
         face_flg    = .true.                                                   
         nrml_flg    = .true.                                                   
c                                                                               
c                        call a separate subroutine piston loading              
c                        parameters at current time. linearly interpolate       
c                        between times listed in the table                      
c                                                                               
         call interp_piston_params(                                             
     &        p3, u3, m3, gam, fdirc, eload_pist, total_model_time )            
c                                                                               
c                        call a separate subroutine to determine the            
c                        face_intens from piston theory                         
c                                                                               
         call piston_face_intens(                                               
     &        face_intens, element, etype, nnode, face,                         
     &        p3, u3, m3, gam, fdirc, tet_elem, elem_nodes,                     
     &        ecoord )                                                          
c                                                                               
         if ( debug ) then                                                      
            write(*,*) '    face load, piston'                                  
            write(*,*) 'face,face_dir,face_intens:'                             
            write(*,*) face,face_dir,face_intens                                
         end if                                                                 
c                                                                               
      else                                                                      
c                                                                               
c                        5) pressure force                                      
c                                                                               
         face        = abs(eload_type)                                          
         face_dir    = eload_dof                                                
         face_intens = elval                                                    
         face_flg    = .true.                                                   
         nrml_flg    = .true.                                                   
         if ( debug ) then                                                      
            write(*,*) '    face load, normal'                                  
            write(*,*) 'face,face_dir,face_intens:'                             
            write(*,*) face,face_dir,face_intens                                
         end if                                                                 
      end if                                                                    
c                                                                               
c                    now process the face and body forces if they have          
c                    been set.                                                  
c                                                                               
c                       1) compute body forces and add to structure             
c                          incremental load vector (by nodal dof).              
c                                                                               
      if ( body_flg ) then                                                      
         if ( debug ) write (*,*) '  processing body forces.'                   
         call body_load( element, etype, nnode, body_dir,                       
     &        body_intens, ecoord, equiv_loads_ele )                            
c                                                                               
         call elem_eqload_allocate( element, nnode )                            
         eforces => elem_eq_loads(element)%forces                               
c                                                                               
         k = body_dir - 1                                                       
         do enode = 1, nnode                                                    
            snode = elem_nodes(enode)                                           
            del_eqload(dstmap(snode)+k) = del_eqload(dstmap(snode)+k) +         
     &           equiv_loads_ele(enode)*mult_fact                               
            eforces(body_dir,enode) = eforces(body_dir,enode) +                 
     &           equiv_loads_ele(enode)*mult_fact                               
         end do                                                                 
      end if                                                                    
c                                                                               
c                       2) face/pressure forces                                 
c                                                                               
      if ( face_flg ) then                                                      
         if ( debug ) write (*,*) '  processing face/pressure loads.'           
         do i = 1, nnode                                                        
            equiv_loads_face(i,1) = zero                                        
            equiv_loads_face(i,2) = zero                                        
            equiv_loads_face(i,3) = zero                                        
         end do                                                                 
         if ( .not. nrml_flg ) then                                             
            call face_load( element, etype, nnode, face,                        
     &           face_intensities(1,1), ecoord,                                 
     &           equiv_loads_face(1,face_dir),                                  
     &           1, face_intens )                                               
         end if                                                                 
         if ( nrml_flg ) then                                                   
c                                                                               
            if ( tet_elem ) then                                                
               call tet_face_nloads( element, etype, nnode, face,               
     &              face_intensities(1,1),                                      
     &              face_intensities(1,2),                                      
     &              face_intensities(1,3),                                      
     &              face_intens, ecoord, failed )                               
            else                                                                
               call hex_face_nloads( element, etype, nnode, face,               
     &              face_intensities(1,1),                                      
     &              face_intensities(1,2),                                      
     &              face_intensities(1,3),                                      
     &              face_intens, ecoord, failed )                               
            end if                                                              
            if ( failed ) return                                                
c                                                                               
            call face_load( element, etype, nnode, face,                        
     &           face_intensities(1,1), ecoord,                                 
     &           equiv_loads_face(1,1),                                         
     &           2, dummy )                                                     
            call face_load( element, etype, nnode, face,                        
     &           face_intensities(1,2), ecoord,                                 
     &           equiv_loads_face(1,2),                                         
     &           2, dummy )                                                     
            call face_load( element, etype, nnode, face,                        
     &           face_intensities(1,3), ecoord,                                 
     &           equiv_loads_face(1,3),                                         
     &           2, dummy )                                                     
         end if                                                                 
c                                                                               
         call elem_eqload_allocate( element, nnode )                            
         eforces => elem_eq_loads(element)%forces                               
c                                                                               
         do enode = 1, nnode                                                    
            snode = elem_nodes(enode)                                           
            del_eqload(dstmap(snode)+0) = del_eqload(dstmap(snode)+0) +         
     &           equiv_loads_face(enode,1)*mult_fact                            
            del_eqload(dstmap(snode)+1) = del_eqload(dstmap(snode)+1) +         
     &           equiv_loads_face(enode,2)*mult_fact                            
            del_eqload(dstmap(snode)+2) = del_eqload(dstmap(snode)+2) +         
     &                        equiv_loads_face(enode,3)*mult_fact               
            eforces(1,enode) = eforces(1,enode) +                               
     &           equiv_loads_face(enode,1)*mult_fact                            
            eforces(2,enode) = eforces(2,enode) +                               
     &           equiv_loads_face(enode,2)*mult_fact                            
            eforces(3,enode) = eforces(3,enode) +                               
     &           equiv_loads_face(enode,3)*mult_fact                            
         end do                                                                 
      end if                                                                    
c                                                                               
      if (debug) write(*,*) '>> Leaving do_elem_load'                           
c                                                                               
 9010 format(4x,'element number:  ',7x,i5)                                      
 9020 format(4x,'load type:       ',7x,i5)                                      
 9030 format(4x,'load dof:        ',7x,i5)                                      
 9040 format(4x,'eload_val:       ',9x,f10.6)                                   
 9050 format(4x,'eload_pist: tab #',i5)                                         
 9070 format(4x,'mult factor:     ',9x,f10.6)                                   
      end                                                                       
c                                                                               
c                                                                               
c ****************************************************************              
c *                                                              *              
c *          normal pressure decomposition for hex elements      *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
      subroutine hex_face_nloads( element, etype, nnode, face,                  
     &                        xload, yload, zload, face_intens,                 
     &                        ecoord, failed )                                  
      implicit double precision (a-h,o-z)                                       
c                                                                               
      dimension  xload(*), yload(*), zload(*), ecoord(3,*)                      
      integer element, etype, face                                              
      logical failed                                                            
c                                                                               
c                local arrays                                                   
c                                                                               
      double precision                                                          
     &    jacob, jacobi, nrmvec                                                 
c                                                                               
      dimension dsf(32,3), jacob(3,3), jacobi(3,3),                             
     &          nrmvec(3), fnodes(10), fcoor(3,10)                              
      logical local_debug, bad                                                  
      integer enode, fnodes                                                     
      data zero, local_debug / 0.0, .false. /                                   
c                                                                               
c             get a list of element node numbers for the loaded                 
c             face and the isoparametric coordinates of nodes                   
c             on the face.                                                      
c                                                                               
      failed = .false.                                                          
      do i = 1, nnode                                                           
        xload(i) = zero                                                         
        yload(i) = zero                                                         
        zload(i) = zero                                                         
      end do                                                                    
      if ( local_debug ) write(*,1000) etype, face                              
      call eqelfn( fnodes, etype, face, nfnode )                                
      call ndpts1( fnodes, nfnode, fcoor, etype, -1, dum1, dum2, dum3 )         
      if ( local_debug ) then                                                   
          write(*,1010) nfnode, (fnodes(i),i=1,nfnode )                         
          write(*,1020) ((fcoor(i,j),i=1,3),j=1,nfnode)                         
      end if                                                                    
c                                                                               
c             loop over all nodes on the loaded face:                           
c                                                                               
c              1) evaluate derivatives of shape functions                       
c                 at the node;                                                  
c              2) compute the jacobian at the node;                             
c              3) select vectors veca,vecb that are in the tangent              
c                 plane to the face at the node.  rows of the                   
c                 jacobian are components of the vectors.                       
c              4) construct an inward normal vector to the face                 
c                 (veca x vecb) and normalize to unit length                    
c                 thus producing direction cosines;                             
c              5) if the normal length <= a small tolerance, we probably        
c                 have a corner node on a collased face. back off from          
c                 the corner node by 10% in parametric space and try to         
c                 compute normal again. if that fails, omit load                
c                 computations for this face of this element and issue          
c                 message                                                       
c              6) x,y,z components of load at the node are the                  
c                 pressure intensity * the corresponding direction              
c                 cosine.                                                       
c                                                                               
      do node = 1, nfnode                                                       
        xsi  = fcoor(1,node)                                                    
        eta  = fcoor(2,node)                                                    
        zeta = fcoor(3,node)                                                    
        call derivs( etype, xsi, eta, zeta, dsf(1,1), dsf(1,2),                 
     &               dsf(1,3) )                                                 
        call eqldjb( dsf, ecoord, nnode, jacob, jacobi, det, ierr )             
        call eqnrmvh( face, jacob, nrmvec, bad, local_debug )                   
        if( bad ) then                                                          
          call eqnrmvh_adjust( face, xsi, eta, zeta )                           
          call derivs( etype, xsi, eta, zeta, dsf(1,1), dsf(1,2),               
     &                 dsf(1,3) )                                               
          call eqldjb( dsf, ecoord, nnode, jacob, jacobi, det, ierr )           
          call eqnrmvh( face, jacob, nrmvec, bad, local_debug )                 
          if( bad ) then                                                        
           call iodevn( idummy, iout, dummy, 1 )                                
           write (iout,1050) face, element                                      
           failed = .true.                                                      
           return                                                               
          end if                                                                
        end if                                                                  
        enode        = fnodes(node)                                             
        xload(enode) = face_intens * nrmvec(1)                                  
        yload(enode) = face_intens * nrmvec(2)                                  
        zload(enode) = face_intens * nrmvec(3)                                  
        if ( local_debug ) write(*,1030) enode, nrmvec                          
      end do                                                                    
c                                                                               
      if ( local_debug ) write(*,1040) (i,xload(i),yload(i),zload(i),           
     &                                  i = 1, nnode )                          
c                                                                               
      return                                                                    
c                                                                               
 1000 format(2x, 'normal load components: ' ,//,                                
     &        5x, 'element type ',i3, 2x,'face ', i3 /)                         
 1010 format(' number of face nodes:',i3,                                       
     &        /,2x,'face nodes:',2(/,2x,10i5) )                                 
 1020 format( 2x,'face node coords:',20(/,2x,3f5.1) )                           
 1030 format('enode:',i3,2x,'nrmvec:',3f15.9 )                                  
 1040 format(2x,'nodal load components: ',                                      
     &        40(/,2x,i4,3f15.8) )                                              
 1050 format('>>> WARNING: pressure load placed on face ',i1,' of ',            
     &       'element ',i8,'.  This face',/,                                    
     &       '             has at least one collapsed edge. Warp3D',            
     & /,    '             cannot resolve the normal direction at',/,           
     &       '             one or more nodes on the loaded face.',/,            
     &       '             This loading ignored...',/)                          
c                                                                               
       end                                                                      
c ****************************************************************              
c *                                                              *              
c *          adjust parametric coordinates for node on collapsed *              
c *          face of hex element                                 *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
      subroutine  eqnrmvh_adjust( face, xsi, eta, zeta )                        
      implicit none                                                             
c                                                                               
      integer face                                                              
      double precision                                                          
     &  xsi, eta, zeta, pt9                                                     
      data pt9 / 0.9 /                                                          
c                                                                               
      if( face .eq. 1 .or. face .eq. 2 ) then                                   
        eta  = eta * pt9                                                        
        zeta = zeta * pt9                                                       
        return                                                                  
      end if                                                                    
      if( face .eq. 3 .or. face .eq. 4 ) then                                   
        xsi  = xsi * pt9                                                        
        zeta = zeta * pt9                                                       
        return                                                                  
      end if                                                                    
      if( face .eq. 5 .or. face .eq. 6 ) then                                   
        xsi  = xsi * pt9                                                        
        eta  = eta * pt9                                                        
        return                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c ****************************************************************              
c *                                                              *              
c *          normal pressure decomposition for tet elements      *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
      subroutine tet_face_nloads( element, etype, nnode, face,                  
     &                        xload, yload, zload, face_intens,                 
     &                        ecoord, failed )                                  
      implicit double precision (a-h,o-z)                                       
c                                                                               
      dimension  xload(*), yload(*), zload(*), ecoord(3,*)                      
      integer element, etype, face                                              
      logical failed                                                            
c                                                                               
c                local arrays                                                   
c                                                                               
      double precision                                                          
     &    nrmvec, mag1, mag2, mag_max, vec1, vec2                               
c                                                                               
      dimension dsf(32,3), temp_dsf(32,3),                                      
     &          sf (32), temp_sf(32),                                           
     &          nrmvec(3), vec1(3), vec2(3),                                    
     &          fnodes(10), fcoor_tri6(3,6),                                    
     &          fcoor_tri3(3,3)                                                 
      logical local_debug, etype_tet10, etype_tet4                              
      integer enode, fnodes, i                                                  
      data   zero, local_debug, tol                                             
     &  / 0.d0, .false., 1.0d-06 /                                              
c                                                                               
c                                                                               
c                    use shape function derivatives for triangular              
c                    interface elements for tractions on the tet                
c                    faces. this is done to avoid taking derivatives            
c                    with respect to the dependent natural                      
c                    coordinate.                                                
c                                                                               
c                    because we're working in 2D, the isoparametric             
c                    coordinates will be the same regardless of the             
c                    face that is loaded.                                       
c                                                                               
c                                                                               
      data fcoor_tri6                                                           
     & / 0.d0,  0.d0,  1.d0,                                                    
     &   0.5d0, 0.d0,  0.5d0,                                                   
     &   1.d0,  0.d0,  0.d0,                                                    
     &   0.5d0, 0.5d0, 0.d0,                                                    
     &   0.d0,  1.d0,  0.d0,                                                    
     &   0.d0,  0.5d0, 0.5d0 /                                                  
c                                                                               
      data fcoor_tri3                                                           
     & / 0.d0, 0.d0, 1.d0,                                                      
     &   1.d0, 0.d0, 0.d0,                                                      
     &   0.d0, 1.d0, 0.d0 /                                                     
c                                                                               
c                                                                               
c                     determine if linear or quadratic tet                      
c                                                                               
c                                                                               
      failed      = .false.                                                     
      etype_tet10 = etype .eq. 6                                                
      etype_tet4  = etype .eq. 13                                               
c                                                                               
c                                                                               
      do i = 1, nnode                                                           
        xload(i) = zero                                                         
        yload(i) = zero                                                         
        zload(i) = zero                                                         
      end do                                                                    
      if ( local_debug ) write(*,1000) etype, face                              
c                                                                               
c                                                                               
c                    get the tet nodes on the face for mapping                  
c                    derivatives to the correct nodal coordinates.              
c                                                                               
c                                                                               
      call tet_get_nodes( fnodes, etype, face, nfnode )                         
c                                                                               
c                                                                               
c                    the index variable causes the 'derivs'                     
c                    routine to go to the correct tri                           
c                    derivatives.                                               
c                                                                               
c                                                                               
      if( etype_tet4 ) then                                                     
        index = 14                                                              
      else if ( etype_tet10 ) then                                              
        index = 15                                                              
      end if                                                                    
c                                                                               
c             loop over all nodes on the loaded face:                           
c                                                                               
c              1) evaluate derivatives of shape functions                       
c                 at the node;                                                  
c              2) build vec1, vec2                                              
c              3) construct an inward normal vector to the face                 
c                 (vec1 x vec2) and normalize to unit length                    
c                 thus producing direction cosines;                             
c              4) x,y,z components of load at the node are the                  
c                 pressure intensity * the corresponding direction              
c                 cosine.                                                       
c                                                                               
c                                                                               
      do node = 1, nfnode                                                       
        if( etype_tet4 ) then                                                   
          index = 14                                                            
          call derivs( index, fcoor_tri3(1,node), fcoor_tri3(2,node),           
     &                 fcoor_tri3(3,node), temp_dsf(1,1), temp_dsf(1,2),        
     &                 temp_dsf(1,3) )                                          
        else if( etype_tet10 ) then                                             
          index = 15                                                            
          call derivs( index, fcoor_tri6(1,node), fcoor_tri6(2,node),           
     &                 fcoor_tri6(3,node), temp_dsf(1,1), temp_dsf(1,2),        
     &                 temp_dsf(1,3) )                                          
        end if                                                                  
        call map_tet_tri( etype, nfnode, fnodes, sf, temp_sf,                   
     &                    dsf, temp_dsf )                                       
        call tet_darea( dsf, nfnode, fnodes, ecoord, darea,                     
     &                  vec1, vec2 )                                            
c                                                                               
c             check to see if either vector has zero length.                    
c                                                                               
        mag1 = sqrt (vec1(1)**2 + vec1(2)**2 + vec1(3)**2 )                     
        mag2 = sqrt (vec2(1)**2 + vec2(2)**2 + vec2(3)**2 )                     
        mag_max = max( mag1, mag2 )                                             
        if( mag_max .eq. zero ) then                                            
           call iodevn( idummy, iout, dummy, 1 )                                
           write (iout,1050) face, element                                      
           failed = .true.                                                      
           return                                                               
        end if                                                                  
c                                                                               
        mag1 = mag1 / mag_max                                                   
        mag2 = mag2 / mag_max                                                   
        if( mag1 .le. tol .or. mag2 .le. tol ) then                             
           call iodevn( idummy, iout, dummy, 1 )                                
           write (iout,1050) face, element                                      
           failed = .true.                                                      
           return                                                               
        end if                                                                  
c                                                                               
c             generate unit normal to surface at node using                     
c             cross product of tangent plane vectors.                           
c             then compute x,y,z components of load.                            
c                                                                               
        nrmvec(1) = vec1(2)*vec2(3) - vec2(2)*vec1(3)                           
        nrmvec(2) = vec2(1)*vec1(3) - vec1(1)*vec2(3)                           
        nrmvec(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)                           
        rlen = sqrt( nrmvec(1)*nrmvec(1) + nrmvec(2)*nrmvec(2) +                
     &               nrmvec(3)*nrmvec(3) )                                      
        nrmvec(1) = nrmvec(1) / rlen                                            
        nrmvec(2) = nrmvec(2) / rlen                                            
        nrmvec(3) = nrmvec(3) / rlen                                            
c                                                                               
        enode        = fnodes(node)                                             
        xload(enode) = -face_intens * nrmvec(1)                                 
        yload(enode) = -face_intens * nrmvec(2)                                 
        zload(enode) = -face_intens * nrmvec(3)                                 
        if ( local_debug ) write(*,1030) enode, nrmvec                          
      end do                                                                    
c                                                                               
      if ( local_debug ) write(*,1040) (i,xload(i),yload(i),zload(i),           
     &                                  i = 1, nnode )                          
c                                                                               
      return                                                                    
c                                                                               
 1000 format(2x, 'normal load components: ' ,//,                                
     &        5x, 'element type ',i3, 2x,'face ', i3 /)                         
 1010 format(' number of face nodes:',i3,                                       
     &        /,2x,'face nodes:',2(/,2x,10i5) )                                 
 1020 format( 2x,'face node coords:',20(/,2x,3f5.1) )                           
 1030 format('enode:',i3,2x,'nrmvec:',3f7.4 )                                   
 1040 format(2x,'nodal load components: ',                                      
     &        40(/,2x,i4,3f15.5) )                                              
 1050 format('>>> WARNING: pressure load placed on face ',i1,' of ',            
     &       'element ',i8,'.  This face',/,                                    
     &       '             has at least one collapsed edge. Warp3D',            
     & /,    '             cannot resolve the normal direction at',/,           
     &       '             one or more nodes on the loaded face.',/,            
     &       '             This loading ignored...',/)                          
c                                                                               
       end                                                                      
                                                                                
c                                                                               
c ****************************************************************              
c *                                                              *              
c *          elem_load_allocate                                  *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
      subroutine elem_eqload_allocate( element, nnode )                         
      use main_data                                                             
c                                                                               
      implicit integer (a-z)                                                    
      double precision                                                          
     & zero                                                                     
                                                                                
      data zero / 0.0 /                                                         
c                                                                               
      elem_equiv_loads_now = .true.                                             
c                                                                               
      if( elem_eq_loads(element)%ncols .eq. 0 ) then                            
        allocate( elem_eq_loads(element)%forces(3,nnode) )                      
        elem_eq_loads(element)%ncols = nnode                                    
        do enode = 1, nnode                                                     
           elem_eq_loads(element)%forces(1:3,enode) = zero                      
        end do                                                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ****************************************************************              
c *                                                              *              
c *   normal vector for pressure decomposition - hex elements    *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
      subroutine eqnrmvh( face, jacob, nrmvec, bad, debug )                     
      implicit none                                                             
c                                                                               
      double precision                                                          
     &  jacob(3,3), nrmvec(3)                                                   
      integer face                                                              
      logical bad, debug                                                        
c                                                                               
c                local arrays                                                   
c                                                                               
      double precision                                                          
     &     maga, magb, mag_max, veca(3), vecb(3), rlen,                         
     &     zero, tol, mag_min                                                   
      integer irowa(6), irowb(6), rowa, rowb                                    
c                                                                               
      data zero, tol / 0.0, 0.0001 /                                            
      data irowa /  2,3,3,1,1,2 /                                               
      data irowb /  3,2,1,3,2,1 /                                               
c                                                                               
      bad     = .false.                                                         
      rowa    = irowa(face)                                                     
      rowb    = irowb(face)                                                     
      veca(1) = jacob(rowa,1)                                                   
      veca(2) = jacob(rowa,2)                                                   
      veca(3) = jacob(rowa,3)                                                   
      vecb(1) = jacob(rowb,1)                                                   
      vecb(2) = jacob(rowb,2)                                                   
      vecb(3) = jacob(rowb,3)                                                   
c      write(*,*) '... eqnrmvh ...'                                             
c      write(*,9900) veca, vecb                                                 
                                                                                
c                                                                               
c             Error conditions:                                                 
c               (1) either vector has exactly zero length                       
c               (2) shorter vector length is < tol2 * longer                    
c                   vector length                                               
c                                                                               
        maga = sqrt (veca(1)**2 + veca(2)**2 + veca(3)**2 )                     
        magb = sqrt (vecb(1)**2 + vecb(2)**2 + vecb(3)**2 )                     
        mag_max = maga                                                          
        mag_min = magb                                                          
c        write(*,9910) maga, magb                                               
        if( magb .gt. maga ) then                                               
          mag_max = magb                                                        
          mag_min = maga                                                        
        end if                                                                  
        if( (mag_min .eq. zero) .or.                                            
     &      (mag_min .lt. tol*mag_max) ) then                                   
           bad = .true.                                                         
           return                                                               
        end if                                                                  
        if( debug ) then                                                        
          write(*,*) '..normal vector computation..'                            
          write(*,*) '   rowa, rowb: ', rowa, rowb                              
          write(*,9000) jacob(1,1:3), jacob(2,1:3),                             
     &           jacob(3,1:3)                                                   
 9000     format(3(5x,3f15.9,/) )                                               
          write(*,*) '    maga, magb: ', maga, magb                             
        end if                                                                  
c                                                                               
c             generate unit normal to surface at node using                     
c             cross product of tangent plane vectors.                           
c             then compute x,y,z components of load.                            
c                                                                               
c             colinear or nearly colinear vectors lead                          
c             to zero or nearly zero length of cross                            
c             product vector                                                    
c                                                                               
c             separately, veca and vecb are ok. normalize them                  
c             before taking cross product.                                      
c             length of the cross-product (normal) vector is then               
c             a good measure of the colinearity of veca and                     
c             vecb.                                                             
c                                                                               
        veca(1:3) = veca(1:3) / maga                                            
        vecb(1:3) = vecb(1:3) / magb                                            
        nrmvec(1) = veca(2)*vecb(3) - vecb(2)*veca(3)                           
        nrmvec(2) = vecb(1)*veca(3) - veca(1)*vecb(3)                           
        nrmvec(3) = veca(1)*vecb(2) - vecb(1)*veca(2)                           
        rlen = sqrt( nrmvec(1)*nrmvec(1) + nrmvec(2)*nrmvec(2) +                
     &               nrmvec(3)*nrmvec(3) )                                      
c        if( debug ) then                                                       
c         write(*,9920) nrmvec, rlen                                            
c        end if                                                                 
        if( rlen .lt. tol*10.0 ) then                                           
           bad = .true.                                                         
           return                                                               
        end if                                                                  
c                                                                               
        nrmvec(1) = nrmvec(1) / rlen                                            
        nrmvec(2) = nrmvec(2) / rlen                                            
        nrmvec(3) = nrmvec(3) / rlen                                            
c                                                                               
        return                                                                  
 9900   format(3x, 'veca: ',3e14.6,/,/,3x, 'vecb: ',3e14.6)                     
 9910   format(3x,'maga, magb: ', 2e14.6)                                       
 9920   format(3x,'nrmvec: ',3e14.6,' rlen: ',e14.6)                            
                                                                                
                                                                                
        end                                                                     
c                                                                               
c                                                                               
c ****************************************************************              
c *                                                              *              
c *          interp_piston_params                                *              
c *                                                              *              
c ****************************************************************              
c                                                                               
c                                                                               
      subroutine interp_piston_params(                                          
     &     p3, u3, m3, gam, fdirc, tabn, cmptime )                              
      use main_data, only: tables                                               
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c     parameter declarations                                                    
c     ----------------------                                                    
c                                                                               
      double precision                                                          
     &     p3, u3, m3, gam, fdirc(3), cmptime                                   
c                                                                               
c     declare local variables                                                   
c     -----------------------                                                   
c                                                                               
      double precision                                                          
     &     t0, tf, tn, ti, pn, un, mn, gn, fn(3),                               
     &     tn1, pn1, un1, mn1, gn1, fn1(3), one, fac1, fac2,                    
     &     mx, my, mz, norm                                                     
      logical debug                                                             
      data zero, one, debug / 0.0d0, 1.0d0, .false. /                           
c                                                                               
      if (debug) write (*,*) ' >>>>> inside of interp_piston_params'            
c                                                                               
c                 initialize number of rows in the table                        
c                                                                               
      nrows = tables(tabn)%num_rows                                             
c                                                                               
c                 check if cmptime is either at the beginning                   
c                 of the end of the loading history, if so                      
c                 no need for interpolatation. set parameters and               
c                 return to calling function                                    
c                                                                               
      t0  = tables(tabn)%table_values_sgl(1,1)                                  
      tf  = tables(tabn)%table_values_sgl(nrows,1)                              
      if (debug) write(*,9000) cmptime, tabn, t0, tf                            
c                                                                               
      if ( cmptime .eq. t0 ) then                                               
         if (debug) write(*,9010)                                               
         p3       = tables(tabn)%table_values_sgl(1,2)                          
         u3       = tables(tabn)%table_values_sgl(1,3)                          
         m3       = tables(tabn)%table_values_sgl(1,4)                          
         gam      = tables(tabn)%table_values_sgl(1,5)                          
         fdirc(1) = tables(tabn)%table_values_sgl(1,6)                          
         fdirc(2) = tables(tabn)%table_values_sgl(1,7)                          
         fdirc(3) = tables(tabn)%table_values_sgl(1,8)                          
         go to 9999                                                             
      end if                                                                    
c                                                                               
      if ( cmptime .gt. tf ) then                                               
         if (debug) write(*,9020)                                               
         p3       = tables(tabn)%table_values_sgl(nrows,2)                      
         u3       = tables(tabn)%table_values_sgl(nrows,3)                      
         m3       = tables(tabn)%table_values_sgl(nrows,4)                      
         gam      = tables(tabn)%table_values_sgl(nrows,5)                      
         fdirc(1) = tables(tabn)%table_values_sgl(nrows,6)                      
         fdirc(2) = tables(tabn)%table_values_sgl(nrows,7)                      
         fdirc(3) = tables(tabn)%table_values_sgl(nrows,8)                      
         go to 9999                                                             
      end if                                                                    
c                                                                               
c                 need to linearly interpolate to determine                     
c                 piston parameters. search for time interval                   
c                 that cmptime falls into then perform interpolation            
c                                                                               
      do row = 2,nrows                                                          
         tn1 = tables(tabn)%table_values_sgl(row,1)                             
c                                                                               
c                 time increase monotonically. search for                       
c                 row where defined time is greater than current                
c                 time                                                          
c                                                                               
         if ( .not. ( tn1 .gt. cmptime ) ) cycle                                
         tn  = tables(tabn)%table_values_sgl(row-1,1)                           
         ti = cmptime - tn                                                      
         fac1 = one - ti/(tn1-tn)                                               
         fac2 = ti/(tn1-tn)                                                     
c                                                                               
         if (debug)                                                             
     &        write(*,9030) tn, tn1, ti, fac1, fac2                             
c                                                                               
c                 get piston paramters at tn and tn1                            
c                                                                               
         pn     = tables(tabn)%table_values_sgl(row-1,2)                        
         un     = tables(tabn)%table_values_sgl(row-1,3)                        
         mn     = tables(tabn)%table_values_sgl(row-1,4)                        
         gn     = tables(tabn)%table_values_sgl(row-1,5)                        
         fn(1)  = tables(tabn)%table_values_sgl(row-1,6)                        
         fn(2)  = tables(tabn)%table_values_sgl(row-1,7)                        
         fn(3)  = tables(tabn)%table_values_sgl(row-1,8)                        
c                                                                               
         pn1    = tables(tabn)%table_values_sgl(row,2)                          
         un1    = tables(tabn)%table_values_sgl(row,3)                          
         mn1    = tables(tabn)%table_values_sgl(row,4)                          
         gn1    = tables(tabn)%table_values_sgl(row,5)                          
         fn1(1) = tables(tabn)%table_values_sgl(row,6)                          
         fn1(2) = tables(tabn)%table_values_sgl(row,7)                          
         fn1(3) = tables(tabn)%table_values_sgl(row,8)                          
c                                                                               
c                 perform interpolation, ensure that flow                       
c                 direction vectors have a unit normal                          
c                                                                               
         p3  = fac1*pn + fac2*pn1                                               
         u3  = fac1*un + fac2*un1                                               
         m3  = fac1*mn + fac2*mn1                                               
         gam = fac1*gn + fac2*gn1                                               
         mx  = fac1*fn(1) + fac2*fn1(1)                                         
         my  = fac1*fn(2) + fac2*fn1(2)                                         
         mz  = fac1*fn(3) + fac2*fn1(3)                                         
c                                                                               
         norm = sqrt( mx*mx + my*my + mz*mz )                                   
         if ( norm .eq. zero ) then                                             
            call errmsg(340,dum,dums,dumr,dumd)                                 
         end if                                                                 
c                                                                               
         fdirc(1) = mx/norm                                                     
         fdirc(2) = my/norm                                                     
         fdirc(3) = mz/norm                                                     
c                                                                               
         go to 9999                                                             
      end do                                                                    
c                                                                               
 9999 continue                                                                  
      if (debug) write(*,9040) p3, u3, m3, gam, (fdirc(k),k=1,3)                
      if (debug) write(*,*) ' >>>>> leaving of interp_piston_params'            
c                                                                               
      return                                                                    
 9000 format(3x,'Current time         ',f13.6,/,                                
     &       3x,'Table number:        ',i6,/,                                   
     &       3x,'Initial time:        ',f13.6,/,                                
     &       3x,'Final time:          ',f13.6 )                                 
 9010 format(3x,'Using piston parameters at initial time.'/)                    
 9020 format(3x,'Using piston parameters at final time.'/)                      
 9030 format(3x,'Current time between ',f13.6,' and ',f13.6,/,                  
     &       3x,'Interpolation factors',3f13.6,/)                               
 9040 format(3x,'p3,   u3:   ',4x,2f13.6,/,                                     
     &       3x,'m3,   gam:  ',4x,2f13.6,/,                                     
     &       3x,'flow direction: ',3f13.6,/)                                    
c                                                                               
      end                                                                       
