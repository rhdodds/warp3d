c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine trnvec                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/01/90                   *          
c     *                                                              *          
c     *     this subroutine transforms element vectors to and from   *          
c     *     constraint compatable global coordinates and uniform     *          
c     *     global coordinates. the element vectors are for one      *          
c     *     element in a block of similar, non-conflicting elements. *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine trnvec(vec,trnmte,trne,ndof,nnode,bele,switch)                 
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision                                                          
     & vec(mxvl,*), trnmte(mxvl,mxedof,*), sums(mxndof)                         
      logical trne(mxvl,*)                                                      
c                                                                               
c                                                                               
c                       note: ndof is an element dependent variable,            
c                       but since ndof= 3 for every element in the              
c                       program, this is implicitly assumed in this             
c                       code. if ever a conflicting element is added            
c                       to the prgram, the code will need to be mod-            
c                       ified.                                                  
c                                                                               
c                                                                               
c                       if the transformation matrix for a node is the          
c                       identity matrix, skip the transformation. the           
c                       the values of switch and type indicate the              
c                       following:                                              
c                                                                               
c                          switch= 1: transformation from uniform global        
c                                     coords. to constraint compatable          
c                                     global coords.                            
c                                                                               
c                          switch= 2: transformation from constraint            
c                                     compatable global coords. to              
c                                     uniform global coords.                    
c                                                                               
c                       the transformation matrices are defined as              
c                       taking displacements in uniform global coords.          
c                       to those in constraint compatable global coords.        
c                       if one transforms from constraint compatable            
c                       global coords. to uniform global coords., i.e.          
c                       the reverse transformation, one uses the                
c                       transpose of the matrices used in the original          
c                       transformation. thus the following matrices are         
c                       used for the two values of switch:                      
c                                                                               
c                          switch= 1 : trans. matrices as stored                
c                                                                               
c                          switch= 2 : transpose                                
c                                                                               
      do nod = 1,nnode                                                          
         if( trne(bele,nod) ) then                                              
            ts = (nod-1)*ndof                                                   
            if( switch.eq.1 ) then                                              
               sums(1)= trnmte(bele,ts+1,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+1,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+1,3)*vec(bele,2*nnode+nod)               
               sums(2)= trnmte(bele,ts+2,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+2,3)*vec(bele,2*nnode+nod)               
               sums(3)= trnmte(bele,ts+3,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+3,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,3)*vec(bele,2*nnode+nod)               
            else                                                                
               sums(1)= trnmte(bele,ts+1,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,1)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,1)*vec(bele,2*nnode+nod)               
               sums(2)= trnmte(bele,ts+1,2)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,2)*vec(bele,2*nnode+nod)               
               sums(3)= trnmte(bele,ts+1,3)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,3)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,3)*vec(bele,2*nnode+nod)               
            end if                                                              
            vec(bele,nod)         = sums(1)                                     
            vec(bele,nnode+nod)   = sums(2)                                     
            vec(bele,2*nnode+nod) = sums(3)                                     
         end if                                                                 
      end do                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine trnvec                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/01/90                   *          
c     *                                                              *          
c     *     this subroutine transforms element vectors to and from   *          
c     *     constraint compatable global coordinates and uniform     *          
c     *     global coordinates. the element vectors are for one      *          
c     *     element in a block of similar, non-conflicting elements. *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine trnvecs(vec,trnmte,trne,ndof,nnode,bele,switch,                
     &                   nrow_vec)                                              
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision                                                          
     & vec(nrow_vec,*), trnmte(mxvl,mxedof,*), sums(mxndof)                     
      logical trne(mxvl,*)                                                      
c                                                                               
c                                                                               
c                       note: ndof is an element dependent variable,            
c                       but since ndof= 3 for every element in the              
c                       program, this is implicitly assumed in this             
c                       code. if ever a conflicting element is added            
c                       to the prgram, the code will need to be mod-            
c                       ified.                                                  
c                                                                               
c                                                                               
c                       if the transformation matrix for a node is the          
c                       identity matrix, skip the transformation. the           
c                       the values of switch and type indicate the              
c                       following:                                              
c                                                                               
c                          switch= 1: transformation from uniform global        
c                                     coords. to constraint compatable          
c                                     global coords.                            
c                                                                               
c                          switch= 2: transformation from constraint            
c                                     compatable global coords. to              
c                                     uniform global coords.                    
c                                                                               
c                       the transformation matrices are defined as              
c                       taking displacements in uniform global coords.          
c                       to those in constraint compatable global coords.        
c                       if one transforms from constraint compatable            
c                       global coords. to uniform global coords., i.e.          
c                       the reverse transformation, one uses the                
c                       transpose of the matrices used in the original          
c                       transformation. thus the following matrices are         
c                       used for the two values of switch:                      
c                                                                               
c                          switch= 1 : trans. matrices as stored                
c                                                                               
c                          switch= 2 : transpose                                
c                                                                               
      do nod = 1,nnode                                                          
         if( trne(bele,nod) ) then                                              
            ts = (nod-1)*ndof                                                   
            if( switch.eq.1 ) then                                              
               sums(1)= trnmte(bele,ts+1,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+1,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+1,3)*vec(bele,2*nnode+nod)               
               sums(2)= trnmte(bele,ts+2,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+2,3)*vec(bele,2*nnode+nod)               
               sums(3)= trnmte(bele,ts+3,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+3,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,3)*vec(bele,2*nnode+nod)               
            else                                                                
               sums(1)= trnmte(bele,ts+1,1)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,1)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,1)*vec(bele,2*nnode+nod)               
               sums(2)= trnmte(bele,ts+1,2)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,2)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,2)*vec(bele,2*nnode+nod)               
               sums(3)= trnmte(bele,ts+1,3)*vec(bele,nod)+                      
     &                  trnmte(bele,ts+2,3)*vec(bele,nnode+nod)+                
     &                  trnmte(bele,ts+3,3)*vec(bele,2*nnode+nod)               
            end if                                                              
            vec(bele,nod)         = sums(1)                                     
            vec(bele,nnode+nod)   = sums(2)                                     
            vec(bele,2*nnode+nod) = sums(3)                                     
         end if                                                                 
      end do                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
