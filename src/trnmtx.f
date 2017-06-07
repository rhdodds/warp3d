c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine trnmtx                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/02/91                   *          
c     *                                                              *          
c     *     this subroutine transforms an element stiffness or       *          
c     *     mass matrix from uniform global coordinates to           *          
c     *     constraint compatable global coordinates by virtue of    *          
c     *     the equation                             t               *          
c     *                     [mat] =  [t] * [mat]  * [t]              *          
c     *                         ccg            ug                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine trnmtx( mat, cp, trnmte, trne, ndof, nnode, totdof,            
     &                   bele, nsize )                                          
      use main_data, only: asymmetric_assembly                                  
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      double precision                                                          
     & mat(nsize,*), trnmte(mxvl,mxedof,*), sums(mxedof,mxndof)                 
        double precision,                                                       
     &      dimension(mxedof,mxedof) :: tmat                                    
      dimension cp(*)                                                           
      logical trne(mxvl,*)                                                      
c                                                                               
      if (.not. asymmetric_assembly) utsz = nsize                               
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
c                       expand the matrix from global upper triangular          
c                       vector form to full matrix form.                        
c                                                                               
c                                                                               
      if (.not. asymmetric_assembly) then                                       
      do j = 1, totdof                                                          
        do i = 1, j                                                             
            tmat(i,j) = mat(cp(j)+i,bele)                                       
         end do                                                                 
         do i = 1, j-1                                                          
            tmat(j,i) = tmat(i,j)                                               
         end do                                                                 
      end do                                                                    
      else                                                                      
            tmat(1:totdof,1:totdof) = reshape(mat(1:totdof*totdof,bele),        
     &            (/totdof, totdof/))                                           
      end if                                                                    
                                                                                
c                                                                               
c                                                                               
c                       transform. first multiply the matrix by the             
c                       transformation matrix. if the transformation            
c                       matrix for a node is identity, skip the                 
c                       multiplication. for each node, store the product        
c                       in temporary storage. place the results in the          
c                       matrix after the node is processed.                     
c                                                                               
c                                                                               
      do nod = 1, nnode                                                         
         if( trne(bele,nod) ) then                                              
            ts = (nod-1)*ndof                                                   
            do j = 1, totdof                                                    
               sums(j,1)= trnmte(bele,ts+1,1)*tmat(nod,j)+                      
     &                    trnmte(bele,ts+1,2)*tmat(nnode+nod,j)+                
     &                    trnmte(bele,ts+1,3)*tmat(2*nnode+nod,j)               
c                                                                               
               sums(j,2)= trnmte(bele,ts+2,1)*tmat(nod,j)+                      
     &                    trnmte(bele,ts+2,2)*tmat(nnode+nod,j)+                
     &                    trnmte(bele,ts+2,3)*tmat(2*nnode+nod,j)               
c                                                                               
               sums(j,3)= trnmte(bele,ts+3,1)*tmat(nod,j)+                      
     &                    trnmte(bele,ts+3,2)*tmat(nnode+nod,j)+                
     &                    trnmte(bele,ts+3,3)*tmat(2*nnode+nod,j)               
            end do                                                              
            do j = 1, totdof                                                    
               tmat(nod,j)         = sums(j,1)                                  
               tmat(nnode+nod,j)   = sums(j,2)                                  
               tmat(2*nnode+nod,j) = sums(j,3)                                  
            end do                                                              
         end if                                                                 
      end do                                                                    
c                                                                               
c                                                                               
c                       multiply the partially transformed matrix by the        
c                       transpose of the transformation matrix. if the          
c                       transformation matrix for a node is identity,           
c                       skip the multiplication. for each node, store           
c                       the product in temporary storage and place the          
c                       results in the matrix after the node is processed.      
c                                                                               
c                                                                               
      do nod = 1, nnode                                                         
         if( trne(bele,nod) ) then                                              
            ts = (nod-1)*ndof                                                   
            do j = 1, totdof                                                    
               sums(j,1)= trnmte(bele,ts+1,1)*tmat(j,nod)+                      
     &                    trnmte(bele,ts+1,2)*tmat(j,nnode+nod)+                
     &                    trnmte(bele,ts+1,3)*tmat(j,2*nnode+nod)               
               sums(j,2)= trnmte(bele,ts+2,1)*tmat(j,nod)+                      
     &                    trnmte(bele,ts+2,2)*tmat(j,nnode+nod)+                
     &                    trnmte(bele,ts+2,3)*tmat(j,2*nnode+nod)               
               sums(j,3)= trnmte(bele,ts+3,1)*tmat(j,nod)+                      
     &                    trnmte(bele,ts+3,2)*tmat(j,nnode+nod)+                
     &                    trnmte(bele,ts+3,3)*tmat(j,2*nnode+nod)               
            end do                                                              
            do j = 1, totdof                                                    
               tmat(j,nod)         = sums(j,1)                                  
               tmat(j,nnode+nod)   = sums(j,2)                                  
               tmat(j,2*nnode+nod) = sums(j,3)                                  
            end do                                                              
         end if                                                                 
      end do                                                                    
c                                                                               
c                       return the matrix to global upper triangular            
c                       vector form. (or full form unrolled form for            
c                       asymmetric assembly)                                    
c                                                                               
      if (.not. asymmetric_assembly) then                                       
      do j = 1, totdof                                                          
         do i = 1, j                                                            
            mat(cp(j)+i,bele) = tmat(i,j)                                       
         end do                                                                 
      end do                                                                    
      else                                                                      
            mat(:,bele) = reshape(tmat, (/totdof*totdof/))                      
      end if                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
