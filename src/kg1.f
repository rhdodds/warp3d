c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine kgstiff                      *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified :  3/15/2017 rhd             *          
c     *                                                              *          
c     *     for this element block, this integration point,          *          
c     *     compute gemetric stiffness and add to [Ke]               *          
c     *     handles [ke] in symmetric and asymmetric forms           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine kgstiff( span, cp, icp, gama, nxi, neta, nzeta, nnode,         
     &                sig, dj, w, ek_full, ek_symm, vol, bbar,                  
     &                totdof )                                                  
      use main_data, only: asymmetric_assembly                                  
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, cp(*), icp(mxutsz,*), nnode, nsz, totdof                 
      double precision ::                                                       
     &    gama(mxvl,ndim,*), nxi(*), neta(*), nzeta(*), sig(mxvl,*),            
     &    dj(*), ek_full(span,*), ek_symm(span,*), w, vol(mxvl,8,*)             
      logical :: bbar                                                           
c                                                                               
c                       locally arrays - on stack                               
c                                                                               
      integer :: j, i, cp1, cp2, cp3, r, c, k, enode,                           
     &           r1, c1, k1, r2, c2, k2, r3, c3, k3                             
      logical :: symmetric_assembly                                             
      double precision                                                          
     &    gtg(mxvl,mxnusz), gxi(mxvl,mxndel),                                   
     &    geta(mxvl,mxndel), gzeta(mxvl,mxndel)                                 
c                                                                               
c              the geometric stiffness, trans([G[) [M] [G] is                   
c              symmetric. calculate building blocks then use separate           
c              paths to update [Ke] based on its symmetry/asymmetry.            
c                                                                               
      symmetric_assembly = .not. asymmetric_assembly                            
c                                                                               
      if( bbar ) then                                                           
        do j = 1, nnode                                                         
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
          do i = 1, span                                                        
           gxi(i,j)   = vol(i,j,1)                                              
           geta(i,j)  = vol(i,j,2)                                              
           gzeta(i,j) = vol(i,j,3)                                              
          end do                                                                
        end do                                                                  
      else                                                                      
        do j = 1, nnode                                                         
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
          do i = 1, span                                                        
            gxi(i,j)  =  gama(i,1,1)*nxi(j)+gama(i,1,2)*neta(j)+                
     &                   gama(i,1,3)*nzeta(j)                                   
            geta(i,j) =  gama(i,2,1)*nxi(j)+gama(i,2,2)*neta(j)+                
     &                   gama(i,2,3)*nzeta(j)                                   
            gzeta(i,j) = gama(i,3,1)*nxi(j)+gama(i,3,2)*neta(j)+                
     &                   gama(i,3,3)*nzeta(j)                                   
          end do                                                                
        end do                                                                  
      end if                                                                    
c                                                                               
c              the gemetric stiffness has an unusual, repeated                  
c              structure with many of the smae submatrices added                
c              multiple times into [Ke]                                         
c                                                                               
c              get geometric stiffness for one (x,y,z) direction.               
c              others are identical.                                            
c              Example. 8-node hex. a row of gtg contains the upper-            
c              triangle of the 8 x 8 geometric stiffnessfor a dof. the          
c              row is the same for other two dof. Thus, cols j=1->8             
c              have row 1 of gtg, cols 9-15 have row 2, etc of the              
c              upper-triangle.                                                  
c                                                                               
      do j = 1, cp(nnode)+nnode                                                 
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
            gtg(i,j)= (gxi(i,icp(j,1))*gxi(i,icp(j,2))*sig(i,1)+                
     &                 geta(i,icp(j,1))*geta(i,icp(j,2))*sig(i,2)+              
     &                 gzeta(i,icp(j,1))*gzeta(i,icp(j,2))*sig(i,3)+            
     &                (gxi(i,icp(j,1))*geta(i,icp(j,2))+                        
     &                 gxi(i,icp(j,2))*geta(i,icp(j,1)))*sig(i,4)+              
     &                (geta(i,icp(j,1))*gzeta(i,icp(j,2))+                      
     &                 geta(i,icp(j,2))*gzeta(i,icp(j,1)))*sig(i,5)+            
     &                (gxi(i,icp(j,1))*gzeta(i,icp(j,2))+                       
     &                 gxi(i,icp(j,2))*gzeta(i,icp(j,1)))*sig(i,6))*            
     &                 dj(i)*w                                                  
         end do                                                                 
      end do                                                                    
c                                                                               
c              add to the symmetric [Ke]                                        
c                                                                               
      if( symmetric_assembly) then                                              
        do enode = 1, nnode ! all element nodes                                 
         cp1 = cp(enode)                                                        
         cp2 = cp(nnode+enode)+nnode                                            
         cp3 = cp(2*nnode+enode)+2*nnode                                        
c                                                                               
         do j = 1, enode                                                        
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
            do  i = 1, span                                                     
               ek_symm(i,cp1+j) = ek_symm(i,cp1+j) + gtg(i,cp1+j)               
               ek_symm(i,cp2+j) = ek_symm(i,cp2+j) + gtg(i,cp1+j)               
               ek_symm(i,cp3+j) = ek_symm(i,cp3+j) + gtg(i,cp1+j)               
            end do                                                              
         end do ! j over element nodes                                          
       end do  ! over element nodes                                             
       return                                                                   
      end if                                                                    
c                                                                               
c               add to the asymmetric assembly. this is more difficult          
c               coding since only the upper-triangle of gtg is                  
c               computed above. so we have loops to add upper-triangle          
c               gtg terms again into lower trangle of the asymmetric            
c               [Ke]. this algortihm keeps ifs out of inner loops and           
c               has stride 1 access all arrays.                                 
c                                                                               
      do enode = 1, nnode                                                       
c                                                                               
         cp1 = cp(enode)                                                        
         c1  = enode                                                            
         c2  = nnode + enode                                                    
         c3  = 2*nnode + enode                                                  
c                                                                               
         do j = 1, enode                                                        
c                                                                               
           r1 = j                                                               
           k1 = (c1-1)*totdof + r1                                              
c                                                                               
           r2 = nnode + j                                                       
           k2 = (c2-1)*totdof + r2                                              
c                                                                               
           r3 = 2*nnode + j                                                     
           k3 = (c3-1)*totdof + r3                                              
c                                                                               
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
           do  i = 1, span! col enode, row j                                    
             ek_full(i,k1) = ek_full(i,k1) + gtg(i,cp1+j)                       
             ek_full(i,k2) = ek_full(i,k2) + gtg(i,cp1+j)                       
             ek_full(i,k3) = ek_full(i,k3) + gtg(i,cp1+j)                       
           end do                                                               
c                                                                               
           if( r1 .ne. c1 ) then  ! cleanup loop 1                              
             k = (r1-1)*totdof + c1                                             
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
             do i = 1, span                                                     
               ek_full(i,k) = ek_full(i,k) + gtg(i,cp1+j)                       
             end do                                                             
           end if                                                               
c                                                                               
           if( r2 .ne. c2 ) then  ! cleanup loop 2                              
             k = (r2-1)*totdof + c2                                             
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
             do i = 1, span                                                     
               ek_full(i,k) = ek_full(i,k) + gtg(i,cp1+j)                       
             end do                                                             
           end if                                                               
c                                                                               
           if( r3 .ne. c3 ) then  ! cleanup loop 3                              
              k = (r3-1)*totdof + c3                                            
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
              do i = 1, span                                                    
                ek_full(i,k) = ek_full(i,k) + gtg(i,cp1+j)                      
              end do                                                            
           end if                                                               
c                                                                               
         end do ! over j                                                        
      end do ! over enode                                                       
c                                                                               
      return                                                                    
      end                                                                       
