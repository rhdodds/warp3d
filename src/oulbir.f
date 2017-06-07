c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulbir                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 05/05/89                   *          
c     *                                                              *          
c     *     this subroutine sets the label type and corresponding    *          
c     *     dof labels for a node for internal forces or residual    *          
c     *     loads output if the above has not already been set, and  *          
c     *     checks the current node against the above if not.        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine oulbir( ir, hedtyp, lbltyp, type, elem, doflbl )               
      implicit integer (a-z)                                                    
      character(len=8) :: doflbl(*)                                             
      character(len=20) :: hedtyp                                               
      real dumr                                                                 
      double precision                                                          
     &     dumd                                                                 
      character :: dums                                                         
c                                                                               
c                       branch on element type.                                 
c                                                                               
      go to (100,100,100,100,100), type                                         
c                                                                               
c                       element types 1-5: 3-d isoparametrics                   
c                                                                               
 100  continue                                                                  
c                                                                               
c                       check to see if a label type exists. if not,            
c                       set one and its corresponding dof labels.               
c                                                                               
      if(lbltyp.eq.0) then                                                      
         lbltyp= 1                                                              
c                                                                               
         if(ir.eq.1) then                                                       
            hedtyp= 'reaction forces     '                                      
         else if(ir.eq.2) then                                                  
            hedtyp= 'reaction forces     '                                      
         else                                                                   
            hedtyp= 'linear residual     '                                      
         end if                                                                 
c                                                                               
         doflbl(1)= '   u    '                                                  
         doflbl(2)= '   v    '                                                  
         doflbl(3)= '   w    '                                                  
         go to 9999                                                             
      end if                                                                    
c                                                                               
c                       if a label type exists, check to see if the             
c                       label type for the 3dqisop matches it. if               
c                       not, print a warning and change the label               
c                       type and dof labels.                                    
c                                                                               
      if(lbltyp.ne.1) then                                                      
         call errmsg(139,elem,dums,dumr,dumd)                                   
         lbltyp= 1                                                              
c                                                                               
         if(ir.eq.1) then                                                       
            hedtyp= 'reaction forces     '                                      
         else if(ir.eq.2) then                                                  
            hedtyp= 'reaction forces     '                                      
         else                                                                   
            hedtyp= 'linear residual     '                                      
         end if                                                                 
c                                                                               
         doflbl(1)= '   u    '                                                  
         doflbl(2)= '   v    '                                                  
         doflbl(3)= '   w    '                                                  
      end if                                                                    
c                                                                               
      go to 9999                                                                
c                                                                               
c                                                                               
 300  continue                                                                  
      go to 9999                                                                
c                                                                               
c                                                                               
 400  continue                                                                  
      go to 9999                                                                
c                                                                               
c                                                                               
 500  continue                                                                  
      go to 9999                                                                
c                                                                               
 9999 return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
