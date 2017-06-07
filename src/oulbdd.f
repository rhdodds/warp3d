c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oulbdd                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 10/30/2012                 *          
c     *                                                              *          
c     *     set the label type and corresponding                     *          
c     *     column labels for a node for displ/vel/acc/temp output   *          
c     *     if the not already been set, and checks the current      *          
c     *     node against the above if not.                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine oulbdd( dva, hedtyp, lbltyp, type, elem, doflbl )              
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
c                       element types 1-5: 3-d isoparametric                    
c                       solid elements.                                         
c                                                                               
 100  continue                                                                  
c                                                                               
c                       check to see if a label type exists. if not,            
c                       set one and its corresponding dof labels.               
c                                                                               
      if( lbltyp .eq. 0 ) then                                                  
         lbltyp = 1                                                             
         if( dva .eq. 1 ) then                                                  
            hedtyp = 'displacements       '                                     
         else if( dva .eq. 2 ) then                                             
            hedtyp = 'velocities          '                                     
         else if( dva .eq. 3 ) then                                             
            hedtyp = 'accelerations       '                                     
         else if( dva .eq. 4 ) then                                             
            hedtyp = 'reaction forces     '                                     
         else if( dva .eq. 5 ) then                                             
            hedtyp = 'temperatures        '                                     
         end if                                                                 
c                                                                               
         doflbl(1) = '   u    '                                                 
         doflbl(2) = '   v    '                                                 
         doflbl(3) = '   w    '                                                 
         if( dva .eq. 4 ) then                                                  
            doflbl(1) = 'Force X '                                              
            doflbl(2) = 'Force Y'                                               
            doflbl(3) = 'Force Z'                                               
         end if                                                                 
         if( dva .eq. 5 ) then                                                  
            doflbl(1) = '   T    '                                              
            doflbl(2) = 'T - Tref'                                              
            doflbl(3) = '        '                                              
         end if                                                                 
         go to 9999                                                             
      end if                                                                    
c                                                                               
c                       if a label type exists, check to see if the             
c                       new label type  matches. if                             
c                       not, print a warning and change the label               
c                       type and new labels.                                    
c                                                                               
      if( lbltyp .ne. 1 ) then                                                  
         call errmsg( 139, elem, dums, dumr, dumd )                             
         lbltyp = 1                                                             
         if( dva .eq. 1 ) then                                                  
            hedtyp = 'displacements       '                                     
         else if( dva .eq. 2 ) then                                             
            hedtyp = 'velocities          '                                     
         else if( dva .eq. 3 ) then                                             
            hedtyp = 'accelerations       '                                     
         else if( dva .eq. 4 ) then                                             
            hedtyp = 'reaction forces     '                                     
         else if( dva .eq. 5 ) then                                             
            hedtyp = 'temperatures        '                                     
         end if                                                                 
         doflbl(1) = '   u    '                                                 
         doflbl(2) = '   v    '                                                 
         doflbl(3) = '   w    '                                                 
         if( dva .eq. 4 ) then                                                  
            doflbl(1) = 'Force X '                                              
            doflbl(2) = 'Force Y'                                               
            doflbl(3) = 'Force Z'                                               
         end if                                                                 
         if( dva .eq. 5 ) then                                                  
            doflbl(1) = '   T    '                                              
            doflbl(2) = 'T - Tref'                                              
            doflbl(3) = '        '                                              
         end if                                                                 
      end if                                                                    
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
 500  continue                                                                  
      go to 9999                                                                
c                                                                               
 9999 return                                                                    
      end                                                                       
