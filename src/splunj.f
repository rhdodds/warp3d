      subroutine splunj                                                         
      implicit integer (a-z)                                                    
c                                                                               
c                     this is a dummy subroutine that does nothing.             
c                     we use it most often in scanning of input                 
c                     to skip over optional words in commands.                  
c                                                                               
c                     suppose the word 'for' is optional in some                
c                     command. then we would normally use:                      
c                                                                               
c                       if ( matchs('for',3) ) continue                         
c                                                                               
c                     but some compilers will just delete the                   
c                     matchs execution since it thinks the                      
c                     outcome is immaterial. we trick the compilers             
c                     with:                                                     
c                                                                               
c                       if ( matchs('for',3) ) call splunj                      
c                                                                               
      return                                                                    
      end                                                                       
