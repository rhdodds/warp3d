c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine license                      *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 01/17/95 by asg            *          
c     *                                                              *          
c     *        This subroutine prints out the license agreement      *          
c     *        for WARP3D.                                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine license                                                        
      use global_data ! old common.main
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
      write(out,*)                                                              
     &"================================================================         
     &===============",                                                         
     &"WARP3D Release License",                                                 
     &"================================================================         
     &===============",                                                         
     &"University of Illinois/NCSA",                                            
     &"Open Source License",                                                    
     &"",                                                                       
     &"Copyright (c) 2011 University of Illinois at Urbana-Champaign",          
     &"All rights reserved.",                                                   
     &"",                                                                       
     &"Developed by:",                                                          
     &"      Professor Robert H. Dodds, Jr. & Research Group",                  
     &"      University of Illinois at Urbana-Champaign",                       
     &"      http://cee.illinois.edu/faculty/robertdodds",                      
     &"",                                                                       
     &"Permission is hereby granted, free of charge, to any person obta         
     &ining a copy of",                                                         
     &"this software and associated documentation files (the ""Software         
     &""), to deal with",                                                       
     &"the Software without restriction, including without limitation t         
     &he rights to",                                                            
     &"use, copy, modify, merge, publish, distribute, sublicense, and/o         
     &r sell copies",                                                           
     &"of the Software, and to permit persons to whom the Software is f         
     &urnished to do",                                                          
     &"so, subject to the following conditions:",                               
     &"",                                                                       
     &"    * Redistributions of source code must retain the above copy          
     &right notice,",                                                           
     &"      this list of conditions and the following disclaimers.",           
     &"    * Redistributions in binary form must reproduce the above co         
     &pyright notice,",                                                         
     &"      this list of conditions and the following disclaimers in t         
     &he",                                                                      
     &"      documentation and/or other materials provided with the di          
     &stribution.",                                                             
     &"    * Neither the names of the WARP3D Team, University of Illin          
     &ois at",                                                                  
     &"      Urbana-Champaign, nor the names of its contributors may b          
     &e used to",                                                               
     &"      endorse or promote products derived from this Software wit         
     &hout specific",                                                           
     &"      prior written permission. ",                                       
     &"",                                                                       
     &"THE SOFTWARE IS PROVIDED ""AS IS"", WITHOUT WARRANTY OF ANY KIND,        
     &EXPRESS OR",                                                              
     &"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHAN          
     &TABILITY, FITNESS",                                                       
     &"FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL          
     &THE",                                                                     
     &"CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMA          
     &GES OR OTHER",                                                            
     &"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,          
     &ARISING FROM,",                                                           
     &"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER            
     &DEALINGS WITH THE",                                                       
     &"SOFTWARE.",                                                              
     &"",                                                                       
     &"===============================================================          
     &================",                                                        
     &"Copyrights and Licenses for Third Party Software Distributed wi          
     &th WARP3D:",                                                              
     &"===============================================================          
     &================",                                                        
     &"The WARP3D software contains code written by third parties.  S           
     &uch software will",                                                       
     &"have its own individual license file in the directory in which           
     &it appears.",                                                             
     &"This file will describe the copyrights, license, and restrictio          
     &ns which apply",                                                          
     &"to that code.",                                                          
     &"",                                                                       
     &"The disclaimer of warranty in the University of Illinois Open            
     &Source License",                                                          
     &"applies to all code in the WARP3D Distribution, and nothing in           
     &any of the",                                                              
     &"other licenses gives permission to use the names of the WARP3            
     &D Team or the",                                                           
     &"University of Illinois to endorse or promote products derived            
     &from this",                                                               
     &"Software.",                                                              
     &"",                                                                       
     &"The following pieces of software have additional or alternate            
     &copyrights,",                                                             
     &"licenses, and/or restrictions:"                                          
                                                                                
      write (out,*)                                                             
     &"Program                 Directory"                                       
      write (out,*)                                                             
     &"-------                 ---------"                                       
      write (out,*)                                                             
     &"hypre                   linux_packages/source/hypre-2.7.0b"              
      write (out,*)                                                             
     &"metis                   linux_packages/source/metis-4.0"                 
      write (out,*)                                                             
     &"Intel MKL libraries     linux_packages/lib"                              
      write (out,*)                                                             
     &""                                                                        
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
                                                                                
