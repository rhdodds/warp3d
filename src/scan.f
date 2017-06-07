c---------------------------------------------------------------------          
c                                                                               
c    Scan --  a Problem Oriented Language interpreter.                          
c             Written by: Robert H. Dodds                                       
c                                                                               
c    Copyright (C) 1994  <Robert H. Dodds>                                      
c                                                                               
c    This program is free software; you can redistribute it and/or modify       
c    it under the terms of the GNU General Public License as published by       
c    the Free Software Foundation; either version 1, or (at your option)        
c    any later version.                                                         
c                                                                               
c    This program is distributed in the hope that it will be useful,            
c    but WITHOUT ANY WARRANTY; without even the implied warranty of             
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
c    GNU General Public License for more details.                               
c                                                                               
c    You should have received a copy of the GNU General Public License          
c    along with this program; if not, write to the Free Software                
c    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                  
c                                                                               
c    questions may be addressed to :                                            
c      Prof. Robert Dodds -- email r-dodds@uiuc.edu                             
c                                                                               
c-----------------------------------------------------------------              
c                                                                               
c                                                                               
C *******************************************************************           
c *******************************************************************           
c *******************************************************************           
c **                                                               **           
c **                                                               **           
c **                          s c a n                              **           
c **                          =======                              **           
c **                                                               **           
c **           a free form input fortran subprogram system         **           
c **                    Unix (ascii) version                       **           
c **                upper and lower case version                   **           
c **                                                               **           
c **              civil engineering systems laboratory             **           
c **                department of civil engineering                **           
c **           university of illinois at urbana-champaign          **           
c **                                                               **           
c **                                                               **           
c *******************************************************************           
c *******************************************************************           
c *******************************************************************           
c **********************************************************************        
c *                                                                    *        
c * scan                                                               *        
c *                                                                    *        
c **********************************************************************        
      subroutine scan                                                           
c                                                                               
c                                                                               
c         scan                                                                  
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      common /scantb/ nseptb, nclass, itab(256), iclass(256)                    
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      dimension istate(110)                                                     
      integer istate, state, exp, ocol, scanmn                                  
      double precision dvalue,flt,f,ten,one,zero,tenth                          
      data ocol/0/, ndot/0/                                                     
      data isig/1/,nsig/1/,int/0/,f/.1/,exp/0/                                  
      data iptlc1/1/, iptlc2/6/                                                 
      data ten/10.0/,one/1.0/,zero/0.0/,tenth/0.1/                              
c                                                                               
c                              d   +   -  e    a   q   .   s   b   $            
c                              i               l   u       e   l                
c                              g               p   o       p   a                
c                              i               h   t           n                
c                              t               a   e           k                
c                                                                               
c        0     get            18, 17, 16, 26, 26, 30, 32, 29,  6,  7,           
c        1     getstring      18, 17, 16, 26, 26, 13, 32, 29,  6,  7,           
c        2     +-             18, 11, 11, 11, 11, 11, 35, 11, 11, 11,           
c        3     integer        19, 11, 11, 28, 28, 11, 20, 11, 11, 11,           
c        4     real           21, 10, 10, 22, 28, 10, 28, 10, 10, 10,           
c        5     start exp      25, 24, 23, 28, 28,  9, 28,  9,  9,  9,           
c        6     exp            25,  9,  9, 28, 28,  9, 28,  9,  9,  9,           
c        7     label           4,  8,  8,  4,  4,  8, 27,  8,  8,  8,           
c        8     name            4, 14, 14,  4,  4, 14, 27, 14, 14, 14,           
c        9     string          4,  4,  4,  4,  4, 12,  4,  4,  4, 33,           
c        a     anytext         4,  8,  8,  4,  4,  8,  4,  8,  8,  8 /          
c                                                                               
c                                                                               
      data istate / 18, 17, 16, 26, 26, 30, 32, 29,  6,  7,                     
     1              18, 17, 16, 26, 26, 13, 32, 29,  6,  7,                     
     2              18, 11, 11, 11, 11, 11, 35, 11, 11, 11,                     
     3              19, 11, 11, 28, 28, 11, 20, 11, 11, 11,                     
     4              21, 10, 10, 22, 28, 10, 28, 10, 10, 10,                     
     5              25, 24, 23, 28, 28,  9, 28,  9,  9,  9,                     
     6              25,  9,  9, 28, 28,  9, 28,  9,  9,  9,                     
     7               4,  8,  8,  4,  4,  8, 27,  8,  8,  8,                     
     8               4, 14, 14,  4,  4, 14, 27, 14, 14, 14,                     
     9               4,  4,  4,  4,  4, 12,  4,  4,  4, 33,                     
     a               4,  8,  8,  4,  4,  8,  4,  8,  8,  8  /                   
c                                                                               
      dvalue = zero                                                             
      flt = zero                                                                
c                                                                               
      nent=nent+1                                                               
      nchar = 0                                                                 
      state = pstate                                                            
      nblank = 0                                                                
      entity(1) = blank                                                         
      jstart(nent)=0                                                            
      ocol=col                                                                  
      lim=limit                                                                 
c                                                                               
c            go to action                                                       
c                                                                               
    1 continue                                                                  
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,            
     1        170,180,190,200,210,220,230,240,250,260,270,280,                  
     2        290,300,310,320,330,340,350), jump                                
c                                                                               
c            get a line                                                         
c                                                                               
    5 continue                                                                  
      doread = .true.                                                           
   10 continue                                                                  
      if(doread)call scanrd                                                     
      col = 0                                                                   
      nent = 1                                                                  
      jstart(1) = 0                                                             
      jump=2                                                                    
      lim=80                                                                    
      if(mode.eq.10)go to 12                                                    
      if(mode.eq.11)go to 14                                                    
      go to 20                                                                  
c                                                                               
c            end of file                                                        
c                                                                               
   12 continue                                                                  
      filpt = filpt-1                                                           
      if(filpt.gt.0)inunit = files(filpt)                                       
      if(filpt.eq.0)mode=-10                                                    
      if(leof)go to 13                                                          
      if(filpt.gt.0)go to 5                                                     
      call scanm(1)                                                             
      call die_gracefully                                                       
      stop                                                                      
   13 jump = 1                                                                  
      doread = .true.                                                           
      go to 150                                                                 
c                                                                               
c             point                                                             
c                                                                               
   14 continue                                                                  
      ival(1) = 1000*idigit(iptlc1+0)+100*idigit(iptlc1+1)+                     
     1          10*idigit(iptlc1+2)+    idigit(iptlc1+3)                        
      ival(2) = 1000*idigit(iptlc2+0)+100*idigit(iptlc2+1)+                     
     1          10*idigit(iptlc2+2)+    idigit(iptlc2+3)                        
      if(.not.menu)go to 16                                                     
      if(scanmn(ival(1),ival(2)).lt.0)go to 5                                   
      jump = 3                                                                  
      icolmn = 1                                                                
      col = 10                                                                  
      jstart(1) = 1                                                             
      go to 150                                                                 
   16 jstart(1) = 1                                                             
      nchar = 10                                                                
      col = 10                                                                  
      go to 80                                                                  
c                                                                               
c            next character                                                     
c                                                                               
   20 continue                                                                  
      col = col+1                                                               
c                                                                               
c            get character                                                      
c                                                                               
   30 jump = istate(10*state+ibuff(col))                                        
      go to 1                                                                   
c                                                                               
c            set jump to next character                                         
c                                                                               
   40 jump = 2                                                                  
c                                                                               
c            pack character                                                     
c                                                                               
   50 nchar = nchar + 1                                                         
      go to 1                                                                   
c                                                                               
c            count blanks                                                       
c                                                                               
   60 continue                                                                  
      nblank = nblank + 1                                                       
      if(nblank.le.lim)go to 20                                                 
c                                                                               
c        check if end of string                                                 
c                                                                               
      if(pstate.eq.1) go to 130                                                 
c                                                                               
c            end line                                                           
c                                                                               
c   70 if(.not.eol)go to 75                                                     
c      jump = 3                                                                 
c      if(autord)jump = 1                                                       
c      if(autord)doread = .true.                                                
c      entity(1)=blank                                                          
c      mode = 9                                                                 
c      nchar=0                                                                  
c      nwd = 0                                                                  
c                                                                               
c                                                                               
c                       this stuff is sort of a                                 
c                       backsp to handle the case of repeated                   
c                       calls to scan after an eol has been                     
c                       reached. looking for a class 8  is for                  
c                       separators. when you leave scan on a                    
c                       separator you are on that column. when                  
c                       you are looking at anything else, you                   
c                       are looking at the first blank after the                
c                       item. setting jump to 2 causes you to                   
c                       go to the  column after the separator                   
c                       on re-entry to scan.                                    
c                                                                               
c                                                                               
c      col=max0(ocol,1)                                                         
c      jstart(nent)=col                                                         
c      nent = nent - 1                                                          
c      if ( ibuff( col ) .eq. 8 ) then                                          
c         doread = .true.                                                       
c         jump = 1                                                              
c         mode = 0                                                              
c         goto 1                                                                
c      endif                                                                    
c      if ( nent .ne. 0 ) go to 150                                             
c      col  = 2                                                                 
c      jump = 2                                                                 
c      go to 150                                                                
c                                                                               
c            end line -- if autoread or auto continue is on                     
c                        set up to go to the next line without letting          
c                        the grammar know it has happened                       
c                                                                               
c                       this stuff is sort of a                                 
c                       backsp to handle the case of repeated                   
c                       calls to scan after an eol has been                     
c                       reached. looking for a class 8  is for                  
c                       separators. when you leave scan on a                    
c                       separator you are on that column. when                  
c                       you are looking at anything else, you                   
c                       are looking at the first blank after the                
c                       item. setting jump to 2 causes you to                   
c                       go to the  column after the separator                   
c                       on re-entry to scan. the nonsense about                 
c                       nent =3 etc... is to prevent an overflow on             
c                       repeated calls to scan after the eol is reached.        
c                                                                               
   70 continue                                                                  
      if(eol)then                                                               
         if(.not.isct) then                                                     
            jump = 3                                                            
            if(autord)jump = 1                                                  
            if(autord)doread = .true.                                           
            entity(1)=blank                                                     
            mode = 9                                                            
            nchar=0                                                             
            nwd = 0                                                             
            col=max0(ocol,1)                                                    
            jstart(nent)=col                                                    
            if ( nent .ge. 4 ) then                                             
               if(jstart(nent-3 ) .eq. jstart(nent) ) nent = nent - 1           
            endif                                                               
            if ( ibuff( col ) .eq. 8 ) jump = 2                                 
            if ( nent .eq. 0 ) then                                             
               col  = 2                                                         
               jump = 2                                                         
            endif                                                               
            return                                                              
         endif                                                                  
      endif                                                                     
c                                                                               
   75 nchar = 0                                                                 
      nwd = 0                                                                   
      nblank = 0                                                                
      state = pstate                                                            
      go to 5                                                                   
c                                                                               
c            pack blank exit                                                    
c                                                                               
   80 jump = 3                                                                  
   85 icolmn = 0                                                                
      if(nent.ne.0)icolmn = jstart(nent)                                        
      call scanpk                                                               
      go to 150                                                                 
c                                                                               
c            end exponent                                                       
c                                                                               
   90 continue                                                                  
      flt = flt * ten ** exp                                                    
      exp = 0                                                                   
      nsig = 1                                                                  
c                                                                               
c            end real                                                           
c                                                                               
  100 continue                                                                  
      f = tenth                                                                 
      dvalue = flt                                                              
      value = flt                                                               
      isig = 1                                                                  
      int = 0                                                                   
      flt = zero                                                                
      goto 80                                                                   
c                                                                               
c            end integer                                                        
c                                                                               
  110 continue                                                                  
      ivalue = int                                                              
      isig = 1                                                                  
      int = 0                                                                   
      go to 80                                                                  
c                                                                               
c            end string                                                         
c                                                                               
  120 skip = col                                                                
  340 jump = 2                                                                  
      go to 85                                                                  
c                                                                               
c            end get string                                                     
c                                                                               
  130 mode = 8                                                                  
      pstate = 0                                                                
      jump = 2                                                                  
      go to 150                                                                 
c                                                                               
c            end name                                                           
c                                                                               
  140 nchar = ndot + 1                                                          
      go to 80                                                                  
c                                                                               
c            return                                                             
c                                                                               
  150 continue                                                                  
      return                                                                    
c                                                                               
c            minus found                                                        
c                                                                               
  160 isig = -1                                                                 
c                                                                               
c            plus found                                                         
c                                                                               
  170 state = 2                                                                 
      mode = 6                                                                  
      int = itab(jbuff(col)+1)                                                  
      ivalue = int                                                              
      go to 310                                                                 
c                                                                               
c            start integer                                                      
c                                                                               
  180 state = 3                                                                 
      int = 0                                                                   
      mode = 1                                                                  
      if(jstart(nent).eq.0)jstart(nent)=col                                     
c                                                                               
c            continue integer                                                   
c                                                                               
  190 int = int*10 + idigit(col) * isig                                         
      go to 40                                                                  
c                                                                               
c            change to real                                                     
c                                                                               
  200 state = 4                                                                 
      mode = 2                                                                  
      flt = real (int)                                                          
      go to 40                                                                  
c                                                                               
c            continue real                                                      
c                                                                               
  210 continue                                                                  
      flt = flt + dble(idigit(col)*isig)*f                                      
      f = f * tenth                                                             
      go to 40                                                                  
c                                                                               
c            start exponent                                                     
c                                                                               
  220 state = 5                                                                 
      go to 40                                                                  
c                                                                               
c            minus exp                                                          
c                                                                               
  230 nsig = -1                                                                 
c                                                                               
c            plus exp                                                           
c                                                                               
  240 state = 6                                                                 
      go to 40                                                                  
c                                                                               
c            continue exp                                                       
c                                                                               
  250 exp = exp * 10 + idigit(col)*nsig                                         
      state = 6                                                                 
      go to 40                                                                  
c                                                                               
c            start label                                                        
c                                                                               
  260 state = 7                                                                 
      ndot = 0                                                                  
      mode = 3                                                                  
      go to 310                                                                 
c                                                                               
c            start name                                                         
c                                                                               
  270 state = 8                                                                 
      ndot = ndot + 1                                                           
      nchar = ndot                                                              
      mode = 4                                                                  
      jstart(nent+ndot) = col+1                                                 
      go to 20                                                                  
c                                                                               
c            any text                                                           
c                                                                               
  280 mode = 5                                                                  
      state = 10                                                                
      isig = 1                                                                  
      int = 0                                                                   
      exp = 0                                                                   
      f = tenth                                                                 
      nsig = 1                                                                  
      go to 40                                                                  
c                                                                               
c            separater                                                          
c                                                                               
  290 int = itab(jbuff(col)+1)                                                  
      ivalue = int                                                              
      mode = 6                                                                  
      jump = 34                                                                 
      jstart(nent) = col                                                        
      isig = 1                                                                  
      go to 50                                                                  
c                                                                               
c            start string                                                       
c                                                                               
  300 state = 9                                                                 
      mode = 7                                                                  
  305 jstart(nent) = col                                                        
      go to 20                                                                  
c                                                                               
c            set beginning col number                                           
c                                                                               
  310 jstart(nent) = col                                                        
      if(mode.eq.6.and.signed)go to 290                                         
      go to 40                                                                  
c                                                                               
c            start real                                                         
c                                                                               
  320 jstart(nent) = col                                                        
      int = 0                                                                   
      flt = zero                                                                
      go to 200                                                                 
c                                                                               
c     end line and end string                                                   
c                                                                               
  330 skip = col                                                                
      go to 80                                                                  
c                                                                               
c      start signed real                                                        
c                                                                               
  350 int = 0                                                                   
      flt = zero                                                                
      go to 200                                                                 
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * rdline                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine rdline                                                         
c                                                                               
c         read a card                                                           
c                                                                               
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      jump = 1                                                                  
      call scanrd                                                               
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scanmn                                                             *        
c *                                                                    *        
c **********************************************************************        
      integer function scanmn( point1, point2 )                                 
c                                                                               
c         menuing routine interface                                             
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      integer point1, point2                                                    
      double precision dvalue                                                   
c                                                                               
c         dummy routine                                                         
c                                                                               
      scanmn = 0                                                                
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scanm                                                              *        
c *                                                                    *        
c **********************************************************************        
      subroutine scanm( errno )                                                 
c                                                                               
c         issue error messages                                                  
c                                                                               
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      integer errno                                                             
c                                                                               
c         send error message and return                                         
c                                                                               
      go to (  10,  20,  30, 40 ), errno                                        
   10 write(iout,1001)                                                          
      go to 9999                                                                
   20 write(iout,1002)                                                          
      go to 9999                                                                
   30 write(iout,1003)                                                          
      go to 9999                                                                
   40 write(iout,1004)                                                          
      call die_gracefully                                                       
      stop                                                                      
 9999 return                                                                    
c                                                                               
 1001 format(40h0 .....end of file - program terminated  /)                     
 1002 format(40h0 .....input error - program terminated  /)                     
 1003 format(40h0 .....file list overflow                /)                     
 1004 format(/,                                                                 
     & '>> ERROR:  text string mismatch in scanms....',/,                       
     & '           lena, lenb, nchar: ',3i6,/,                                  
     & '           program terminated.....',//)                                 
c                                                                               
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * setin                                                              *        
c *                                                                    *        
c **********************************************************************        
      subroutine setin(in)                                                      
c                                                                               
c         set scan input unit                                                   
c                                                                               
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      inunit = in                                                               
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * setout                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine setout(out)                                                    
c                                                                               
c         set scan echo output unit                                             
c                                                                               
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      integer out                                                               
      iout = out                                                                
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * setrem                                                          *           
c *                                                                 *           
c *******************************************************************           
      subroutine setrem( iunit, ounit )                                         
c                                                                               
c          set the remote unit for prompting                                    
c                                                                               
      integer ounit                                                             
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      inremo = iunit                                                            
      iotrem = ounit                                                            
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * setfil                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine setfil(filist,nfiles)                                          
c                                                                               
c         stack a list of input files                                           
c                                                                               
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      dimension filist(10)                                                      
      integer filist                                                            
      nfil = nfiles                                                             
      if(filpt+nfil.le.fillim) go to 10                                         
      nfil=fillim-filpt                                                         
      call scanm(3)                                                             
   10 j=filpt+nfil+1                                                            
      if(nfil.eq.0)go to 30                                                     
      do 20 k=1,nfil                                                            
        files(j-k)=filist(k)                                                    
   20 continue                                                                  
      filpt=filpt+nfil                                                          
      inunit=files(filpt)                                                       
   30 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scinit                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine scinit(nblank, reclen, endchr, promsw, echosw, comsw,          
     1                  atrdsw, eolsw, eofsw, menusw, ptsw, signsw )            
c                                                                               
c         initialize                                                            
c                                                                               
c                                                                               
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      logical promsw, echosw, comsw, atrdsw, eolsw, eofsw, ptsw,                
     1        menusw, signsw                                                    
      integer reclen                                                            
      init = .true.                                                             
      echar = endchr                                                            
      echo = echosw                                                             
      leof = eofsw                                                              
      eol = eolsw                                                               
      menu = menusw                                                             
      autord = atrdsw                                                           
      commnt = .not.comsw                                                       
      promt = promsw                                                            
      point = ptsw                                                              
      signed = signsw                                                           
      autoct = .true.                                                           
      limit = min0(recsiz,nblank)                                               
      mark = min0(recsiz,reclen)                                                
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * addlab                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine addlab(label)                                                  
c                                                                               
c         add label for scan echo                                               
c                                                                               
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      ilabel = label                                                            
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * backsp                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine backsp(nument)                                                 
c                                                                               
c         backspace the scanner                                                 
c                                                                               
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      if(nent.lt.nument) go to 10                                               
      nent = nent - nument                                                      
      col = jstart(nent + 1)                                                    
      jump = 3                                                                  
      go to 20                                                                  
   10 col = 0                                                                   
      jump = 2                                                                  
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * getstr                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine getstr                                                         
c                                                                               
c         parse a string                                                        
c                                                                               
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      col = jstart(nent)                                                        
      pstate = 1                                                                
      jump = 2                                                                  
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * reset                                                              *        
c *                                                                    *        
c **********************************************************************        
      subroutine reset                                                          
c                                                                               
c         start line over                                                       
c                                                                               
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      col = 0                                                                   
      jump = 2                                                                  
      nent = 0                                                                  
      entity(1) = blank                                                         
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * skpstr                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine skpstr                                                         
c                                                                               
c         skip to end of string                                                 
c                                                                               
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      col = skip                                                                
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * septab                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine septab(newtab)                                                 
c                                                                               
c         change seperator table                                                
c                                                                               
      dimension newtab(1)                                                       
      common /scantb/ nseptb, nclass, itab(256), iclass(256)                    
      do i = 1, nseptb                                                          
        itab(i) = newtab(i)                                                     
      end do                                                                    
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * doscan                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function doscan(dummy)                                            
c                                                                               
c          advance the scanner now                                              
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      if(next)call scan                                                         
      next = .false.                                                            
      doscan = .true.                                                           
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * endcrd                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function endcrd(dummy)                                            
c                                                                               
c         is the current scan entity an end of line                             
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modeol/9/                                                            
c                                                                               
      if( mode .eq. modeol ) then                                               
        endcrd = .true.                                                         
        next = .false.                                                          
        return                                                                  
      end if                                                                    
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 endcrd = .false.                                                          
      if(mode.ne.modeol)go to 20                                                
      endcrd = .true.                                                           
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * endfil                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function endfil(last)                                             
c                                                                               
c        is the current scan entity an end of file                              
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      logical last                                                              
      data modeof/10/                                                           
      if(.not.next)go to 10                                                     
      if ( mode .ne. iabs ( modeof) )call scan                                  
      next = .false.                                                            
   10 endfil = .false.                                                          
      last = .false.                                                            
      if(mode.ne.iabs(modeof))go to 20                                          
      endfil = .true.                                                           
      if(mode.eq.-modeof)last = .true.                                          
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * entit                                                              *        
c *                                                                    *        
c **********************************************************************        
      subroutine entit(text,nc,nw)                                              
c                                                                               
c          return the current characters from entity                            
c                                                                               
      dimension text(1)                                                         
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      nc = nchar                                                                
      nw = nwd                                                                  
      do i = 1, nw                                                              
       text(i) = entity(i)                                                      
      end do                                                                    
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * entits                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine entits(cstr,nc)                                                
c                                                                               
c         return the contents of the scanner in a string                        
c                                                                               
      character(len=*) :: cstr                                                  
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      double precision dvalue                                                   
      character(len=80) :: string                                               
      equivalence ( string,entity )                                             
c                                                                               
c                       get the hollerith from scan and convert                 
c                       to a string.  simply use equivalence of hollerith       
c                       and characters.                                         
c                                                                               
      nc = min0 (len(cstr),nchar)                                               
      do i = 1,nc                                                               
      cstr(i:i) = string(i:i)                                                   
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * integr                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function integr(integ)                                            
c                                                                               
c         is the current scan entity an integer                                 
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modint /1/                                                           
c                                                                               
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 integr = .false.                                                          
      if(mode.ne.modint)go to 20                                                
      integ = ivalue                                                            
      integr = .true.                                                           
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * label                                                              *        
c *                                                                    *        
c **********************************************************************        
      logical function label(dummy)                                             
c                                                                               
c         is the current scan entity a label                                    
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modlab/3/                                                            
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 label = .false.                                                           
      if(mode.ne.modlab)go to 20                                                
      label = .true.                                                            
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * match                                                              *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      logical function match(string,nc)                                         
c                                                                               
c         match a hollerith array against the current scan entity               
c                                                                               
      integer string(*)                                                         
      logical scanmc                                                            
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
c                                                                               
      match = .false.                                                           
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 if(nchar.lt.nc)go to 20                                                   
      if(.not.scanmc(string(1),entity(1),nc))go to 20                           
      match = .true.                                                            
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * match_exact                                                        *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      logical function match_exact( string, nc )                                
c                                                                               
c         exact match a character string against the current scan               
c         entity. input entity must have exactly the same                       
c         number of characters as the input string cstr.                        
c                                                                               
      integer string(*)                                                         
      logical scanmc                                                            
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
c                                                                               
      match_exact = .false.                                                     
      if( .not. next ) go to 10                                                 
      call scan                                                                 
      next = .false.                                                            
   10 continue                                                                  
      if( nchar .ne. nc ) go to 20                                              
      if(.not.scanmc(string(1),entity(1),nc))go to 20                           
      match_exact = .true.                                                      
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * matchs                                                             *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      logical function matchs(cstr,nc)                                          
c                                                                               
c         match a character array against the current scan entity               
c                                                                               
      character(len=*) :: cstr, text*80                                         
      integer  string(20)                                                       
      logical match                                                             
      equivalence ( string, text )                                              
c                                                                               
c                                                                               
c                       convert the string to hollerith                         
c                                                                               
      text(1:80) = ' '                                                          
      text(1:nc) = cstr(1:nc)                                                   
c                                                                               
c                       call the hollerith match routine to do                  
c                       the work                                                
c                                                                               
      matchs = match( string, nc )                                              
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * matchs_exact                                                       *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      logical function matchs_exact( cstr )                                     
c                                                                               
c         exact match a character array against the current                     
c         scan entity                                                           
c                                                                               
      character(len=*) :: cstr, text*80                                         
      integer  string(20)                                                       
      logical match_exact                                                       
      equivalence ( string, text )                                              
c                                                                               
c                                                                               
c                       convert the string to hollerith                         
c                                                                               
      text(1:80) = ' '                                                          
      nc = len(cstr)                                                            
      do 10 i = 1, nc                                                           
      text(i:i) = cstr(i:i)                                                     
   10 continue                                                                  
c                                                                               
c                       call the hollerith match routine to do                  
c                       the work                                                
c                                                                               
      matchs_exact = match_exact ( string, nc )                                 
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * name                                                               *        
c *                                                                    *        
c **********************************************************************        
      logical function name(dummy)                                              
c                                                                               
c         is the current scan entity a name                                     
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modnam/4/                                                            
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 name = .false.                                                            
      if(mode.ne.modnam)go to 20                                                
      name = .true.                                                             
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * noscan                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function noscan(dummy)                                            
c                                                                               
c         do not advance the scanner                                            
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      next = .false.                                                            
      noscan = .true.                                                           
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * numi                                                               *        
c *                                                                    *        
c **********************************************************************        
      logical function numi(integr)                                             
c                                                                               
c         is the current scan entity a number, convert to integer               
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modrel/2/                                                            
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 numi = .false.                                                            
      if (mode .gt. modrel)  go to 20                                           
      if (mode .eq. modrel ) integr = value                                     
      if (mode .ne. modrel ) integr = ivalue                                    
      numi = .true.                                                             
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * numd                                                               *        
c *                                                                    *        
c **********************************************************************        
      logical function numd(double)                                             
c                                                                               
c         is the current scan entity a number, convert to double                
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue, double                                           
      data modrel/2/, modint/1/                                                 
c                                                                               
      iout = 6                                                                  
      if(.not.next)go to 10                                                     
      call scan                                                                 
   10 next = .false.                                                            
      numd = .false.                                                            
      if (mode .gt. modrel) go to 20                                            
      if (mode .eq. modint ) double = real(ivalue)                              
      if (mode .ne. modint ) double = dvalue                                    
      numd = .true.                                                             
      next = .true.                                                             
   20 continue                                                                  
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * numr                                                               *        
c *                                                                    *        
c **********************************************************************        
      logical function numr(real)                                               
c                                                                               
c         is the current scan entity a number, convert to real                  
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modrel/2/, modint/1/                                                 
c                                                                               
      if(.not.next)go to 10                                                     
      call scan                                                                 
   10 next = .false.                                                            
      numr = .false.                                                            
      if (mode .gt. modrel) go to 20                                            
      if (mode .eq. modint ) real = ivalue                                      
      if (mode .ne. modint ) real = value                                       
      numr = .true.                                                             
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * point                                                              *        
c *                                                                    *        
c **********************************************************************        
      logical function point(point1,point2)                                     
c                                                                               
c         is the scan entity a point                                            
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      integer point1, point2                                                    
      data modpt/11/                                                            
c                                                                               
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 point = .false.                                                           
      if(mode.ne.modpt)go to 20                                                 
      point1 = ival(1)                                                          
      point2 = ival(2)                                                          
      point = .true.                                                            
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * readsc                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine readsc                                                         
c                                                                               
c         advance the scanner now                                               
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      data modeol/9/                                                            
c                                                                               
      if(.not.autord.or.(autord.and.mode.ne.modeol))call rdline                 
      next = .true.                                                             
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * realn                                                              *        
c *                                                                    *        
c **********************************************************************        
      logical function realn(real)                                              
c                                                                               
c         is the current scan entity a real                                     
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modrel/2/                                                            
c                                                                               
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 realn = .false.                                                           
      if(mode.ne.modrel)go to 20                                                
      real = value                                                              
      realn = .true.                                                            
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * sep                                                                *        
c *                                                                    *        
c **********************************************************************        
      logical function sep(septyp)                                              
c                                                                               
c         is the current scan entity a seperator                                
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      integer septyp                                                            
      data modsep/6/                                                            
c                                                                               
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 sep = .false.                                                             
      if(mode.ne.modsep)go to 20                                                
      septyp = ivalue                                                           
      sep = .true.                                                              
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * isstring - do not touch next variable                              *        
c *                                                                    *        
c **********************************************************************        
      logical function isstring(dummy)                                          
c                                                                               
c         is the current scan entity a string                                   
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modstr/7/                                                            
c                                                                               
      isstring = .false.                                                        
      if(mode.eq.modstr) isstring = .true.                                      
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * string                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function string(dummy)                                            
c                                                                               
c         is the current scan entity a string                                   
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      data modstr/7/                                                            
c                                                                               
      if(.not.next)go to 10                                                     
      call scan                                                                 
      next = .false.                                                            
   10 string = .false.                                                          
      if(mode.ne.modstr)go to 20                                                
      string = .true.                                                           
      next = .true.                                                             
   20 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * true                                                               *        
c *                                                                    *        
c **********************************************************************        
      logical function true(dummy)                                              
c                                                                               
c         advance the scanner before the next test                              
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      next = .true.                                                             
      true = .true.                                                             
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * integerlist                                                     *           
c *                                                                 *           
c *******************************************************************           
      subroutine trscan_list ( list, mlist, iall, nlist, ierr )                 
c                                                                               
c          scan action to input a list of integer terms                         
c              each action may be delimitted by a comma                         
c              a comma preceeding an eol indicates continuation                 
c              terms may be:                                                    
c                   a) <integer>                                                
c                   b) <integer1> - <integer2>                                  
c                   c) <integer1> to <integer2>                                 
c                   d) <integer1> - <integer2> by <integer3>                    
c                   e) <integer1> to <integer2> by <integer3>                   
c                   f) all                                                      
c              type a stores <integer>  in list                                 
c              type b and c stores <integer1>, -<integer2>, 1 in list           
c              type d and e stores <integer1>, -<integer2>, <integer3>          
c              in list                                                          
c              type f stores 1, -iall, 1 in list                                
c                                                                               
c         dummy arguments                                                       
c              list      (output) - list of parsed input - as described         
c              mlist     (input)  - allowable size of list                      
c              iall      (input)  - value of 'all'                              
c                                   = 0 - 'all' is not acceptable               
c              nlist     (output) - number of terms stored in list              
c              ierr      (output) - error code                                  
c                                   = 1 - no error                              
c                                   = 2 - parse rules failded                   
c                                   = 3 - list overflow                         
c                                   = 4 - list not found                        
c                                                                               
c         common blocks                                                         
c              labelled common /scaner/                                         
c                                                                               
c         called subprograms                                                    
c              scan                                                             
c              rdline                                                           
c              scanmc                                                           
c                                                                               
c         local variables                                                       
c              istate            - fsa states                                   
c              nstate            - next state table                             
c              iclass            - class of input                               
c              ifsa              - action table                                 
c              iact              - action to do                                 
c              iby               - hollerth 'by'                                
c              ito               - hollerth 'to'                                
c              jall              - hollerth 'all'                               
c              istart            - flag set to true on first scan               
c              iovfl             - flag set to true on overflow                 
c                                                                               
c         algorithm terms                                                       
c              states   - 1 = start                                             
c                         2 = item (integer)                                    
c                         3 = delimiter (,)                                     
c                         4 = iteration (-, to)                                 
c                         5 = increment (by)                                    
c                         6 = delta (increment integer)                         
c              classes  - 1 = integer                                           
c                         2 = alpha 'to' or seperator '-'                       
c                         3 = alpha 'by'                                        
c                         4 = seperator ','                                     
c                         5 = end of line                                       
c                         6 = else                                              
c              actions  - 0 = switch state                                      
c                         1 = done                                              
c                         2 = error                                             
c                         3 = save inger, - value implies iteration te          
c                         4 = read line                                         
c                         5 = save -1*integer                                   
c                         6 = save 1, save integer                              
c                         7 = save 1                                            
c                         8 = save 1, done                                      
c                         9 = test overflow                                     
c                        10 = save iteration increment                          
c                                                                               
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      dimension list(mlist), ifsa(6,6), nstate(6,6)                             
      integer entity                                                            
      logical istart, iovfl                                                     
      logical scanmc                                                            
      dimension iby(1), ito(1), jall(1)                                         
      data iby/2hby/, ito/2hto/, jall/3hall/                                    
c                                                                               
c         transfer table                                                        
c                                                                               
      data ifsa/                                                                
     1            3, 1, 1, 0, 1, 1,                                             
     2            3, 9, 2, 0, 1, 1,                                             
     3            3, 2, 2, 2, 4, 2,                                             
     4            5, 2, 2, 2, 2, 2,                                             
     5            6, 2, 0, 7, 8, 8,                                             
     6           10, 2, 2, 0, 4, 2/                                             
c                                                                               
c         next state table                                                      
c                                                                               
      data nstate/                                                              
     1            2, 0, 0, 1, 0, 0,                                             
     2            2, 4, 0, 3, 0, 0,                                             
     3            2, 0, 0, 0, 1, 0,                                             
     4            5, 0, 0, 0, 0, 0,                                             
     5            2, 0, 6, 3, 0, 0,                                             
     6            2, 0, 0, 6, 6, 0/                                             
c                                                                               
c         initialy in start state                                               
c                                                                               
      istart = .true.                                                           
      iovfl = .false.                                                           
      nlist = 1                                                                 
      istate = 1                                                                
c                                                                               
c         scan and determine class                                              
c                                                                               
 100  iclass = 6                                                                
      if(.not.istart)call scan                                                  
      if(mode.ne.9)go to 110                                                    
      iclass = 5                                                                
      go to 140                                                                 
 110  if(mode.ne.6)go to 120                                                    
      if(ivalue.eq.3)iclass = 2                                                 
      if(ivalue.eq.16)iclass = 4                                                
      go to 140                                                                 
 120  if( mode.ne.3 )go to 130                                                  
      if( scanmc(entity(1),iby,2) .and. nchar.eq.2 ) iclass = 3                 
      if( scanmc(entity(1),ito,2) .and. nchar.eq.2 ) iclass = 2                 
      if( scanmc(entity(1),jall,3) .and. istart .and. iall.ne.0                 
     &    .and. nchar.eq.3 ) go to 350                                          
      go to 140                                                                 
 130  if(mode.ne.1)go to 140                                                    
      iclass = 1                                                                
c                                                                               
c         determine next state and action                                       
c                                                                               
 140  iact = ifsa(iclass,istate)+1                                              
      istate = nstate(iclass,istate)                                            
      if(iact.ne.2)istart=.false.                                               
c                                                                               
c         action transfer, skip store on overflow                               
c                                                                               
      if(iovfl.and.(iact.eq.4.or.iact.eq.6.or.iact.eq.7.or.                     
     1              iact.eq.8.or.iact.eq.10))go to 100                          
      go to (100, 210, 220, 230, 250, 260, 270, 280, 290, 310, 320              
     1       ), iact                                                            
c                                                                               
c         done                                                                  
c                                                                               
 210  ierr = 1                                                                  
      nlist = nlist-1                                                           
      if(istart)ierr = 4                                                        
      if(iovfl)ierr = 3                                                         
      if(list(1).lt.0.and.nlist.gt.1)ierr = 2                                   
      return                                                                    
c                                                                               
c         error - parse                                                         
c                                                                               
 220  ierr = 2                                                                  
      nlist = nlist-1                                                           
      return                                                                    
c                                                                               
c         save integer, -value becomes iteration bound                          
c                                                                               
 230  if(nlist.le.mlist)go to 240                                               
      iovfl = .true.                                                            
      go to 100                                                                 
 240  list(nlist) = ivalue                                                      
      nlist = nlist+1                                                           
      if(ivalue.lt.0.and.nlist.eq.2)go to 100                                   
      if(ivalue.lt.0)istate = 5                                                 
      go to 100                                                                 
c                                                                               
c         read line                                                             
c                                                                               
 250  call rdline                                                               
      go to 100                                                                 
c                                                                               
c         save -1*integer                                                       
c                                                                               
 260  list(nlist) = -ivalue                                                     
      if(ivalue.lt.0)go to 220                                                  
      nlist = nlist+1                                                           
      go to 100                                                                 
c                                                                               
c         save 1, save integer                                                  
c                                                                               
 270  list(nlist) = 1                                                           
      list(nlist+1) = ivalue                                                    
      nlist = nlist+2                                                           
      go to 100                                                                 
c                                                                               
c         save 1                                                                
c                                                                               
 280  list(nlist) = 1                                                           
      nlist = nlist+1                                                           
      go to 100                                                                 
c                                                                               
c         save 1, done                                                          
c                                                                               
 290  if(iovfl)go to 300                                                        
      list(nlist) = 1                                                           
 300  ierr = 1                                                                  
      if(iovfl)ierr = 3                                                         
      if(list(1).lt.0.and.nlist.ne.1)ierr = 2                                   
      return                                                                    
c                                                                               
c         test for overflow in iteration term                                   
c                                                                               
 310  if(nlist+1.gt.mlist)iovfl = .true.                                        
      go to 100                                                                 
c                                                                               
c         save integer for increment                                            
c                                                                               
 320  if(nlist.le.mlist)go to 330                                               
      iovfl = .true.                                                            
      go to 100                                                                 
 330  list(nlist) = ivalue                                                      
      nlist = nlist+1                                                           
      go to 100                                                                 
c                                                                               
c        all, done                                                              
c                                                                               
 350  if(mlist.lt.3)go to 360                                                   
      list(1) = 1                                                               
      list(2) = -iall                                                           
      list(3) = 1                                                               
      nlist = 3                                                                 
      ierr = 1                                                                  
      call scan                                                                 
      if(mode.ne.6)return                                                       
      if(ivalue.ne.16)return                                                    
      call scan                                                                 
      if(mode.ne.9)return                                                       
      call rdline                                                               
      call scan                                                                 
      return                                                                    
 360  ierr = 3                                                                  
      nlist = 0                                                                 
      return                                                                    
c                                                                               
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * alpha match                                                     *           
c *                                                                 *           
c *******************************************************************           
      subroutine tramat ( word, mchar, matchl, loc, scnact )                    
c                                                                               
c         alphamatch                                                            
c         match a word against a list of words                                  
c                                                                               
c         dummy arguments                                                       
c               word     (input)  - word to match against list                  
c               mchar    (input)  - length of word in characters                
c               matchl   (input)  - list of words to match against              
c                                   first entry is the number of words          
c                                   the word entries follow, each begins        
c                                   with one word which is the length           
c               loc      (output) - the location of word in matchl              
c                                   = 0 if not found                            
c               scnact   (input)  - scan action if match is found               
c                                   = 1 - scan                                  
c                                   = 2 - scan & skip ,                         
c                                   = 3 - scan & skip , & end                   
c                                         do the readline but don't scan        
c                                   = 4 - scan & skip , & end                   
c                                         readline & scan                       
c                                   = else - noscan                             
c                                                                               
c                                                                               
      dimension word(1), matchl(1)                                              
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      integer scnact                                                            
      logical scanmc                                                            
c                                                                               
c         go through the list to find the match                                 
c                                                                               
      n = matchl(1)                                                             
      loc = 0                                                                   
      j = 2                                                                     
      do 10 i = 1, n                                                            
      l = matchl(j)                                                             
      if(l.gt.nchar)go to 10                                                    
      if(scanmc(word,matchl(j+1),l))go to 20                                    
      j = j+(l+ncpw-1)/ncpw+1                                                   
   10 continue                                                                  
      return                                                                    
   20 loc = i                                                                   
      if(scnact.lt.1.or.scnact.gt.4)return                                      
      call scan                                                                 
      if(scnact.eq.1)return                                                     
      if(mode.ne.6.or.ivalue.ne.16)return                                       
      call scan                                                                 
      if(scnact.eq.2)return                                                     
      if(mode.ne.9)return                                                       
      call rdline                                                               
      if(scnact.eq.3)return                                                     
      call scan                                                                 
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * alpha list                                                      *           
c *                                                                 *           
c *******************************************************************           
      subroutine tralst ( matchl, ifound, nmax, nfound, ierror )                
c                                                                               
c         alpha list                                                            
c         match alpha input against a list of words                             
c                                                                               
c         dummy arguments                                                       
c               matchl   (input)  - list of words to match against              
c                                   first entry is the number of words          
c                                   the word entries follow, each begins        
c                                   with one word which is the length           
c               ifound   (output) - a vector of what matches were found         
c               nmax     (input)  - the length of ifound                        
c               nfound   (output) - the number of matches, entries in           
c                                   ifound                                      
c               ierror   (output) - error flag                                  
c                                   = 1 - no errors                             
c                                   = 2 = syntax                                
c                                   = 3 - overflow of ifound                    
c                                                                               
c                                                                               
      dimension matchl(1), ifound(1)                                            
      logical scanmc                                                            
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
c                                                                               
c         go through the list to look for a match                               
c                                                                               
      ierror = 1                                                                
      nfound = 0                                                                
      n = matchl(1)                                                             
      if(mode.eq.3)go to 10                                                     
      return                                                                    
   10 j = 2                                                                     
      do 20 i = 1, n                                                            
      l = matchl(j)                                                             
      if(l.gt.nchar)go to 20                                                    
      if(scanmc(entity(1),matchl(j+1),l))go to 30                               
      j = j+(l+ncpw-1)/ncpw+1                                                   
  20  continue                                                                  
      if(nfound.gt.nmax)ierror = 3                                              
      return                                                                    
c                                                                               
c         save it and get the next item                                         
c         skip , or autoread on ,&eol                                           
c                                                                               
   30 nfound = nfound+1                                                         
      if(nfound.le.nmax)ifound(nfound) = i                                      
      call scan                                                                 
      if(mode.eq.3)go to 10                                                     
      if(mode.eq.6.and.ivalue.eq.16)go to 50                                    
      if(nfound.gt.nmax)ierror = 3                                              
      return                                                                    
   50 call scan                                                                 
      if(mode.eq.9)go to 60                                                     
      if(mode.eq.3)go to 10                                                     
      ierror = 2                                                                
      return                                                                    
   60 call rdline                                                               
      call scan                                                                 
      if(mode.eq.3)go to 10                                                     
      if(nfound.gt.nmax)ierror = 3                                              
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * next on integer list                                            *           
c *                                                                 *           
c *******************************************************************           
      subroutine trxlst ( list, nl, iplist, icn, next )                         
c                                                                               
c         polo action to determine next entry on an integerlist                 
c                                                                               
c         dummy arguments                                                       
c              list      (input)  - list of data                                
c              nl        (input)  - length of list                              
c              iplist    (update) - input  - current list element               
c                                          - = 1 on first call                  
c                                   output - next list element                  
c                                            - = 0 on end of list               
c                                                next is returned               
c              icn       (update) - input  - current iteration count            
c                                            - = 0 on first call                
c                                   output - next iteration count               
c              next      (output) - next item on list to be proccessed          
c                                                                               
c                                                                               
      dimension list(nl)                                                        
c                                                                               
c         get result                                                            
c                                                                               
      if(iplist.le.0)return                                                     
      next = list(iplist)                                                       
      iplist = iplist+1                                                         
      if(iplist.gt.nl)go to 120                                                 
c                                                                               
c         next is plus, single term                                             
c                                                                               
      if(list(iplist).ge.0)return                                               
c                                                                               
c         iteration term (always done once), loop exhausted                     
c                                                                               
      next = next+icn*list(iplist+1)                                            
      ir = next+list(iplist+1)                                                  
      if(list(iplist+1).gt.0.and.                                               
     1  (ir.le.iabs(list(iplist))))                                             
     2  go to 110                                                               
      if(list(iplist+1).lt.0.and.                                               
     1  (ir.ge.iabs(list(iplist))))                                             
     2  go to 110                                                               
      icn =  0                                                                  
      iplist = iplist+2                                                         
      if(iplist.gt.nl)go to 120                                                 
      return                                                                    
c                                                                               
c         iterate                                                               
c                                                                               
 110  iplist = iplist-1                                                         
      icn = icn+1                                                               
      return                                                                    
c                                                                               
c         end of list                                                           
c                                                                               
 120  iplist = 0                                                                
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * iszlst - size of a list                                         *           
c *                                                                 *           
c *******************************************************************           
      integer function iszlst ( list, nl )                                      
c                                                                               
c                                                                               
c          determine the number of terms in an integerlist                      
c                                                                               
c                                                                               
      dimension list(nl)                                                        
      ip = 0                                                                    
      ir = 0                                                                    
   10 if ( ip.ge.nl ) go to 20                                                  
      ip = ip + 1                                                               
      ir = ir + 1                                                               
      if ( list(ip).ge.0 ) go to 10                                             
      ir = ir + ((iabs(list(ip))-list(ip-1)+list(ip+1))/list(ip+1))-2           
      ip = ip + 1                                                               
      go to 10                                                                  
   20 iszlst = ir                                                               
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c * scandb                                                          *           
c *                                                                 *           
c *******************************************************************           
      subroutine scandb                                                         
c                                                                               
c         debug - all character output at 4 per word                            
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      write(iout,1001)                                                          
      write(iout,1002)recsiz, limit, mark, echar, echo, promt, ilabel,          
     1                init, point, menu, leof, eol, autord, commnt,             
     2                signed                                                    
      write(iout,1003)inunit, iout, inremo, filpt, fillim, files                
      write(iout,1004)col, jump, nent, pstate, skip, doread, incol,             
     3                card, jbuff, ibuff, idigit                                
      write(iout,1005)entity, ival, mode,nchar, nwd,                            
     1                next, icolmn                                              
      return                                                                    
 1001 format(1h1/24h ***********************/                                   
     1           24h * s c a n   d e b u g */                                   
     2           24h ***********************/)                                  
 1002 format( 5x, 25hc o n t r o l   s t a t e/                                 
     1       10x, 14hrecord size - , i2/                                        
     2       10x, 15hrecord limit - , i2/                                       
     3       10x, 13hrecord end - , i2/                                         
     4       10x, 23hrecord end character - , a1/                               
     5       10x, 14hecho switch - , l1/                                        
     6       10x, 16hprompt switch - , l1/                                      
     7       10x, 13hline label - , i5/                                         
     8       10x, 24hinitialization switch - , l1/                              
     9       10x, 15hpoint switch - , l1/                                       
     a       10x, 14hmenu switch - , l1/                                        
     b       10x, 21hend of file switch - , l1/                                 
     c       10x, 21hend of line switch - , l1/                                 
     d       10x, 19hauto read switch - ,l1/                                    
     e       10x, 17hcomment switch - , l1/                                     
     f       10x, 14hsign switch - , l1/)                                       
 1003 format( 5x, 17hi / o   s t a t e/                                         
     1       10x, 13hinput unit - , i2/                                         
     2       10x, 14houtput unit - , i2/                                        
     3       10x, 14hremote unit - , i2/                                        
     4       10x, 21hfile stack pointer - , i2/                                 
     5       10x, 13hfile limit - , i2/                                         
     6       10x, 13hfile stack - , 10(i2, 2x)/)                                
 1004 format( 5x, 19hs c a n   s t a t e/                                       
     1       10x, 9hcolumn - , i2/                                              
     2       10x, 7hjump - , i2/                                                
     3       10x, 16hentity number - , i2/                                      
     4       10x, 15hstring state - , i2/                                       
     5       10x, 13hskip state - , i2/                                         
     6       10x, 14hread switch - , l1/                                        
     7       10x, 14hrecord size - , i2/                                        
     8       10x, 13hinput line - , 81a1/                                       
     9       10x, 18hinternal record - , 20(i3,1x)/28x,20(i3,1x)/               
     a            28x,20(i3,1x)/28x,21(i3,1x)/                                  
     b       10x, 8hclass - , 20(i3,1x)/18x,20(i3,1x)/18x,20(i3,1x)/            
     c            18x,21(i3,1x)/                                                
     d       10x, 9hdigits - , 20(i3,1x)/19x,20(i3,1x)/19x,20(i3,1x)/           
     e            19x,21(i3,1x)/)                                               
c                                                                               
c                the a4 format is machine dependent                             
c                                                                               
 1005 format( 5x, 11he n t i t y/                                               
     1       10x, 9hentity - , 20a4/                                            
     4       10x, 7hival - , z8,1x,z8/                                          
     5       10x, 7hmode - , i2/                                                
     6       10x, 12hcharacter - , i2/                                          
     7       10x, 8hwords - , i2/                                               
     8       10x, 7hnext - , l1/                                                
     9       10x, 9hcolumn - , i2///)                                           
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scancl                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine scancl                                                         
c                                                                               
c       classify the input characters                                           
c                                                                               
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      common /scantb/ nseptb, nclass, itab(256), iclass(256)                    
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      character(len=1) :: xcard(640), xbuff(640)                                
      equivalence ( xcard(1), card(1) ), ( xbuff(1), jbuff(1) )                 
      j = 1                                                                     
      k = ncpw                                                                  
      nn = min0(incol,mark)                                                     
      nn = nn+1                                                                 
      card(nn)   = echar                                                        
c                                                                               
c                                                                               
c         the loop first constructs xbuff, where each character is              
c         moved to the last byte in the word.  Then if the character is         
c         a digit in ascii (0-9), it is converted to the actual number.         
c         if the character is found to be a tab, then it is converted to        
c         a space.  The classification of the character is then found in        
c         the iclass array.                                                     
c                                                                               
c         note: for most word orientations, : xbuff(k) = xcard(j)               
c                                                                               
      do 10 i = 1, nn                                                           
         xbuff(j) = xcard(j)                                                    
         idigit(i) = -1                                                         
         if(jbuff(i).ge.intzer.and.jbuff(i).lt.intnin)                          
     &        idigit(i) = jbuff(i)-intzer                                       
         if(jbuff(i).eq.inttab) jbuff(i) = intblk                               
         ibuff(i) = iclass(jbuff(i)+1)                                          
         if(card(i).eq.echar)ibuff(i) = 10                                      
         j = j+ncpw                                                             
         k = k+ncpw                                                             
 10   continue                                                                  
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scanin                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine scanin(in,buff,reclen,errcod)                                  
c                                                                               
c         input a record from a device                                          
c                                                                               
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      integer reclen, errcod, recsiz                                            
      dimension buff(80)                                                        
      character(len=1) :: chartab(4)                                            
      equivalence (holtab, chartab(1))                                          
c                                                                               
c             create a holerith tab for checking for tabs                       
c                                                                               
      holtab = blank                                                            
      chartab(1) = char(9)                                                      
c                                                                               
c             do the read                                                       
c                                                                               
      recsiz = reclen                                                           
      errcod = 0                                                                
      read(in,1001,err=20,end=30)(buff(i),i=1,recsiz)                           
c                                                                               
c             find length of line. count tabs at end of line as                 
c             blanks                                                            
c                                                                               
      do i = 1, recsiz                                                          
         if(.not.(buff(reclen).eq.blank) .and.                                  
     &      .not.(buff(reclen).eq.holtab)     ) then                            
            go to 40                                                            
         endif                                                                  
         reclen = reclen - 1                                                    
      enddo                                                                     
      reclen = 0                                                                
      return                                                                    
c                                                                               
c             process errors during read. 20 = error, 30 = eof                  
c                                                                               
 20   errcod = 1                                                                
      go to 40                                                                  
 30   errcod = -1                                                               
c                                                                               
c             ready to return.                                                  
c                                                                               
   40 reclen = max0(reclen,2)                                                   
      return                                                                    
 1001 format(80a1)                                                              
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scanpk                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine scanpk                                                         
c                                                                               
c           pack the scaner result in entity at 4 characters per word           
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      character(len=1) :: xentit(80), xcard(320)                                
      equivalence ( xentit(1), entity(1) ), ( xcard(1), card(1) )               
c                                                                               
      nwd = (nchar+ncpw-1)/ncpw                                                 
      if(nwd.le.0)nwd = 1                                                       
      entity(nwd) = blank                                                       
      if(nwd.lt.20)entity(nwd+1) = blank                                        
      istart = (icolmn-1)*ncpw+1                                                
      if(nchar.eq.0)nwd = 0                                                     
      if(nchar.eq.0)go to 20                                                    
      if(mode.eq.7)istart = istart+ncpw                                         
      do 10 i = 1, nchar                                                        
      xentit(i) = xcard(istart)                                                 
      istart = istart+ncpw                                                      
      if(mode.eq.4)istart = (jstart(nent+i)-1)*ncpw+1                           
   10 continue                                                                  
   20 return                                                                    
      end                                                                       
                                                                                
c **********************************************************************        
c *                                                                    *        
c * scanrd                                                             *        
c *                                                                    *        
c **********************************************************************        
      subroutine scanrd                                                         
      use file_info                                                             
c                                                                               
c         input for scan                                                        
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
      logical reread, comlin                                                    
      logical matchs                                                            
      character(len=1) :: dums                                                  
c                                                                               
      logical comlin2                                                           
                                                                                
      logical nowopen                                                           
      character(len=80) :: nowname                                              
      data inthash / 33 /, intexclam / 35 /                                     
c                                                                               
c         read a record.  do a re-read for blank line.                          
c                                                                               
      doread = .false.                                                          
   10 continue                                                                  
      reread = .false.                                                          
      incol = recsiz                                                            
      if (  promt .and. ( inunit .eq. inremo ) )                                
     1                call wrnocr( 3h > , iotrem )                              
      call scanin( inunit, card, incol, ierr )                                  
      if ( incol .eq. 0 ) go to 10                                              
      if( ierr < 0 ) go to 20                                                   
      if( ierr == 0 ) go to 40                                                  
      if( ierr > 0 ) go to 30                                                   
c                                                                               
c         end of file -- pop the file stack                                     
c                                                                               
   20 continue                                                                  
      close (inlun(filcnt))                                                     
      filcnt=filcnt-1                                                           
      if (filcnt.eq.0) then                                                     
         call errmsg(204,dum,dums,dumr,dumd)                                    
         call die_gracefully                                                    
         stop                                                                   
      endif                                                                     
      call setin(inlun(filcnt))                                                 
      in = inlun(filcnt)                                                        
      call errmsg(179,dum,dums,dumr, dumd)                                      
      go to 10                                                                  
c                                                                               
c        read error                                                             
c                                                                               
   30 call scanm(2)                                                             
      inquire(unit=inunit,opened=nowopen,name=nowname)                          
      write(*,*) '..unit: ',inunit                                              
      write(*,*) '..opened: ',nowopen                                           
      write(*,*) '..name: ',nowname                                             
      call die_gracefully                                                       
      stop                                                                      
c                                                                               
c         classify                                                              
c                                                                               
   40 call scancl                                                               
      mode = 0                                                                  
c                                                                               
c         comment  'c ' in columns 1 & 2, # in column 1 or                      
c                   ! in column 1                                               
c                                                                               
      comlin  = jbuff(1) .eq. intc .or. jbuff(1) .eq. intlcc                    
      comlin  = comlin .and. jbuff(2) .eq. intblk                               
      comlin2 = (jbuff(1) .eq. inthash) .or. (jbuff(1) .eq. intexclam)          
      comlin = comlin .or. comlin2                                              
      if ( .not. comlin ) go to 50                                              
      if( commnt ) reread = .true.                                              
      go to 70                                                                  
c                                                                               
c         point 'dddd,dddd' in cc1-9 = (4@zf0-f9),z6b,(4@zf0-f9)                
c                                      = (4@240-249),107,(4@240-249)            
c                                                                               
 50   if(.not.point) go to 70                                                   
      if(incol.ne.9)go to 70                                                    
      if(jbuff(5).ne.intcom)go to 70                                            
      do 60 i = 1, 4                                                            
         if(idigit(i).eq.-1.or.idigit(i+5).eq.-1)go to 70                       
 60   continue                                                                  
      mode = 11                                                                 
      go to 70                                                                  
c                                                                               
c         echo (not a point with menuing)                                       
c                                                                               
   70 if(mode.eq.11.and.menu)go to 75                                           
      if(echo) then                                                             
         nn = min0(incol,mark)                                                  
         if(ilabel.eq.0) then                                                   
            write(iout,1002)(card(i),i=1,nn)                                    
         else                                                                   
            write(iout,1001)ilabel,(card(i),i=1,nn)                             
            ilabel = 0                                                          
         endif                                                                  
      endif                                                                     
      if(reread)go to 10                                                        
c                                                                               
c         check for non-blank characters in columns 80+ of                      
c         this non-comment line. if found, stop program.                        
c                                                                               
 75   continue                                                                  
      if ( incol .gt. 80 ) then                                                 
         write(iout,1005)                                                       
         call die_gracefully                                                    
         stop                                                                   
      end if                                                                    
C                                                                               
c         check for an asterisk ('*') before the command.                       
c         must set next=true so that scan will be called.                       
c         if there is an asterisk, then call star_com to                        
c         process the star command, then read the next line.                    
c         if there isn't an asterisk, then rest scan and go on.                 
c                                                                               
 80   continue                                                                  
      if (isct) then                                                            
         goto 85                                                                
      endif                                                                     
      next = .true.                                                             
      if (matchs('*',1)) then                                                   
         call star_com                                                          
         jump = 1                                                               
         goto 10                                                                
      else                                                                      
         call reset                                                             
      endif                                                                     
c                                                                               
c                                                                               
c                     patch to  look at the end of the line. if it              
c                     has a comma just before the end of line it                
c                     has a continuation line following it. throw away          
c                     the comma but set the variable 'isct' true.               
c                     this will cause scan to do an automatic read of           
c                     the next line when the end of this line is                
c                     encountered. the effect is to treat multilpe              
c                     lines as one statement without the calling                
c                     program seeing the intermediate end of lines.             
c                                                                               
c                                                                               
 85   continue                                                                  
      isct=.false.                                                              
C                                                                               
      if( autoct )then                                                          
         nn=min(incol,mark)  + 1                                                
         do 91 kk=1,nn                                                          
            if(ibuff(kk).eq.10) then                                            
               nn=kk-1                                                          
               if(nn .ne.  0 )then                                              
                  ii = kk                                                       
                  do 93 jj=1,nn                                                 
                     ii=ii-1                                                    
                     if(jbuff(ii).ne.intblk) then                               
                        isct=jbuff(ii).eq.intcom                                
                        if(isct) then                                           
                           ibuff(ii)=10                                         
                        endif                                                   
                        goto 100                                                
                     endif                                                      
   93             continue                                                      
               endif                                                            
               goto 100                                                         
            endif                                                               
   91    continue                                                               
      endif                                                                     
  100 continue                                                                  
      return                                                                    
 1001 format(1x,i5,1x,80a1)                                                     
 1002 format(7x,      80a1)                                                     
 1003 format(3h ? )                                                             
 1005 format('>>>> Fatal Error: non-blank data in columns beyond 80.',          
     &   /,  '                  job terminated to prevent subsequent',          
     &   /,  '                  errors.',//)                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scanmc                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function scanmc(texta,textb,nchar)                                
c                                                                               
c         character match of a and b for nchars.  ignores                       
c         case of characters during compare.  this routine                      
c         depends on use of the ascii character set in scan.                    
c                                                                               
      dimension texta(1), textb(1)                                              
      real ja, jb, jc, jd                                                       
      integer char1, char2                                                      
      character(len=1) :: aa(4), bb(4), cc(4), dd(4)                            
      equivalence (aa(1),ja), (bb(1),jb), (cc(1),jc), (dd(1),jd)                
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
c                                                                               
      scanmc = .false.                                                          
      nw = nchar/ncpw                                                           
      nr = nchar-nw*ncpw                                                        
      if(nw.eq.0)go to 20                                                       
      do 10 i = 1, nw                                                           
c                                                                               
c             get each character in the word, convert to upper                  
c             case and do the compare.                                          
c                                                                               
      ja = texta(i)                                                             
      jb = textb(i)                                                             
      do ii = 1, ncpw                                                           
      char1 = ichar(aa(ii))                                                     
      char2 = ichar(bb(ii))                                                     
      if ( char1.ge.97 .and. char1.le.122 ) char1 = char1 - 32                  
      if ( char2.ge.97 .and. char2.le.122 ) char2 = char2 - 32                  
      if ( char1 .ne. char2 ) go to 50                                          
      end do                                                                    
   10 continue                                                                  
c                                                                               
c             do characters in partial last word                                
c                                                                               
   20 if(nr.eq.0)go to 40                                                       
      ja = texta(nw+1)                                                          
      jb = textb(nw+1)                                                          
      do 30 i = 1, nr                                                           
      char1 = ichar(aa(i))                                                      
      char2 = ichar(bb(i))                                                      
      if ( char1.ge.97 .and. char1.le.122 ) char1 = char1 - 32                  
      if ( char2.ge.97 .and. char2.le.122 ) char2 = char2 - 32                  
      if ( char1 .ne. char2 ) go to 50                                          
   30 continue                                                                  
   40 scanmc = .true.                                                           
   50 return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * scanms                                                             *        
c *                                                                    *        
c **********************************************************************        
      logical function scanms( texta, textb, nchar )                            
c                                                                               
c         character match of a and b for nchar.  ignores                        
c         case of characters during compare.  this routine                      
c         depends on use of the ascii character set in scan.                    
c                                                                               
c         texta and textb are GUARANTEED to be declared                         
c         in the calling program as type character so that                      
c         their lengths and memory alignment are correct.                       
c                                                                               
      implicit integer (a-z)                                                    
      character(len=*) :: texta, textb                                          
c                                                                               
      scanms = .false.                                                          
      lena   = len( texta )                                                     
      lenb   = len( textb )                                                     
      if ( (nchar .gt. lena) .or. (nchar .gt. lenb) ) then                      
       call scanm( 4 )                                                          
      end if                                                                    
c                                                                               
c         extract each pair of characters, convert to upper                     
c         case and do the comparsison                                           
c                                                                               
      do i = 1, nchar                                                           
       char1 = ichar( texta(i:i) )                                              
       char2 = ichar( textb(i:i) )                                              
       if ( char1.ge.97 .and. char1.le.122 ) char1 = char1 - 32                 
       if ( char2.ge.97 .and. char2.le.122 ) char2 = char2 - 32                 
       if ( char1 .ne. char2 ) return                                           
      end do                                                                    
c                                                                               
      scanms = .true.                                                           
      return                                                                    
c                                                                               
 9000 format(/,                                                                 
     & '>> ERROR:  text string mismatch in scanms....',/,                       
     & '           lena, lenb, nchar: ',3i6,/,                                  
     & '           program terminated.....',//)                                 
c                                                                               
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * issue prompt character to terminal without carriage return         *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      subroutine wrnocr( prompt, outrem )                                       
c                                                                               
c                write a prompt to a terminal without a carriage                
c                return.                                                        
c                                                                               
      integer   outrem, prompt                                                  
      write(outrem,fmt='(a3$)' ) prompt                                         
      return                                                                    
      end                                                                       
c **********************************************************************        
c *                                                                    *        
c * block data                                                         *        
c *                                                                    *        
c **********************************************************************        
       block data                                                               
c                                                                               
c          initialze the scanner                                                
c                                                                               
      common /scaner/ entity(20), ival(2), mode, nchar, nwd, next,              
     1                icolmn, ispace, dvalue                                    
      equivalence (ival(1),value,ivalue)                                        
      logical next                                                              
      double precision dvalue                                                   
      common /scanct/ echar, echo, ilabel, limit, mark, promt,                  
     1                point, recsiz, init, leof, eol, menu,                     
     2                autord, commnt, signed, autoct                            
      integer recsiz                                                            
      logical echo, promt, point, init, leof, eol, menu,                        
     1        commnt, autord, signed, autoct                                    
      common /scanio/ inunit, iout, files(10), filpt, fillim, inremo,           
     1                iotrem                                                    
      integer files, filpt, fillim                                              
      common /scanln/ col, jump, nent, pstate, skip, doread,                    
     1                incol, jstart(80), ibuff(81), jbuff(81),                  
     2                idigit(81), card(81), isct                                
      integer col, skip, pstate                                                 
      logical doread, isct                                                      
      common /scantb/ nseptb, nclass, itab(256), iclass(256)                    
      common /scanim/ ncpw, blank, intzer, intnin, intc, intcom,                
     &                intblk, intlcc, inttab, intstar                           
c                                                                               
c         /scaner/                                                              
c                                                                               
      data entity/20*4h     /, ival/2*0/, mode/0/, nchar/0/,                    
     1     nwd/0/, next/.false./, icolmn/0/                                     
c                                                                               
c         /scanct/                                                              
c                                                                               
      data echar/1h$/, echo/.true./, ilabel/0/, limit/80/, mark/80/,            
     1     point/.false./, recsiz/80/, init/.false./, leof/.false./,            
     2     eol/.true./, menu/.false./, autord/.false./,                         
     3     commnt/.true./, signed/.false./, promt/.true./                       
c                                                                               
c         /scanio/                                                              
c                                                                               
      data inunit/5/, iout/6/, files/10*0/, filpt/0/, fillim/10/,               
     1     inremo/5/, iotrem/6/                                                 
c                                                                               
c         /scanln/                                                              
c                                                                               
      data col/0/, jump/1/, nent/1/, pstate/0/, skip/1/,                        
     1     doread/.true./, incol/80/, jstart/80*0/,                             
     2     ibuff/81*9/, jbuff/81*9/, idigit/81*0/,                              
     3     card/81*4h     /                                                     
c                                                                               
c         /scantb/                                                              
c                                                                               
      data nseptb/256/, nclass/256/                                             
c                 0,1,2,3,4,5,6,7    0,1,2,3,4,5,6,7                            
      data iclass/8*8,               8*8,                                       
     2            8*8,               8*8,                                       
     4            9,8,6,8,8,8,8,6,   8,8,8,2,8,3,7,8,                           
     6            8*1,               1,1,8,8,8,8,8,8,                           
     8            8,5,5,5,5,4,5,5,   8*5,                                       
     a            8*5,               5,5,5,8,8,8,8,5,                           
     c            8,5,5,5,5,4,5,5,   8*5,                                       
     e            8*5,               5,5,5,8,8,8,8,8,                           
     1            8*8,               8*8,                                       
     3            8*8,               8*8,                                       
     5            9,8,6,8,8,8,8,6,   8,8,8,2,8,3,7,8,                           
     7            8*1,               1,1,8,8,8,8,8,8,                           
     9            8,5,5,5,5,4,5,5,   8*5,                                       
     b            8*5,               5,5,5,8,8,8,8,5,                           
     d            8,5,5,5,5,4,5,5,   8*5,                                       
     f            8*5,               5,5,5,8,8,8,8,8/                           
c                 0,1,2,3,4,5,6,7    0,1,2,3,4,5,6,7                            
c                                                                               
c                                                                               
c               00,01,02,03,04,05,06,07   00,01,02,03,04,05,06,07               
      data itab/8*0,                      8*0,                                  
     2          8*0,                      8*0,                                  
     4          00,00,26,22,09,17,07,24,  04,11,10,05,16,03,02,15,              
     6          8*0,                      00,00,21,12,00,25,19,20,              
     8          23,00,00,00,00,00,00,00,  8*0,                                  
     a          8*0,                      00,00,00,01,00,08,13,00,              
     c          8*0,                      8*0,                                  
     e          8*0,                      00,00,00,00,06,00,00,00,              
     1          8*0,                      8*0,                                  
     3          8*0,                      8*0,                                  
     5          00,00,26,22,09,17,07,24,  04,11,10,05,16,03,02,15,              
     7          8*0,                      00,00,21,12,00,25,19,20,              
     9          23,00,00,00,00,00,00,00,  8*0,                                  
     b          8*0,                      00,00,00,01,00,08,13,00,              
     d          8*0,                      8*0,                                  
     f          8*0,                      00,00,00,00,06,00,00,00/              
c                                                                               
c         /scanim/                                                              
c                                                                               
      data blank/4h     /, intzer/48/, intnin/58/, intc/67/,                    
     1     intcom/44/, intblk/32/, intlcc/99/, inttab/9/,                       
     2     intstar/42/                                                          
      data ncpw /4/                                                             
c                                                                               
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *            subroutine check_to_prompt                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/20/99                   *          
c     *     this determines if scan should issue a prompt before     *          
c     *     a read.                                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine check_to_prompt( prompt )                                      
      logical prompt                                                            
      prompt =  .false.                                                         
                                                                                
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine scan_flushline                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/25/2011                 *          
c     *                                                              *          
c     *     this subroutine skips a logical line of the input file   *          
c     *     when the autoct flag is set true (scan automatically     *          
c     *     reads next physical line when current line ends with ,   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine scan_flushline                                                 
      integer  dum                                                              
      logical  endcrd, true                                                     
c                                                                               
      do                                                                        
       if( endcrd(dum) ) return                                                 
       if( true(dum) ) call splunj; cycle                                       
      end do                                                                    
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
