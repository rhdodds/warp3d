c **********************************************************************        
c *                                                                    *        
c * scan  modules                                                      *        
c *                                                                    *        
c **********************************************************************        
c
c
      module scaner
c     -------------
      implicit none
c
      save 
      real :: value
      real :: entity(20)
      integer :: ival(2), mode, nchar, nwd, icolmn, ispace,
     &           ientity(20), ivalue
      logical :: next
      double precision :: dvalue
      character(len=80) :: scaner_string = " "    
      character(len=1) :: xentit(80)                                          
      equivalence ( ival, value, ivalue ), 
     &            ( entity, ientity, scaner_string, xentit )
c
      data ival/2*0/, mode/0/, nchar/0/,                    
     1     nwd/0/, next/.false./, icolmn/0/                                     
c
      end module scaner
c
c
      module scanim
c     -------------
      implicit none
c
      integer :: ncpw, intzer, intnin, intc, intcom, intblk,
     &           intlcc, inttab, intstar
      real :: blank
      character(len=4) :: scanim_blank = " "
      equivalence ( blank, scanim_blank ) 
c                                                                               
      data ncpw /4/, intzer/48/, intnin/58/, intc/67/,                    
     1     intcom/44/, intblk/32/, intlcc/99/, inttab/9/,                       
     2     intstar/42/                                                          
c
      end module scanim
c
c
      module scanln
c     -------------
      implicit none
c
      save
      integer :: col, jump, nent, pstate, skip, incol, jstart(80),
     &           ibuff(81), jbuff(81), idigit(81), card(81)
      logical :: doread, isct
      real :: rcard(81)
      character(len=81*4) :: string_card = " "
      character(len=1) :: xcard(320)  
      character(len=1) :: xbuff(640) !  xcard(640),                              
      equivalence  ( xbuff(1), jbuff(1) ),                 
     &             ( card, string_card, xcard, rcard )
c
      data col/0/, jump/1/, nent/1/, pstate/0/, skip/1/,                        
     1     doread/.true./, incol/80/, jstart/80*0/,                             
     2     ibuff/81*9/, jbuff/81*9/, idigit/81*0/
c
      end module scanln
c
c
      module scanio
c     -------------
      implicit none
c
      save
      integer :: inunit = 5, iout = 6, files(1:10) = 0, filpt = 0, 
     &           fillim = 10, inremo = 5, iotrem = 6          
c
      end module scanio
c 
c
      module scantb
c     -------------
      implicit none
c
      integer, save ::  nseptb, nclass, itab(256), iclass(256) 
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
      end module scantb
c
c
      module scanct
c     -------------
      implicit none
c
      save
      integer:: ilabel, limit, mark, recsiz
      real :: echar = real( Z'20202024', kind=kind(echar) ) ! $
      logical :: echo, promt, point, init, leof, eol, menu,                        
     1           commnt, autord, signed, autoct                                    
      data  echo/.true./, ilabel/0/, limit/80/, mark/80/,            
     1     point/.false./, recsiz/80/, init/.false./, leof/.false./,            
     2     eol/.true./, menu/.false./, autord/.false./,                         
     3     commnt/.true./, signed/.false./, promt/.true./, 
     4     autoct/.true./                       
c
      end module scanct
c
c
      module scan_macros
c     ------------------
      implicit none
c
      save
      integer:: num_macros = 0
      integer, parameter :: max_macros = 300
c
      type :: smacros
        integer :: nchars_id
        integer :: nchars_value
        character(len=80) :: id
        character(len=80) :: value
      end type
      type( smacros ), allocatable :: macros(:)
c
      end module scan_macros

