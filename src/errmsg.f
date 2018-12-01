c     ****************************************************************
c     *                      suboutine errmsg                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/26/2018 rhd             *
c     *                                                              *
c     *     this subroutine prints assorted error messages in re-    *
c     *     ponse to calls from all over the program. virtually all  *
c     *     error messages in the program are generated here.        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errmsg( errnum, param, sparam, rparam, dparam )
      use global_data ! old common.main
      use main_data, only : max_step_limit
      implicit integer (a-z)
      character(len=50) :: string
      character(len=35) :: strng1
      character(len=8) :: shrtst,shrtst1,shrtst2
      character(len=24) :: lngstr
      character(len=21) :: strng21
      character(len=29) :: strng29
      character(len=*) :: sparam
c
      double precision :: dparam, hundred
      real :: rparam
      integer, save :: high_lvl_count = 0
      logical, save :: mess_61, write_msg_255, write_msg_321
      integer, save :: count_150=0
      integer, save :: count_136=0
      data hundred /100.0/
      data mess_61, write_msg_255, write_msg_321
     &     / .true., .true., .true. /
c
c
c                       print the appropriate error message as
c                       indicated by the input variable errnum
c
c                       current holes:
c
c
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,
     &       160,170,180,190,200,210,220,230,240,250,260,270,280,
     &       290,300,310,320,330,340,350,360,370,380,390,400,410,
     &       420,430,440,450,460,470,480,490,500,510,520,530,540,
     &       550,560,570,580,590,600,610,620,630,640,650,660,670,
     &       680,690,700,710,720,730,740,750,760,770,780,790,800,
     &       810,820,830,840,850,860,870,880,890,900,910,920,930,
     &       940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,
     &       1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,1160,
     &       1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,
     &       1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,
     &       1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,
     &       1500,1510,1520,1530,1540,1550,1560,1570,1580,1590,1600,
     &       1610,1620,1630,1640,1650,1660,1670,1680,1690,1700,1710,
     &       1720,1730,1740,1750,1760,1770,1780,1790,1800,1810,1820,
     &       1830,1840,1850,1860,1870,1880,1890,1900,1910,1920,1930,
     &       1940,1950,1960,1970,1980,1990,2000,2010,2020,2030,2040,
     &       2050,2060,2070,2080,2090,2100,2110,2120,2130,2140,2150,
     &       2160,2170,2180,2190,2200,2210,2220,2230,2240,2250,2260,
     &       2270,2280,2290,2300,2310,2320,2330,2340,2350,2360,2370,
     &       2380,2390,2400,2410,2420,2430,2440,2450,2460,2470,2480,
     &       2490,2500,2510,2520,2530,2540,2550,2560,2570,2580,2590,
     &       2600,2610,2620,2630,2640,2650,2660,2670,2680,2690,2700,
     &       2710,2720,2730,2740,2750,2760,2770,2780,2790,2800,2810,
     &       2820,2830,2840,2850,2860,2870,2880,2890,2900,2910,2920,
     &       2930,2940,2950,2960,2970,2980,2990,3000,3010,3020,3030,
     &       3040,3050,3060,3070,3080,3090,3100,3110,3120,3130,3140,
     &       3150,3160,3170,3180,3190,3200,3210,3220,3230,3240,3250,
     &       3260,3270,3280,3290,3300,3310,3320,3330,3340,3350,3360,
     &       3370,3380,3390,3400,3410,3420,3430,3440,3450,3460,3470,
     &       3480,3490,3500,3510),
     &         errnum
c
c
 10   continue
      num_error = num_error + 1
      high_lvl_count = high_lvl_count + 1
      if( high_lvl_count == 20 ) then
        write(out,8999)
        go to 9999
      end if
      if( high_lvl_count > 20 ) go to 9999
      write(out,9000)
 8999 format(/1x,'>>>>> no further messages about high level command',
     &     /14x,'will be displayed.'/)

 9000 format(/1x,'>>>>> error: a high level command is expected. ',
     &       'continue'/14x,'scanning until one is encountered.'/)
      go to 9999
c
c
 20   continue
      num_error = num_error + 1
      write(out,9001) sparam
 9001 format(/1x,'>>>>> error: the number of ',a5,' has already been',
     &           ' input. any attempt',/,14x,'to alter this number ',
     &           'will be denied.'/)
      go to 9999
c
c
 30   continue
      num_error = num_error + 1
      write(out,9002)
 9002 format(/,1x,'>>>>> error: the name of the material is ',
     &         'expected. abort the input'/16x,'of this material ',
     &        'and scan for another high level',/,16x,'command.',/)
      go to 9999
c
c
 40   continue
      num_error = num_error + 1
      call entits( string, strlng )
      write(out,9003)string(1:strlng)
 9003 format(/1x,'>>>>> error: there is no match with an exist',
     &           'ing property. ignoring'/16x,'the current entity ',
     &           'and scan for another property.',/,16x,
     &           'scanning: ',a,/)
      go to 9999
c
c
 50   continue
      num_error = num_error + 1
      write(out,9004) sparam
 9004 format(/1x,'>>>>> error: a real or integer number is ',
     &           'expected for ',a4/)
      go to 9999
c
c
 60   continue
      num_warn = num_warn + 1
      write(out,9005)
 9005 format(/1x,'>>>>> warning: the label nodes or elements is ',
     &           'expected. continue'/16x,'scanning until such is ',
     &           'encountered.'/)
      go to 9999
c
c
 70   continue
      num_error = num_error + 1
      write(out,9006)
 9006 format(/1x,'>>>>> error: an integer indicating the number of ',
     &           'nodes in the'/14x,'structure is expected. a solu',
     &           'tion cannot be'/14x,'performed without such info',
     &           'rmation.'/)
      input_ok = .false.
      go to 9999
c
c
 80   continue
      go to 9999
c
c
 90   continue
      num_error = num_error + 1
      write(out,9008)
 9008 format(/1x,'>>>>> error: an integer indicating the number of ',
     &           'elements in the'/14x,'current structure is expect',
     &           'ed. a solution cannot be'/14x,'performed without ',
     &           'such information.'/)
      input_ok = .false.
      go to 9999
c
c
 100  continue
c                 available     
      go to 9999
c
c
 110  continue
      num_error = num_error + 1
      go to (111,112,113,114,115,116) param
 111  continue
      string= 'number of nodes                                 '
      write(out,9010) string
      go to 9999
 112  continue
      string= 'number of elements                              '
      write(out,9010) string
      go to 9999
 113  continue
      string= 'elements, number of nodes                       '
      write(out,9010) string
      go to 9999
c
 114  continue
      string= 'element properties & incidences, number of nodes'
      write(out,9010) string
      go to 9999
c
 115  continue
      string= 'element props, incidences, & blocking, coords   '
      write(out,9010) string
      go to 9999
c
 116  continue
      string= '  constraints                                   '
      write(out,9010) string
      go to 9999
c
 9010 format(/1x,'>>>>> error: a dependent hl command has been',
     &           ' encountered before the'/14x,'input of its pre-',
     &           'requisite command(s). these command(s)'/14x,'are: ',
     &           a50)
c
c
 120  continue
      num_warn = num_warn + 1
      write(out,9011)
 9011 format(/1x,'>>>>> warning: not enough default assignments ',
     &           'have been made. there'/16x,'is a probable error ',
     &           'in the coordinate data.'/)
      go to 9999
c
c
 130  continue
      num_warn = num_warn + 1
      write(out,9012)
 9012 format(/1x,'>>>>> warning: the current entity is not a valid ',
     &           'coordinate. no'/16x,'assignment is made.'/)
      go to 9999
c
c
 140  continue
      num_warn = num_warn + 1
      write(out,9013)
 9013 format(/1x,'>>>>> warning: repeated default assignment. ',
     &           'there is a probable'/16x,'error in the coordinate',
     &           ' data.'/)
      go to 9999
c
c
 150  continue
      num_fatal = num_fatal + 1
      count_150 = count_150 + 1
      if( count_150 <= 20 ) then
        write(out,9014) param
 9014 format(/1x,'>>>>> error: coordinates for node ',i7,' have',
     & ' not been input.',
     &           '  displacement computation cannot'/14x,
     &     'proceed without node coordiantes.')
      end if
      if( count_150 == 20 ) then
        write(out,*)
     &     '      >>> additional messages this type suppressed ...'
      end if
      input_ok = .false.
      go to 9999
c
c
 160  continue
      num_error = num_error + 1
      write(out,9015) param,nonode
 9015 format(/1x,'>>>>> error: ',i7,': the current node number input ',
     &           'exceeds ',i7,': the'/14x,'number of nodes in the ',
     &           'structure. it will be ignored.'/)
      go to 9999
c
c
 170  continue
      num_warn = num_warn + 1
      write(out,9016)
 9016 format(/1x,'>>>>> warning: a dof label is expected. ignore ',
     &           'the current entity'/16x,'and search for one.'/)
      go to 9999
c
c
 180  continue
      num_warn = num_warn + 1
      write(out,9017)
 9017 format(/1x,'>>>>> warning: a real or integer number is expe',
     &           'cted. the current'/16x,'entity will be ignored ',
     &           'and a number will be sought.'/)
      go to 9999
c
c
 190  continue
      num_warn = num_warn + 1
      write(out,9018)
 9018 format(/1x,'>>>>> warning: an attempt has been made to input ',
     &           'more coordinate'/16x,'data than there are dof for',
     &         ' this node. input for this'/16x,'node is terminated.'/)
      go to 9999
c
c
 200  continue
      num_warn = num_warn + 1
      write(out,9019) sparam
 9019 format(/1x,'>>>>> warning: a real or integer number is expe',
     &           'cted. an end of card'/16x,'has been encountered. ',
     &          'there is missing data for the'/16x,'dof ',a4,' and ',
     &           'the previous node"s value will be used.'/)
      go to 9999
c
c
 210  continue
      num_warn = num_warn + 1
      write(out,9020)
 9020 format(/1x,'>>>>> warning: a real or integer number is expe',
     &           'cted. the current'/16x,'entity will be ignored ',
     &           'and a new dof will be sought.'/)
      go to 9999
c
c
 220  continue
      num_fatal = num_fatal + 1
      write(out,9021)
 9021 format(/1x,'>>>>> fatal error: expecting the number of nodes ',
     &           'or elements to be input'/14x,'here. execution ',
     &           'cannot proceed without these numbers.'/)
      input_ok = .false.
      go to 9999
c
c
 230  continue
      num_warn = num_warn + 1
      write(out,9022)
 9022 format(/1x,'>>>>> warning: an integer indicating the current ',
     &           'node being input'/16x,'is expected. the remainder',
     &           ' of this line will be'/16x,'ignored and another ',
     &           'node will be sought.'/)
      go to 9999
c
c
 240  continue
      num_error = num_error + 1
      go to (241,242,243,244) param
 241  continue
      write(out,9023)
 9023 format(/1x,'>>>>> error: while attempting to translate a list',
     &           ' the parse rules'/14x,'failed. the current list ',
     &           ' will be discarded.'/)
      go to 9999
 242  continue
      write(out,1023)
 1023 format(/1x,'>>>>> error: while attempting to translate a list',
     &           ' a list overflow was'/14x,'encountered. the current ',
     &           ' list will be discarded.'/)
      go to 9999
 243  continue
      write(out,2023)
 2023 format(/1x,'>>>>> error: while attempting to translate a list',
     &           ' an unknown error'/14x,'was encountered. the ',
     &           'current list will be discarded.'/)
      go to 9999
 244  continue
      write(out,3023)
 3023 format(/1x,'>>>>> error: while attempting to translate a list',
     &           ' the list was not'/14x,'found. the current entity ',
     &           'will be discarded.'/)
      go to 9999
c
c
 250  continue
      num_error = num_error + 1
      write(out,9024) mxmat
 9024 format(/1x,'>>>>> error: there is no room for the current mat',
     &           'erial being input.'/14x,'the number of materials ',
     &           'already stored equals ',i6,' :'/14x,'the maximum ',
     &           'number of materials allowed in the problem.'/)
      go to 9999
c
c
 260  continue
      num_error = num_error + 1
      write(out,9025)
 9025 format(/1x,'>>>>> error: the element type specified cannot be',
     &           ' found in the element'/14x,'library. the type must',
     &           ' exist before it can be used. the'/14x,'current ',
     &           'element list will be discarded and a new card ',
     &           'sought.'/)
      go to 9999
c
c
 270  continue
      num_error = num_error + 1
      write(out,9026)
 9026 format(/1x,'>>>>> error: an element type for the current elem',
     &           'ent list is expected.'/14x,'discard the list and ',
     &           'the default type and seek a new card.'/)
      go to 9999
c
c
 280  continue
      num_warn = num_warn + 1
      write(out,9027)
 9027 format(/1x,'>>>>> warning: the order of integration for the ',
     &           'current element list'/16x,'is expected. the defa',
     &           'ult value for the element type'/16x,'of the list ',
     &           'will be used.'/)
      go to 9999
c
c
 290  continue
      num_error = num_error + 1
      write(out,9028)
 9028 format(/1x,'>>>>> error: the material specified cannot be fou',
     &           'nd in the material'/14x,'library. the material ',
     &           'must be defined before it can be'/14x,'used. the',
     &           'current element list will be discarded and a'/14x,
     &           'new card sought.'/)
      go to 9999
c
c
 300  continue
      num_error = num_error + 1
      write(out,9029)
 9029 format(/1x,'>>>>> error: the name of the material is expected.',
     &           ' discard the current'/14x,'element list and the ',
     &           'default material and seek a new card.'/)
      go to 9999
c
c
 310  continue
      num_error = num_error + 1
      write(out,9030)
 9030 format(/1x,'>>>>> error: the material type for the current ',
     &           'element list has been'/14x,'omitted. the current ',
     &           'element list will be discarded and'/14x,'a new ',
     &           'card sought.'/)
      go to 9999
c
c
 320  continue
      num_error = num_error + 1
      write(out,9031)
 9031 format(/1x,'>>>>> error: an element type for the current elem',
     &           'ent list is expected.'/14x,'discard the list and ',
     &           'seek a new card.'/)
      go to 9999
c
c
 330  continue
      num_warn = num_warn + 1
      write(out,9032) param
 9032 format(/1x,'>>>>> warning: the order of integration input has ',
     &           'no meaning. the'/16x,'default value for element ',
     &           i7,' will be used.'/)
      go to 9999
c
c
 340  continue
      num_fatal = num_fatal + 1
      write(out,9033) param
 9033 format(/1x,'>>>>> fatal error: properties for element ',i7,
     &           ' have not been input.'/14x,'execution cannot ',
     &           'proceed without this information.'/)
      input_ok = .false.
      go to 9999
c
c
 350  continue
      num_error = num_error + 1
      write(out,9034) param,noelem
 9034 format(/1x,'>>>>> error: ',i7,': the current element number ',
     &           ' exceeds ',i7,':'/14x,'the number of elements in ',
     &           'the structure. it will be ignored.'/)
      go to 9999
c
c
 360  continue
      num_warn = num_warn + 1
      write(out,9035)
 9035 format(/1x,'>>>>> warning: the new name of the material is ',
     &           'expected. input will'/16x,' continue using the ',
     &           'old material name.'/)
      go to 9999
c
c
 370  continue
      num_warn = num_warn + 1
      count= param/two16
      nnode= param-count*two16
      write(out,9036) count,nnode
 9036 format(/1x,'>>>>> warning: ',i7,': the current number of nodes',
     &           ' input exceeds ',i7/16x,': the number of element ',
     &           'nodes. incidence input will'/16x,'terminate.'/)
      go to 9999
c
c
 380  continue
      num_error = num_error + 1
      count= param/two16
      nnode= param-count*two16
      write(out,9037) count,nnode
 9037 format(/1x,'>>>>> error: ',i7,': the number of nodes input is ',
     &           'less than ',i7,': the'/16x,'number of element nodes.',
     &           ' there are holes in the incidences'/16x,'for this',
     &           ' element'/)
      input_ok = .false.
      go to 9999
c
c
 390  continue
      num_warn = num_warn + 1
      write(out,9038)
 9038 format(/1x,'>>>>> warning: an integer indicating the current ',
     &           'element being input'/16x,'is expected. the remainder',
     &           ' of this line will be'/16x,'ignored and another ',
     &           'element will be sought.'/)
      go to 9999
c
c
 400  continue
      num_fatal = num_fatal + 1
      write(out,9039) param
 9039 format(/1x,'>>>>> fatal error: incidences have not been input ',
     &           'for element ',i7,'.'/)
      input_ok = .false.
      go to 9999
c
c
 410  continue
      num_fatal = num_fatal + 1
      write(out,9040) param
 9040 format(/1x,'>>>>> fatal error: the incidences for element ',i7,
     &           ' contain holes(zeros).'/)
      input_ok = .false.
      go to 9999
c
c
 420  continue
      num_error = num_error + 1
      write(out,9041)
 9041 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'checking the incidences',
     & /,16x,'job terminated....'/)
      go to 9999
c
c
 430  continue
      num_warn = num_warn + 1
      write(out,9042)
 9042 format(/1x,'>>>>> warning: an integer indicating the current ',
     &           'block being input'/16x,'is expected. the remainder',
     &           ' of this line will be'/16x,'ignored and another ',
     &           'block will be sought.'/)
      go to 9999
c
c
 440  continue
      num_warn = num_warn + 1
      write(out,9043) stname
 9043 format(/1x,'>>>>> warning: the name of the structure has not ',
     &           'been input. the name'/16x,'of the structure is ',
     &           'currently ', a8,'.'/)
      go to 9999
c
c
 450  continue
      num_error = num_error + 1
      write(out,9044)
 9044 format(/1x,'>>>>> error: expecting type of information to be ',
     &           'displayed. display'/14x,'command will be ignored.'/)
      go to 9999
c
c
 460  continue
      num_error = num_error + 1
      write(out,9045)
 9045 format(/1x,'>>>>> error: expecting valid information for displ',
     &           'ay. keep searching'/14x,'for such.'/)
      go to 9999
c
c
 470  continue
      num_error = num_error + 1
      write(out,9046)
 9046 format(/1x,'>>>>> error: element properties and/or incidences ',
     &           'have not been'/14x,'input correctly. this must be ',
     &           'done before display'/14x,'of the elements.'/)
      go to 9999
c
c
 480  continue
      num_error = num_error + 1
      write(out,9047)
 9047 format(/1x,'>>>>> error: element incidences and/or coordinates ',
     &           'have not'/14x,'been input correctly. this must be ',
     &           'done before display'/14x,'of the nodes.'/)
      go to 9999
c
c
 490  continue
      num_error = num_error + 1
      write(out,9048) param,nonode
 9048 format(/1x,'>>>>> error: ',i7,': the current node number ',
     &           'exceeds ',i7,': the'/16x,'number of struct',
     &           'ure nodes. this node will not be dis-'/16x,'played.'/)
      go to 9999
c
c
 500  continue
      num_error = num_error + 1
      write(out,9049)
 9049 format(/1x,'>>>>> error: the material specified cannot be fou',
     &           'nd in the material'/14x,'library. the material ',
     &           'must be defined before it can be'/14x,'displayed'/)
      go to 9999
c
c
 510  continue
      num_error = num_error + 1
      write(out,9050)
 9050 format(/1x,'>>>>> error: the name of the material is expected.',
     &           ' there are no generic'/14x,'material names for dis-',
     &           'play.'/)
      go to 9999
c
c
 520  continue
      num_warn = num_warn + 1
      write(out,9051) stname
 9051 format(/1x,'>>>>> warning: the name of the structure is ',a8,
     &           '. this structure'/16x,'will be displayed.'/)
      go to 9999
c
c
 530  continue
      num_error = num_error + 1
      write(out,9052) mxlc
 9052 format(/1x,'>>>>> error: there is no room for the current load',
     &           'ing being input.'/14x,'the number of loadings alre',
     &           'ady stored equals ',i6,':'/14x,'the maximum number',
     &           ' of loadings allowed in the problem.'/)
      go to 9999
c
c
 540  continue
      num_error = num_error + 1
      write(out,9053)
 9053 format(/1x,'>>>>> error: the name of the loading is expected. ',
     &           'disregard the present'/,14x,'command or request.'/)
      go to 9999
c
c
 550  continue
      num_error = num_error + 1
      go to (551,552) param
c
 551  shrtst= 'nonlinr '
      go to 553
c
 552  shrtst= 'combnatn'
      go to 553
c
 553  write(out,9054) shrtst
 9054 format(/1x,'>>>>> error: attempt to input ',a8,' loading when ',
     &           'have already input'/14x,'nodal and/or element load',
     &           'ing. look for new loading mode'/14x,'command.'/)
      go to 9999
c
c
 560  continue
      num_error = num_error + 1
      write(out,9055)
 9055 format(/1x,'>>>>> error: expecting command for loading mode. ',
     &           'keep searching for one.'/)
      go to 9999
c
c
 570  continue
      num_error = num_error + 1
      write(out,9056)
 9056 format(/1x,'>>>>> error: expecting real number.',
     &           ' read a '/14x,'new line.'/)
      go to 9999
c
c
 580  continue
      num_error = num_error + 1
      write(out,9057) param
 9057 format(/1x,'>>>>> error: ',i7,': the current node is less than',
     &           ' zero. it will be'/14x,'ignored.'/)
      go to 9999
c
c
 590  continue
      num_error = num_error + 1
      write(out,9058)
 9058 format(/1x,'>>>>> error: expecting loading name . keep search',
     &           'ing.'/)
      go to 9999
c
c
 600  continue
      num_fatal = num_fatal + 1
      write(out,9059)
 9059 format(/1x,'>>>>> fatal: the loading specified cannot be found',
     &           ' in the loadings'/14x,'library. the loading must ',
     &           'be defined before it can be'/14x,'used.'/)
      go to 9999
c
c
 610  continue
      if ( .not. mess_61 ) return
      num_error = num_error + 1
      mess_61 = .false.
      write(out,9060) param
 9060 format(/1x,'>>>>> error: the maximum number of steps has been',
     & ' exceeded.',
     & /14x,'Step limit: ',i8,
     & /14x,'Subsequent messages of this type omitted...'/)
      go to 9999
c
c
 620  continue
      num_error = num_error + 1
      write(out,9061)
 9061 format(/1x,'>>>>> error: name of specified loading is the same',
     &           ' as the name of the'/14x,'combination loading. sea',
     &           'rch for new loading.'/)
      go to 9999
c
c
 630  continue
      write(out,9062)
 9062 format(/1x,
     & '>>>> WARNING: turning off crack growth suspends the force',
     & /,'               reduction process for elements being killed or',
     & /,'               nodes being released. This may not be what you',
     & /,'               want to happen....',/)
      go to 9999
c
c
 640  continue
      num_error = num_error + 1
      write(out,9063)
 9063 format(/1x,'>>>>> error: there is an error in the re-entry int',
     &           'o  subroutine inlod.'/14x,'inlod will not be re-',
     &           'entered.'/)
      go to 9999
c
c
 650  continue
      num_error = num_error + 1
      write(out,9064)
 9064 format(/1x,'>>>>> error: an end of card has been reached with',
     &           'out encountering a'/14x,'loading mode command. dis',
     &           'regard current loading name and'/14x,'search for a',
     &           ' high level command.'/)
      go to 9999
c
c
 660  continue
      num_error = num_error + 1
      write(out,9065)
 9065 format(/1x,'>>>>> error: type of element load expected.',
     &           ' keep searching.'/)
      go to 9999
c
c
 670  continue
      num_error = num_error + 1
      write(out,9066)
 9066 format(/1x,'>>>>> error: expecting real number denoting multip',
     &           'lication factor. find'/14x,'new loading.'/)
      go to 9999
c
c
 680  continue
      num_error = num_error + 1
      write(out,9067)
 9067 format(/1x,'>>>>> error: the name of the loading is expected.',
     &           ' there are no generic'/14x,'loading names for dis-',
     &           'play.'/)
      go to 9999
c
c
 690  continue
      num_error = num_error + 1
      write(out,9068) param
 9068 format(/1x,'>>>>> error: ',i7,': the current step is less than',
     &           ' zero. it will be'/14x,'ignored.'/)
      go to 9999
c
c
 700  continue
      num_error = num_error + 1
      write(out,9069) param
 9069 format(/1x,'>>>>> error: there is no transformation matrix for',
     &           ' node ',i7,'.'/)
      go to 9999
c
c
 710  continue
      num_error = num_error + 1
      write(out,9070)
 9070 format(/1x,'>>>>> error: there is a syntax error in constraints',
     &           ' input. read a'/14x,'new line'/)
      go to 9999
c
c
 720  continue
      num_error = num_error + 1
      write(out,9071)
 9071 format(/1x,'>>>>> error: expecting type of dof for constraint ',
     &           'or initial cond.'/14x,'input. ignore entity and ',
     &           'keep searching.'/)
      go to 9999
c
c
 730  continue
      num_error = num_error + 1
      write(out,9072)
 9072 format(/1x,'>>>>> error: expecting number constraint or initi',
     &           'al cond. value.'/14x,'ignore entity and search for',
     &           ' new dof type.'/)
      go to 9999
c
c
 740  continue
c               available
      go to 9999
c
c
 750  continue
      num_error = num_error + 1
      write(out,9073)
 9073 format(/1x,'>>>>> error: area must be >0 for bar element',/)
      go to 9999
c
c
 760  continue
      num_error = num_error + 1
      write(out,9074) param
 9074 format(/1x,'>>>>> error: ',i6,': the current block is less ',
     &           'than zero. it will be'/14x,'ignored.'/)
      go to 9999
c
c
 770  continue
      num_error = num_error + 1
      write(out,9175)
 9175 format(/1x,'>>>>> error: expecting the span of the block ',
     &           'being input.',/14x,' skip this line of in',
     &           'put and proceed to the next.'/)
       go to 9999
c
c
 780  continue
      num_error = num_error + 1
      write(out,9075)
 9075 format(/1x,'>>>>> error: the number of degrees of freedom ',
     &           'input for node ',i7/14x,' exceeds the maximum ',
     &           'number of degrees of freedom allowed'/14x,'any ',
     &           'node. the number of dof stored will be the above',
     &           /14x,'maximum.'/)
      go to 9999
c
c
 790  continue
      num_error = num_error + 1
      write(out,9076)
 9076 format(/1x,'>>>>> error: expecting key word matrix. this input',
     &           ' line will be ig-'/14x,'nored and a new line ',
     &           'sought.'/)
      go to 9999
c
c
 800  continue
      num_error = num_error + 1
      write(out,9077) param
 9077 format(/1x,'>>>>> error: there is a zero pointer to the trans-',
     &           'formation matrix'/14x,'for node ',i7,'. the given',
     &           ' constraints are invalid.'/)
      go to 9999
c
c
 810  continue
      num_error = num_error + 1
      write(out,9078) param
 9078 format(/1x,'>>>>> error: the transformation matrix for node ',i7,
     &           ' contains'/14x,'illegal data. the given constraints',
     &           ' are invalid.'/)
      go to 9999
c
c
 820  continue
      num_fatal = num_fatal + 1
      write(out,9079)
 9079 format(/1x,'>>>>> Fatal error: no constraints have been input.',
     &       'the solution cannot'/14x,'proceed without constraints.'/)
      input_ok = .false.
      go to 9999
c
c
 830  continue
      num_error = num_error + 1
      write(out,9080) param
 9080 format(/1x,'>>>>> error: the constraint value for dof ',i6,
     &           ' is illegal. the'/14x,'given constraints are ',
     &           'invalid.'/)
      go to 9999
c
c
 840  continue
      num_error = num_error + 1
      write(out,9081)
 9081 format(/1x,'>>>>> error: expecting keyword indicating which row',
     &           ' of the matrix is'/14x,'being input. keep search',
     &           'ing.'/)
      go to 9999
c
c
 850  continue
      num_error = num_error + 1
      write(out,9082)
 9082 format(/1x,'>>>>> error: expecting keyword area',/)
      go to 9999
c
c
 860  continue
      num_error = num_error + 1
      write(out,9083) param
 9083 format(/1x,'>>>>> error: ',i7,': the current element is less ',
     &           'than zero. it will be'/14x,'ignored.'/)
      go to 9999
c
c
 870  continue
      num_error = num_error + 1
      write(out,9084)
 9084 format(/1x,'>>>>> error: each link stiffness value must',
     & ' be >= 0',/)
      go to 9999
c
c
 880  continue
      num_error = num_error + 1
      shrtst= lodnam(param)
      write(out,9085) shrtst
 9085 format(/1x,'>>>>> error: loading ',a8,' does not consist of ',
     &           'actual loads. it'/14x,'is not applicable in this ',
     &           'context.'/)
      go to 9999
c
c
 890  continue
      num_error = num_error + 1
      shrtst= lodnam(param)
      write(out,9086) shrtst
 9086 format(/1x,'>>>>> error: the loading ',a8,' is a time step ',
     &           'definition. to display this'/7x,'loading the steps',
     &           'to be displayed must be specified.'/)
      go to 9999
c
c
 900  continue
      num_error = num_error + 1
      shrtst= lodnam(param)
      write(out,9087) shrtst
 9087 format(/1x,'>>>>> error: the loading ',a8,' is of unknown type',
     &           '. it will not be'/14x,'processed.'/)
      go to 9999
c
c
 910  continue
      num_error = num_error + 1
      write(out,9088)
 9088 format(/1x,'>>>>> error: incorrect syntax for dynamic analysis',
     &           ' parameters command.'/14x,'a new high level comman',
     &           'd will be sought.'/)
      go to 9999
c
c
 920  continue
      num_error = num_error + 1
      write(out,9089)
 9089 format(/1x,'>>>>> error: incorrect syntax for dynamic paramete',
     &           'rs commands. a new'/14x,'such command will be ',
     &           'sought.'/)
      go to 9999
c
c
 930  continue
      num_error = num_error + 1
      write(out,9090)
 9090 format(/1x,'>>>>> error: expecting switch on or off, or keyword',
     &           ' iteration. ignore'/14x,'command and search for a ',
     &           'new one.'/)
      go to 9999
c
c
 940  continue
      num_error = num_error + 1
      write(out,9091)
 9091 format(/1x,'>>>>> error: expecting keyword residual. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 950  continue
      num_error = num_error + 1
      write(out,9092)
 9092 format(/1x,'>>>>> error: expecting keyword updates. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 960  continue
      num_error = num_error + 1
      write(out,9093)
 9093 format(/1x,'>>>>> error: expecting keyword stiffness or mass.',
     &           ' ignore command'/14x,'and search for a new one.'/)
      go to 9999
c
c
 970  continue
      num_error = num_error + 1
      write(out,9094)
 9094 format(/1x,'>>>>> error: expecting real number. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 980  continue
      num_error = num_error + 1
      write(out,9095)
 9095 format(/1x,'>>>>> error: expecting keyword beta. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 990  continue
      num_error = num_error + 1
      write(out,9096)
 9096 format(/1x,'>>>>> error: expecting keyword step. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1000 continue
      go to 9999
c
c
 1010 continue
      num_error = num_error + 1
      write(out,9098)
 9098 format(/1x,'>>>>> error: expecting keyword stiffness. ignore ',
     &           'command and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1020  continue
      num_error = num_error + 1
      write(out,9099)
 9099 format(/1x,'>>>>> error: expecting stop or continue. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1030  continue
      num_error = num_error + 1
      write(out,9100)
 9100 format(/1x,'>>>>> error: expecting integer number. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1040  continue
      num_error = num_error + 1
      write(out,9101)
 9101 format(/1x,'>>>>> error: expecting keyword iterations. ignore ',
     &           'command and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1050  continue
      num_error = num_error + 1
      write(out,9102)
 9102 format(/1x,'>>>>> error: expecting residual or displacement. ',
     &           'ignore test type and'/14x,'search for a new one.'/)
      go to 9999
c
c
 1060  continue
      num_error = num_error + 1
      write(out,9103)
 9103 format(/1x,'>>>>> error: expecting norm or maximum. ignore ',
     &           'test type and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1070  continue
      num_error = num_error + 1
      write(out,9104)
 9104 format(/1x,'>>>>> error: expecting keyword solutions. ignore ',
     &           'command and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1080 continue
      num_error = num_error + 1
      write(out,9105)
 9105 format(/1x,'>>>>> error: the element properties for the struc',
     &           'ture have already'/14x,'been input. they cannot be',
     &           ' changed unless the problem is'/14x,'restarted.'/)
      go to 9999
c
c
 1090 continue
      num_warn = num_warn + 1
      write(out,9106) sparam
 9106 format(/1x,'>>>>> warning: dof ',a4,' has already been constr',
     &        'ained. this constraint'/14x,'request will be ignored'/)
      go to 9999
c
c
 1100 continue
      num_error = num_error + 1
      shrtst= lodnam(param)
      write(out,9107) shrtst
 9107 format(/1x,'>>>>> error: loading ',a8,' has already been defin',
     &           'ed. a different'/14x,'name must be used.'/)
      go to 9999
c
c
 1110 continue
      num_error = num_error + 1
c
      if(param.eq.1) then
         shrtst= ' nodes  '
      else if(param.eq.2) then
         shrtst= 'elements'
      else
         shrtst= 'cr nodes'
      end if
c
      write(out,9108) shrtst
 9108 format(/1x,'>>>>> error: the number of ',a8,' input is less th',
     &           'an or equal to zero.'/14x,'it will be ignored.'/)
      go to 9999
c
c
 1120 continue
      num_error = num_error + 1
      write(out,9109)
 9109 format(/1x,'>>>>> error: the element incidences for the struc',
     &           'ture have already'/14x,'been input. they cannot be',
     &           ' changed unless the problem is'/14x,'restarted.'/)
      go to 9999
c
c
 1130 continue
      num_error = num_error + 1
      write(out,9110)
 9110 format(/1x,'>>>>> error: the nodal coordinates for the struc',
     &           'ture have already'/14x,'been input. they cannot be',
     &           ' changed unless the problem is'/14x,'restarted.'/)
      go to 9999
c
c
 1140 continue
      num_error = num_error + 1
c
      if(param.eq.1) then
         shrtst= ' vector '
         shrtst1= ' matnam '
         shrtst2= 'material'
      else if(param.eq.2) then
         shrtst= ' vector '
         shrtst1= ' lodnam '
         shrtst2= ' loading'
      else if(param.eq.3) then
         shrtst= ' vector '
         shrtst1=' tabnam '
         shrtst2=' table '
      end if
c
      write(out,9111) shrtst,shrtst1,shrtst2
 9111 format(/1x,'>>>>> error: the opening in ',2a8,' slated ',
     &           'for the current'/14x,a8,'is out of bounds. it cann',
     &           'not be processed.'/)
      go to 9999
c
c
 1150 continue
      num_error = num_error + 1
      lngstr = matnam(param)
      write(out,9112) lngstr
 9112 format(/1x,'>>>>> error: material ',a24,'has already been defin',
     &           'ed, and cannot'/14x,'be overridden. it can be de-',
     &           'leted and redefined.'/)
      go to 9999
c
c
 1160 continue
      num_warn = num_warn + 1
      write(out,9113)
 9113 format(/1x,'>>>>> warning: the element properties for the stru',
     &       'cture have already '/16x,'been set. the material curre',
     &       'ntly being defined will not'/16x,'be used.'/)
      go to 9999
c
c
 1170 continue
      num_error = num_error + 1
      write(out,9114) sparam
 9114 format(/1x,'>>>>> error: dof ',a4,' has already been given an',
     &           'initial value.'/14x,'the current initial conditio',
     &           'n will be ignored.'/)
      go to 9999
c
c
 1180 continue
      num_error = num_error + 1
      write(out,9115)
 9115 format(/1x,'>>>>> error: the number of dof given exceeds the ',
     &           'number of dof for'/14x,'the node listed. only the ',
     &           'correct number of dof will be'/14x,'stored'/)
      go to 9999
c
c
 1190 continue
      num_error = num_error + 1
      write(out,9116)
 9116 format(/1x,'>>>>> error: there is a syntax error in initial ',
     &           'conditions input.'/14x,'read a new line'/)
      go to 9999
c
c
 1200 continue
      go to 9999
c
c
 1210 continue
      num_error = num_error + 1
      write(out,9118)
 9118 format(/1x,'>>>>> error: there is a syntax error in the compu',
     &           'tation request.'/14x,'find a new request.'/)
      go to 9999
c
c
 1220 continue
      num_error = num_error + 1
      shrtst= lodnam(param)
      write(out,9119) shrtst
 9119 format(/1x,'>>>>> error: loading ',a8,' is not a time step def',
     &           'inition. it'/14x,'is not applicable in this ',
     &           'context.'/)
      go to 9999
c
c
 1230 continue
      num_error = num_error + 1
      write(out,9120) param
 9120 format(/1x,'>>>>> error: step number ',i7,' has already been ',
     &           'computed. the next'/14x,'step in the list will be',
     &           ' processed.'/)
      go to 9999
c
c
 1240 continue
      num_error = num_error + 1
      write(out,9121) sparam
 9121 format(/1x,'>>>>> error: expecting a link stiffness value after',
     & /,     1x,'             link property: ',a,/)
      go to 9999
c
c
 1250 continue
      num_error = num_error + 1
      write(out,9122) param, max_step_limit
 9122 format(/1x,'>>>>> error: ',i76,',',' the step to be computed',
     &           ',',' is greater than'/14x,i7,',',' the maximum nu',
     &           'mber of steps defined. it will'/14x,'be ignored',
     &           ',',' and the next step in the list will be ',
     &           'processed.'/)
      go to 9999
c
c
 1260 continue
      num_error = num_error + 1
      iter= param/two16
      type= param-iter*two16
      if(type.eq.1) then
         shrtst= ' stiff  '
      else
         shrtst= '  mass  '
      end if
      write(out,9123) iter,shrtst
 9123 format(/1x,'>>>>> error: ',i6,',',' the iteration specified',
     &           ' for the next ',a8/14x,'update is greater than ',
     &           'the maximum number of iterations'/14x,'permitted ',
     &           'without convergence. it will be ignored.'/)
      go to 9999
c
c
 1270 continue
      num_error = num_error + 1
      iter= param/two16
      type= param-iter*two16
      if(type.eq.1) then
         shrtst= ' stiff  '
      else
         shrtst= '  mass  '
      end if
      write(out,9124) iter,shrtst
 9124 format(/1x,'>>>>> error: ',i6,',',' the iteration specified',
     &           ' for the next ',a8/14x,'update is less than zero.',
     &           ' it will be ignored.'/)
      go to 9999
c
c
 1280  continue
      num_error = num_error + 1
      write(out,9125)
 9125 format(/1x,'>>>>> error: expecting keyword gausspts. ignore ',
     &           'command and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1290 continue
      num_error = num_error + 1
      write(out,9126) param
 9126 format(/1x,'>>>>> error: ',i6,': the current gauss point is ',
     &           'less than zero. it'/14x,'will be ignored.'/)
      go to 9999
c
c
 1300 continue
      num_error = num_error + 1
      write(out,9127) param
 9127 format(/1x,'>>>>> error: data type ',a4,' is illegal for the ',
     &           'current element.'/14x,'the current element or node',
     &           ' will not be processed.'/)
      go to 9999
c
c
 1310 continue
      num_error = num_error + 1
      write(out,9128) param
 9128 format(/1x,'>>>>> error: data type ',a4,' was referenced for ',
     &           'the current element.'/14x,'it is not legal and ',
     &           'will not be processed.'/)
      go to 9999
c
c
 1320 continue
      num_error = num_error + 1
      write(out,9129)
 9129 format(/1x,'>>>>> error: expecting type of stress for initial ',
     &           'stress input.'/14x,'ignore entity and keep search',
     &           'ing.'/)
      go to 9999
c
c
 1330  continue
      num_error = num_error + 1
      write(out,9130)
 9130 format(/1x,'>>>>> error: expecting real number initial stress ',
     &           'value. ignore'/14x,'entity and search for new ',
     &           'stress type.'/)
      go to 9999
c
c
 1340 continue
      num_error = num_error + 1
      write(out,9131) param
 9131 format(/1x,'>>>>> error: stress ',a4,' has already been given ',
     &           'an initial value.'/14x,'the current initial condit',
     &           'ion will be ignored.'/)
      go to 9999
c
c
 1350 continue
      num_error = num_error + 1
      write(out,9132) param
 9132 format(/1x,'>>>>> error: row ',i6,' of the current transformat',
     &           'ion matrix has'/14x,'already been input. it will ',
     &           'not be overridden.'/)
      go to 9999
c
c
 1360 continue      
      num_warn = num_warn + 1
      if( count_136 == 0 ) write(out,8136)
      count_136 = 1
      write(out,9133) int(rparam), param
 9133 format(10x,2i10)
 8136 format(/1x,'>>>>> warning: these nodes and dof numbers already',
     &       ' have constraints defined:')
      go to 9999
c
c
 1370 continue
      num_error = num_error + 1
      write(out,9134) param
 9134 format(/1x,'>>>>> error: load type ',a4,' has already been ',
     &           'given a magnitiude.'/14x,'the current nodal load ',
     &           'definition will be ignored.'/)
      go to 9999
c
c
 1380 continue
      call entits( string, strlng )
      write(out,9135)string(1:strlng)
 9135 format(/1x,'>>>>> warning: a type of output is expected. ',/,
     & 14x,'scanning: ',a,
     & /,14x,'remainder of line flushed...'/)
      go to 9999
c
c
 1390 continue
      num_warn = num_warn + 1
      write(out,9136) param
 9136 format(/1x,'>>>>> warning: the label type for element ',i7,
     &           'does not match the'/16x,'existing label type for ',
     &           'output headers. the existing'/16x,'label type will',
     &           'changed.'/)
      go to 9999
c
 1400 continue
      num_error = num_error + 1
      write(out,9137)
 9137 format(/1x,'>>>>> error: mass for link element must be >=0 ',/)
      go to 9999
c
c
 1410 continue
      num_error = num_error + 1
c
      if(param.eq.1) then
         strng1= 'a stiffness matrix of any kind     '
      else if(param.eq.2) then
         strng1= 'a mass matrix of any kind          '
      else if(param.eq.3) then
         strng1= 'a stiffness/mass matrix of any kind'
      end if
c
      write(out,9138) strng1
 9138 format(/1x,'>>>>> error: the ',a35,' has not been '/14x,
     &           'computed yet and cannot be output.'/)
      go to 9999
c
c
 1420 continue
c
      write(out,9139)
 9139 format(/1x,'>>>>> error: the type of trace must be followed by ',
     &           'the keywords on or off.'/14x,'the trace type will ',
     &           'be ignored.'/)
      go to 9999
c
c
 1430 continue
      num_error = num_error + 1
c
      write(out,9140)
 9140 format(/1x,'>>>>> error: constraints and or transformation ',
     &           'matrices have not yet'/14x,'been specified. neith',
     &           'er can therefore be displayed.'/)
      go to 9999
c
c
 1440 continue
      num_error = num_error + 1
c
      if(param.eq.4hpcm ) then
         shrtst= 'dyn stfn'
      else
         shrtst= '  mass  '
      end if
c
      write(out,9141) shrtst
 9141 format(/1x,'>>>>> error: the ',a8,' matrix used as a pcm in the ',
     &           'calculations for the'/14x,' current step is not pos',
     &           'itive definite. calculations for'/14x,'this step',
     &           ' cannot proceed, and thus no further steps can'/14x,
     &           'be analyzed.'/)
      go to 9999
c
c
 1450 continue
      num_error = num_error + 1
      param= ltmstp+1
      write(out,9142) param
 9142 format(/1x,'>>>>> error: there has been an error in the solutio',
     &           'n of step ',i7/14x,'the remaining steps of the ',
     &           'current request cannot be'/14x,'performed, and ',
     &           'further requests will be ignored.'/)
      go to 9999
c
c
 1460 continue
      num_error = num_error + 1
      write(out,9143)
 9143 format(/1x,'>>>>> error: the name of the structure must be give',
     &           'n. request'/14x,'denied.')
      go to 9999
c
c
 1470 continue
      num_error = num_error + 1
      dbnum= param+10
      write(out,9144) dbnum,snames(param)
 9144 format(/1x,'>>>>> error: data base nld_db',i2,' does not exist ',
     &           'on disk.'/14x,'structure ',a8,' cannot be deleted.')
      go to 9999
c
c
 1480 continue
      num_error = num_error + 1
      write(out,9145) sparam
 9145 format(/1x,'>>>>> error: database: ',/14x,a,/7x,' does not ',
     &           'exist on disk.  structure cannot be reopened.',/14x,
     &           'request denied.')
      go to 9999
c
c
 1490 continue
      num_error = num_error + 1
      write(out,9146)
 9146 format(/1x,'>>>>> error: the specified structure cannot be ',
     &           'found in the structure'/14x,'library. request ',
     &           'denied.')
      go to 9999
c
c
 1500 continue
      num_error = num_error + 1
      write(out,9147) stname
 9147 format(/1x,'>>>>> error: the structure to be saved is not the',
     &           ' current structure.'/14x,'the current structure is ',
     &           a8,'.  request denied.')
      go to 9999
c
c
 1510 continue
      num_warn = num_warn + 1
      write(out,9148)
 9148 format(/1x,'>>>>> warning: no filename was specified, and the',
     &           ' structure',/14x,'has no name.  saving as',
     &           ' default_db in executing directory.')
      go to 9999
c
c
 1520 continue
      num_error = num_error + 1
      write(out,9149)
 9149 format(/1x,'>>>>> error in specifying preconditioner type.',
     &           'assuming diagonal.')
      go to 9999
c
c
 1530 continue
      num_error = num_error + 1
      write(out,9150)
 9150 format(/1x,'>>>>> error: expecting keyword one. ignore co',
     &           'mmand and search'/14x,'for a new one.'/)
      go to 9999
c
c
 1540 continue
      num_error = num_error + 1
      write(out,9151)
 9151 format(/1x,'>>>>> error: the element blocking for the struc',
     &           'ture has already'/14x,'been input. it cannot be',
     &           ' changed unless the problem is'/14x,'restarted.'/)
      go to 9999
c
c
 1550 continue
      num_error = num_error + 1
      write(out,9152)
 9152 format(/1x,'>>>>> error: accelerate command must be followed ',
     &           'by keywords on or off.'/14x,'accelerate command ',
     &           'will be ignored.'/)
      go to 9999
c
c
 1560 continue
      num_warn = num_warn + 1
      write(out,9153) sparam
 9153 format(/1x,'>>>>> warning: a real number is expected. ',
     &           'the default value for'/16x,a4,' will be used.'/)
      go to 9999
c
c
 1570 continue
      num_error = num_error + 1
      write(out,9154)
 9154 format(/1x,'>>>>> error: there has been an error in the solutio',
     &           'n of the'/14x,'accelerations for unconstrained ',
     &           'dof at time zero.'/14x,'accelerations at time ',
     &           'zero cannot be processed.'/)
      go to 9999
c
c
 1580 continue
      num_warn = num_warn + 1
      write(out,9158)
 9158 format(/1x,'>>>>> warning: the dot product of the current sear',
     &           'ch direction and'/16x,'the previous residual for ',
     &           'the preconditioned conjugate'/16x,'gradient algori',
     &           'thm is negative. the algorithm may not'/16x,'conv',
     &           'erge.'/)
      go to 9999
c
c
 1590 continue
      num_error = num_error + 1
      write(out,9159)
 9159 format(/1x,'>>>>> error: expecting the first element of ',
     &           'the block being',/14x,' input. skip this line of ',
     &           'input and proceed to the next.'/)
      go to 9999
c
c
 1600 continue
      num_warn = num_warn + 1
      write(out,9160) param
 9160 format(/1x,'>>>>> warning: the line search has failed to brac',
     &           'ket the root within'/16x,i6,' bracketing attempts.',
     &           ' the last estimate of the'/16x,'step length will ',
     &           'be used. the preconditioned conjugate'/16x,'gradi',
     &           'ent algorithm will probably not converge.'/)
      go to 9999
c
c
 1610 continue
      num_error = num_error + 1
      write(out,9173)
 9173 format(/1x,'>>>>> error: the element block data structure ',
     &           'is incomplete.'/)
      go to 9999
c
c
 1620 continue
      num_error = num_error + 1
      write(out,9161)
 9161 format(/1x,'>>>>> error: expecting keyword convergence. ',
     &           'the current '/14x,'command will be ',
     &           'ignored and a new dynamic'/14x,'analysis command ',
     &           'sought.'/)
      go to 9999
c
c
 1630 continue
      num_error = num_error + 1
      write(out,9162)
 9162 format(/1x,'>>>>> error: expecting line search parameter ',
     &           'keyword eta_#. the'/14x,'current entity will be ',
     &           'ignored and a new parameter'/14x,'keyword ',
     &           'sought.'/)
      go to 9999
c
c
 1640 continue
      num_error = num_error + 1
      write(out,9163)
 9163 format(/1x,'>>>>> error: expecting real number denoting ',
     &           'line search parameter.'/14x,'ignore current ',
     &           'entity and search for new parameters'/14x,
     &           'keyword.'/)
      go to 9999
c
c
 1650 continue
      num_error = num_error + 1
      write(out,9164)
 9164 format(/1x,'>>>>> error: expecting line search convergence ',
     &           'test keywords'/14x,'step_length or dot_product. ',
     &           'the current entity will'/14x,'be ignored and a ',
     &           'new convergence test keyword sought.'/)
      go to 9999
c
c
 1660 continue
      num_error = num_error + 1
      write(out,9165)
 9165 format(/1x,'>>>>> error: expecting real number denoting ',
     &           'line search convergence'/14x,'test tolerance. ',
     &           'ignore the current entity and search'/14x,'for ',
     &           'a new convergence test.'/)
      go to 9999
c
c
 1670 continue
      num_error = num_error + 1
      write(out,9166)
 9166 format(/1x,'>>>>> error: expecting keywords mnr or pcg. the',
     &           ' default solution'/14x,'algorithm will be used ',
     &           'and a new'/14x,'dynamic analysis command sought.'/)
      go to 9999
c
c
 1680 continue
      num_error = num_error + 1
      write(out,9167)
 9167 format(/1x,'>>>>> error: expecting keywords polak_ribiere or ',
     &           'sorenson. the'/14x,'default pcg beta parameter wi',
     &           'll be used and a new'/14x,'dynamic analysis comma',
     &           'nd sought.'/)
      go to 9999
c
c
 1690 continue
      num_error = num_error + 1
      write(out,9168)
 9168 format(/1x,'>>>>> error: expecting keywords dynamic or diagon',
     &           'al. the default'/14x,'preconditioner type will be',
     &           ' used and a new dynamic'/14x,'analysis command ',
     &           'sought.'/)
      go to 9999
c
c
 1700 continue
      num_error = num_error + 1
      write(out,9169)
 9169 format(/1x,'>>>>> error: expecting keywords: next step .',
     &       /14x,'command ignored.'/)
      go to 9999
c
c
 1710 continue
      go to 9999
c
c
 1720 continue
c
      write(out,9171)
 9171 format(/1x,'>>>>> error: expecting the type of trace to be per',
     &           'formed. continue'/14x,'searching until such a type',
     &           ' is found.'/)
      go to 9999
c
c
 1730 continue
      num_error = num_error + 1
c
      write(out,9174)
 9174 format(/1x,'>>>>> error: the sum of the block spans does not',
     &           ' equal the number'/14x,'of elements in the ',
     &           'structure.'/)
      go to 9999
c
c
 1740 continue
      num_error = num_error + 1
c
      write(out,9176)
 9176 format(/1x,'>>>>> error: extrapolate command must be followed ',
     &           'by keywords on or off.'/14x,'command ',
     &           'will be ignored.'/)
      go to 9999
c
c
 1750 continue
      num_error = num_error + 1
      write (out,9185)
 9185 format(/1x,'>>>>> error: name of file expected.  ignoring.',/,
     &           '             command.')
      goto 9999
c
c
 1760 continue
      num_error = num_error + 1
      write(out,9186) sparam
 9186 format(/1x,'>>>>> error: error in opening file.',a)
      goto 9999
c
c
 1770 continue
      num_error = num_error + 1
      write(out,9187) sparam
 9187 format(/1x,'>>>>> error: the file specified: ',a,/,
     &           '             does not exist.')
      goto 9999
c
c
 1780 continue
      num_error = num_error + 1
      write(out,9188)
 9188 format(/1x,'>>>>> error: too many input files open at once.',/,
     &           '             aborting attempt to read input file.')
c
c
 1790 continue
      write(out,9189)
 9189 format(/1x,'>>>>> reached end of file.')
      goto 9999
c
c
 1800 continue
      num_error = num_error + 1
      write(out,9190)
 9190 format(/1x,'>>>>> error: already writing output to file.',/,
     &           '             ignoring new outfile request.')
      goto 9999
c
c
 1810 continue
      num_error = num_error + 1
      write(out,9191)
 9191 format(/1x,'>>>>> error: the home environment variable',/,
     &           '             must be set to use a ~ in a filename')
      goto 9999
c
c
 1820 continue
      write(out,9192) rparam
 9192 format(1x,'>>>> wall time so far (secs): ',f20.7)
      goto 9999
c
c
 1830 continue
      num_error = num_error + 1
      write(out,9193)
 9193 format(/1x,'>>>>> error: fatal error in input was never fixed.',
     &         /,'              unable to continue.')
      goto 9999
c
c
 1840 continue
      write(out,9194)
 9194 format(/1x,' >> previous fatal errors are now being ignored.')
      goto 9999
c
c
 1850 continue
      num_error = num_error + 1
      write(out,9195) param, mxvl
 9195 format(/1x,'>>>>> error: the blocking span of ',i3,' is higher',
     &           'than the maximum, ',i3,'.',/,7x,' execution cannot',
     &           ' continue with this error.'/)
      input_ok = .false.
      goto 9999
c
c
 1860 continue
      num_warn = num_warn + 1
      write(out,9196)
 9196 format(/1x,'>>>>> warning: one or more of the principle axes',
     &           ' may not be constrained',/,7x,'in which case',
     &           ' rigid boby rotation or translation could occur.',
     &        /,7x,'This check does not account for possible effect',
     &        ' of skewed constraints.'/)
      goto 9999
c
c
 1870 continue
      num_error = num_error + 1
      write(out,9197) sparam
 9197 format(/1x,'>>>>> error: cannot access files from another users',
     &           ' directory.  filename invalid:',/,7x,a80)
      goto 9999
c
c
 1880 continue
      num_error = num_error + 1
      write(out,9198) sparam
 9198 format(/1x,'>>>>> error: filename invalid:',/,7x,a80)
      goto 9999
c
c
 1890 continue
      num_error = num_error + 1
      write(out,9199)
 9199 format(/1x,'>>>>> error: the home environment variable is not',
     &       ' set, thus could',/,7x,'not find the home directory.',
     &       '  filename invalid.')
      goto 9999
c
c
 1900 continue
      num_error = num_error + 1
      write(out,9200)
 9200 format(/1x,'>>>>> error: the wrkdir environment variable is not',
     &       ' set, thus could',/,7x,'not save the database as named.')
      goto 9999
c
c
 1910 continue
      num_warn = num_warn + 1
      write(out,9201)
 9201 format(/1x,'>>>>> warning: given filename was unacceptible. ',
     &       'saving database as',/7x,'default_db in executing ',
     &       'directory.')
      goto 9999
c
c
 1920 continue
      num_error = num_error + 1
      write(out,9201)
 9202 format(/1x,'>>>>> error: given filename was unacceptible.',
     &       '  cannot',/7x,'retrieve database.  command ignored.')
      goto 9999
c
c
 1930  continue
      write(out,9203) sparam
 9203 format(/1x,'>>>>> saved structure as ',a)
      goto 9999
c
c
 1940  continue
      write(out,9204) sparam,param,dparam,rparam
 9204 format(/1x,'>>>>> retrieved structure from ',a100,
     &      /1x,'>>> step recovered from restart file:',i7,
     &      /1x,'>>> analysis time so far            : ',e14.6,' sec',
     &      /1x,'>>> wall time now :                 ',f10.1,' sec',
     &      /)
      goto 9999
c
c
 1950   continue
      num_error = num_error + 1
      write(out,9205) param, sparam
 9205 format(/1x,'>>>>> error: doing another step would cause the',
     &           ' time to go beyond',/7x,'the alloted time. saving',
     &           ' structure at step ',i7,/7x,' as the following',
     &           ' name:',/7x,a80)
      goto 9999
c
c
 1960 continue
      num_error = num_error + 1
      write(out,9206)
 9206 format(/1x,'>>>>> error: time limit quantity missing.')
      goto 9999
c
c
 1970 continue
      write(out,9207) rparam
 9207 format(/1x,'>>>>> time limit set to ',f10.3,' seconds.')
      goto 9999
c
c
 1980 continue
      num_error = num_error + 1
      write(out,9208)
 9208 format(/1x,'>>>>> error: porosity value must be a real number.')
      goto 9999
c
c
 1990 continue
      num_error = num_error + 1
      write(out,9209)
 9209 format(/1x,'>>>>> error: porosity value must be greater than ',
     &     'zero.',/,14x,'defaulting to zero for critical porosity',
     &     ' value.')
      goto 9999
c
c
 2000 continue
      num_error = num_error + 1
      write(out,9210)
 9210 format(/1x,'>>>>> error: crack growth type unrecognizable.',
     &     /,14x,'defaulting to no crack growth.')
      goto 9999
c
c
 2010 continue
      num_error = num_error + 1
      write(out,9211)
 9211 format(/1x,'>>>>> error: incorrect syntax for crack growth',
     &           'parameter command. a new'/14x,'such command will be ',
     &           'sought.'/)
      go to 9999
c
 2020 continue
      num_error = num_error + 1
      write(out,9212) param
 9212 format(/1x,'>>>>> error: number of release steps must be ',
     &           'greater than zero. the'/14x,'default value of ',
     &           i2,' will be used.')
      go to 9999
c
c
 2030 continue
      num_error = num_error + 1
      write(out,9213) param
 9213 format(/1x,'>>>>> error: the number of release steps must be ',
     &           'an integer. the'/14x,'default value of ',
     &           i2,' will be used.')
      goto 9999
c
c
 2040 continue
      num_error = num_error + 1
      write(out,9214)
 9214 format(/1x,'>>>>> error: the end of the input file has ',
     &           'been reached. stopping warp.')
      goto 9999
c
c
 2050 continue
      num_error = num_error + 1
      write(out,9215)
 9215 format(/1x,'>>>>> error:  incorrect syntax for domain',
     &           ' definition ')
      goto 9999
c
c
 2060 continue
      num_error = num_error + 1
      write(out,9216)
 9216 format(/1x,'>>>>> error:  star command not found.')
      goto 9999
c
c
 2070 continue
      write(out,9217)
 9217 format(/1x,'>>>>> warning:  changing solver type not permitted.',
     &           ' command ignored ...',//)
      go to 9999
c
c
 2080 continue
      num_warn = num_warn + 1
      write(out,9218)
 9218 format(/1x,'>>>>> warning:  none of the material types were ',
     &           'specified to be killable.',/14x,'crack growth ',
     &           'cannot occur and is now disabled.')
      goto 9999
c
 2090 continue
      num_warn = num_warn + 1
      write(out,9219)
 9219 format(/1x,'>>>>> warning:  the current loading has been ',
     &           'previously defined.',
     &       /1x,'                this definition supercedes all',
     &                          ' information from',
     &       /1x,'                the previous definition.',/)
      goto 9999
c
 2100 continue
      num_error = num_error + 1
      write(out,9220)
 9220 format(/1x,'>>>>> error: dof must be u, v, w in mpc',/ )
      goto 9999
c
 2110 continue
      num_warn = num_warn + 1
      call entits( string, strlng )
      write(out,9221) string(1:strlng)
 9221 format(/1x,'>>>>> warning:  invalid thermal expansion data item:',
     &       1x, a )
      go to 9999
c
 2120 continue
      num_error = num_error + 1
      write(out,9222)
 9222 format(/1x,'>>>>> error: the number of release steps cannot',
     &          ' be changed after',/,
     &          '              elements have been killed.  ignoring',
     &          ' change request.',/)
      goto 9999
c
 2130 continue
      write(out,9223)
 9223 format(/1x,'>>>>> note: due to the absence of an element list,',
     &          ' all killable elements',/,
     &          '              are assumed.',/)
      goto 9999
c
 2140 continue
      num_error = num_error + 1
      write(out,9224)
 9224 format(/1x,'>>>>> error: gurson crack growth cannot be turned',
     &          ' back on after some',/,
     &          '              elements have been killed previously.',
     &          ' leaving crack growth off.',/)
      goto 9999
c
 2150 continue
      num_error = num_error + 1
      write(out,9225) sparam
 9225 format(/1x,'>>>>> error: command has already been specified',
     &          ' once -- cannot be ',/,
     &          '              specifed again in this run.',
     &          ' variable ',a20,/,' cannot be reallocated.',/)
      goto 9999
c
 2160 continue
      num_fatal = num_fatal + 1
      write(out,9226) param
 9226 format(/1x,'>>>>> FATAL ERROR: there is an error in the',
     &          ' blocking for this model.',
     &   /,      '                    block:',i4,
     &  ' for the following reason(s):' )
      goto 9999
c
 2170 continue
      num_warn = num_warn + 1
      write(out,9227) sparam
 9227 format(/1x,'>>>>> warning: an  integer number is ',
     &           'expected. the default'/16x,'value for ',a5,
     &           ' will be used.' )
      goto 9999
c
 2180 continue
      num_warn = num_warn + 1
      write(out,9228)
 9228 format(/1x,'>>>>> warning: expected the keyword  : curve ',
     &          '            skipping rest of curve definition.')
      goto 9999
c
 2190 continue
      num_warn = num_warn + 1
      write(out,9229) param
 9229 format(/1x,'>>>>> warning: unpaired strain-stress values.',
     &        /,'                ignoring current entry and scanning',
     &        /,'                for next stress-strain point.',
     &        /,'                have processed ',i3,' points so far.',
     &        /)
      goto 9999
c
 2200 continue
      num_warn = num_warn + 1
      write(out,9230) param
 9230 format(/1x,'>>>>> warning: exceeded max number of stress-strain',
     &          ' values : ',i5,/,
     &          '                stopping curve definition at this ,',
     &          'point.',/)
        goto 9999
c
 2210 continue
      num_warn = num_warn + 1
      write(out,9231)
 9231 format(/1x,'>>>>> warning: error in stress-strain input.',
     &        /,'                ignoring curve definition.',/)
       goto 9999
c
 2220 continue
      num_error = num_error + 1
      write (out,9232)
 9232 format(/1x,'>>>>> error in stress-strain input. scanning for ',
     &          'next line.',/)
      goto 9999
c
 2230 continue
      num_fatal = num_fatal + 1
      write (out,9233)
 9233 format(/1x,'>>>>> FATAL Error: set of stress-strain curves does',
     &/,         '                   not exist. set: ',i3,
     &        /)
      input_ok = .false.
      goto 9999
c
 2240 continue
      write(out,9234)
 9234 format(/1x,'>>>>> WARNING: all existing constraint data and local',
     &       /1x,'               coordinate systems at nodes now',
     &       ' deleted...'/ )
       goto 9999
c
 2250 continue
      num_warn = num_warn + 1
      write(out,9235)
 9235 format(/1x,'>>>>>  Warning: none of the specified elements',
     &          ' are killable,',/,
     &          '                 so no information about these',
     &          ' elements will be',/,
     &          '                    printed.',/)
      goto 9999
c
 2260 continue
      num_error = num_error + 1
      write (out,9236)
 9236 format(/1x,'>>>>> Error: expecting face number. Skipping ',
     &          ' rest of line.',/)
      goto 9999
c
 2270 continue
      num_error = num_error + 1
      write (out,9237) param
 9237 format(/1x,'>>>>> Error: face ',i2,', is invalid. Skipping',
     &          ' rest of line.',/)
      goto 9999
c
 2280 continue
      num_error = num_error + 1
      write (out,9238) param
 9238 format(/1x,'>>>>> Error: the pressure for face ',i1,' has',
     &          ' already been given on',/,
     &          '                the same line.  The new pressure',
     &          ' value will be ignored.',/)
      goto 9999
c
 2290 continue
      num_fatal = num_fatal + 1
      element = param/10
      face = param - element*10
      write (out,9239) face,element
 9239 format(/1x,'>>>>> Fatal Error: in processing the pressure',
     &          ' loading for face ',i1,/,
     &          '                of element ',i7,', at least one of',
     &          ' the face normals was found to',/,
     &          '               have a zero length.',
     &          '  Execution cannot continue.',/)
      goto 9999
c
 2300 continue
      num_warn = num_warn + 1
      write (out,9240) param
 9240 format (/1x,' >>>>> Warning: no loads, no non-zero ',
     &          ' constraints and no',/,
     &          '         inertial forces were applied to the',
     &          ' structure for',/,
     &          '         step ',i7,'.  Skipping to next step.',/)
      goto 9999
c
 2310 continue
      num_fatal = num_fatal + 1
      write (out,9241)
 9241 format (/1x,' >>>>> Fatal Error:  the restart file is either not',
     &          ' compatible with this',/,
     &          '        version of WARP3D or it has been corrupted.',
     &          '  Cannot restart.',/,
     &          '        Execution terminated.',/)
      goto 9999
c
 2320 continue
      num_warn = num_warn + 1
      write (out,9250)
 9250 format (/1x,' >>>>>  Warning:  Specified element output of a ',
     &          'nodal quantity',/,
     &          '        (e.g. displacement, velocity ). ',
     &          '  This quantity shall be output as',/,
     &          '        a nodal value.',/)
      goto 9999
c
 2330 continue
      num_warn = num_warn + 1
      write (out,9251)
 9251 format (/1x,' >>>>>  Warning:  first strain point of a specified',
     &       /,1x,'        stress-strain curve must not be zero.',
     &       /,1x,'        curve definition deleted'/)
      goto 9999
c
 2340 continue
      num_warn = num_warn + 1
      write (out,9252)
 9252 format (/1x,' >>>>>  Warning:  unrecognized keyword.',
     &        /1x,'        steps or traction-separation expected.',
     &        /1x,'        command ignored.',/)
      goto 9999
c
 2350 continue
      num_warn = num_warn + 1
      write (out,9253)
 9253 format (/1x,' >>>>>  Warning:  no valid crack plane',
     &        /1x,'                  normal given.',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2360 continue
      num_warn = num_warn + 1
      write (out,9254)
 9254 format (/1x,' >>>>>  Warning:  expecting coordinate of',
     &        /1x,'                  crack plane.',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2370 continue
      num_warn = num_warn + 1
      write (out,9255)
 9255 format (/1x,' >>>>>  Warning:  expecting crack plane',
     &        /1x,'                  normal description.',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2380 continue
      num_warn = num_warn + 1
      write (out,9256)
 9256 format (/1x,' >>>>>  Warning:  expecting value of cell height',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2390 continue
      num_warn = num_warn + 1
      write (out,9257)
 9257 format (/1x,' >>>>>  Warning:  expecting value of release',
     &        /1x,'                  fraction.',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2400 continue
      num_warn = num_warn + 1
      write (out,9258)
 9258 format (/1x,' >>>>>  Warning:  expecting keyword steps',
     &        /1x,'                  or fraction.',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2410 continue
      num_fatal = num_fatal + 1
      write (out,9259)
      input_ok = .false.
 9259 format (/1x,' >>>>> Fatal Error:  no valid direction defined',
     &        /1x,'                     for crack plane normal.'/)
      goto 9999
c
 2420 continue
      num_fatal = num_fatal + 1
      write (out,9260)
      input_ok = .false.
 9260 format (/1x,' >>>>> Fatal Error:  no valid cell size',
     &        /1x,'                     specified.'/)
      goto 9999
c
 2430 continue
      num_fatal = num_fatal + 1
      if (param.eq.1) then
         strng21 = 'release fraction     '
      else if (param.eq.2) then
         strng21 = 'characteristic length'
      endif
      write (out,9261) strng21
      input_ok = .false.
 9261 format (/1x,' >>>>> Fatal Error:  no valid ',a21,
     &        /1x,'                     specified.'/)
      goto 9999
c
 2440 continue
      num_warn = num_warn + 1
      write (out,9262)
 9262 format (/1x,' >>>>>  Warning:  no list found, all',
     &        /1x,'                  nodes/elements will be output')
      goto 9999
c
 2450 continue
      num_warn = num_warn + 1
      write (out,9263)
 9263 format (/1x,' >>>>>  Warning:  expecting value of cell height',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
c
 2460 continue
      num_fatal = num_fatal + 1
      if ( param .eq. 1) then
         strng29 = 'critical angle for release   '
      else if ( param .eq. 2) then
         strng29 = 'critical angle for initiation'
      endif
      write (out,9264) strng29
      input_ok = .false.
 9264 format (/1x,' >>>>> Fatal Error:  no ',a29,
     &        /1x,'                     specified.'/)
      goto 9999
c
 2470 continue
      num_fatal = num_fatal + 1
      write (out,9265)
      input_ok = .false.
 9265 format (/1x,' >>>>> Fatal Error: no crack plane nodes were ',
     &        /1x,'                    found. Node release crack',
     &        /1x,'                    growth cannot execute.',/)
      goto 9999
c
 2480 continue
      num_warn = num_warn + 1
      write (out,9266) param
 9266 format (/1x,' >>>>> Warning: a new constraint was specified',
     &        /1x,'                that constrained node ',i7,',',
     &        /1x,'                which had been previously released',
     &        /1x,'                by the crack growth routines.',
     &        /1x,'                The new constraint will be',
     &        /1x,'                removed.',/)
      goto 9999
c
 2490 continue
      num_error = num_error + 1
      write (out,9267)
 9267 format (/1x,' >>>>> Error: Element Extinction and Node Release',
     &        /1x,'               crack growth cannot both operate in',
     &        /1x,'               a single analysis. The new request',
     &        /1x,'               for crack growth will be ignored.',
     &        /)
      goto 9999
c
 2500 continue
      num_fatal = num_fatal + 1
      write (out,9268)
      input_ok = .false.
 9268 format (/1x,' >>>>> Fatal Error: no crack front was found on',
     &        /1x,'                 the crack plane. Either no ',
     &        /1x,'                 initial crack was provided, or ',
     &        /1x,'                 the crack plane is not a',
     &        /1x,'                 symmetry plane.',/)
      goto 9999
c
 2510 continue
      num_error = num_error + 1
      write (out,9269)
 9269 format (/1x,' >>>>> Error: The crack plane cannot be redefined',
     &        /1x,'               after node release or element',
     &        /1x,'               death has occurred.  Ignoring new',
     &        /1x,'               crack plane definition.',
     &        /)
      goto 9999
c
 2520 continue
      num_fatal = num_fatal + 1
      write (out,9270) param, mxconn
      input_ok = .false.
 9270 format (/1x,' >>>>> Fatal Error: the number of elements',
     &        /1x,'                 connected to node ',i7,' exceeds',
     &        /1x,'                 the maximum number of elements',
     &        /1x,'                 that may be connected to a single',
     &        /1x,'                 node, which is ',i3,'.',
     &        /)
      goto 9999
c
 2530 continue
      num_warn = num_warn + 1
      write (out,9271) param, dparam
 9271 format (/1x,' >>>>>  Warning:  The released node ',i7,' has',
     &        /1x,'             travelled opposite to the direction',
     &        /1x,'             of the normal to the crack plane,',
     &        /1x,'             either closing the crack or travelling',
     &        /1x,'             past the crack plane in the wrong',
     &        /1x,'             direction. Current displacement of',
     &        /1x,'             the node in the direction of the',
     &        /1x,'             normal to the crack plane:',e13.6,
     &        /)
      goto 9999
c
 2540 continue
      num_warn = num_warn + 1
      write (out,9272)
 9272 format (/1x,' >>>>>  Warning: expecting value of characteristic',
     &        /1x,'                 length.',
     &        /1x,'                 command ignored.',/)
      goto 9999
c
c
 2550 continue
      if( write_msg_255 ) then
        write(out,9273)
        write_msg_255 = .false.
      else
        write(out,9900)
      end if
 9273 format (/1x,' >>>>>  Warning (type 100):',
     & /,12x,'Loading patterns + constraints (w/ their multipliers)',
     & /,12x,'defined for this step create an incremental loading ',
     & /,12x,'vector not a simple multiple of the prior step vector. ',
     & /,12x,'Extrapolation of displacements to start a step works ',
     & /,12x,'best with proportional, incremental loads between two ',
     & /,12x,'steps. Non-proportional loads w/ extrapolation may lead',
     & /,12x,'to convergence issues.',
     & ' To improve convergence: ', /,
     & 14x,'-> Extrapolation cancelled for only this load step',/,
     & 14x,'-> This also forces use of the linear [D]s at',
     &  ' start of step',/)
 9900 format(7x,'>> Warning: see previous type 100 message')
      goto 9999
c
 2560 continue
      num_error = num_error + 1
      write(out,9274)
 9274 format(/1x,'>>>>> error: the data for the traction-separation',
     &          ' law cannot be ',
     &       /1x,'             changed after nodes have been',
     &          ' released. ignoring',
     &       /1x,'             change request.',/)
      goto 9999
c
 2570 continue
      goto 9999
c
 2580 continue
      num_fatal = num_fatal + 1
      lngstr = matnam(param)
      write(out,9276) lngstr
 9276 format(/1x,'>>>>> Fatal error in material ',a24,':',/
     &           '         either tan_e or n_power must be specified',/,
     &           '         and be greater than zero for this material.',
     &           /)
      input_ok = .false.
      goto 9999
c
 2590 continue
      num_fatal = num_fatal + 1
      write(out,9277)
 9277 format(/1x,'>>>>> Fatal error in crack tip definition. ',
     &       'Node number exceeds number of nodes in structure ',
     &           /)
      input_ok = .false.
      goto 9999
c
 2600 continue
      num_fatal = num_fatal + 1
      write(out,9278)
 9278 format(/1x,'>>>>> Fatal error in crack tip definition. ',
     &       'Non-positive node number given')
      input_ok = .false.
      goto 9999
c
 2610 continue
      num_fatal = num_fatal + 1
      write(out,9279)
 9279 format(/1x,'>>>>> Fatal error in crack tip definition. ',
     &       'Expected an integer node set id.')
      input_ok = .false.
      goto 9999
c
 2620 continue
      num_fatal = num_fatal + 1
      write(out,9280) int(rparam)
 9280 format(/1x,'>>>>> Fatal error in crack tip definition. ',
     &       'Exceeded maximum allowable crack-tip nodes'/
     &       'per node set. Maximum equals ',i6 )
      input_ok = .false.
      goto 9999
c
 2630 continue
      num_fatal = num_fatal + 1
      write(out,9281) int(rparam)
 9281 format(/1x,'>>>>> Fatal error in crack tip definition. ',
     &       'Exceeded maximum allowable node sets.'/
     &       'Maximum equals ',i6 )
      input_ok = .false.
      goto 9999
c
 2640 continue
      num_warn = num_warn + 1
      write(out,9282)
 9282 format(/1x,'>>>>> Warning: node sets contain different ',
     &       'numbers of nodes.')
      goto 9999
c
c
 2650 continue
      num_fatal = num_fatal + 1
      write(out,9283)
 9283 format(/1x,'>>>>> FATAL ERROR: overflow of table that stores',
     &       /1x,'                   transformation matrices for',
     &       /1x,'                   local constraint coordinate',
     &       /1x,'                   systems...'/)
      goto 9999
c
 2660 continue
      num_fatal = num_fatal + 1
      write(out,9284)
 9284 format(/1x,'>>>>> FATAL ERROR: the specified transformation',
     &       /1x,'                   matrix is illegal. it is not',
     &       /1x,'                   orthogonal and/or rows are',
     &       /1x,'                   all not of unit length'/)
      goto 9999
c
 2670 continue
      num_warn = num_warn + 1
      write(out,9285)
 9285 format(/1x,'>>>>> Warning: value for alpha or beta parameter',
     &       /1x,'               for smcs growth criterion expected',
     &       /1x,'                                     '/)
      goto 9999
c
 2680 continue
      num_warn = num_warn + 1
      write(out,9286) rparam
 9286 format(/1x,'>>>>> Warning: value for percent overshoot',
     &       /1x,'               expected. Using default value',
     &       /1x,'               of ',f6.2,' percent.'/)
      goto 9999
c
 2690 continue
      num_warn = num_warn + 1
      write(out,9287) dparam * hundred
 9287 format(/1x,'>>>>> Warning: minimum load reduction cannot',
     &       /1x,'               be zero.  Setting value to',
     &       /1x,'               ',f6.2,' percent.'/)
      goto 9999
c
 2700 continue
      write(out,9288) rparam * hundred, dparam * hundred
 9288 format(/1x,'>>>>> Note: the overshoot control algorithm ',
     &           'predicts that',
     &       /1x,'            the CTOA will overshoot the critical ',
     &           'CTOA by',
     &       /1x,'            greater than ',f6.3,' percent in the ',
     &           'next load step.',
     &       /1x,'            Next load step will be reduced to ',
     &           f6.3,' percent',
     &       /1x,'            of its original value.',/)
      goto 9999
c
 2710 continue
      num_warn = num_warn + 1
      write(out,9289)
 9289 format(/1x,'>>>>> Warning: expecting keyword on or off'/)
      goto 9999
c
 2720 continue
      num_warn = num_warn + 1
      write(out,9290) param
 9290 format(/1x,'>>>>> Warning: the minimum number of steps ',
     &           'between node releases',
     &           '               must be greater than zero. Using ',
     &           'default of ',i3,' steps.',/)
      goto 9999
c
 2730 continue
      write(out,9291) dparam * hundred
 9291 format(/1x,'>>>>> Note: the current load step size is too ',
     &           'large for',
     &       /1x,'            accurate crack growth. The size of ',
     &           'future load',
     &       /1x,'            steps will be cut to ',f5.2,' percent ',
     &           'of the',
     &       /1x,'            original size.'/)
      goto 9999
c
c
 2740 continue
      write(out,9292)
      num_fatal = num_fatal + 1
      lngstr = matnam(param)
      write(out,9292) lngstr
 9292 format(/1x,'>>>>> Fatal error in material ',a24,':',/
     &           '         cannot specify both isotropic & ',/,
     &           '         anisotropic thermal properties.',
     &           /)
      input_ok = .false.
      goto 9999
c
 2750 continue
      num_error = num_error + 1
      write(out,9293)
 9293 format(/1x,'>>>>> Error: the maximum percent allowed change in ',
     &           'porosity between steps',/,
     &           '               must be greater than zero. Turning',
     &           ' off automatic load reduction.',/)
      goto 9999
c
 2760 continue
      num_error = num_error + 1
      write(out,9294)
 9294 format(/1x,'>>>>> Error: the maximum percent allowed change in ',
     &           'plastic strain between steps',/,
     &           '               must be greater than zero. Turning',
     &           ' off automatic load reduction.',/)
      goto 9999
c
 2770 continue
      num_warn = num_warn + 1
      write(out,9295)
 9295 format(/1x,'>>>>> Warning: memory (in mb) expected for solver use',
     &           /)
      goto 9999
c
 2780 continue
      num_warn = num_warn + 1
      write(out,9296)
 9296 format(/1x,'>>>>> Warning: string for solver scratch directory ',
     &           'expected',
     &           /)
      goto 9999
c
 2790 continue
      num_warn = num_warn + 1
      write(out,9297)
 9297 format(/1x,'>>>>> Warning: invalid option after solver command',
     &           /)
      goto 9999
c
 2800 continue
      num_warn = num_warn + 1
      write(out,9298)
 9298 format(/1x,'>>>>> Warning: must use on or off',
     &           /)
      goto 9999
c
 2810 continue
      num_fatal = num_fatal + 1
      write(out,9299)
 9299 format(/1x,'>>>>> Fatal Error: No acceptable master nodes',
     &           ' specified. Stopping execution.',
     &           /)
      goto 9999
c
 2820 continue
      num_error = num_error + 1
      write(out,9300)
 9300 format(/1x,'>>>>> Error: no acceptable master nodes in list.',
     &           ' Ignoring list.',
     &           /)
      goto 9999
c
 2830 continue
      num_warn = num_warn + 1
      write(out,9301) param
 9301 format(/1x,'>>>>> Warning: node ',i7,' is not a crack front',
     &           ' node, and thus cannot',/,
     &           '        be a master node.  Skipping node.',
     &           /)
      goto 9999
c
 2840 continue
      num_fatal = num_fatal + 1
      write(out,9302)
 9302 format(/1x,'>>>>> Error: Number of nodes through thickness',
     &           ' must be greater than zero.',/,
     &           '         Stopping execution.',
     &           /)
      goto 9999
c
 2850 continue
      num_fatal = num_fatal + 1
      write(out,9303)
 9303 format(/1x,'>>>>> Fatal Error: Programmer Error in ',
     &           'const_front routines. Stopping.',/)
      goto 9999
c
 2860 continue
      num_error = num_error + 1
      write(out,9304)
 9304 format(/1x,'>>>>> Error: Crack plane normal must be specified',
     &           ' before',/,
     &           '            master nodes can be input.  Ignoring',
     &           ' list.',
     &       /)
      goto 9999
c
 2870 continue
      num_error = num_error + 1
      write (out,9305)
 9305 format (/1x,' >>>>>  Error:  expecting ctoa distance value ',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2880 continue
      num_fatal = num_fatal + 1
      write(out,9306)
 9306 format(/1x,'>>>>> Fatal Error:  constant front growth requires',
     &           ' a valid distance'/14x,'for measurement of CTOA.',
     &           ' Stopping Execution.',/)
      go to 9999
c
 2890 continue
      num_error = num_error + 1
      write (out,9307)
 9307 format (/1x,' >>>>>  Error:  expecting ctoa angle value ',
     &        /1x,'                  command ignored.',/)
      goto 9999
c
 2900 continue
      num_fatal = num_fatal + 1
      write(out,9308) param
 9308 format(/1x,'>>>>> Fatal Error: Warp3d stores up to ',i6,' nodes',
     &           'for each master line.',/,
     &           '          The distance specified is further than',
     &           ' this many nodes.  Warp3d must',/,
     &           '          be recompiled with this value set higher',
     &           ' (set in initst.f',/,
     &           '          as num_back_nodes). Stopping Execution.',/)
      go to 9999
c
 2910 continue
      num_fatal = num_fatal + 1
      write(out,9309) param
 9309 format(/1x,'>>>>> Fatal Error: distance to measure CTOA',
     &           ' exceeds the master line length',/,
     &           '          for master node ',i7,'. Specify a shorter',
     &           ' length or check the mesh.',/,
     &           '          Stopping Execution.',/)
      go to 9999
c
 2920 continue
      num_warn = num_warn + 1
      write(out,9310)
 9310 format(/1x,'>>>>> Note: constant front data structures have',
     &           ' already been allocated.',/,
     &           '          Any changes in the number of nodes through',
     &           ' the thickness, or',/,
     &           '          in the master list will be ignored.',/)
      go to 9999
c
 2930 continue
      num_error = num_error + 1
      write(out,9311)
 9311 format(/1x,'>>>>> Error: a restart file cannot be created',
     &           ' until at least',/,
     &           '          one load step has been computed. Ignoring',
     &           ' request.',/)
      go to 9999
c
 2940 continue
      num_error = num_error + 1
      write(out,9312)
 9312 format(/1x,'>>>>> Error: missing surface number. Ignoring ',
     &        'surface definition.',/)
      go to 9999
c
 2950 continue
      num_error = num_error + 1
      write(out,9313) param
 9313 format(/1x,'>>>>> Error: shape number must be between 0 and ',i5,/,
     &           '           Ignoring surface definition.',/)
      go to 9999
c
 2960 continue
      num_error = num_error + 1
      write(out,9314)
 9314 format(/1x,'>>>>> Error: missing surface type. Ignoring ',
     &                'surface definition.',/)
      go to 9999
c
 2970 continue
      num_error = num_error + 1
      write(out,9315)
 9315 format(/1x,'>>>>> Error: more than three points given to',
     &                 ' define plane. Ignoring new point.',/)
      go to 9999
c
 2980 continue
      num_error = num_error + 1
      write(out,9316)
 9316 format(/1x,'>>>>> Error: mistake in coordinate values.',
     &           ' Ignoring point.',/)
      go to 9999
c
 2990 continue
      num_error = num_error + 1
      write(out,9317)
 9317 format(/1x,'>>>>> Error: missing spring stiffness value.',/)
      go to 9999
c
 3000 continue
      num_error = num_error + 1
      write(out,9318)
 9318 format(/1x,'>>>>> Error: missing value for frictional ',
     &      'coefficient.',/)
      go to 9999
c
 3010 continue
      num_error = num_error + 1
      write(out,9319)
 9319 format(/1x,'>>>>> Error: mistake in velocity value. ',
     &           ' Ignoring velocity.',/)
      go to 9999
c
 3020 continue
      num_error = num_error + 1
      write(out,9320)
 9320 format(/1x,'>>>>> Error: missing value for depth. Using default',
     &           ' value.'/)
      go to 9999
c
 3030 continue
      num_error = num_error + 1
      write(out,9321)
 9321 format(/1x,'>>>>> Error: at least three points must be input to',
     &           ' define a plane. Ignoring contact surface.',/)
      go to 9999
c
 3040 continue
      num_error = num_error + 1
      write(out,9322)
 9322 format(/,1x,'>>>>> Error: stiffness must be greater than zero.',
     &           ' Ignoring contact surface definition.',/)
      go to 9999
c
 3050 continue
      num_error = num_error + 1
      write(out,9323)
 9323 format(/1x,'>>>>> Error: points do not define orthogonal',
     &           ' vectors. Ignoring contact surface.',/)
      go to 9999
c
 3060 continue
      num_error = num_error + 1
      write(out,9324)
 9324 format(/1x,'>>>>> Error: radius of cylinder must be greater than',
     &           ' zero.',/)
      go to 9999
c
 3070 continue
      num_error = num_error + 1
      write(out,9325)
 9325 format(/1x,'>>>>> Error: length of cylinder must be greater than',
     &           ' zero.',/)
      go to 9999
c
 3080 continue
      num_error = num_error + 1
      write(out,9326)
 9326 format(/1x,'>>>>> Error: mistake in coordinate value. ',
     &           ' Ignoring direction.',/)
      go to 9999
c
 3090 continue
      num_error = num_error + 1
      write(out,9327)
 9327 format(/1x,'>>>>> Error: mistake in direction or point. ',
     &           ' Ignoring contact surface.',/)
      go to 9999
c
 3320  continue
      num_error = num_error + 1
      write(out,8115)
 8115 format(/1x,'>>>>> error: expecting a HYPRE keyword'
     &                        /)
      go to 9999
c
 3100 continue
      num_fatal = num_fatal + 1
      write(out,9328) param
 9328 format(/1x,'>>>>> FATAL ERROR: contact node ',i7,
     &           ' has a coordinate transformation from',/,
     &           '      user input.  Contact nodes cannot',
     &           ' have user-defined coordinate transformations.',/,
     &           '      Stopping exectution.',/)
      go to 9999
c
 3110 continue
      num_error = num_error + 1
      write(out,9329)
 9329 format(/1x,'>>>>> Error: processor number cannot be less than',
     &           ' zero.',/,
     &           '         setting ownership to processor 0.',/)
      go to 9999
c
c
 3120 continue
      num_fatal = num_fatal + 1
      write(out,9330)
 9330 format(/1x,'>>>>> Fatal Error: This version of WARP3D only',
     &           ' supports',/,
     &           '         up to ',i5,' processors. The specified',
     &           ' number of processors',
     &           '         exceeds this value.',/)
      go to 9999
c
 3130 continue
      num_warn = num_warn + 1
      write(out,9331)
 9331 format(/,1x,'>>>>> Warning: Processor assignment of blocks not',
     &           ' specified in blocking',/,
     &           '         definition.  Assigning blocks to processors',
     &           ' in round-robin',
     &           '         fashion.',/)
      go to 9999
c
 3140 continue
      num_fatal = num_fatal + 1
      write(out,9332)
 9332 format(/,1x,
     &  '>>>>> Fatal Error: Blocks assigned to more processors',
     &  ' than are',/,
     &  '         available.  Either redefine the processor',
     &  '-block assignment,',/,
     &  '         or run with more processors.',/)
      go to 9999
c
 3150 continue
      num_fatal = num_fatal + 1
      write(out,9333) param
 9333 format(/1x,'>>>>> Fatal Error: Curve number must be between 0',
     &           ' and ',i6,'.',/,
     &           '         Execution terminating.',/)
      go to 9999
c
 3160 continue
      num_warn = num_warn + 1
      string(1:) = ' '
      call entits( string, strlng )
      write(out,9334) string(1:10)
 9334 format(/1x,'>>>>> Warning: unrecognized reset command keyword: ',
     & a10,/,25x, 'Command ignored....',/)
      go to 9999
c
 3170 continue
      num_warn = num_warn + 1
      write(out,9335) rparam
 9335 format(/1x,'>>>>> Warning: invalid factor value: ',
     & f14.6,/,25x, 'Command ignored....',/)
      go to 9999
c
 3180 continue
      num_warn = num_warn + 1
      write(out,9336) rparam
 9336 format(/1x,'>>>>> Warning: invalid reference surface selection',
     & f14.6,/,25x, 'Command ignored....',/)
      go to 9999
c
 3190 continue
      num_error = num_error + 1
      write(out,9337)
 9337 format(/1x,'>>>>> Error: <real> required after node as ',
     &       'coefficient in mpc',/)
      go to 9999
c
 3200 continue
      num_error = num_error + 1
      write(out,9338)
 9338 format(/1x,'>>>>> Error: rhs const must be 0.0 in mpc',
     &  /,7x,   'in this version of WARP3D',/ )
      go to 9999
c
 3210 continue
      if( write_msg_321 ) then
        write(out,9339)
        write_msg_321 = .false.
      else
        write(out,9902)
      end if
 9339 format (/1x,' >>>>>  Warning (type 101):',
     & /,12x,'Loading patterns + constraints (w/ their multipliers)',
     & ' defined',
     & /,12x,'for this step create an incremental loading vector ',
     &       'not a ',
     & /,12x,'simple multiple of the prior step vector. ',
     &       'Non-proportional ',
     & /,12x,'loads may lead to convergence issues. ',
     & /,12x,'To improve convergence: ', /,
     & 14x,'-> Linear stiffness used at start of step',/)
 9902 format(7x,'>> Warning: see previous type 101 message')
c
      go to 9999
c
 3220 continue
      write(out,9340)
 9340 format(/1x,'>>>>> A <label> or <string> expected to provide file',
     &           ' name...',/)
      go to 9999
c
 3230 continue
      elem = param
      write(out,9341) sparam, elem
 9341 format(/1x,'>>>>> Fatal Error: ',A6,' surface of the cohesive ',
     &'element', I7,/,7x,'is not a  quadralateral')
      go to 9999
c
 3240 continue
      elem = param
      write(out,9342) elem
 9342 format(/,1x,'>>>>> warning: top and bottom nodes of the',
     &' cohesive  element',I7,/,7x,' are not coincident')
      go to 9999
c
 3250 continue
      write(out,9343)
 9343 format(/1x,'>>>>> A search tolerance value is expected ...',/)
      go to 9999
c
 3260 continue
      num_error = num_error + 1
      write(out,9344)
 9344 format(
     & /1x,'>>>>> Error: On/Off expected',/)
      go to 9999
c
 3270 continue
      write(out,9345) dparam
 9345 format(
     &/,1x,'>>>>> Note:  increment of (effective) jump displacement'
     & ,' over last 3 steps is smaller than',
     &            /,'       target value by more than 20%. Step ',
     &           'size increased to',/,
     &           '       to ',f7.3,' of original loading.',/)
      go to 9999
c
c
 3280 continue
      num_error = num_error + 1
      write(out,9346)
 9346 format(/1x,'>>>>> Error: the maximum percent allowed change in ',
     &           '(effective) relative displacement between steps',/,
     &           '               must be greater than zero. Turning',
     &           ' off automatic load reduction.',/)
      goto 9999
c
 3290 continue
      num_error = num_error + 1
      write(out,9347)
 9347 format(/1x,'>>>>> error: critical cohesive displacement',
     &      /,14x,'(fraction) value must be a real number.')
      goto 9999
c
 3300 continue
      num_error = num_error + 1
      write(out,9348)
 9348 format(/1x,'>>>>> error: critical cohesive displacement',
     &   /,14x,'(fraction) must be .gt. 1. Value set to 5.0')
      goto 9999
c
 3310 continue
      num_error = num_error + 1
      write(out,9349) param
 9349 format(/1x,'>>>>> error: number of front nodes or sets exceeds',
     &   /,14x,'program limit of: ',i5,' input ignored...',//)
      goto 9999
c
 3330 continue
      num_error = num_error + 1
      write (out,9350)
 9350 format(/1x,
     & '>>>>> error: cannot use direct solvers with mpi',
     & /,16x,'version of WARP3D.  Use hypre instead.'/)
      go to 9999
c
 3340 continue
      num_error = num_error + 1
      write(out,9351)
 9351 format(/1x,
     &'>>>>> error: piston table loading parameters not defined.'/,
     &'      searching for new loading command.'/)
      go to 9999
c
 3350 continue
      num_warn = num_warn + 1
      write(out,9352)
 9352 format(/1x,'>>>>> warning:  the current table has been ',
     &           'previously defined.',
     &       /1x,'                this definition supercedes all',
     &                          ' information from',
     &       /1x,'                the previous definition.',/)
      goto 9999
c
c
 3360 continue
      num_error = num_error + 1
      write(out,9353) max_tables
 9353 format(/1x,'>>>>> error: there is no room for the current table',
     &           'being input.'/14x,'the number of tables alre',
     &           'ady stored equals ',i6,':'/14x,'the maximum number',
     &           ' of tables allowed in the problem.'/)
      go to 9999
c
c
 3370 continue
      num_error = num_error + 1
      write(out,9354)
 9354 format(/1x,'>>>>> error: an end of card has been reached with',
     &           'out fully defining'/14x,'table parameters. dis',
     &           'regard current table name and'/14x,'search for a',
     &           ' high level command.'/)
      go to 9999
c
c
 3380  continue
      num_error = num_error + 1
      write(out,9355)
 9355 format(/1x,'>>>>> error: expecting command for table type. ',
     &           'keep searching for one.'/)
      go to 9999
c
c
 3390 continue
      num_error = num_error + 1
      write(out,9356)
 9356 format(/1x,'>>>>> error: expecting more than one row for ',
     &           'table definition. keep searching.'/)
      go to 9999
c
c
 3400 continue
      num_fatal = num_fatal + 1
      write(out,9357)
 9357 format(/1x,'>>>>> Fatal Error: flow direction in piston ',
     &           'table has a zero valued norm.'/)
      go to 9999
c
c
 3410 continue
      num_fatal = num_fatal + 1
      write(out,9358)
 9358 format(/1x,'>>>>> Fatal Error: time must increase ',
     &           'monotonically in pistontype tables.'/)
      go to 9999
c
c
 3420 continue
      num_fatal = num_fatal + 1
      write(out,9359)
 9359 format(/1x,'>>>>> Fatal Error: time must start at zero ',
     &           'in piston type tables.'/)
      go to 9999
c
c
 3430 continue
      num_error = num_error + 1
      write(out,9360)
 9360 format(/1x,'>>>>> error: table definition expected.',
     &           ' keep searching.'/)
      go to 9999
c
c
 3440 continue
      num_error = num_error + 1
      write(out,9361)
 9361 format(/1x,'>>>>> error: table column name expected.',
     &           ' keep searching.'/)
      go to 9999
c
c
 3450 continue
      num_error = num_error + 1
      write(out,9362)
 9362 format(/1x,'>>>>> error: number of table column names ',
     &           'exceeds maximum number for table type.'/,
     &           '       proceeding to read table data if possible.' )
      go to 9999
c
c
 3460 continue
      num_error = num_error + 1
      write(out,9363)
 9363 format(/1x,'>>>>> error: all column names have not been ',
     &           'defined for table. disregard table definition.'/)
      go to 9999
c
c
 3470 continue
      num_warn = num_warn + 1
      write(out,9364)
 9364 format(/1x,'>>>>> warning: an attempt has been made to input ',
     &           'more table'/16x,'entries than columns for',
     &         ' this table type. input for this'/16x,
     &         'row is terminated.'/)
      go to 9999
c
c
 3480 continue
      num_fatal = num_fatal + 1
      write(out,9365) param
 9365 format(/1x,'>>>>> Fatal Error: row ', i5, ' in table does ',
     &           'not contain all column data.'/)
      go to 9999
c
c
 3490 continue
      num_fatal = num_fatal + 1
      write(out,9366)
 9366 format(/1x,'>>>>> Fatal Error: table not completely defined.'/)
      go to 9999
c
c
 3500 continue
      num_error = num_error + 1
c
      write(out,9367)
 9367 format(/1x,'>>>>> error: unloading command must be followed ',
     &           'by keywords on or off.'/14x,'unloading command ',
     &           'will be ignored.'/)
      go to 9999
c
 3510 continue
      num_error = num_error + 1
      write(out,9368)
 9368 format(/1x,'>>>>> error: the material specified cannot be fou',
     &           'nd in the material'/14x,'library. the material ',
     &           'must be defined before it can be'/14x,'used.' /)
      go to 9999
c
 3520 continue
      num_error = num_error + 1
      write(out,9369)
 9369 format(/1x,'>>>>> error: expecting value for link_mass' /)
      go to 9999
c
c
 9999 return
      end

c     ****************************************************************
c     *                      subroutine errmsg2                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/13/2014 rhd              *
c     *                                                              *
c     *     this subroutine prints assorted error messages in re-    *
c     *     ponse to calls from all over the program. virtually all  *
c     *     error messages in the program are generated here.        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine errmsg2( errnum, param, sparam, rparam, dparam )
      use global_data ! old common.main
      implicit none
c
      integer :: errnum, param
      real :: rparam
      double precision :: dparam
      character(len=*) :: sparam
c
      integer :: strlng
      character(len=50) :: string
      character(len=35) :: strng1
      character(len=8) :: shrtst,shrtst1,shrtst2
      character(len=24) :: lngstr
      character(len=21) :: strng21
      character(len=29) :: strng29
c
      double precision, parameter :: hundred = 100.0d0
      logical, parameter :: mess_61 = .true.
c
c
c                       print the appropriate error message as
c                       indicated by the input variable errnum
c
c
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,
     &  150, 160, 170, 180, 190, 200, 210, 220, 230,
     &  240, 250, 260, 270, 280, 290, 300, 310, 320,
     &  330, 340, 350, 360, 370, 380, 390, 400, 410,
     &  420, 430, 440, 450, 460, 470, 480, 490, 500,
     &  510, 520, 530, 540, 550, 560, 570, 580, 590,
     &  600, 610, 620, 630, 640, 650, 660, 670, 680,
     &  690, 700, 710, 720, 730, 740, 750, 760, 770,
     &  780, 790, 800, 810, 820, 830, 840, 850, 860,
     &  870, 880, 890, 900 ), errnum
c
c
 10   continue
      num_warn = num_warn + 1
      write(out,9001)
 9001 format(/1x,'>>>>> warning: expecting temperature value.',
     & /,16x,'curve treated as temperature invariant...'/)
      go to 9999
c
 20   continue
      num_warn = num_warn + 1
      write(out,9002)
 9002 format(/1x,'>>>>> warning: expecting plastic strain rate value.',
     & /,16x,'curve treated as rate independent...'/)
      go to 9999
c
 30   continue
      num_warn = num_warn + 1
      write(out,9003)
 9003 format(/1x,'>>>>> warning: unknown option on stress-strain curve.',
     & /,16x,'curve treated as rate and temperature independent...'/)
      go to 9999
c
 40   continue
      num_error = num_error + 1
      write(out,9004)
 9004 format(/1x,'>>>>> error: unknown plane id for constraints.',
     & /,14x,'this command ignored...'/)
      go to 9999
c
 50   continue
      num_error = num_error + 1
      write(out,9005)
 9005 format(/1x,'>>>>> error: plane coordinate value expected.',
     & /,14x,'this command ignored...'/)
      go to 9999
c
 60   continue
      num_error = num_error + 1
      write(out,9006)
 9006 format(/1x,'>>>>> error: invalid displacement component.',
     & /,14x, 'this command ignored...'/)
      go to 9999
c
 70   continue
      num_error = num_error + 1
      write(out,9007)
 9007 format(/1x,'>>>>> error: displacement value expected.',
     & /,14x,'this command ignored...'/)
      go to 9999
c
 80   continue
      num_error = num_error + 1
      write(out,9008)
 9008 format(1x,'>>>>> error: list of segmental curves ignored...',/)
      go to 9999
c
 90   continue
      num_error = num_error + 1
      write(out,9009) param
 9009 format(/1x,'>>>>> error: too many curve sets have been specified',
     & /,14x,'max no. curve sets allowed: ',i3,
     & /,14x,'this list of segmental curves ignored..',/)
      go to 9999
c
 100  continue
      num_error = num_error + 1
      write(out,9010) param
 9010 format(/1x,'>>>>> error: invalid curve number: ',i5,
     & /,14x,'this list of segmental curves ignored..',/)
      go to 9999
c
 110  continue
      num_error = num_error + 1
      write(out,9011) param
 9011 format(/1x,'>>>>> error: all stress-strain curves of the set',
     & /,14x,'must have the same type of dependency. curve: ',i3,
     & /,14x,'this list of segmental curves ignored..',/)
      go to 9999
c
 120  continue
      num_warn = num_warn + 1
      write(out,9012)
 9012 format(/1x,'>>>>> warning: list of curves has no dependency',
     & /,17x,'on temperature or flow rate...',/)
      go to 9999
c
 130  continue
      num_error = num_error + 1
      write(out,9013) param
 9013 format(/1x,'>>>>> error: the stress-strain curves in the set',
     & /,14x,'do not have the same number of points. curve: ',i3,
     & /,14x,'this list of segmental curves ignored..',/)
      go to 9999
c
 140  continue
      num_error = num_error + 1
      write(out,9014) param
 9014 format(/1x,'>>>>> error: the stress-strain curves in the set',
     & /,14x,'do not have the same strain values. curve: ',i3,
     & /,14x,'this list of segmental curves ignored..',/)
      go to 9999
c
 150  continue
      num_error = num_error + 1
      write(out,9015)
 9015 format(/1x,'>>>>> error: the strain values must increase',
     & /,14x,'monotonically. curve ignored..',/)
      go to 9999
c
 160  continue
      num_error = num_error + 1
      call entits( string, strlng )
      write(out,9016) string(1:strlng)
 9016 format(/1x,'>>>>> error: data value expected...',
     & /,14x,'scanning: ',a, /)
      go to 9999
c
 170  continue
      num_error = num_error + 1
      call entits( string, strlng )
      write(out,9017) string(1:strlng)
 9017 format(/1x,
     &'>>>>> error: unrecognized value for temperature dependent curve',
     & /,14x,'scanning: ',a,
     & /,14x,'curve marked as temperature independent',/)
      go to 9999
c
 180  continue
      num_error = num_error + 1
      write(out,9018)
 9018 format(/1x,
     &'>>>>> error: modulus, alpha, poisson ratio must be given',
     & /,14x,'curve marked as temperature independent',/)
      go to 9999
c
 190  continue
      call errmsg3( out, param )
      go to 9999
c
 200  continue
      write(out,9020)
 9020 format(/1x,
     &'>>>>> Job terminated due to fatal blocking errors...',/)
      go to 9999
c
 210  continue
      num_fatal = num_fatal + 1
      write(out,9021)
 9021 format(/1x,
     &    '>>>>> Error: the first point on this type of',
     &  /,'              stress-strain curve must have a zero',
     &  /,'              plastic strain value...',
     & /, '              curve definition deleted'/)
      go to 9999
c
 220  continue
      num_error = num_error + 1
      call entits( string, strlng )
      write(out,9022) string(1:strlng)
 9022 format(/1x,
     &'>>>>> error: unrecognized value for strain-rate dependent curve',
     & /,14x,'scanning: ',a,
     & /,14x,'curve marked as strain-rate independent',/)
c
 230  continue
      num_error = num_error + 1
      write(out,9023)
 9023 format(/1x,
     &'>>>>> error: modulus and poisson ratio must be given',
     & /,14x,'curve marked as strain-rate independent',/)
      go to 9999
c
 240  continue
      num_warn = num_warn + 1
      call entits( string, strlng )
      write(out,9024) string(1:strlng)
 9024 format(/1x,
     &'>>>>> warning: recognized keyword. scanning: ',a,
     &/,16x,'command ignored',/)
      go to 9999
c
 250  continue
      num_warn = num_warn + 1
      call entits( string, strlng )
      write(out,9025) string(1:strlng)
 9025 format(/1x,
     &'>>>>> warning: packets file name expected',
     & /,16x,'scanning: ',a,
     &/,16x,'command ignored',/)
      go to 9999
c
 260  continue
      num_warn = num_warn + 1
      write(out,9026) sparam(1:40)
 9026 format(/1x,
     &'>>>>> warning: binary packets file: ',a,
     & /,16x,'exists. existing file opened in append mode',/)
      go to 9999
c
 270  continue
      num_warn = num_warn + 1
      write(out,9027) sparam(1:40)
 9027 format(/1x,
     &'>>>>> warning: binary packets file: ',a,
     & /,16x,'now closed for additional output.',/)
      go to 9999
c
 280  continue
c     num_warn = num_warn + 1
      write(out,9028)
 9028 format(/1x,
     &'>>>>> warning: binary packets file was not specified in',
     & /,16x,'the solution parameters. packet output unavailable.',/)
      go to 9999
c
 290  continue
      num_error = num_error + 1
      write(out,9029)
 9029 format(/1x,
     &'>>>>> warning: unrecognized material property value...',/)
      go to 9999
c
 300  continue
      num_error = num_error + 1
      call entits( string, strlng )
      write(out,9030) sparam, string(1:strlng)
 9030 format(/1x,
     &'>>>>> error: unrecognized value for material property: ',a,
     & /,14x,'scanning: ',a,/)
      go to 9999
c
 310  continue
      num_error = num_error + 1
      write(out,9031)
 9031 format(/1x,
     &'>>>>> fatal error: material properties cannot be both',
     & /,14x,'functionally graded and use segmental curves',
     & /,14x,'job terminated....')
      call die_abort
      go to 9999
c
 320  continue
      num_warn = num_warn + 1
      write(out,9032)
 9032 format(/1x,
     &'>>>>> warning: only 1 of packet, patran, or flat ',
     & /,16x,'allowed in output command. This command ignored.',/)
      go to 9999
c
 330  continue
      num_warn = num_warn + 1
      write(out,9033) param,sparam,rparam
 9033 format(/1x,
     &'>>>>> warning: elem ',i7,' has reference surface ',a6,'.',
     & /,16x,'non-zero beta_coh (input=',e13.6,') not acceptable.',
     & /,16x,'beta_coh set to 0.0.',
     & /,16x,'Also in this case all nodes on symmetry plane and ',
     & /,16x,'connected only to interface elements must have u=v=w=0.'/)
      go to 9999
c
 340  continue
      num_warn = num_warn + 1
      write(out,9034) param
 9034 format(/1x,
     &'>>>>> warning: node ',i7,' is out of bounds.',
     & /,16x,'This mpc equation will be skipped....'/)
      go to 9999
c
 350  continue
      num_error = num_error + 1
      write(out,9035)
 9035 format(/1x,
     &'>>>>> error: missing one or more direction cosines'/)
      go to 9999
c
 360  continue
      num_warn = num_warn + 1
      write(out,9036)
 9036 format(/1x,
     &'>>>>> warning: the surace must have a name to be used',
     & /,16x,'to define a tied contact region'
     & /,16x,'This surface will be skipped....'/)
      go to 9999
c
 370  continue
      num_warn = num_warn + 1
      write(out,9037)
 9037 format(/1x,
     &'>>>>> warning: a corresponding face number is expected',
     & /,16x,'for the current surface list'
     & /,16x,'This surface list will be skipped....'/)
      go to 9999
c
 380  continue
      num_error = num_error + 1
      write(out,9038) param
 9038 format(/1x,
     &'>>>>> error: element ',i7,' is out of bounds.',
     & /,16x,'This surface entry will be ignored....'/)
      go to 9999
c
 390  continue
      num_error = num_error + 1
      write(out,9039) param
 9039 format(/1x,
     &'>>>>> error: the face number cannot be negative.',
     & /,16x,'This surface entry will be ignored....'/)
      go to 9999
c
 400  continue
      num_error = num_error + 1
      write(out,9040) param
 9040 format(/1x,
     &'>>>>> fatal error: too many surfaces have been defined',
     & /,16x,'only ',i8,' surfaces are allowed',
     & /,16x,'job terminated....'/)
      go to 9999
c
 410  continue
      num_warn = num_warn + 1
      write(out,9041)
 9041 format(/1x,
     &'>>>>> warning: the tied contact set must have a name to be used',
     & /,16x,'to generate multipoint constraint equations'
     & /,16x,'This set will be skipped....'/)
      go to 9999
c
 420  continue
      num_error = num_error + 1
      write(out,9042) param
 9042 format(/1x,
     &'>>>>> fatal error: too many tied contact sets have been defined',
     & /,16x,'only ',i8,' sets are allowed',
     & /,16x,'job terminated....'/)
      go to 9999
c
 430  continue
      num_error = num_error + 1
      write(out,9043) param
 9043 format(/1x,
     &'>>>>> fatal error: too many terms are in the mpc equation',
     & /,16x,'only ',i8,' terms are allowed',
     & /,16x,'job terminated....'/)
      go to 9999
c
 440  continue
      num_error = num_error + 1
      write(out,9044) param
 9044 format(/1x,
     &'>>>>> fatal error: too many elements are in the surface',
     & /,16x,'only ',i8,' elements are allowed',
     & /,16x,'job terminated....'/)
      go to 9999
c
 450  continue
      num_error = num_error + 1
      write(out,9045)
 9045 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'reading tied contact surface data',
     & /,16x,'job terminated....'/)
      go to 9999
c
 460  continue
      num_error = num_error + 1
      write(out,9046)
 9046 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'reading tied contact data',
     & /,16x,'job terminated....'/)
      go to 9999
c
 470  continue
      num_error = num_error + 1
      write(out,9047)
 9047 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'reading mpc equations',
     & /,16x,'job terminated....'/)
      go to 9999
c
 480  continue
      num_error = num_error + 1
      write(out,9048)
 9048 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'implementing the sparse solver',
     & /,16x,'job terminated....'/)
      go to 9999
c
 490  continue
      num_warn = num_warn + 1
      write(out,9049) rparam
 9049 format(/1x,
     &'>>>>> warning: a number is expected for the tolerance value but',
     & /,16x,'none was found.'
     & /,16x,'The default value of ', f6.3, ' will be used....'/)
      go to 9999
c
 500  continue
      num_warn = num_warn + 1
      write(out,9050)
 9050 format(/1x,
     &'>>>>> warning: a name is expected for the surface but',
     & /,16x,'none was found.'
     & /,16x,'This tied contact pair will be skipped....'/)
      go to 9999
c
 510  continue
      num_warn = num_warn + 1
      write(out,9051)
 9051 format(/1x,
     &'>>>>> warning: the slave surface was expected but',
     & /,16x,'none was found.'
     & /,16x,'This tied contact pair will be skipped....'/)
      go to 9999
c
 520  continue
      num_warn = num_warn + 1
      write(out,9052)
 9052 format(/1x,
     &'>>>>> warning: the master surface was expected but',
     & /,16x,'none was found.'
     & /,16x,'This tied contact pair will be skipped....'/)
      go to 9999
c
 530  continue
      num_warn = num_warn + 1
      write(out,9053) sparam
 9053 format(/1x,
     &'>>>>> warning: no tied contact pairs have been entered',
     & /,16x,'The tied contact set, ', a16,' will be skipped....'/)
      go to 9999
c
 540  continue
      num_warn = num_warn + 1
      write(out,9054)
 9054 format(/1x,
     &'>>>>> warning: the master and slave cannot be the same surface',
     & /,16x,'This tied contact pair will be skipped....'/)
      go to 9999
c
 550  continue
      num_warn = num_warn + 1
      write(out,9055)
 9055 format(/1x,
     &'>>>>> warning: the master surface given does not exist',
     & /,16x,'This tied contact pair will be skipped....'/)
      go to 9999
c
 560  continue
      num_warn = num_warn + 1
      write(out,9056)
 9056 format(/1x,
     &'>>>>> warning: the slave surface given does not exist',
     & /,16x,'This tied contact pair will be skipped....'/)
      go to 9999
c
 570  continue
      num_warn = num_warn + 1
      write(out,9057) sparam
 9057 format(/1x,
     &'>>>>> warning: no elements have been entered',
     & /,16x,'The surface, ',a16,', will be skipped....'/)
      go to 9999
c
 580  continue
      num_error = num_error + 1
      write(out,9058)
 9058 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'processing the tied contact data',
     & /,16x,'job terminated....'/)
      go to 9999
c
 590  continue
      num_warn = num_warn + 1
      write(out,9059)
 9059 format(/,1x,'>>>>> WARNING: all existing multipoint constraint',
     &       /,1x,'               data now deleted...',/ )
       goto 9999
c
 600  continue
      num_warn = num_warn + 1
      write(out,9060) sparam
 9060 format(/1x,
     &'>>>>> warning: the tied contact set, ',a16,', has already been',
     & /,16x,'defined.  The previous entry will be deleted....'/)
      go to 9999
c
 610  continue
      num_warn = num_warn + 1
      write(out,9061) sparam
 9061 format(/1x,
     &'>>>>> warning: the surface, ',a16,', has already been defined.',
     & /,16x,'The previous entry will be deleted....'/)
      go to 9999
c
 620  continue
      num_error = num_error + 1
      write(out,9062) param
 9062 format(/1x,
     &'>>>>> fatal error: too many surface pairs have been defined',
     & /,16x,'only ',i8,' sets are allowed',
     & /,16x,'job terminated....'/)
      go to 9999
c
 630  continue
      write(out,9063)
 9063 format(/1x,
     &'>>>>> warning: expecting either on or off, but neither found.',
     & /,16x,'Display of tied mesh mpcs will be left off by default.'/)
      go to 9999
c
 640  continue
      num_error = num_error + 1
      write(out,9064)
 9064 format(/1x,
     &'>>>>> fatal error: a memory allocation error occurred while',
     & /,16x,'processing the multipoint constraint data',
     & /,16x,'job terminated....'/)
      go to 9999
c
 650  continue
      num_warn = num_warn + 1
      write(out,9065) param
 9065 format(/1x,
     &'>>>>> warning: the same dof is used as both a dependent and',
     & /,16x,'independent term in MPC equation number',i8,
     & /,16x,'This equation will not be applied....'/)
      go to 9999
c
 660  continue
      num_error = num_error + 1
      write(out,9066)
 9066 format(/1x,
     &'>>>>> fatal error: the stiffness matrix has either a zero',
     & /,16x,'or negative number on the diagonal',
     & /,16x,'job terminated....'/)
      go to 9999
c
 670  continue
      num_warn = num_warn + 1
      write(out,9067)
 9067 format(/1x,
     &'>>>>> Warning: Binary packet output is currently turned off.',
     & /,16x,'Conversion to ascii command will be ignored.')
      go to 9999
c
 680  continue
      num_warn = num_warn + 1
      write(out,9068)
 9068 format(/,1x,
     &'>>>>> Warning: No ascii packet output file name found.',
     & /,16x,'The default name will be used.')
      go to 9999
c
 690  continue
      num_warn = num_warn + 1
      write(out,9069)
 9069 format(/1x,
     &'>>>>> Warning: All packets in the binary packet file will be',
     & /,16x,'converted by default.')
      go to 9999
c
 700  continue
      num_error = num_error + 1
      write(out,9070)
 9070 format(/,1x,'>>>>> Error: Could not open binary packet file.',
     & /,14x,'Returning to look for the next high level command.')
      go to 9999
c
 710  continue
      num_error = num_error + 1
      write(out,9071)
 9071 format(/,1x,'>>>>> Error: Could not read from binary packet file.',
     & /,14x,'Returning to look for next high level command.')
      go to 9999
c
 720  continue
      num_warn = num_warn + 1
      write(out,9072)
 9072 format(/1x)

      go to 9999
c
 730  continue
      num_warn = num_warn + 1
      write(out,9073) param
 9073 format(/,1x)

      go to 9999
c
 740  continue
      num_error = num_error + 1
      write(out,9074)
 9074 format(/,1x)

      go to 9999
c
 750  continue
      num_error = num_error + 1
      write(out,9075)
 9075 format(/,1x,'>>>>> Error: Ascii packet file could not be opened.',
     & /,16x,'Returning to look for next high level command.')
      go to 9999
c
 760  continue
      num_error = num_error + 1
      write(out,9076)
 9076 format(/,1x,'>>>>> Error: ring number exceeds program limit',
     & /,16x,'of: ',i5,' input ignored...')
      go to 9999
c
 770  continue
      num_error = num_error + 1
      write(out,9077)
 9077 format(/,1x,'>>>>> Error: ring number < 0',
     & /,16x,'input ignored...')
      go to 9999
c
 780  continue
      num_warn = num_warn + 1
      write(out,9078) param,sparam
 9078 format(/1x,
     &'>>>>> warning: elem ',i7,' has reference surface ',a6,'.',
     & /,16x,'for PPR cohesive option, specify only_normal_mode',
     & /,16x,'or only_shear_mode in material definition.',
     & /,16x,'Also in this case all nodes on symmetry plane and ',
     & /,16x,'connected only to interface elements must have u=v=w=0.'/)
      go to 9999
c
 790  continue
      num_error = num_error + 1
      write(out,9079)
 9079 format(/1x,
     &'>>>>> Error: all 4 properties of the GP model must be given.',
     & /,14x,'see write-up on generalized_plasticity option of the',
     & ' cyclic model.',
     & /,14x,'This curve ignored ...'/)
      go to 9999
c
 800  continue
      num_warn = num_warn + 1
      write(out,9080)
 9080 format(/1x,
     &'>>>>> Warning: curve pts not allowed for generalized_plasticity',
     & /,16x,'option of cyclic model. Given points ignored'/)
      go to 9999
c
 810  continue
      num_warn = num_warn + 1
      write(out,9081)
 9081 format(/1x,
     &'>>>>> Warning: must be on or off. on assumed'/)
      go to 9999
c
 820  continue
      num_warn = num_warn + 1
      write(out,9082)
 9082 format(/1x,
     &'>>>>> Warning: cannot request both stream and text file'/)
      go to 9999
c
 830  continue
      num_warn = num_warn + 1
      write(out,9083)
 9083 format(/1x,
     &'>>>>> Warning: must request either stream or flat'/)
      go to 9999
c
 840  continue
      write(out,9084)
 9084 format(/1x,
     &'>>>>> Warning: compressed not supported with stream.'
     & /,16x,'compressed option ignored'/)
      go to 9999
c
 850  continue
      write(out,9085)
 9085 format(/1x,
     &'>>>>> Warning: material states output supported only'
     & /,16x,'for flat and Patran files. Command ignored'/)
      go to 9999
c
 860  continue
      write(out,9086)
 9086 format(/1x,
     &'>>>>> Warning: compressed not supported on Windows.'
     & /,16x,'compressed option ignored'/)
      go to 9999
c
 870  continue
      write(out,9087)
 9087 format(/1x,
     &'>>>>> Warning: remainder of this output command ignored.',
     & / )
      go to 9999
c
 880  continue
      num_warn = num_warn + 1
      write(out,9088)
 9088 format(/1x,
     &'>>>>> Warning: must specify Patran or WARP3D convention.',
     & /,16x,'Command ignored'/)
      go to 9999
c
 890  continue
      num_warn = num_warn + 1
      write(out,9089)
 9089 format(/1x,
     &'>>>>> Warning: must specify text or stream format.',
     & /,16x,'Command ignored'/)
      go to 9999
c
 900  continue
      num_warn = num_warn + 1
      write(out,9090)
 9090 format(/1x,
     &'>>>>> Warning: must specify file <name>.',
     & /,16x,'Command ignored'/)
      go to 9999
c
 9999 return
      end



c     ****************************************************************
c     *                      subroutine errmsg3                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/29/2016 rhd             *
c     *                                                              *
c     *            supporting, low-level error messgae routine       *
c     *                                                              *
c     ****************************************************************
c
      subroutine errmsg3( out, param )
      implicit integer (a-z)
c
      go to ( 191, 192, 193, 194, 195, 196, 197, 198, 199,
     &        200, 210, 220, 230, 240, 250 ), param
 191  continue
      write(out,fmt='(10x,a)') '>> material model mismatch'
      go to 9999
 192  continue
      write(out,fmt='(10x,a)') '>> element type mismatch'
      go to 9999
 193  continue
      write(out,fmt='(10x,a)') '>> integration order mismatch'
      go to 9999
 194  continue
      write(out,fmt='(10x,a)') '>> integration points mismatch'
      go to 9999
 195  continue
      write(out,fmt='(10x,a)') '>> geometric nonlinear/linear mismatch'
      go to 9999
 196  continue
      write(out,fmt='(10x,a)') '>> b-bar option mismatch'
      go to 9999
 197  continue
      write(out,fmt='(10x,a)')
     & '>> some elements with or without segmental stress-strain curve'
      go to 9999
 198  continue
      write(out,fmt='(10x,a)')
     & '>> segmental curve sets are not identical'
 199  continue
      write(out,9019)
 9019 format(
     & 10x,'>> for exponential cohesive option, all elements in block ',
     &/10x,'   must have same beta_coh value. For PPR cohesive, all ',
     &/10x,'   elements in block must have same only_normal_mode or ',
     &/10x,'   only_shear_mode flag if these options are used.' )
      go to 9999
c
 200  continue
      write(out,9020)
 9020 format(
     & 10x,'>> for cyclic model, all elements in block must be the ',
     &/10x,'   same material. this restriction to be removed ',
     &/10x,'   in future.',/)
      go to 9999
c
 210  continue
      write(out,9021)
 9021 format(
     & 10x,'>> fatal error: invalid dva value in ouocdd',
     &/10x,'   job terminated.')
      go to 9999
c
 220  continue
      write(out,9022)
 9022 format(
     & 10x,'>> fatal error: invalid state = 1 in ouocdd',
     &/10x,'   job terminated.')
      go to 9999
c
 230  continue
      write(out,9023)
 9023 format(
     & 10x,'>> fatal error: at least one convergence test must be',
     &/10x,'                defined for global Newton iterations',
     &/10x,'   job terminated.')
      call die_gracefully
c
 240  continue
      write(out,*) " "
      write(out,9024)
 9024 format(/1x,
     &    '>>>>> Error: no computational results are available',
     &/1x,'             yet for output.',
     &/1x,'             output command ignored.')
      write(out,*) " "
           return
      go to 9999
c
 250  continue
      write(out,9025)
 9025 format(
     & 10x,'>> ***** available  ******',
     &/10x,'            default value = 0 assumed')
      go to 9999

c
 9999 continue
      return
      end



c     ****************************************************************
c     *                      subroutine errmsg4                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/18/2017 rhd             *
c     *                                                              *
c     *          specific error message routine for material input   *
c     *                                                              *
c     ****************************************************************
c
      subroutine errmsg4( matnum, sparam, mess, num_fatal, out, matnam,
     &                    input_ok )
      implicit none
c
      integer :: matnum, out, num_fatal
      logical :: input_ok
      character(len=*) :: sparam, matnam(*)
      character(len=50) :: mess
      character(len=24) :: lngstr
c
      num_fatal = num_fatal + 1
      lngstr = matnam(matnum)
      write(out,9275) lngstr, sparam, mess
 9275 format(/1x,'>>>>> Fatal error in material ',a24,':',/
     &          '          ',a7,' must be ',a50,/)
      input_ok = .false.
      return
      end

