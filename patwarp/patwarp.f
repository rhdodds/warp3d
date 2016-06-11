c ****************************************************************************
c *                                                                          *
c *     patran-to-warp3d input translation system - Linux version            *
c *     includes Metis for domain decomposition                              *
c *                                                                          *
c ****************************************************************************
c *                                                                          *
c *    Modifications  8/1/2014 rhd                                           *
c *                                                                          *
c ****************************************************************************
c
c
      program main
c
      use patwarp_data
      implicit integer (a-z)
      logical   status
c
c            do some initialization. set neutral file and output
c            file ids to those files.
c
      call init
c
c            write the program id message. scan for command line
c            arguments.
c            get neutral file name, translated file name.  strip
c            leading blanks from name in case user goofs up.
c
      write(termot,9000) maxnod, maxele
      write(termot,9010)
      read(termin,9020) neunam
      call stripf( neunam )
c
      write(termot,9040)
      read(termin,9020) outnam
      call stripf( outnam )
c
c            try to open the files.
c
 100  continue
      if ( neunam(1:1) .eq. ' ' ) then
          neunam(1:) = 'patran.out.1'
      end if
      if ( outnam(1:1) .eq. ' ' ) then
         outnam(1:) = 'warp3d_input'
      end if
      open( unit=ofile, file=outnam, status='unknown' )
      open( unit=nfile, file=neunam, status='old', err=110 )
      go to 200
c
c            neutral file could not be opened.  get name again.
c
 110  continue
      write(termot,9050) neunam
      write(termot,9010 )
      read(termin,9020) neunam
      call stripf( neunam )
      open( unit=nfile, file=neunam, status='old', err=110 )
c
 200  continue
c
c            call questions routine to obtain the formating information
c            for the warp input file
c
      call questions
c
c            begin an open loop to read the neutral file.  a header
c            data card is read and the appropriate routine called
c            to process remaining data cards of that packet type.
c            the routine returns here ready for the next header card to
c            be processed.  quit reading on a packet type 99.
c
 300  continue
      read(nfile,9070,end=5000) it, id, iv, kc, n1, n2, n3, n4, n5
      if ( debug ) write(termot,9090) it, id, iv, kc, n1, n2, n3,
     &                                n4, n5
      if ( it .eq. 99 ) go to 2000
      if ( it .le. 0  .or.  it .gt. 30 ) go to 4000
c
      go to ( 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009,
     &        1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018,
     &        1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027,
     &        1028, 1029, 1030 ), it
c
c            packet type 1 -- node data
c
1001  continue
      call ptyp01( id, status )
      if ( status ) go to 300
      call send( 1 )
      go to 3000
c
c            packet type 2 -- element data
c
 1002 continue
      call ptyp02( id, iv, kc, n1, status )
      if ( status ) go to 300
      call send( 2 )
      go to 3000
c
c            packet type 3 -- material properties
c
 1003 continue
      call pflush( nfile, kc, status )
      call ptyp03( status )
      if ( status ) go to 300
      call send( 3 )
      go to 3000
c
c            packet type 4 -- physical properties
c
 1004 continue
      call pflush( nfile, kc, status )
      call ptyp04( status )
      if ( status ) go to 300
      call send( 4 )
      go to 3000
c
c            packet type 5 -- coordinate frames
c
 1005 continue
      call pflush( nfile, kc, status )
      call ptyp05( status )
      if ( status ) go to 300
      call send( 5 )
      go to 3000
c
c            packet type 6 -- distributed loads
c
 1006 continue
      call ptyp06( id, iv, kc, status )
      if ( status ) go to 300
      call send( 6 )
      go to 3000
c
c            packet type 7 -- node forces
c
 1007 continue
      call ptyp07( id, iv, status )
      if ( status ) go to 300
      call send( 7 )
      go to 3000
c
c            packet type 8 -- node displacements
c
 1008 continue
      call ptyp08( id, iv, status )
      if ( status ) go to 300
      call send( 8 )
      go to 3000
c
c            packet type 9 -- bar element initial displacements
c
 1009 continue
      call pflush( nfile, kc, status )
      call ptyp09( status )
      if ( status ) go to 300
      call send( 9 )
      go to 3000
c
c            packet type 10 -- node temperatures
c
 1010 continue
      call ptyp10( id, iv, status )
      if ( status ) go to 300
      call send( 10 )
      go to 3000
c
c            packet type 11 -- element temperatures
c
 1011 continue
      call pflush( nfile, kc, status )
      call ptyp11( status )
      if ( status ) go to 300
      call send( 11 )
      go to 3000
c
c            packet type 12-13 not currently used by patran
c
 1012 continue
 1013 continue
c
c            packet type 14 -- multipoint constraints
c
 1014 continue
      call ptyp14( kc, n1, n2, status )
      if ( status) go to 300
      call send( 14 )
      go to 3000
c
c            packet type 15-20 not currently used by patran
c
 1015 continue
 1016 continue
 1017 continue
 1018 continue
 1019 continue
 1020 continue
      go to 3000
c
c            packet type 21 -- named components.  translate those named
c                              components that define coupled sets.
c                              a coupled set is a bunch of nodes which are
c                              linked for absolute or relative constraints.
c                              the nodes generally lie on the outer faces
c                              of a model.
c
 1021 continue
      call ptyp21( id, iv, kc, status )
      if ( status ) go to 300
      call send( 21 )
      go to 3000
c
c            packet type 22 - 24 not currently used by patran.
c
 1022 continue
 1023 continue
 1024 continue
      go to 3000
c
c            packet type 25 -- neutral file title
c
 1025 continue
      call ptyp25( status )
      if ( status ) go to 300
      call send( 25 )
      go to 3000
c
c            packet type 26 -- summary data
c
 1026 continue
      call ptyp26( n1, n2, n3, n4, status )
      call size_data_structures
      if ( status ) go to 300
      call send( 26 )
      go to 3000
c
c            packet type 27-30 not currently used by patran
c
 1027 continue
 1028 continue
 1029 continue
 1030 continue
      go to 3000
c
c            finished reading neutral file with normal condition.
c            ask which is wanted: finite format or warp3d format
c
 2000 continue
      call trnl
      close ( unit = ofile )
      write(termot,9080)
      call exit
c
c            fatal error detected by one of the packet processing
c            routines.
c
 3000 continue
      write(termot,9110)
      stop
c
c            invalid packet type detected (assume we have encountered
c            a phase 1 packet).  flush it.
c
 4000 continue
      call pflush( nfile, kc, status )
      if ( status ) go to 300
      go to 3000
c
c            unexpected end of file on neutral file.
c
 5000 continue
      write(termot,9120)
      write(termot,9110)
      stop
c
c
 9000 format(/,
     &       ' ****************************************************',
     &  /,   ' *                                                  *',
     &  /,   ' *     PATRAN to WARP3D Neutral File Translator     *',
     &  /,   ' *                    (Linux)                       *',
     &  /,   ' *                                                  *',
     &  /,   ' *   Processes Patran 2 (formatted) Neutral File    *',
     &  /,   ' *       (',i7,' nodes - ',i7,' elements)         *',
     &  /,   ' *            Build Date:  6-11-2016                *',
     &  /,   ' *                                                  *',
     &  /,   ' * includes:                                        *',
     &  /,   ' *  o support for 8, 9, 12, 15, 20-node hexs        *',
     &  /,   ' *  o support for 4, 10 node tets                   *',
     &  /,   ' *  o support for 6, 15 node wedges                 *',
     &  /,   ' *  o output of blocking and partitioning info      *',
     &  /,   ' *       in Patran-readable (element) results files *',
     &  /,   ' *  o tet4 and tet10 elements now supported         *',
     &  /,   ' *  o MPCs defined in Patran model now supported    *',
     &  /,   ' *  o Domain decomposition using METIS for MPI      *',
     &  /,   ' *       based, parallel analyses in WARP3D.        *',
     &  /,   ' *       Hex, tet & mixed hex-tet meshes supported. *',
     &  /,   ' *                                                  *',
     &  /,   ' ****************************************************'
     &   )
 9010 format(2(/),1x, '>> patran neutral file name ',
     &   '(default: patran.out.1) ? ',$)
 9020 format(a80)
 9040 format(/,1x, '>> warp3d input file name ',
     &   '(default: warp3d_input) ? ',$)
 9050 format(1x,'>>>> error -- neutral file ',a30,' could not ',
     &          'be opened' )
 9070 format( i2, 8i8 )
 9080 format(/,8x,'>> analysis file generation completed.',
     &        /,8x,'>> job terminated normally. ')
 9090 format(1x,'>> packet header information: ',
     &   /,5x,'it, id, iv, kc = ',4i8,
     &   /,5x,'n1, n2, n3, n4 = ',4i8,
     &   /,5x,'            n5 = ', i8 )
 9110 format(/,1x,'>> job terminated due to fatal error.')
 9120 format(/,1x,'>> fatal error -- unexpected eof on neutral file')
 9140 format('>> Illegal option specified. Try again.')
c
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- s t r i p f                      *
c *                                                               *
c *****************************************************************
c
      subroutine stripf( string )
      implicit integer (a-z)
      character *(*) string, copy*256
c
c             strip leading blanks from string and return it.
c
      last = len(string)
      copy(1:last) = string(1:last)
      do 100 cpos = 1, last
       if ( string(cpos:cpos) .ne. ' ' ) go to 200
 100  continue
      return
 200  continue
      string(1:) = copy(cpos:last)
      return
      end
c
c ************************************************************************
c *                                                                      *
c *   routine  init --  initialize variables for the translator.        *
c *                                                                      *
c ************************************************************************
c
c
      subroutine init
      use patwarp_data
      implicit integer (a-z)
c
      two16 = 65536
      debug = .false.
      nfile = 20
      ofile = 21
      ofile_store = ofile
      termin = 5
      termot = 6
      eincpt = 0
      eadvpt = 0
      cd_file = 10
      inbl_file = 11
      ct_file  = 12
      cdf_len  = 0
      inbl_len = 0
      cstf_len = 0
      prnode = .false.
      princi = .false.
      prcons = .false.
      prdistl= .false.
      nxtnls = 1
      nxtndt = 1
      nxtnvl = 1
      nxteld = 1
      nxtdld = 1
      prnlod = .false.
      prntmp = .false.
      prnmpc = .false.
      prcset = .false.
      scalar = .false.
      sep_file_control = .false.
      ncset(1:6) = 0
      nmpc = 0
      nptfle(1:80) = ' '
c
      return
      end
c
c ************************************************************************
c *                                                                      *
c *   routine  size_data_structures -- allocates space from model size   *
c *                                                                      *
c ************************************************************************
c
c
      subroutine size_data_structures
      use patwarp_data
      implicit integer (a-z)
c
      allocate ( nodcon(numnod),
     &           nodndf(numnod),
     &           nodcfg(numnod),
     &           nodcid(numnod),
     &           nodpcn(numnod),
     &           nod_trn_list(numnod),
     &           cset(numnod,6),
     &           coord(3,numnod),
     &           nodtyp(numnod) )
c
      allocate ( elepid(numele),
     &           elecid(numele),
     &           eletyp(numele),
     &           elennd(numele),
     &           elecfg(numele),
     &           eleipt(numele),
     &           elenad(numele),
     &           eleadp(numele),
     &           eletran(numele),
     &           eletrpoin(numele),
     &           elreno(numele),
     &           eleord(numele),
     &           elegnl(numele),
     &           elem_access(numele),
     &           elreon(numele),
     &           eleproc(numele),
     &           grplst(numele),
     &           etypls(numele,2),
     &           eleang(3,numele) )
c
      maxinc = numele*30
      allocate ( eleinc(maxinc),
     &           eletrinc(maxinc) )
c
      maxadv = numele*3
      allocate ( eleadv(maxadv) )
c
      maxnls = numnod
      allocate ( ndcnls(maxnls,2) )
c
      maxndt = numnod*3
      allocate ( ndcndt(maxndt),
     &           rndcndt(maxndt) )
c
      maxnvl = numnod*3
      allocate ( nodval(maxnvl),
     &           nodptr(maxnvl),
     &           rnodval(maxnvl) )
c
      maxeld = numele
      allocate ( eloads(maxeld,3) )
c
      maxdld = 4*(numele+1)
      allocate ( deload(maxdld,2) )
c
      maxmpc = numnod/20
      if (maxmpc .lt. 50)  maxmpc = 50
      allocate ( mpctrm(maxmpc),
     &           mpcnod(maxmpc,10),
     &           mpcdof(maxmpc,10),
     &           mpcmlt(maxmpc,10) )
c
      return
      end
c
c ************************************************************************
c *                                                                      *
c *   routine  send --  send a message without parameters to terminal    *
c *                                                                      *
c ************************************************************************
c
c
      subroutine send( errno )
      use patwarp_data
      implicit integer (a-z)
c
      go to ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
     &       130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230,
     &       240, 250, 260, 270, 280, 290, 300 ), errno
c
 10   continue
      write(termot,9001)
      return
c
 20   continue
      write(termot,9002)
      return
c
 30   continue
      write(termot,9003)
      return
c
 40   continue
      write(termot,9004)
      return
c
 50   continue
      write(termot,9005)
      return
c
 60   continue
      write(termot,9006)
      return
c
 70   continue
      write(termot,9007)
      return
c
 80   continue
      write(termot,9008)
      return
c
 90   continue
      write(termot,9009)
      return
c
 100  continue
      write(termot,9010)
      return
c
 110  continue
      write(termot,9011)
      return
c
 120  continue
      write(termot,9012)
      return
c
 130  continue
      write(termot,9013)
      return
c
 140  continue
      write(termot,9014)
      return
c
 150  continue
      write(termot,9015)
      return
c
 160  continue
      write(termot,9016)
      return
c
 170  continue
      write(termot,9017)
      return
c
 180  continue
      write(termot,9018)
      return
c
 190  continue
      write(termot,9019)
      return
c
 200  continue
      write(termot,9020)
      return
c
 210  continue
      write(termot,9021)
      return
c
 220  continue
      write(termot,9022)
      return
c
 230  continue
      write(termot,9023)
      return
c
 240  continue
      write(termot,9024)
      return
c
 250  continue
      write(termot,9025)
      return
c
 260  continue
      write(termot,9026)
      return
c
 270  continue
      write(termot,9027)
      return
c
 280  continue
      write(termot,9028)
      return
c
 290  continue
      write(termot,9029)
      return
c
 300  continue
      write(termot,9030)
      return
c
 9001 format(1x,' ')
 9002 format(1x,' ')
 9003 format(1x,'fatal error in material properties packet')
 9004 format(1x,'fatal error in physical properties packet')
 9005 format(1x,'fatal error in coordinate frames packet')
 9006 format(1x,'fatal error in distributed loads packet')
 9007 format(1x,'fatal error --                           ')
 9008 format(1x,'fatal error --                                 ')
 9009 format(1x,'fatal error in bar init. displacement packet')
 9010 format(1x,'fatal error --                               ')
 9011 format(1x,'fatal error in element temperatures packet')
 9012 format(1x,'fatal error -- packet type 12 not defined by patran')
 9013 format(1x,'fatal error -- packet type 13 not defined by patran')
 9014 format(1x,'fatal error in multipoint constraint packet')
 9015 format(1x,'fatal error -- packet type 15 not defined by patran')
 9016 format(1x,'fatal error -- packet type 16 not defined by patran')
 9017 format(1x,'fatal error -- packet type 17 not defined by patran')
 9018 format(1x,'fatal error -- packet type 18 not defined by patran')
 9019 format(1x,'fatal error -- packet type 19 not defined by patran')
 9020 format(1x,'fatal error -- packet type 20 not defined by patran')
 9021 format(1x,'fatal error in named component packet ' )
 9022 format(1x,'fatal error -- packet type 22 not defined by patran')
 9023 format(1x,'fatal error -- packet type 23 not defined by patran')
 9024 format(1x,'fatal error -- packet type 24 not defined by patran')
 9025 format(1x,'fatal error -- packet type 25 not defined by patran')
 9026 format(1x,'fatal error -- packet type 26 not defined by patran')
 9027 format(1x,'fatal error -- packet type 27 not defined by patran')
 9028 format(1x,'fatal error -- packet type 28 not defined by patran')
 9029 format(1x,'fatal error -- packet type 29 not defined by patran')
 9030 format(1x,'fatal error -- packet type 30 not defined by patran')
c
      end
c
c ************************************************************************
c *                                                                      *
c *    routine  ptyp25 -- read type 25 data packet                       *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp25( status )
      use patwarp_data
      implicit integer (a-z)
      logical  status
c
c
c            read the user title from data line type 25.
c
      if ( debug ) write(termot,1000)
      read(nfile,1010) ustitl
      write(termot,1020) ustitl(1:80)
      status = .true.
      return
c
 1000 format(1x,'>> entered routine ptyp25')
 1010 format(a80)
 1020 format(/,8x,'>> user title as read from neutral file is:',
     &  /,/,1x,a80 )
c
      end
c
c ************************************************************************
c *                                                                      *
c *   routine ptyp26  --  read type 26 data packet                       *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp26( n1, n2, n3, n4, status )
      use patwarp_data
      implicit integer (a-z)
      logical status
c
c
c            set the type 26 data (model sizes) passed from
c            control routine.  then read
c            next line which contains date, time, version.
c
c
      if ( debug ) write(termot,1080)
      numnod = n1
      numele = n2
      nummtl = n3
      numprp = n4
      msgnod = 100
      if ( numnod .gt. 500 ) msgnod = 200
      if ( numnod .gt. 1500) msgnod = 1000
      if ( numnod .gt. 5000) msgnod = 5000
      msgele = 100
      if ( numele .gt. 500 ) msgele = 200
      if ( numele .gt. 1500) msgele = 500
      if ( numele .gt. 5000) msgele = 5000
c
c            read date, time, and patran version number.
c
      read(nfile,1000) date, time, versn
c
c            check that the sizes for this model are within
c            maximum limits assigned for common areas.
c
      status = .true.
      if ( numnod .gt. maxnod ) then
             status = .false.
             write(termot,1010) numnod, maxnod
      end if
c
      if ( numele .gt. maxele ) then
             status = .false.
             write(termot,1020) numele, maxele
      end if
c
      if ( nummtl .gt. maxmtl ) then
             status = .false.
             write(termot,1030) nummtl, maxmtl
      end if
c
      if ( numprp .gt. maxprp ) then
             status = .false.
             write(termot,1040) numprp, maxprp
      end if
c
      if ( .not. status ) write(termot,1050)
      if ( .not. status ) return
c
      write(termot,1060) date, time, versn(6:12)
      write(termot,1070) numnod, numele, nummtl, numprp
      return
c
1000  format( a12, a8, a12 )
1010  format(1x,'>>>> fatal error -- number of nodes: ',i8,
     &        /,'                    exceeds current limit of: ',i8)
1020  format(1x,'>>>> fatal error -- number of elements: ',i8,
     &        /,'                    exceeds current limit of: ',i8)
1030  format(1x,'>>>> fatal error -- number of material properties',
     &       ': ',i5,
     &        /,'                    exceeds current limit of: ',i8)
1040  format(1x,'>>>> fatal error -- number of physical properties',
     &       ': ',i5,
     &        /,'                    exceeds current limit of: ',i8)
1050  format(1x,'>>>> fatal errors have occurred.  neutral file',
     &        /,'     processing aborted.')
1060  format(/,8x,'>> neutral file created on: ', a12,' time: ',a8,
     &       //,8x,  '>> patran version:', a )
1070  format(/,8x,'>> model size parameters:'
     &  /,14x,'>> number of nodes ...............',i8,
     &  /,14x,'>> number of elements ............',i8,
     &  /,14x,'>> number of materials ...........',i8,
     &  /,14x,'>> number of physical properties .',i8 )
1080  format(1x,'>> entering subroutine ptyp26')
c
      end
c
c ************************************************************************
c *                                                                      *
c *    routine ptyp01 -- read data packet 01                             *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp01( n1, status )
      use patwarp_data
      implicit integer (a-z)
      logical  status
      double precision  x, y, z
      character * 1 gtype
      dimension pspc(6)
c
c            complete reading and saving of this type 01 data packet.
c
      if ( debug ) write(termot,1000)
c
      nodeid = n1
      if ( .not. prnode ) write(termot,1050)
      prnode = .true.
      if ( mod( nodeid,msgnod ) .eq. 0 )
     &                  write(termot,1040) nodeid
c
c            read 2 lines of data for this packet type.
c
      read(nfile,1010) x, y, z
      read(nfile,1020) icf, gtype, ndf, config, cid, pspc
c
c            store coordinates and other data.  constraint flags are
c            compressed into the rightmost 6 bits of the nodpcn
c            variable for this nodeid.
c
      coord(1,nodeid) = x
      coord(2,nodeid) = y
      coord(3,nodeid) = z
      nodcon(nodeid) = icf
      nodtyp(nodeid) = gtype
      nodndf(nodeid) = ndf
      nodcfg(nodeid) = config
      nodcid(nodeid) = cid
c
      word = 0
      do 100 i = 1, 6
        if ( pspc(i) .eq. 1 ) word = biton( word, i)
 100  continue
      nodpcn(nodeid) = word
c
      status = .true.
      if ( .not. debug ) return
c
      write(termot,1030) nodeid, x, y, z, icf, gtype, ndf, config,
     &                 pspc
      return
c
 1000 format( 1x,'>> enter routine ptyp01' )
 1010 format( 3e16.9 )
 1020 format( i1, a1, 3i8, 2x, 6i1 )
 1030 format( 1x,'>> nodeid: ',i5,' x,y,z: ',3e19.6,
     &         /,'   icf: ',i1,' gtype: ',a1,' ndf: ',i8,
     &         /,'   config: ',i8,' cid: ',i8,' pspc: ',6i1)
 1040 format( 14x,'>> processing data for node.......',i8)
 1050 format( /,8x,'>> begin processing nodal data')
c
      end
c
c ************************************************************************
c *                                                                      *
c *    routine ptyp02 -- read and store data packet type 02.             *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp02( id, iv, kc, n1, status )
      use patwarp_data
      implicit integer (a-z)
      real   dvals, angle1, angle2, angle3
      logical  status
      dimension   nodlst(10), dvals(5)
      data  nnpc, nadvpc / 10, 5 /
c
c
c            read and store all data for this element data packet.
c
      if ( debug ) write(termot,1000)
      elemid = id
      if ( elemid .gt. numele ) then
        write(termot,1300)
        stop
      end if
      if ( .not. princi ) write(termot,1070)
      princi = .true.
      if ( mod( elemid,msgele ) .eq. 0 )
     &                write(termot,1060) elemid
c
c            read next line for element.
c
      read(nfile,1010) nodes, config, pid, ceid, angle1, angle2,
     &                 angle3
c
      elennd(elemid) = nodes
      elecfg(elemid) = config
      elepid(elemid) = pid
      elecid(elemid) = ceid
      eleang(1,elemid) = angle1
      eleang(2,elemid) = angle2
      eleang(3,elemid) = angle3
      if ( debug ) write(termot,1100) id, iv, nodes, config
c
c            read element incidences.  compute number of data lines
c            to read.  values are stored nnpc per line into a local
c            vector.  they are then placed into a pointer and vector
c            data structure.
c
      eleipt(elemid) = eincpt + 1
      nlines = ( nodes + nnpc - 1 ) / nnpc
      count = 0
      do 200 line = 1, nlines
             read(nfile,1020) nodlst
             do 150 j = 1, nnpc
                    count = count + 1
                    if ( count .gt. nodes ) go to 300
                    eincpt = eincpt + 1
                    if ( eincpt .gt. maxinc ) go to 700
                    eleinc(eincpt) = nodlst(j)
 150         continue
 200  continue
c
c            completed reading and storing element incidences without
c            list overflow.  repeat same process for associated
c            data values if any are present.
c
 300  continue
c
c            Types of elements: in warp all elements are solids
c             iv = 5   tets (patwarp element types 11,12)
c             iv = 7   wedges (patwarp element types 21, 22)
c             iv = 8   hexes (patwarp element types 1-5)
c
      select case ( iv )
      case ( 5 )
c
c                  tet elements. 4 or 10 node versions
c
         if ( nodes .eq.  4 ) then
           eletyp(elemid) = 11
         elseif ( nodes .eq. 10 ) then
           eletyp(elemid) = 12
         end if
c
      case ( 7 )
c
c                  wedge elements. 6 or 15 node versions
c
         if ( nodes .eq.  6 ) then
           eletyp(elemid) = 21
         elseif ( nodes .eq. 15 ) then
           eletyp(elemid) = 22
         end if
c
      case ( 8 )
c
c                  hex elements. from patran these can be only
c                  8 or 20 node elements. look for zeros in
c                  nodes of 20-node element. It could really be an
c                  8-node element. set element type for processing
c                  later in patwarp.
c
         nnode_elem = nodes
         if ( nnode_elem .eq. 20 ) then
            call order1215( eleinc(eleipt(elemid)), nnode_elem, 1 )
            if ( nnode_elem .eq. 8 ) elennd(elemid) = 8
         end if
         if ( nnode_elem .eq.  8 ) then
           eletyp(elemid) = 1
         elseif ( nnode_elem .eq. 12 ) then
           eletyp(elemid) = 2
         elseif ( nnode_elem .eq. 15 ) then
           eletyp(elemid) = 3
         elseif ( nnode_elem .eq. 20 ) then
           eletyp(elemid) = 4
         elseif ( nnode_elem .eq.  9 ) then
           eletyp(elemid) = 5
         end if
c
      case default
         write(termot,1200)  elemid
         stop
      end select
c
c
      if ( debug ) write(termot,1110) eleipt(elemid), (eleinc(j),
     &             j = eleipt(elemid), eincpt )
      elenad(elemid) = 0
      if ( n1 .eq. 0 ) go to 600
c
      elenad(elemid) = n1
      eleadp(elemid) = eadvpt + 1
      nlines = ( n1 + nadvpc - 1 ) / nadvpc
      count = 0
      do 500 line = 1, nlines
             read(nfile,1030) dvals
             do 400 j = 1, nadvpc
                    count = count + 1
                    if ( count .gt. n1 ) go to 600
                    eadvpt = eadvpt + 1
                    if ( eadvpt .gt. maxadv ) go to 800
                    eleadv(eadvpt) = dvals(j)
 400         continue
 500  continue
c
c            completed reading associated data values.
c            completed reading this packet of element data.
c
 600  continue
      status = .true.
      return
c
c            overflow of available space in the vector of element
c            incidences.
c
 700  continue
      write(termot,1040)
      status = .false.
      return
c
c            overflow of available space in the vector of element
c            associated data values.
c
 800  continue
      write(termot,1050)
      status = .false.
      return
c
 1000 format(1x,'>> enter routine ptyp02')
 1010 format( 4i8, 3e16.9 )
 1020 format( 10i8 )
 1030 format( 5e16.9 )
 1040 format( 1x,'>>>> fatal error -- overflow of system incidence',
     &          ' vector' )
 1050 format( 1x,'>>>> fatal error -- overflow of system associated'
     &        /,'                    element data vector.' )
 1060 format( 14x,'>> processing data for element....',i8 )
 1070 format(/,8x,'>> begin processing element data' )
 1100 format(1x,'>> element id: ',i8,' shape: ',i5,' nodes: ',i8,
     &         ' config: ',i5 )
 1110 format(1x,'>> element incidence pointer: ', i8,
     &        /,'>> element incidences: ',10(/,10x,10i8) )
 1200 format(1x,'>> element id: ',i8,' has unsupported type..',
     &       /, '   translation aborted' )
 1300 format(1x,'>> element id > number of elements. job aborted' )
c
      end
c ************************************************************************
c *                                                                      *
c *   routine  ptyp03 -- read packet type 03 -- material properties      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp03( status )
      use patwarp_data
      implicit integer (a-z)
      logical status
c
c            this function not yet implemented.
c
      return
      end
c ************************************************************************
c *                                                                      *
c *   routine  ptyp04 -- read packet type 04 -- physical properties      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp04( status )
      use patwarp_data
      implicit integer (a-z)
      logical status
c
c            this function not yet implemented.
c
      return
      end
c ************************************************************************
c *                                                                      *
c *   routine  ptyp05 -- read packet type 05 -- coordinate frames        *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp05( status )
      use patwarp_data
      implicit integer (a-z)
      logical status
c
c            this function not yet implemented.
c
      return
      end
c ************************************************************************
c *                                                                      *
c *   routine  ptyp06 -- read packet type 06 -- distributed loads        *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp06( id, iv, kc, status )
      use patwarp_data
      implicit integer (a-z)
      logical status
      dimension icomp(6), node(8), ivals(5), rvals(5)
      real      rvals
      equivalence ( rvals, ivals )
c
c             read and store all data for distributed element loads.
c
      if ( debug ) write(termot,1000)
      elemid = id
      if ( .not. prdistl ) write(termot,1070)
      prdistl = .true.
      if ( mod( elemid,msgele ) .eq. 0 )
     &                write(termot,1060) elemid
c
c             check for overflow of vectors.
c
      if ( nxteld .gt. maxeld ) then
         write(termot,1040)
         status = .false.
         return
      end if
      if ( nxtdld + (kc-1)*5 + 2 .gt. maxdld ) then
         write(termot,1050)
         status = .false.
         return
      end if
c
c             data for load on this element will fit into vectors.
c             store pointer into the data vector.
c
      lodset = iv
      eloads(nxteld,1) = lodset
      eloads(nxteld,2) = nxtdld
      eloads(nxteld,3) = elemid
c
c             read and process main description line.
c
      read(nfile,1100) ltype, eflag, gflag, icomp, node, nfe
c
c             convert integer flags to a bit map. count the number
c             of non-zero element centroid and corner node loads.
c
      map = 0
      if ( ltype .eq. 1 ) map = biton( map, 1 )
      if ( eflag .eq. 1 ) map = biton( map, 2 )
      if ( gflag .eq. 1 ) map = biton( map, 3 )
      nc = 0
      bitno = 3
      do 10 i = 1, 6
        bitno = bitno + 1
        if ( icomp(i) .eq. 1 ) then
           nc = nc + 1
           map = biton( map, bitno )
        end if
 10   continue
      nn = 0
      do 20 i = 1, 8
        bitno = bitno + 1
        if ( node(i) .eq. 1 ) then
           nn = nn + 1
           map = biton( map, bitno )
        end if
 20   continue
c
c             store the loaded edge/face and bit map into data vector
c
      deload(nxtdld,1) = map
      deload(nxtdld,2) = nfe
c
c             read and store the non-zero edge/face load components.
c
      npv = nc * ( eflag + nn*gflag )
      nxtdld = nxtdld + 1
      deload(nxtdld,1) = npv
      nline = (npv+4)/5
      count = 0
      do 50 line = 1, nline
        read(nfile,1200) rvals
        do 40 j = 1, 5
          nxtdld = nxtdld + 1
          deload(nxtdld,1) = ivals(j)
          count = count + 1
          if ( count .eq. npv ) go to 60
 40     continue
 50   continue
c
c             all done with this distributed load packet. bump counters
c             for use by next packet of this type.
c
 60   continue
      nxtdld = nxtdld + 1
      nxteld = nxteld + 1
      status = .true.
      return
c
 1000 format(1x,'>> enter routine ptyp06')
 1040 format( 1x,'>>>> fatal error -- overflow of distributed load',
     &          ' list' )
 1050 format( 1x,'>>>> fatal error -- overflow of distributed load'
     &        /,'                      data vector.' )
 1060 format( 10x,'>> processing data for element: ',i8 )
 1070 format(/,8x,'>> begin processing distributed loads')
 1100 format( i1,i1,i1,6i1,8i1,i2 )
 1200 format( 5e16.9 )
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp07 -- read packet type 07 -- node forces              *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp07( id, iv, status )
      use patwarp_data
      implicit integer (a-z)
      dimension  icomp(6), lodbuf(6)
      real       lodbuf
      logical   status
      save counter
      data counter /0/
c
c            read data packet type 7 which contains prescribed
c            nodal loads at a node.
c
      status = .true.
      node = id
      lodset = iv
      counter = counter + 1
      if ( .not. prnlod ) write(termot,9070)
      prnlod = .true.
      if (mod(counter,msgnod).eq.1) write(termot,9060) node, lodset
c
c            read the alternate coord. frame id  and the dof flags.
c
      read(nfile,9010) cid, icomp
c
c            count the number of dof with loads applied.
c            build a bit mask for the loaded dof.
c
      ldof = 0
      mask = 0
      do dof = 1, 6
        if ( icomp(dof) .eq. 0 ) cycle
        ldof = ldof + 1
        mask = biton( mask, dof )
      end do
c
c            read the data cards containing values in compressed format.
c            either one or two cards.
c
      read(nfile,9020) ( lodbuf(dof), dof = 1, 5 )
      if ( ldof .eq. 6 ) read(nfile,9020) lodbuf(6)
      if ( debug ) write(termot,9050) node, cid, icomp, lodbuf
c
c            save node number, load set and loads.  status returns
c            false if there is no room in the data vector.
c
      call ptypns( 'loads', node, lodset, ldof+1, datloc, status )
      if ( status ) then
                    nodval(datloc) = mask
                    datloc = datloc + 1
                    ldof = 1
                    do 20 dof = 1, 6
                       if ( icomp(dof) .ne. 0 ) then
                           rnodval(datloc) = lodbuf(ldof)
                           datloc = datloc + 1
                           ldof = ldof + 1
                       end if
 20                 continue
             else
                    write(termot,9030) node
      end if
c
      status = .true.
c
c            all done with this node.
c
      return
 9010 format( i8,6i1 )
 9020 format( 5e16.9 )
 9030 format(1x,'>>>> warning -- insufficient array space to store',
     &   /,  1x,'                loads for node ',i5 )
 9050 format(1x,'>> enter ptyp07.',
     & /,5x,'node, cid, icomp = ',i5,i8,1x,6i1,
     & /,5x,'lodbuf = ',3e16.9,/,14x,3e16.9 )
 9060 format( 14x,'>> processing nodal loads for node',i8,
     &            ' load set: ',i5 )
 9070 format( /,7x,' >> begin processing nodal loads' )
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp08 -- read packet type 08 -- node displacements.      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp08( id, iv, status )
      use patwarp_data
      implicit integer (a-z)
      dimension  icomp(6), dspbuf(6)
      real       dspbuf
      logical   status
      save counter
      data counter /0/
c
c            read data packet type 8 which contains prescribed
c            displacements at a node.
c
      status = .true.
      node = id
      conset = iv
      if ( .not. prcons ) write(termot,9070)
      prcons = .true.
      counter = counter + 1
      if (mod(counter,msgnod).eq.0 )write(termot,9060) node
c
c            read the alternate coord. frame id  and the dof flags.
c
      read(nfile,9010) cid, icomp
c
c            count the number of dof with constraints applied.
c            build a bit mask for the constrained dof.
c
      cdof = 0
      mask = 0
      do dof = 1, 6
        if ( icomp(dof) .eq. 0 ) cycle
        cdof = cdof + 1
        mask = biton( mask, dof )
      end do
c
c            read the data cards containing values in compressed format.
c            either one or two cards.
c
      read(nfile,9020) ( dspbuf(dof), dof = 1, 5 )
      if ( cdof .eq. 6 ) read(nfile,9020) dspbuf(6)
      if ( debug ) write(termot,9050) node, cid, icomp, dspbuf
c
c            determine if there is sufficient room in the data arrays
c            to save constraints for this node.  we combine all constraints
c            into one set.
c
c      if ( conset .ne. 1 ) then
c             write(termot,9040)
c             return
c      end if
      if ( nxtnls .gt. maxnls ) then
             write(termot,9030) node
             return
      end if
      if ( nxtndt + cdof .gt. maxndt ) then
             write(termot,9030) node
             return
      end if
c
c            save node number and constraints values
c
      ndcnls(nxtnls,1) = nxtndt
      ndcnls(nxtnls,2) = node
      nxtnls = nxtnls + 1
      ndcndt(nxtndt) = mask
      nxtndt = nxtndt + 1
      cdof = 1
      do dof = 1, 6
        if ( icomp(dof) .eq. 0 ) cycle
        rndcndt(nxtndt) = dspbuf(cdof)
        nxtndt = nxtndt + 1
        cdof = cdof + 1
      end do
c
c            all done with this node.
c
      return
 9010 format( i8,6i1 )
 9020 format( 5e16.9 )
 9030 format(1x,'>>>> warning -- insufficient array space to store',
     &   /,  1x,'                constraints for node ',i8 )
 9040 format(1x,'>>>> warning -- constraint on set other than',
     &          ' 1 ignored.')
 9050 format(1x,'>> enter ptyp08.',
     & /,5x,'node, cid, icomp = ',i5,i8,1x,6i1,
     & /,5x,'dspbuf = ',3e16.9,/,14x,3e16.9 )
 9060 format( 14x,'>> processing displacements for node:    ',i8)
 9070 format( /,8x,'>> begin processing nodal displacement data')
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp09 -- read packet type 09 -- bar init. displ          *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp09( status )
      use patwarp_data
      implicit integer (a-z)
      logical status
c
c            this function not yet implemented.
c
      return
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp10 -- read packet type 10 -- node temperatures        *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp10 ( id, iv, status )
      use patwarp_data
      implicit integer (a-z)
      real      temp
      logical   status
c
c            read data packet type 10 which contains prescribed
c            temperatures at a node.
c
      status = .true.
      node = id
      lodset = iv
      if ( .not. prntmp ) write(termot,9070)
      prntmp = .true.
      if ( mod( node,msgnod ) .eq. 0 )
     &                 write(termot,9060) node, lodset
c
c            read the data card containing temperature.
c
      read(nfile,9020) temp
      if ( debug ) write(termot,9050) node, temp, lodset
c
c            save node number, load set, and temperature in the nodal
c            values list.
c
      call ptypns( 'temperatures', node, lodset, 1, datloc, status )
      if ( status ) then
                    rnodval(datloc) = temp
             else
                    write(termot,9030) node
      end if
      status = .true.
c
c            all done with this node.
c
      return
 9020 format( e16.9 )
 9030 format(1x,'>>>> warning -- insufficient array space to store',
     &   /,  1x,'                temperature for node ',i8 )
 9050 format(1x,'>> enter ptyp10.',
     & /,5x,'node, temp, lodset = ',i8,e16.9,i5 )
 9060 format( 14x,'>> processing nodal temperatures for node..',i8,
     &            ' load set: ',i5 )
 9070 format( /,8x,'>> begin processing nodal temperatures' )
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp11 -- read packet type 11 -- element temperatures     *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp11( status )
      use patwarp_data
      implicit integer (a-z)
      logical status
c
c            this function not yet implemented.
c
      return
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp14 -- read packet type 14 -- multipoint contraints    *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp14( kc, n1, n2, status )
      use patwarp_data
      implicit integer (a-z)
      character*12 junk1, junk2, type
      integer, allocatable, dimension (:) :: tmpnod, tmpdof
      double precision, allocatable, dimension (:) :: tmpmlt
      double precision mlt, mlt2
      logical status
c
      if ( .not. prnmpc ) write(termot,9090)
      prnmpc = .true.
      read(nfile,9030) junk1, junk2, type
c
c        if n2 = 0, then this is an explicit mpc and coefs must be read
c
      trmcnt = 0
      if (n2 .eq. 0) then
         allocate ( tmpnod(n1+1),
     &              tmpdof(n1+1),
     &              tmpmlt(n1+1), stat=err )
         if (err .ne. 0) then
            write(termot,9040)
            stop
         end if
         read(nfile,9000) nod, dof, mlt
         trmcnt = trmcnt + 1
         tmpnod(trmcnt) = nod
         tmpdof(trmcnt) = dof
         tmpmlt(trmcnt) = -1.0
         if ( mod(n1,2) .eq. 0 ) then
            do row = 1, kc-2
               read(nfile,9010) nod, dof, mlt, nod2, dof2, mlt2
               trmcnt = trmcnt + 1
               tmpnod(trmcnt) = nod
               tmpdof(trmcnt) = dof
               tmpmlt(trmcnt) = mlt
               trmcnt = trmcnt + 1
               tmpnod(trmcnt) = nod2
               tmpdof(trmcnt) = dof2
               tmpmlt(trmcnt) = mlt2
            end do
         else
            do row = 1, kc-3
               read(nfile,9010) nod, dof, mlt, nod2, dof2, mlt2
               trmcnt = trmcnt + 1
               tmpnod(trmcnt) = nod
               tmpdof(trmcnt) = dof
               tmpmlt(trmcnt) = mlt
               trmcnt = trmcnt + 1
               tmpnod(trmcnt) = nod2
               tmpdof(trmcnt) = dof2
               tmpmlt(trmcnt) = mlt2
            end do
            read(nfile,9000) nod, dof, mlt
            trmcnt = trmcnt + 1
            tmpnod(trmcnt) = nod
            tmpdof(trmcnt) = dof
            tmpmlt(trmcnt) = mlt
         end if
         ntrm = 0
         do trm = 1, trmcnt
            dof = tmpdof(trm)
            if (dof .gt. 3) then
               write(termot,9050)
               if (trm .gt. 1) nmpc = nmpc - 1
               return
            end if
            if (trm .eq. 1)  nmpc = nmpc + 1
            ntrm = ntrm + 1
            mpcnod(nmpc,ntrm) = tmpnod(trm)
            mpcdof(nmpc,ntrm) = dof
            mpcmlt(nmpc,ntrm) = tmpmlt(trm)
         end do
         if (ntrm .gt. 0)  mpctrm(nmpc) = ntrm
         deallocate (tmpnod,tmpdof,tmpmlt)
      end if
c
c        if n2 > 0, this is a 'rigid' type--all coefs are 1.0 or -1.0,
c           there are only 2 terms per mpc, only node numbers needed
c        also note that the last node read is the independent node for
c           each mpc
c
      if (n2 .gt. 0) then
         if (type .ne. 'RIGID') then
            write(termot,9060) type
            do row = 1, kc-1
               read(nfile,*)
            end do
            return
         end if
         allocate( tmpnod(kc), stat=err )
         if (err .ne. 0) then
            write(termot,9040)
            stop
         end if
         nodcnt = 0
         do row = 1, kc-1
            read(nfile,9020) newnod
            if (newnod .ne. oldnod) then
               nodcnt = nodcnt + 1
               tmpnod(nodcnt) = newnod
               oldnod = newnod
            end if
         end do
         indnod = tmpnod(nodcnt)
         do nod = 1, nodcnt-1
            nmpc = nmpc + 3
            depnod = tmpnod(nod)
            do mpc = 1, 3
               ptr = mpc - 3
               mpcnod(nmpc+ptr,1) = depnod
               mpcdof(nmpc+ptr,1) = mpc
               mpcmlt(nmpc+ptr,1) = -1.0
               mpcnod(nmpc+ptr,2) = indnod
               mpcdof(nmpc+ptr,2) = mpc
               mpcmlt(nmpc+ptr,2) = 1.0
               mpctrm(nmpc+ptr)   = 2
            end do
         end do
         deallocate(tmpnod)
      end if
c
      status = .true.
c
      return
 9000 format( i8, i8, e16.9 )
 9010 format( i8, i8, e16.9, i8, i8, e16.9 )
 9020 format( i8 )
 9030 format( 3a12 )
 9040 format( /,8x,'>>>> error -- a memory allocation error occured',
     &        /,8x,'              during mpc processing')
 9050 format( /,8x,'>>>> warning -- rotational dofs are not supported',
     &        /,8x,'                this mpc equation will be skipped')
 9060 format( /,8x,'>>>> warning -- ', a12,'is not a supported mpc type',
     &        /,8x,'                this mpc equation will be skipped')
 9090 format( /,8x,'>> begin processing multipoint constraints' )
      end

c ************************************************************************
c *                                                                      *
c *   routine  ptyp21 -- read packet type 21 -- named components         *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptyp21( id, iv, kc, status )
      use patwarp_data
      implicit integer (a-z)
      character  compna*12, ucname(6)*6, lcname(6)*6
      logical status
c
      data ucname/'X+FACE','X-FACE','Y+FACE','Y-FACE','Z+FACE','Z-FACE'/
      data lcname/'x+face','x-face','y+face','y-face','z+face','z-face'/
c
c            read data packet type 21 which contains named components.
c            the only named components of interest are those which
c            define coupled sets of nodes for either absolute or
c            relative constraints.  these are recognized by special
c            component names.
c
      if ( .not. prcset ) write( termot, 1000 )
      prcset = .true.
c
c              read the name of the component.  if it isn't one that
c              we want, flush the packet and return.
c
      read( nfile, 1010, err = 100 ) compna
      if ( debug ) write( termot, 1030 ) compna
      do 10 i = 1, 6
      if ( index( compna, ucname(i) ) .ne. 0 .or.
     &    index( compna, lcname(i) ) .ne. 0      ) go to 30
   10 continue
c
      call pflush( nfile, kc-1, status )
      return
c
c              found a coupled set. finish reading the packet.
c
   30 nnode = iv / 2
      ncset( i ) = nnode
      read( nfile, 1020, err=100 ) (ntype, cset(inode,i),inode=1,nnode)
      status = .true.
      return
c
c              error during read.  set up for termination.
c
  100 status = .false.
      return
c
 1000 format(/,8x,'>> begin processing named components' )

 1010 format( a12 )
 1020 format( 10i8 )
 1030 format( ' processing named component:  ', a12 )
c
      end
c
c ************************************************************************
c *                                                                      *
c *   routine bitoff -- turn a bit off                                   *
c *                                                                      *
c ************************************************************************
c
c
      integer function bitoff( iword, bitno )
      implicit integer (a-z)
c
c            turn off bit bitno in word iword if it is on.
c
      logical  bitchk
      bitoff = iword
      if ( bitchk( iword, bitno ) ) bitoff = iword - 2**(bitno-1)
      return
      end

c ************************************************************************
c *                                                                      *
c *    routine  biton  --  turn on a bit in a word                       *
c *                                                                      *
c ************************************************************************
c
c
      integer function biton( iword, bitno )
      implicit integer (a-z)
c
c            turn on bit bitno in word iword if it is not already on.
c
      logical  bitchk
      biton = iword
      if ( bitchk( iword, bitno ) ) return
      biton = iword + 2**(bitno-1)
c
      return
      end

c ************************************************************************
c *                                                                      *
c *  routine bitchk -- determine if a bit in a word is on                *
c *                                                                      *
c ************************************************************************
c
c
      logical function bitchk( iword, bitno )
      implicit integer (a-z)
c
c            determine if bit bitno in word iword is on.
c
      pow1 = 2 ** bitno
      pow2 = 2 ** (bitno-1)
      left = ( iword / pow1 ) * pow1
      right = iword - ( iword / pow2 ) * pow2
      if ( right .lt. 0 ) right = 0
      bitchk = iword - left - right .ne. 0
      return
      end

c ************************************************************************
c *                                                                      *
c *      routine ptypns -- update linked list of nodal loads, temps      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine ptypns( lodtyp, node, lodset, numval, datloc, status )
      use patwarp_data
      implicit integer (a-z)
      character *(*) lodtyp
      logical  status
c
c
c            search table of load sets to find the passed load set.
c            if not in table add at first available slot.
c
      if ( debug ) write(termot,9000) lodtyp
      avslot = 0
      do 10 slot = 1, maxset
             if ( lodset .eq. nlodtb(slot,1) ) go to 20
             if ( avslot .eq. 0 .and. nlodtb(slot,1) .eq. 0 ) then
                    avslot = slot
             end if
 10   continue
c
c            load set not in table. add it. zero head and tail pointers.
c
      slot = avslot
      nlodtb(slot,1) = lodset
      nlodtb(slot,2) = 0
      nlodtb(slot,3) = 0
      nlodtb(slot,4) = 0
      nlodtb(slot,5) = 0
      if ( debug ) write(termot,9010) lodset, slot
c
c            update link pointer for new data values to be inserted
c            by the calling routine.  get head and tail based on load
c            type.
c
 20   continue
      if ( lodtyp .eq. 'loads' .or. lodtyp .eq. 'LOADS' ) then
                    hcol = 2
                    tcol = 3
             else
                    hcol = 4
                    tcol = 5
      end if
c
      head = nlodtb(slot,hcol)
      tail = nlodtb(slot,tcol)
      if ( debug ) write(termot,9020) hcol, tcol, head, tail
c
c            will new data values fit onto end of data vector.
c
      if ( nxtnvl + numval .gt. maxnvl ) then
                    status = .false.
                    return
             else
                    status = .true.
      end if
c
c            update pointers for data values inserted at tail of list.
c
      if ( head .eq. 0 ) then
                    head = nxtnvl
                    tail = head
             else
                    nodptr(tail) = nxtnvl
                    tail = nxtnvl
      end if
c
c            update next available space location, set row number for
c            calling routine to use for loading in data value(s).
c
      nodval(tail) = node
      nodptr(tail) = 0
      nxtnvl = nxtnvl + 1
      datloc = nxtnvl
      nxtnvl = nxtnvl + numval
      nlodtb(slot,hcol) = head
      nlodtb(slot,tcol) = tail
      if ( debug ) write(termot,9030) node, nxtnvl, datloc,
     &                                head, tail
c
      return
 9000 format(' >>> enter routine ptypns.  load type= ',a )
 9010 format(' >>> open new lodset. lodset, slot = ',2i5 )
 9020 format(' >>> hcol, tcol, head, tail = ',4i5 )
 9030 format(' >>> node, nxtnvl, datloc, head, tail = ',i8,4i5 )
      end
c
c ************************************************************************
c *                                                                      *
c *     routine  pflush -- flush the remainder of the current packet     *
c *                                                                      *
c ************************************************************************
c
c
      subroutine pflush( nfile, kc, status )
      implicit integer (a-z)
      character * 1  dummy
      logical status
c
      do i = 1, kc
         read(nfile, 1000, end = 100, err=100 ) dummy
      end do
      status = .true.
      return
c
  100 status = .false.
      return
c
 1000 format( a1 )
c
      end
c
c ************************************************************************
c ************************************************************************
c ************************************************************************
c ************************************************************************
c *                                                                      *
c *     routine  trnl -- main driver for warp3d input generator          *
c *                                                                      *
c ************************************************************************
c ************************************************************************
c ************************************************************************
c ************************************************************************
c
c
      subroutine trnl
      use patwarp_data
      implicit integer (a-z)
c
c            convert the model description now stored in data arrays
c            into a warp3d input file.
c
c            currently supported input data are:
c
c                   1) model sizes
c                   2) element types
c                   3) nodal coordinates
c                   4) element incidences
c                   5) class a constraints (absolute)
c                   6) multipoint constraints
c                   7) nodal loads (forces and temperatures)
c                   8) element loads (constant pressure over faces)
c
c            initialize the warp3d specific variables.
c
      call trnlin
c
c            if required, convert 8-node hex solids to
c            hex transition elements to
c            eliminate displacement incompatibility
c
      call trntran
c
c            reorder elements in to blocks to warp input file.
c
      call trnler
c
c            check for output format (coordinates, incidences and
c            blocking, and constraints as separate files)
c
      call separate_file
c
c
c            write out header
c
      write(termot,9000)
      write(ofile,9010) ustitl(1:78)
      write(ofile,9020)
c
c            write sample material to file
c
      call sample_materials
c
      write(ofile,9030) numnod, numele
      write(termot,9040)
c
c            write element types and numbers to file.
c            first build element list then write it out.
c
      call trnlec
      call trnlpe
      write(termot,9080)
c
c            write nodal coordinate data to file.
c
      call trnlnd
      write(termot,9050)
c
c            write element incidences to file.
c
      call trnled
      write(termot,9060)
c
c            write blocking command
c
      call trnlbc
      write(termot,9120)
c
c            write nodal and element loads all load sets.
c
      call trnlnl
      write(termot,9100)
c
c            write constraints to file.
c
      call trnlcn
      write(termot,9090)
c
c            nonlinear loading, analysis parameters,
c            example compute and output commands.
c
      write(ofile,9200)
      write(ofile,9300)
c
c            all done.
c
      write(termot,9070)
c
c            write modified patran neutral file
c
      call trnlnf
c
c            write patran elem results file with blocking and processor
c            info
c
      call trnlpo
c
c
      return
c
 9000 format(8x,'>> begin warp3d input file generation')
 9010 format(2('c',/),'c ',a78,/,'c ',/,'c ' )
 9020 format('structure model',/,'c ',/,'c ',/,'material',
     & ' default',/,'    properties bilinear e 30000 nu 0.3',
     & ' yld_pt 1.0e20,'/,
     & 15x,'tan_e 100 alphax 1.0e-06 alphay 1.0e-06 alphaz 1.0e-06',
     & /,'c ',/,'c ')
 9030 format('number of nodes ',i8,'  elements ',i8)
 9040 format(14x,'>> model title and sizes written' )
 9050 format(14x,'>> nodal coordinates written' )
 9060 format(14x,'>> element incidences written' )
 9070 format(14x,'>> warp3d input file completed' )
 9080 format(14x,'>> element types written' )
 9090 format(14x,'>> constraints written' )
 9100 format(14x,'>> nodal and element loads written' )
 9110 format(14x,'>> coupled set equations written' )
 9120 format(14x,'>> blocking command written')
 9200 format('c ',/,'c ',/,' loading test',/,'  nonlinear',/,
     &  '    step 1-10 set_01 1.0 set_02 2.5 constraint 1.4 ',
     &   '$ just an example !',/,
     &  'c ',/,
     &  ' nonlinear analysis parameters',/,
     &  '   solution technique sparse direct ',/,
     &  'c   solution technique sparse iterative ',/,
     &  'c   solution technique hypre',/,
     &  'c   hypre tolerance 0.000001',/,
     &  '   time step 1.0e06 ',/,
     &  '   maximum iterations 5 $  global Newton iterations',/,
     &  '   minimum iterations 1',/,
     &  '   divergence check on',/,
     &  '   convergence test norm res tol 0.01',/,
     &  '   nonconvergent solutions stop',/,
     &  '   adaptive on  $ global Newton iterations',/,
     &  '   linear stiffness for iteration one off',/,
     &  '   batch messages off',/,
     &  '   wall time limit off',/,
     &  '   material messages off',/,
     &  '   bbar stabilization factor 0.0',/,
     &  '   consistent q-matrix on' )
 9300 format('   trace solution on',/,
     &  '   extrapolate on',/,
     &  '   display tied mesh mpcs off',/,
     &  '   user_routine off',/,'c ',/,'c ',/,
     &  ' compute displacements for loading test step 1',/,
     &  ' output displacements 1-8',/,
     &  ' output stresses 1',/,
     &  ' output strains 1',/,'c ',/,
     &  ' output patran binary displacements',/,
     &  ' output patran binary strains',/,
     &  ' output patran binary stresses',/,'c',/,
     &  ' output flat text compressed nodal stresses',/,
     &  ' output flat stream element strains',/,
     &  'c ',/,'stop' )
c
      end
c
c
c ************************************************************************
c *                                                                      *
c *  routine separate_file -- user specifies if coordinates, incidences  *
c *                           and blocking, and constraints are to be    *
c *                           written to separate output files (default  *
c *                           writes them to main output file)           *
c *                                                                      *
c *  written by: jp               8-1-2014 rhd                           *
c *                                                                      *
c ************************************************************************
c
      subroutine separate_file
      use patwarp_data
      implicit integer (a-z)
c
c               write coordinates, incidences and blocking, and
c               constraints to separate file or main output file
c
c
c                         seventh question asked (q_num = 7)
c                         'coord., incid-blocking, const. written
c                         to separate files...'
c                         ynanswers(7) = 1, => yes
c
c                         create files for coordinates, incidences-blocking,
c                         and constraints, open...
c
c
      q_num = 7
c
      if( ynanswers(q_num) .ne. 1 ) return
c      
      sep_file_control = .true.
c      
      prefix_name = trim(adjustl( prefix_name ))
      if( prefix_name(1:1) .eq. ' ' ) prefix_name(1:) = 'default'
      last = len( trim(prefix_name) )
c      
      coord_file(1:) = prefix_name(1:last) // '_coords.inp'
      incid_block_file(1:) = prefix_name(1:last) //  '_incid.inp'
      const_file(1:) = prefix_name(1:last) // '_constraints.inp'
c      
      cdf_len = len(trim(coord_file))
      open( unit=cd_file, file=coord_file, status='unknown' )
c      
      inbl_len = len(trim(incid_block_file))
      open( unit=inbl_file, file=incid_block_file, 
     &         status='unknown' )
c     
      ctf_len = len(trim(const_file))
      open( unit=ct_file, file=const_file, status='unknown' )
c
      return
      end
c
c ************************************************************************
c *                                                                      *
c *   routine sample_materials -- writes out sample materials to output  *
c *                               file                                   *
c *   written by: jp                    8-1-2014 updated RHD             *
c *                                                                      *
c ************************************************************************
c
      subroutine sample_materials
      use patwarp_data
      implicit integer (a-z)
c
c
c           sample material title block
c
      write(ofile,9000)
      write(ofile,9200)
      write(ofile,9300)
      write(ofile,9400)
      write(ofile,9500)
      write(ofile,9600)
      write(ofile,9700)
      write(ofile,9800)
      write(ofile,9900)
c
      return
c
 9000 format('c *******************************************',
     &       '**************',/,
     &       'c *                                          ',
     &       '             *',/,
     &       'c *     sample material definitions for use in',
     &       ' warp3d     *',/,
     &       'c *           un-comment and modify as desired',
     &       '            *',/,
     &       'c *                                           '
     &       '            *',/,
     &       'c ********************************************'
     &       '*************',/,'c')
 9200 format('c   rate independent mises material. bilinear',
     &       ' isotropic',/,
     &       'c   hardening, anisotropic thermal expansion',
     &       ' coefficients',/,
     &       'c   -----------------------------------------',
     &       '------------',/,
     &       'c ',/,
     &       'c material steel',/,
     &       'c     properties bilinear e 30000  nu 0.3 yld_pt',
     &       ' 60.0 tan_e 300,',/,
     &       'c                         alphax 3.5e-05,',/,
     &       'c                         alphay 1.5e-05,',/,
     &       'c                         alphaz 2.0e-05,',/,
     &       'c                         alphaxy 1.5e-05,',/,
     &       'c                         alphayz 1.0e-05,',/,
     &       'c                         alphaxz 1.1e-05,',/,
     &       'c                         rho 7.3e-07',/,
     &       'c ',/,'c ')
 9300 format('c   simple rate independent mises with bilinear,',
     &       ' mixed hardening',/,
     &       'c   25% kinematic and 75% isotropic',/,
     &       'c   --------------------------------------------',
     &       '-----------------',/,
     &       'c ',/,
     &       'c material steel',/,
     &       'c     properties bilinear e 30000  nu 0.3 yld_pt',
     &       ' 60.0 tan_e 300,',/,
     &       'c                         beta 0.25 alpha 1.5e-05',
     &       ' rho 7.3e-07',/,
     &       'c ',/,'c ')
 9400 format('c   deformation plasticity model (i.e., nonlinear',
     &       ' elasticity)',/,
     &       'c   ---------------------------------------------',
     &       '------------',/,
     &       'c ',/,
     &       'c material steel',/,
     &       'c    properties deformation  e 30000 nu 0.3 yld_pt 60,',/,
     &       'c               n_power 10 rho 7.3e-07',/,
     &       'c ',/,'c ')
 9500 format('c   rate dependent mises material',/,
     &       'c   -----------------------------',/,
     &       'c ',/,
     &       'c material steel',/,
     &       'c    properties mises  e 30000 nu 0.3 yld_pt 60',
     &       ' n_power 10,',/,
     &       'c                      m_power 5.0 ref_eps 40',
     &       ' rho 7.3e-07',/,
     &       'c ',/,'c ')
 9600 format('c   rate dependent gurson-tvergaard material for',
     &       ' use with',/,
     &       'c   killable elements including (strain-controlled)',
     &       ' nucleation',/,
     &       'c   -----------------------------------------------',
     &       '------',/,
     &       'c ',/,
     &       'c material void_strip',/,
     &       'c   properties gurson e 29000 nu 0.3 yld_pt 58',
     &       ' n_power 10,',/,
     &       'c                     m_power 5.0 ref_eps 40,',/,
     &       'c                     f_0 0.005 q1 1.25 q2 1.0 q3',
     &       ' 1.5625,',/,
     &       'c                     nucleation e_n 0.4 s_n 0.05',
     &       ' f_n 0.50,',/,
     &       'c                     rho 7.3e-07 killable',/,
     &       'c ',/,'c ')
 9700 format('c    rate independent mises model using a segmental',
     &       ' stress-strain',/,
     &       'c    curve',/,
     &       'c    ----------------------------------------------',
     &       '--------------',/,
     &       'c ',/,
     &       'c stress-strain curve 1   $ A533B base metal',/,
     &       'c    0.002210 64.0877,',/,
     &       'c    0.007414 69.9927,',/,
     &       'c    0.012689 77.984,',/,
     &       'c    0.022948 85.4822,',/,
     &       'c    0.033138 90.9935,',/,
     &       'c    0.043241 93.9811,',/,
     &       'c    0.063448 99.9855,',/,
     &       'c    0.123792 109.978,',/,
     &       'c    0.204069 117.984,',/,
     &       'c      50.0      2000',/,
     &       'c ',/,'c ')
 9800 format('c material a533b',/,
     &       'c    properties mises e 28998.953 nu 0.3 curve',
     &        ' 1 rho 7.3e-07',/,
     &        'c ')
 9900 format('c     ***************************************',/,
     &       'c     *       end of sample materials       *',/,
     &       'c     ***************************************',/,
     &       'c ',/,'c ')
c
c
c
      end
c
c
c
c ************************************************************************
c *                                                                      *
c *  routine questions -- Asks questions to be used in creating the      *
c *                       warp3d input file                              *
c *                                                                      *
c ************************************************************************
c
      subroutine questions_old
      use patwarp_data
      implicit integer (a-z)
c
c
      character * 1 ques_ans, default_ans
c
c
c                   asks questions to determine format of warp input
c                   file and other variables
c
c                   integer answers are stored in array 'ynanswers'
c                   ynanswers has a length of 9, possible values
c
c                   ynanswers(q_num) = 1, => answer was yes
c                   ynanswers(q_num) = 0, => answer was no
c                   ynanswers(q_num) = some integer
c
c                   stored as follows:
c
c                       ynanswers(q_num)
c
c                              q_num = 1   create transition elements?
c                                          * used in subroutine trntran *
c                              q_num = 2   run in mpi parallel?
c                                          * used in subroutine trnler *
c                              q_num = 3   block size (threads or mpi)
c                                          * used in subroutine trnler *
c                              q_num = 4   number of processors? (mpi)
c                                          * used in subroutine trnler *
c                              q_num = 5   setup to use ebe preconditioner?
c                                            (parallel)
c                                              -> obsoleted. 5/2013
c                                          * used in subroutine trnler *
c                              q_num = 6   print the new=>old element
c                                            listing?
c                                          * used in subroutine trnler *
c                              q_num = 7   coord, incid-block, const placed
c                                            in separate files?
c                                          * used in subroutine separate_file *
c                              q_num = 8   make updated patran neut file?
c                                          * used in subroutine trnlnf *
c                              q_num = 9   make a patran-readable results
c                                            file?
c                                          * used in subroutine trnlpo *
c
c                       q_num = 6,7,8,9 request file names if yes was
c                                       answered
c
c
c
c
      default_ans = 'y'
      q_num = 1
      write(termot,110)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
c
      default_ans = 'n'
      q_num = 2
      write(termot,120)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
      if( ynanswers(q_num) .eq. 0 ) then  ! threads only
         q_num = 3
         write(termot,131)
         isize = 128
         read(termin,200) block_size
         if( block_size .eq. 0 ) block_size = isize
         ynanswers(q_num) = block_size
      else    ! MPI + threads
         q_num = 4
         write(termot,140)
         read(termin,200) numprocs
         ynanswers(q_num) = numprocs
         q_num = 5
         default_ans = 'n'
         call process_answer(q_num,ques_ans,default_ans)
         isize = 128
         write(termot,132)
         read(termin,200) block_size
         if( block_size .eq. 0 ) block_size = isize
         q_num = 3
         ynanswers(q_num) = iabs( block_size )
      end if
c
c
      default_ans = 'n'
      q_num = 6
      write(termot,160)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
      if( ynanswers(q_num) .eq. 1 ) then
         write(termot,161)
         read(termin,300)elem_print_file
      end if
c
c
      default_ans = 'n'
      q_num = 7
      write(termot,170)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
      if( ynanswers(q_num) .eq. 1 ) then
         write(termot,171)
         read(termin,300) prefix_name
      end if
c
c
      default_ans = 'n'
      q_num = 8
      write(termot,180)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
      if( ynanswers(q_num) .eq. 1 ) then
         write(termot,181)
         read(termin,300)nptfle
      end if
c
c
      default_ans = 'n'
      q_num = 9
      write(termot,190)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
      if( ynanswers(q_num) .eq. 1 ) then
         write(termot,191)
         read(termin,300)resfile
      end if
c
c
      return
c
c
 100  format(a1)
 110  format(/,' >> create transition elements (y/n, default=y)? ',$)
 120  format(/,' >> run with MPI + threads parallel:'
     &       /,'     (y/n, default is no & means threads only)? ', $)
 131  format(/,' >> threads only mode. input element block size:',/,
     &         ' >> (default = 128)? ',$)
 210  format(/,' >> set up to use conjugate-gradient solver:',
     &       /,' >>  (y/n, default = n)? ',$)
 132  format(/,' >> element block size used in domains: ',
     &           '(default = 128)? ',$)
 140  format(/,' >> input the number of MPI processes. will equal',/,
     &         ' >>    the number of model domains (default = 1)? ',$)
 150  format(/,' >> setup to use conjugate gradient solver:',
     &         ' (y/n, default=n)? ',$)
 160  format(/,' >> due to domain decomposition and/or blocking',/,
     &         ' >>    element numbers may change.',/,
     &         ' >>    print the new=>old element listing ',
     &         '(y/n, default=n)? ' ,$)
 161  format(/,' >> name of file for printing element listing? ',
     &        /' >> (press enter to print to screen, or press dd',
     &        /'    to use filename old_to_new_ele.txt)? ',$)
 170  format(/,' >> coordinates, incidences-blocking, and ',
     &         'constraints',/,
     &         ' >>    input data placed in separate files, '
     &         '(y/n, default = n)? ',$)
 171  format(/,' >> specify file prefix  *.coordinates',/,
     &         ' >>                      *.incid_and_blocking',/,
     &         ' >>                      *.constraints',/,
     &         ' >>   (default file names: default.*)? ',$)
 180  format(/,' >> make an updated patran neutral file ',
     &         '(y/n, default=n)? ',$)
 181  format(/,' >> input new patran neutral file name',/,
     &         ' >>   (default: newpat.out)? ',$)
 190  format(/,' >> make a patran-readable element results file',/,
     &         ' >>   which displays the processor and block',/,
     &         ' >>   assignments for the elements ',
     &         '(y/n, default=n)? ',$)
c    rwh -- modification
 191  format(/,' >> input file name (default: metis_pat_info.els)? ',$)
 200  format(i6)
 300  format(a80)
 9000 format(/,10x,'Note: block assignments will prevent elements',
     &       /,10x,'      within a block from sharing a common node.',
     &       /,10x,'      This is termed vectorized blocking.')
c
      end
c
c ************************************************************************
c *                                                                      *
c *  routine questions -- Asks questions to be used in creating the      *
c *                       warp3d input file                              *
c *                                                                      *
c ************************************************************************
c
      subroutine questions
      use patwarp_data
      implicit integer (a-z)
c
c
      character * 1 ques_ans, default_ans
c
c
c                   asks questions to determine format of warp input
c                   file and other variables
c
c                   integer answers are stored in array 'ynanswers'
c                   ynanswers has a length of 9, possible values
c
c                   ynanswers(q_num) = 1, => answer was yes
c                   ynanswers(q_num) = 0, => answer was no
c                   ynanswers(q_num) = some integer
c
c                   stored as follows:
c
c                       ynanswers(q_num)
c
c                              q_num = 1   create transition elements?
c                                          * used in subroutine trntran *
c                              q_num = 2   run in mpi parallel?
c                                          * used in subroutine trnler *
c                              q_num = 3   block size (threads or mpi)
c                                          * used in subroutine trnler *
c                              q_num = 4   number of processors? (mpi)
c                                          * used in subroutine trnler *
c                              q_num = 5   setup to use ebe preconditioner?
c                                            (parallel)
c                                             -> obsolete 5/2013
c                                          * used in subroutine trnler *
c                              q_num = 6   print the new=>old element
c                                            listing?
c                                          * used in subroutine trnler *
c                              q_num = 7   coord, incid-block, const placed
c                                            in separate files?
c                                          * used in subroutine separate_file *
c                              q_num = 8   make updated patran neut file?
c                                          * used in subroutine trnlnf *
c                              q_num = 9   make a patran-readable results
c                                            file?
c                                          * used in subroutine trnlpo *
c
c                       q_num = 6,7,8,9 request file names if yes was
c                                       answered
c
c
c
c
      default_blk_size = 128
 10   continue
      write(termot,400); write(termot,410)
      read(termin,*) ianswer
      if( ianswer .gt. 0 .and. ianswer .le. 2 ) then
          block_method = ianswer
      else
          write(termot,*) '... invalid choice. try again...'
          go to 10
      end if
c
      if( block_method .eq. 1 ) then
        ynanswers(1) = 0
        ynanswers(2) = 0
        ynanswers(3) = default_blk_size
        ynanswers(4) = 1
        ynanswers(5) = 0
        ynanswers(6) = 0
        ynanswers(8) = 0
        ynanswers(9) = 0
      end if
c
      if( block_method .eq. 2 ) then
        ynanswers(1) = 0
        ynanswers(2) = 1
        ynanswers(3) = 128
        ynanswers(4) = 0
        ynanswers(5) = 0
        ynanswers(6) = 0
        ynanswers(8) = 0
        ynanswers(9) = 0
        q_num = 4
        write(termot,140)
        read(termin,200) numprocs
        if( numprocs .eq. 0 ) numprocs = 1
        ynanswers(q_num) = numprocs
        q_num = 3
        write(termot,132)
        read(termin,200) block_size
        if( block_size .eq. 0 ) block_size = default_blk_size
        ynanswers(q_num) = block_size
        default_ans = 'n'
        q_num = 6
        write(termot,162)
        read(termin,100)ques_ans
        call process_answer(q_num,ques_ans,default_ans)
        if( ynanswers(q_num) .eq. 1 ) then
           write(termot,161)
           read(termin,300)elem_print_file
        end if
        default_ans = 'n'
        q_num = 8
        write(termot,180)
        read(termin,100)ques_ans
        call process_answer(q_num,ques_ans,default_ans)
        if( ynanswers(q_num) .eq. 1 ) then
           write(termot,181)
           read(termin,300)nptfle
        end if
        default_ans = 'n'
        q_num = 9
        write(termot,190)
        read(termin,100)ques_ans
        call process_answer(q_num,ques_ans,default_ans)
        if( ynanswers(q_num) .eq. 1 ) then
           write(termot,191)
           read(termin,300)resfile
        end if
      end if
c

      default_ans = 'n'
      q_num = 7
      write(termot,170)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
      if( ynanswers(q_num) .eq. 1 ) then
         write(termot,171)
         read(termin,300) prefix_name
      end if
c
      default_ans = 'y'
      q_num = 1
      write(termot,110)
      read(termin,100)ques_ans
      call process_answer(q_num,ques_ans,default_ans)
c
      return
c
 100  format(a1)
 110  format(/,' >> create transition elements (y/n, default=y)? ',$)
 120  format(/,' >> run with MPI + threads parallel:'
     &       /,'     (y/n, default is no & means threads only)? ', $)
 131  format(/,' >> threads only mode. input element block size:',/,
     &         ' >> (default = 128)? ',$)
 210  format(/,' >> set up to use conjugate-gradient solver:',
     &       /,' >>  (y/n, default = n)? ',$)
 132  format(/,' >> element block size used in domains: ',
     &           '(default = 128)? ',$)
 140  format(/,' >> number of MPI processes (same as',/,
     &         '    the number of model domains (default = 1)? ',$)
 150  format(/,' >> setup to use conjugate gradient solver:',
     &         ' (y/n, default=n)? ',$)
 160  format(/,' >> vectorized blocking may change ',
     &         'element numbers',/,
     &         ' >> print the new => old element listing ',
     &         '(y/n, default=n)? ' ,$)
 161  format(/,' >> name of file for printing element listing? ',
     &        /'    (press enter to print to screen, or press dd',
     &        /'    to use filename old_to_new_ele.txt)? ',$)
 162  format(/,' >> domain decomposition may change ',
     &         'element numbers',/,
     &         '    print the new => old element listing ',
     &         '(y/n, default=n)? ' ,$)
 170  format(/,' >> coordinates, incidences-blocking, and ',
     &         'constraints',/,
     &         '    input data placed in separate files, '
     &         '(y/n, default = n)? ',$)
 171  format(/,' >> specify file prefix  *.coordinates',/,
     &         ' >>                      *.incid_and_blocking',/,
     &         ' >>                      *.constraints',/,
     &         ' >>   (default file names: default.*)? ',$)
 180  format(/,' >> make an updated patran neutral file ',
     &         '(y/n, default=n)? ',$)
 181  format(/,' >> input new patran neutral file name ',
     &         '(default: newpat.out)? ',$)
 190  format(/,' >> make a patran-readable element results file ',
     &         'to display the block',/,
     &         '    assignments for the elements ',
     &         '(y/n, default=n)? ',$)
 191  format(/,' >> input file name (default: metis_pat_info.els)? ',$)
 200  format(i6)
 300  format(a80)
 400  format(/,' >> execution procedure: ',
     & /,6x,'(1) threads-only, Pardiso sparse solver ',
     &      '(direct/iterative)',
     & /,6x,'(2) MPI + threads, hypre solver ' )
 410  format(/,' >> choice: ',$)
c
      end
c
c ************************************************************************
c *                                                                      *
c *  routine process_answer -- Processes the answer to questions used    *
c *                            to specify warp input file format         *
c *                                                                      *
c ************************************************************************
c
      subroutine process_answer(q_num,ques_ans,default_ans)
      use patwarp_data
      implicit integer (a-z)
      character * 1 ques_ans, default_ans
      logical yn_control
c
c
c                    takes the answer to a yes/no question
c                    and stores the answer to be used later
c                    in the program
c
c                    if answer was not a 'y' or 'n' or ' ',
c                    the answer is asked for again
c
c                    when a valid answer is obtained, the yes
c                    or no is stored in the array ynanswers as
c                    a 0 or 1 (no and yes, respectively)
c

      yn_control = .true.
c
      do while(yn_control)
         if( ques_ans .ne. 'y' .and. ques_ans .ne. 'Y' .and.
     &       ques_ans .ne. 'n' .and. ques_ans .ne. 'N' .and.
     &       ques_ans .ne. ' ' )then
            write(termot,9000)
            read(termin,9010)ques_ans
         else
            yn_control = .false.
         end if
      end do
c
c
      if( ques_ans .eq. ' ' ) ques_ans = default_ans
c
      if( ques_ans .eq. 'y' .or. ques_ans .eq. 'Y' )
     &                                ynanswers(q_num) = 1
c
      if( ques_ans .eq. 'n' .or. ques_ans .eq. 'N' )
     &                                ynanswers(q_num) = 0
c
c
c
      return
c
c
c
 9000 format(' >> please enter y or n ? ',$)
 9010 format(a1)
c
      end
c
c
c ************************************************************************
c *                                                                      *
c *  routine trntran -- Determine what kind of hex transition elements   *
c *                     are needed and then add nodes to original 8-node *
c *                     element. ignores non-hex elements                *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trntran
      use patwarp_data
      implicit integer (a-z)
c
c             local variables
c
      integer elem_count(maxnod), nodal_incid(maxnod, maxcon)
      integer nodes_tr(20),nodes_20(20),ele_trnod_count(maxele)
      integer corner_nodes_in_tran(6)
      integer scr_vector(maxtrn), common_face_elem(maxtrn)
      integer num_rep(maxtrn)
c
      logical found_8, found_20, local_debug, hex_element
c
c
      data etype_8_node, etype_9_node, etype_12_node, etype_15_node,
     &     etype_20_node /1, 5, 2, 3, 4/
c
      local_debug = .false.
c
c
c           data structures used to convert 8-node elements to transition
c           elements
c
c             elem_count(node)  : number of elements that share the current
c                                   node
c             nodal_incid(node,maxcon)  : List of elements that share the current
c                                  node, there are elem_count elements.
c             nod_trn_list(node) : = 1 if node is a corner node connecting
c                                  an 8 to a 20-node element.
c                                  = 0 if not
c             ele_trnod_count(element) : number of nodes for the  element which
c                                        are on the transition boundary
c             eletran(element)   : = 1 transition element
c                                  = 0  normal (8 or 20-node element)
c             scr_vector         : auxiliary vector
c
c             transition zone: nodes that belong to 8-node and 20-node element
c                              simultaneously.
c             transition element: element that have some nodes in the
c                              transition zone.
c
c           global variables
c
c             eletrinc(elem)  :contains the new incidences for transition
c                              elements
c             eletrpoin(elem) :index to eletrinc
c
c
c             first find out if user wants transition elements.  If not,
c             return.
c
      made_transition = .false.
 10   continue
c
c
c              the first question (q_num = 1) asked was
c              'create transition elements...?'
c              ynanswers = 1, => yes, ynanwsers = 0, => no
c
c
      q_num = 1
      if( ynanswers(q_num) .eq. 0 ) return
c
c
      made_transition = .true.
c
      write(termot, 1000)
      write(termot, 1100)
c
c             Build nodal incidences. List of all elements connected to each
c             node we use a fixed data structure for simplicity
c
      do snode = 1, numnod
         elem_count(snode) = 0
      end do
c
      do elem = 1, numele
         etype = eletyp(elem)
         num_enodes = elennd(elem)
         incptr = eleipt(elem) - 1
         do enode = 1, num_enodes
            snode = eleinc(incptr + enode)
c
c               if nodes were not renumbered in patran, then we have
c               problems... do a fatal error.
c
            if (snode .gt. numnod) then
               write(termot,9000)
               write(termot,9020)
               stop
            end if
c
c               if structural node number, "snode," appears only
c               once on the current element, record the element
c               number in nodal_incid().
c
            elem_count(snode) = elem_count(snode) + 1
            nodal_incid(snode,elem_count(snode)) = elem
c
c               if "snode" appears multiple times on this element,
c               the next lines undo the previous two lines of code.
c
            if( elem_count(snode) .gt. 1 ) then
               last_elem = nodal_incid(snode,elem_count(snode) - 1)
               if( last_elem .eq. elem ) then
                  nodal_incid(snode,elem_count(snode)) = 0
                  elem_count(snode) = elem_count(snode) - 1
               end if
            end if
c
c               if maximum allowable number of elements per node
c               has been exceeded, do a fatal error.
c
            if ( elem_count(snode) .gt. maxcon ) then
               write(termot,9000)
               write(termot,9050) snode
               stop
            end if
c
         end do
      end do
c
c             construct nod_trn_list vector.
c             nod_trn_list(node)=1 if the structure node is a corner node
c             between an 8 and a 20-node element
c
      do snode = 1, numnod
         num_elems = elem_count(snode)
         found_8  = .false.
         found_20 = .false.
	 nod_trn_list(snode) = 0
         do j = 1, num_elems
            elem = nodal_incid(snode,j)
            etype = eletyp(elem)
            hex_element = etype .le. 5
            if ( .not. hex_element ) cycle
            if ( etype .eq. etype_8_node ) found_8  = .true.
            if ( etype .eq. etype_20_node) found_20 = .true.
         end do
         if (found_8 .and. found_20 )  nod_trn_list(snode) = 1
      end do
c
c             1. loop over all elements to build count for each element of the
c                number of corner nodes shared between an 8-node and a 20-node
c                element (ele_trnod_count)
c                for each element, check how many nodes are in the transition
c                zone.
c
      write(termot, 1200)
c
      do elem = 1, numele
         ele_trnod_count(elem) = 0
         eletran(elem) = 0
      end do
c
      do elem = 1, numele
         etype = eletyp(elem)
         hex_element = etype .le. 5
         if ( .not. hex_element ) cycle
         num_enodes = elennd(elem)
         incptr = eleipt(elem) - 1
         do enode = 1, num_enodes
            snode = eleinc(incptr + enode)
            if ( nod_trn_list(snode) .eq. 1) then
              ele_trnod_count(elem) = ele_trnod_count(elem) + 1
            end if
         end do
      end do
c
c             2. loop over all elements. based on the number of corner nodes for
c                a transition element, classify the type of transition element.
c                mark 8-node element as a transition element (eletran)
c                Otherwise we do not recognize the type of transition
c
c                If element has 2 corner nodes in transition zone then
c                   it is a  9-node transition element
c                If element has 4 corner nodes in transition zone then
c                   it is a 12-node transition element
c                If element has 6 corner nodes in transition zone then
c                   it is a 15-node transition element
c
      number_of_transition_elem = 0
      do elem = 1, numele
         etype = eletyp(elem)
         hex_element = etype .le. 5
         if ( .not. hex_element ) cycle
         if ( etype .eq. etype_20_node ) cycle
c
         num_corner_nodes = ele_trnod_count(elem)
         incptr = eleipt(elem) - 1
c
         if  (num_corner_nodes .eq. 2) then
             if (local_debug)
     &       write(termot,2000)elem
             eletyp(elem)  = etype_9_node
             eletran(elem) = 1
             number_of_transition_elem = number_of_transition_elem + 1
         elseif (num_corner_nodes .eq. 4) then
             if (local_debug)
     &       write(termot,2100)elem
             eletyp(elem)  = etype_12_node
             eletran(elem) = 1
             number_of_transition_elem = number_of_transition_elem + 1
         elseif (num_corner_nodes .eq. 6) then
             if (local_debug)
     &       write(termot,2200)elem
             eletyp(elem)  = etype_15_node
             eletran(elem) = 1
             number_of_transition_elem = number_of_transition_elem + 1
         else
             if (num_corner_nodes .ne. 0) then
                write(termot,9000)
                write(termot, 9100) elem
                stop
             end if
         end if
      end do
c
c             write the number of 8-node elements that became
c             transition elements.
c
      write(termot,1400)  number_of_transition_elem
c
c             modify incidences to convert 8-node elements to the
c             required transition element. non-transition elements
c             are skipped
c
      write(termot, 1300)
      incid_index = 0
      do elem = 1, numele
         eletrpoin(elem) = incid_index
         num_neighbors = 0
         if (eletran(elem) .eq. 0) goto 300
c
c             for each transition element,a list containing
c             the corner nodes on the transition zone is formed.
c             this vector is corner_nodes_in_tran.
c             check consistency of number of corners nodes. if bad, stop.
c
         num_enodes = elennd(elem)
         incptr = eleipt(elem) - 1
         counter = 0
         do enode = 1, num_enodes
            snode = eleinc(incptr + enode)
            if (nod_trn_list(snode) .eq. 1) then
               counter = counter + 1
               corner_nodes_in_tran(counter) = snode
            end if
         end do
         if (counter .ne. ele_trnod_count(elem)) then
             write(termot,9000)
             write(termot,9500) elem
             stop
         end if
c
c             for each transition element: loop over its corner transition
c             nodes; for each corner transition node, put all elements
c             (nodal_incid) except the current element (elem) in the vector
c             scr_vector.
c             from this list we choose the element appearing the most number
c             of times. this element will be the element that shares a face
c             with the current transition element. for 15-node element we
c             have 2 such elements.
c
c             share_1: element appearing most
c             share_2: second most appearing element (valid for 15-node element)
c             both must be 20-node elements
c             for 9-node element, we can have several elements appearing an
c             equal number of times, in this case we choose the first 20-node.
c
c             according to the corner transition number of nodes (2, 4 or 6) we
c             form the new incidences in routines element_12_15 or element_9.
c
         do j = 1, ele_trnod_count(elem)
            corner_node = corner_nodes_in_tran(j)
            number_of_elements = elem_count(corner_node)
            do k = 1, number_of_elements
               selem = nodal_incid(corner_node,k)
               if (selem .ne. elem) then
                   num_neighbors = num_neighbors + 1
                   if (num_neighbors .gt. maxtrn) then
                      write(termot,9400) elem
                      stop
                   end if
                   scr_vector(num_neighbors)= selem
               end if
            end do
         end do
c
c             tran_freq provides the element numbers share_1 and share_2
c
         call tran_freq(scr_vector, num_neighbors, common_face_elem,
     &                  num_of_elem, num_rep)
         share_1 = common_face_elem(1)
c
c             9-node case: check if share_1 is quadratic
c
         if (ele_trnod_count(elem) .eq. 2) then
            do j = 1, num_of_elem
               if (elennd(common_face_elem(j)) .eq. 20) then
                  share_1 = common_face_elem(j)
                  if ( num_rep(j) .ne. 2) then
                     write(termot,9000)
                     write(termot,9600) elem, num_rep(j)
                  end if
                  goto 200
               end if
            end do
 200        continue
            number_faces = 0
	    call element_9(elem, share_1,nodes_20,nodes_tr,incid_index)
c
c             12-node case : check if number of repetitions is four
c
         elseif (ele_trnod_count(elem) .eq. 4) then
            if ( num_rep(1) .ne. 4) then
                write(termot,9000)
                write(termot,9700) elem, num_rep(1)
            end if
            number_faces = 1
	    call element_12_15(elem,number_faces,share_1,0,incid_index)
c
c             15-node case : check if shared_1 and shared_2 are repeated
c                            exactly four times.
c
         elseif (ele_trnod_count(elem) .eq. 6) then
            if ( num_rep(1) .ne. 4 .or. num_rep(2) .ne. 4) then
                write(termot,9000)
                write(termot,9800) elem,num_rep(1)
            end if
            share_2 = common_face_elem(2)
            number_faces = 2
	    call element_12_15(elem,number_faces,share_1,share_2,
     &                         incid_index)
         else
            write(termot,9000)
            write(termot,9300) elem
	    stop
         end if
  300    continue
      end do
c
      write(termot, 1500)
      return
c
 1000 format(/,8x'>> begin processing hex transition elements    ')
 1100 format(7x,'        > building nodal incidences      ')
 1200 format(7x,'        > locating transition elements   ')
 1300 format(7x,'        > converting 8-node elements into transition',
     &                   ' elements')
 1400 format(7x,'           > ',i8,' transition elements needed')
 1500 format(7x,'        > done building transition elements')
 2000 format(' Element',i7,' is a 9 node transition')
 2100 format(' Element',i7,' is a 12 node transition')
 2200 format(' Element',i7,' is a 15 node transition')
 9000 format(1x,'>>>> Fatal error: routine trtran:')
 9020 format(1x,'>>>> There are nodes with no elements attached.',
     &     'Renumber nodes',/,4x,'and elements in patran and try again')
 9050 format(8x,'overflow in number of elements reaching node', i8)
 9100 format(8x,'in element',i5,'transition type not available in warp')
 9300 format(8x,'in element',i5,'more than 3 faces are quadratic')
 9400 format(8x,'overflow in the number of elements @',i8)
 9500 format(8x,'in element ',i5,'transition not available in warp')
 9600 format(8x,'in  9-node element ',i8,' have',i8,' transition nodes')
 9700 format(8x,'in 12-node element ',i8,' have',i8,' transition nodes')
 9800 format(8x,'in 15-node element ',i8,' have',i8,' transition nodes')
c
      end
c
c ************************************************************************
c *                                                                      *
c *                     routine element_12_15                            *
c *                                                                      *
c *   this subroutine redefines the incidences for a 8-node to convert   *
c *   it in to a 12 or 15 node transition element                        *
c *                                                                      *
c ************************************************************************
c
c
      subroutine element_12_15(elem, number_faces, share_1, share_2,
     &     incid_index)
      use patwarp_data
      implicit integer (a-z)
c
c     local variables
c
      integer scr_vector(maxtrn), scr_value(maxtrn)
      integer nodes_tr(20),nodes_20(20), list_tr(6), list_20(6)
      integer vertice(8,4), list(4,4), list2(3,4)
      logical local_debug
c
      data vertice/1,  2,  3,  4,  5,  6,  7,  8,
     &             2,  3,  4,  1,  6,  7,  8,  5,
     &             3,  4,  1,  2,  7,  8,  5,  6,
     &             4,  1,  2,  3,  8,  5,  6,  7/
c
      data list  / 9, 10, 11, 12,
     &            10, 11, 12,  9,
     &            11, 12,  9, 10,
     &            12,  9, 10, 11/
c
      data list2 / 15, 14, 16,
     &             16, 15, 13,
     &             13, 16, 14,
     &             14, 13, 15/

c
      local_debug = .false.
c
c             global variables
c
c             eletrinc(elem) :contains the new incidences for transition
c                             elements
c             eletrpoin(elem):index to eletrinc
c
c
c                o-----*                        *-----o
c               /     /|                       /     /|
c              o-----* |<< face_tr            *-----o |
c              |     | *            face_20>> |     | o
c              |     |/                       |     |/
c              o-----*                        *-----o
c
c             Transition element              20-node element
c                   (elem)                     (share_1)
c             * nodes in transition zone
c
c
c	      number_faces = 1 : 12-node element
c             number_faces = 2 : 15-node element
c             all nodes of transition element (elem) will be compared with
c             all nodes of its 20-node element(s) share_1 (and share_2).
c             common (local) nodes will be stored in vectors list_tr and
c             list_20.
c             With  list_tr and list_20 we represent two 4-digit numbers face_tr
c             and face_20, which will be reordered according to warp numbering
c             of nodes, i.e.2367 should be 2673 and so forth.
c             for the 15-node element we have to determine in addition the
c             common vertice of the quadratic faces. we form 2 sets of 12
c             incidences each, the first 8 are the same but in different
c             order; we add the four different incidences and have a vector
c             of 16 components, then compare incidences 9 to 16 to find the
c             repeated incidence and find the common vertice; with this
c             vertice, the element is renumbered in such way that we have
c             quadratic faces as faces 1 and 2.
c
c             loop over the  number of faces
c                1. form incidences vector for transition element (elem)
c                2. form incidences vector for quadratic element (share)
c                3. compare the vectors of incidences and determine the
c                     common nodes (common faces) in each element
c                     (elem and share)
c                4. determine the 12 incidences
c                5. if element has one quadratic face then save new
c                   incidences and exit
c                6. if element has two quadratic faces repeat 1,2,3,4
c                   and form a new vector of incidences. from this new
c                   vector extract the last four to add the old vector
c                   and form a vector of 16 incidences.
c
c                   {1,2,3,4,5,6,7,8}  +  {9,10,11,12}  +  {13,14,15,16}
c                  from original 8-node    from share_1     from share_2
c
c                7. from the vector of 16 incidences we look for the
c                   repeated incidence that is on the common vertice of
c                   the two quadratic faces
c                8. remove the repeated index and renumber the nodes in
c                   order to have the first and second faces as quadratic.

      num_enode_trn = elennd(elem)
      do nface = 1, number_faces
c
c            incidences of transition element
c
         incptr = eleipt(elem) - 1
         do enode = 1, num_enode_trn
            nodes_tr(enode) = eleinc(incptr + enode)
            scr_vector(enode) = nodes_tr(enode)
         end do
c
c            incidences of 20-node element (share)
c
         incptr = eleipt(share_1) - 1
         num_enode_20 = elennd(share_1)
	 if (num_enode_20 .ne. 20) then
            write(termot,9000)
	    write(termot,9300) elem,share_1
	 end if
         do enode = 1, num_enode_20
            nodes_20(enode) = eleinc(incptr + enode)
         end do
c
c             renumber the 20-node according to warp
c
         call order1215(nodes_20, num_enode_20, 2)
c
c             searching for common nodes (common face)
c
         common_nodes = 0
         do enode_tr = 1, num_enode_trn
            if (nod_trn_list(nodes_tr(enode_tr)) .eq. 1) then
               do enode_20 = 1, num_enode_20
                  if ( nodes_tr(enode_tr). eq . nodes_20(enode_20)) then
                      common_nodes = common_nodes + 1
                      list_tr(common_nodes) = enode_tr
                      list_20(common_nodes) = enode_20
                  end if
               end do
             end if
         end do
c
c             incidences 1-8 renumbered and new incidences 9-12
c
	 call incidences12(elem,share_1,list_20,list_tr,scr_vector,
     &	                   nodes_20,nodes_tr)
c
c             additional things for 15-node element
c
         if (number_faces .eq. 2) then
            if (nface .eq. 1) then
               do enode = 1, 12
                  scr_value(enode) = nodes_tr(enode)
               end do
               share_1 = share_2
            end if
         end if
      end do
c
c             save new incidences for 12-node element and update index
c
      if (number_faces .eq. 1) then
         eletrpoin(elem) = incid_index
         do enode = 1, 12
            eletrinc(incid_index + enode) = nodes_tr(enode)
         end do
         incid_index = incid_index + 12
c
         if(local_debug)
     &       write(termot,'(12i5)')(nodes_tr(enode),enode = 1,12)
c
         goto 100
      end if
c
c             15-node element: second quadratic face
c
      scr_value(13) = nodes_tr(9)
      scr_value(14) = nodes_tr(10)
      scr_value(15) = nodes_tr(11)
      scr_value(16) = nodes_tr(12)
c
      if(local_debug) write(termot,'(16i5)')(scr_value(k), k = 1, 16)
c
c             comparing incidences {9,10,11,12} against {13,14,15,16}
c             to find common vertice
c
      vert_common = 0
      do enode1 = 9, 12
         do enode2 = 13, 16
            if ( scr_value(enode1) .eq. scr_value(enode2) ) then
                vert_common = enode1 - 8
                vert_2nd    = enode2 - 12
            end if
         end do
      end do
c
c             if vert_common=0 then there is no common vertice and therefore
c             element is not the 15-node described in warp
c
      if (vert_common .eq. 0) then
          write(termot,9000)
          write(termot,9600) elem, vert_common
	  stop
      end if
c
c             first eight nodes are the same but in different order
c
c             incidences from 1 to 8
c
      do enode = 1, 8
         nodes_tr(enode) = scr_value(vertice(enode,vert_common))
      end do
c      if(number_faces.eq.2) then
c        write(termot,9700)elem,vert_common
c        stop
c      end if
c
c             incidences from 9 to 12
c
      nodes_tr(9 ) = scr_value(list(1,vert_common))
      nodes_tr(10) = scr_value(list(2,vert_common))
      nodes_tr(11) = scr_value(list(3,vert_common))
      nodes_tr(12) = scr_value(list(4,vert_common))
c
c             incidences from 13 to 15
c
      nodes_tr(13) = scr_value(list2(1,vert_2nd))
      nodes_tr(14) = scr_value(list2(2,vert_2nd))
      nodes_tr(15) = scr_value(list2(3,vert_2nd))
c
c             save new incidences for 15-node element and update index
c
      eletrpoin(elem) = incid_index
      do enode = 1, 15
         eletrinc(incid_index + enode) = nodes_tr(enode)
      end do
      incid_index = incid_index + 15
c
      if (local_debug) write(termot,'(15i5)')(nodes_tr(k), k = 1,15)
c
  100 continue
c
 9000 format(1x,'>>>>Fatal error: routine element_12_15:')
 9300 format(5x,'in element',i5,' share nodes with element', i5)
 9600 format(5x,'in element',i5,' vertice. ', i5)
 9700 format(5x,'Elem 15 =',i5,1x,i5)
      return
      end
c
c ************************************************************************
c *                                                                      *
c *                      routine element 9                               *
c *                                                                      *
c *           this subroutine redefines the incidences for a 8-node      *
c *           to convert it in to a 9-node transition element            *
c *                                                                      *
c ************************************************************************
c
c
      subroutine element_9(elem,share_1,nodes_20,nodes_tr,incid_index)
c
      use patwarp_data
      implicit integer (a-z)
      integer nodes_20(1), nodes_tr(1)
c
c     local variables
c
      integer list(12), scr_vector(20), list_tr(6), list_20(6)
      integer vert(9,12), vert_typ(12), num_dummy
c
      logical local_debug, flag
c
      data vert  /1, 2, 3, 4, 5, 6, 7, 8, 9,
     &            2, 3, 4, 1, 6, 7, 8, 5, 9,
     &            3, 4, 1, 2, 7, 8, 5, 6, 9,
     &            4, 1, 2, 3, 8, 5, 6, 7, 9,
     &            6, 5, 8, 7, 2, 1, 3, 4, 9,
     &            7, 6, 5, 8, 3, 2, 1, 4, 9,
     &            7, 8, 4, 3, 6, 5, 1, 2, 9,
     &            5, 8, 7, 6, 1, 4, 2, 3, 9,
     &            2, 6, 7, 3, 1, 5, 8, 4, 9,
     &            7, 3, 2, 6, 8, 4, 1, 5, 9,
     &            4, 8, 5, 1, 3, 7, 6, 2, 9,
     &            5, 1, 4, 8, 6, 2, 3, 7, 9/
c
      data list  /9,10,11,12,13,14,15,16,18,19,20,17/
c
      data vert_typ /12,23,34,14,56,67,78,58,26,37,48,15/
c
      local_debug = .false.
c
c             same procedure as element_12_15 but instead
c             of forming the four node common-face we form
c             the 2 node common-vertex
c
c             incidences of transition element
c
      incptr = eleipt(elem) - 1
      num_enode_trn = elennd(elem)
      do enode = 1, num_enode_trn
         nodes_tr(enode) = eleinc(incptr + enode)
         scr_vector(enode) = nodes_tr(enode)
      end do
c
c             incidences 20-node element (share_1)
c
      incptr = eleipt(share_1) - 1
      num_enode_20 = elennd(share_1)
      do enode = 1, num_enode_20
         nodes_20(enode) = eleinc(incptr + enode)
      end do
c
c             renumber the 20-node according to warp
c
      num_dummy = 20
      call order1215(nodes_20, num_dummy, 2)
c
c             looking for the common vertex
c
      common_nodes = 0
      do enode_tr = 1, num_enode_trn
         do enode_20 = 1, num_enode_20
            if (nod_trn_list(nodes_tr(enode_tr)) .eq. 1) then
               if ( nodes_tr(enode_tr). eq . nodes_20(enode_20)) then
                  common_nodes = common_nodes + 1
                  list_tr(common_nodes) = enode_tr
                  list_20(common_nodes) = enode_20
               end if
            end if
         end do
      end do
      if ( list_tr(1) .gt. list_tr(2) ) then
         call swap(list_tr(1),list_tr(2))
      end if
      if ( list_20(1) .gt. list_20(2) ) then
         call swap(list_20(1),list_20(2))
      end if
c
c             if the vertex has no 2 nodes (exactly) something is wrong
c
      if (common_nodes .ne. 2 ) then
         write(termot,9000)
         write(termot,9500) elem
         stop
      end if
c
c             form the 2-digit number that represent the vertices
c             vert_20 : Vertice of 20-node element
c             vert_tr : Vertice of transition element
c
      vert_20 = 10*list_20(1) + list_20(2)
      vert_tr = 10*list_tr(1) + list_tr(2)
c
c             selecting the 20-node element vertex from the 12 available
c
      flag = .false.
      do i = 1, 12
         if (vert_20 .eq. vert_typ(i)) then
            case_20 = i
            flag = .true.
            goto 100
         end if
      end do
  100 continue
      if ( .not. flag ) then
         write(termot,9000)
         write(termot,9400)  elem, share_1, vert_20
         stop
      end if
c
c            9th node is determined
c
      scr_vector(9) = nodes_20(list(case_20))
c
c             selecting the transition element vertex from the 12 available
c
      flag = .false.
      do i = 1, 12
         if (vert_tr .eq. vert_typ(i)) then
            case_tr = i
            flag = .true.
            goto 200
         end if
      end do
  200 continue
c
c                 if that vertice does not exit...
c
      if ( .not.  flag ) then
         write(termot,9000)
         write(termot,9300)  elem, vert_tr, share_1
         stop
      end if
c
c             save new incidences for 9-node element and update index
c
      do enode = 1, 9
         nodes_tr(enode) = scr_vector(vert(enode,case_tr))
      end do
      eletrpoin(elem) = incid_index
      do enode = 1, 9
         eletrinc(incid_index + enode) = nodes_tr(enode)
      end do
      incid_index = incid_index + 9
c
      if(local_debug)
     &       write(termot,'(9i5)') (nodes_tr(enode), enode = 1, 9)
c
 9000 format(1x,'>>>>Fatal error: routine element_9:')
 9300 format(1x,'in element',i5,' vertex ',i5,'with neighbor elem ',i5)
 9400 format(1x,'in element',i5,' neighbor elem ',i5,'vertex ',i5)
 9500 format(5x,'in element',i5,'common vertex must have only 2 nodes')
 9600 format(5x,'in element',i5,' vertice. ',i5)
      return
      end
c ************************************************************************
c *                                                                      *
c *                        routine incidences12                          *
c *                                                                      *
c *     this subroutine returns the 12 first incidences for a 12 or      *
c *     a 15-node element                                                *
c *                                                                      *
c ************************************************************************
c
c
      subroutine incidences12(elem,share,list_20,list_tr,scr_vector,
     &      	        nodes_20,nodes_tr)
      use patwarp_data
      implicit integer (a-z)
c
      integer scr_vector(*), nodes_20(*), nodes_tr(*)
      integer list_tr(*), list_20(*)
c
c             local variables
c
      integer  face(12,6), list(4, 24), face_typ_20(24), face_typ_tr(6)
      logical  flag
c
      data face/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
     &          1, 5, 6, 2, 4, 8, 7, 3, 9, 10, 11, 12,
     &          2, 6, 7, 3, 1, 5, 8, 4, 9, 10, 11, 12,
     &          1, 4, 8, 5, 2, 3, 7, 6, 9, 10, 11, 12,
     &          3, 7, 8, 4, 2, 6, 5, 1, 9, 10, 11, 12,
     &          5, 8, 7, 6, 2, 1, 4, 3, 9, 10, 11, 12/
c
      data list/ 9,12,11,10,12,11,10, 9,11,10, 9,12,10, 9,12,11,
     &          13,17, 9,18,17, 9,18,13, 9,18,13,17,18,13,17, 9,
     &          18,10,19,14,10,19,14,18,19,14,18,10,14,18,10,19,
     &          16,20,12,17,20,12,17,16,12,17,16,20,17,16,20,12,
     &          15,19,11,20,19,11,20,15,11,20,15,19,20,15,19,11,
     &          13,14,15,16,14,15,16,13,15,16,13,14,16,13,14,15/
c
      data face_typ_20/  2143,1432,4321,3214,
     &                   6512,5126,1265,2651,
     &                   6237,2376,3762,7623,
     &                   5841,8415,4158,1584,
     &                   8734,7348,3487,4873,
     &                   5678,6785,7856,8567/
c
      data face_typ_tr/ 1234,1562,2673,1485,3784,5876/
c
c
c
c               forming the 4-digit face number of the transition
c               element (elem) and its 'share' 20-node element.
c
c
      face_tr=1000*list_tr(1)+100*list_tr(2)+10*list_tr(3)+list_tr(4)
      face_20=1000*list_20(1)+100*list_20(2)+10*list_20(3)+list_20(4)
c
c               because faces were obtained in sequential order we need
c               to renumber them according to warp
c               i.e. face 1234 could be 2143,1432,4321 or 3214
c
      if     (face_tr .eq. 1458) then
	      call swap(list_20(4), list_20(3))
	      call swap(list_tr(4), list_tr(3))
      elseif (face_tr .eq. 3478 .or. face_tr .eq. 1256 .or.
     &        face_tr .eq. 2367) then
	      call swap(list_20(3), list_20(2))
	      call swap(list_tr(3), list_tr(2))
	      call swap(list_20(3), list_20(4))
	      call swap(list_tr(3), list_tr(4))
      elseif (face_tr .eq. 5678) then
	      call swap(list_20(4), list_20(2))
	      call swap(list_tr(4), list_tr(2))
      elseif (face_tr .eq. 1234) then
              continue
      else
              write(termot,9000)
              write(termot,9300) elem, face_tr
              stop
      end if
c
      face_tr=1000*list_tr(1)+100*list_tr(2)+10*list_tr(3)+list_tr(4)
      face_20=1000*list_20(1)+100*list_20(2)+10*list_20(3)+list_20(4)
c
c             selecting the 20-node element face from the
c             24 available (face_typ_20)
c
      flag = .false.
      do i = 1, 24
         if (face_20 .eq. face_typ_20(i)) then
            case_20 = i
            flag = .true.
            goto 100
         end if
      end do
  100 continue
c
c             check if 20-node face is one of the 24 available
c
      if (.not. flag ) then
         write(termot,9000)
         write(termot,9400) elem, share, face_20
         stop
      end if
c
c             selecting the transition element face from
c             the 6 available (face_typ_tr)
c
      flag = .false.
      do i = 1, 6
         if (face_tr .eq. face_typ_tr(i)) then
            case_tr = i
            flag = .true.
            goto 200
         end if
      end do
  200 continue
c
c             check if transition element face is one of the 6 available
c
      if ( .not.flag ) then
         write(termot,9000)
         write(termot,9500) elem, face_tr
         stop
      end if
c
c             save intermediate new nodes 9, 10, 11, and 12
c
      do enode = 9, 12
         scr_vector(enode) = nodes_20(list(enode-8,case_20))
      end do
      do enode = 1, 12
         nodes_tr(enode) = scr_vector(face(enode,case_tr))
      end do
c
 9000 format(1x,'>>>>Fatal error: routine incidences_12:')
 9300 format(5x,'in element',i5,' face ',i5)
 9400 format(5x,'in element',i5,' neighbor elem. ',i5,' face ',i5)
 9500 format(5x,'in element',i5,' face with nodes ', i5)
      return
      end
c
c ************************************************************************
c *                                                                      *
c *                      routine tran_freq                               *
c *                                                                      *
c *           frequency subroutine returns the most repeated numbers of  *
c *           a set numbers                                              *
c *                                                                      *
c ************************************************************************
c
c
      subroutine tran_freq(list, num_data, sort_list, number_elem,
     &                     freq)
      use patwarp_data
      implicit integer (a-z)
      integer list(1), num_data, sort_list(1), freq(1)
      logical local_debug

      local_debug = .false.

      if (local_debug) then
         write(termot,*) num_data
         write(termot,*) (list(i), i = 1, num_data)
      end if
c
c            cleaning vector of repetitions
c
      do i = 1, maxtrn
         freq(i) = 0
      end do
c
c            number of different elements
c
      number_elem = 0
      do i = 1, num_data
         if (number_elem .ne. 0) then
            do j = 1, number_elem
               if (list(i) .eq. sort_list(j)) goto 10
            end do
         end if
         number_elem = number_elem + 1
         sort_list(number_elem) = list(i)
         do j = 1, num_data
            if (sort_list(number_elem) .eq. list(j)) then
                freq(number_elem) = freq(number_elem) + 1
            end if
         end do
   10    continue
      end do
c
c     reorder according to most appearing elements first
c
      do i = 1, number_elem
         do j = 1, number_elem
            if ( freq(i) .gt. freq(j) ) then
               call swap(sort_list(j),sort_list(i))
               call swap(freq(j),freq(i))
            end if
         end do
      end do
      if (local_debug) then
         do i = 1, number_elem
            write(termot, *) i, sort_list(i), freq(i)
         end do
      end if
      return
      end
c
c ************************************************************************
c *                                                                      *
c *          routine swap                                                *
c *                                                                      *
c *          this subroutine interchange two values a and b              *
c *                                                                      *
c ************************************************************************
c
c
      subroutine swap(a,b)
      implicit integer (a-z)
c
      temp = a
      a    = b
      b    = temp
c
      return
      end
c
c ************************************************************************
c *                                                                      *
c *   routine  trnler -- originally subroutine elrodr                    *
c *                                                                      *
c *     this subroutine reorders the elements into blocks following a    *
c *     set of rules (same type, nonlinearity, material model, etc.)     *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnler
      use patwarp_data
      implicit integer (a-z)
      dimension iblock(numele)
      logical print_to_file, start_blk, parallel_ebe
c
      data numprocs / 1 /
c
      if (debug) write (termot,*) '>>>> inside trnler'
      print_file_num = 13
c
c
c                    block_method = 1 => just use simple
c                                        WARP3D built-in automatic
c                                        blocking
c
      if( block_method .eq. 1 ) then
        write(termot,9500)
        do i = 1, numele
          elreno(i) = i
        end do
        blksz = ynanswers(3)
        return
      end if

      write(termot,9400)
c
c
 10   continue
c
c                     second question asked (q_num = 2)
c                     'run in MPI + threads parallel...?'
c                     ynanswers(2) = 1, => yes
c                     ynanswers(2) = 0, => no
c
      q_num = 2
      if( ynanswers(q_num) .eq. 0 ) then
         parallel = .false.
       else
         parallel = .true.
c
c                     fourth question asked (q_num = 4)
c                     'number of MPI processors...?'
c                     ynanswers(4) = number of processors
c
         q_num = 4
         numprocs = ynanswers(q_num)
         if ( numprocs .eq. 0 ) numprocs = 1
c
c                     fifth question asked (q_num = 5)
c                     'setup for conjugate gradient solver...?'
c                     we now assume ebe (HW) pre-conditioner will be used
c                     ynanswers(5) = 1, => yes
c                     ynanswers(5) = 0, => no
c
         q_num = 5
         if ( ynanswers(q_num) .eq. 1 )then
            parallel_ebe = .true.
         else
            parallel_ebe = .false.
         end if
      end if
c
c                      when in threads only:
c                      third question asked (q_num = 3)
c                      'threads only mode, input block size...?'
c                      (negative means set up for conjugate gradient
c                      solver -> vectorized blocking
c                      ynanswers(3) = block size
c
c                      when in MPI_threads parallel:
c                      fifth question asked, but stored (q_num = 3)
c                      questions routine sets block size negative if
c                      user also wnts set up for conjugate gradient solver.
c                      we assume HW preconditioner will be used ->
c                      vectorized blocking.
c
c                      If blocking is given as postive, then use
c                      scalar blocking; otherwise use vectorized blocking.
c
c
      q_num = 3
      blksz = ynanswers(q_num)
c
      if ( blksz .eq. 0 ) blksz = 128
      if ( blksz .gt. 0 ) then
         write (termot,9080)
         scalar = .true.
      else
         blksz = - blksz
         scalar = .false.
      end if
c
c
c                       determine blocks of non-conflicting, similar
c                       elements.
c
c                       if this is MPI, first determine the element
c                       processor assignment.  If the ebe preconditioner
c                       is to be used, determine the processors that
c                       own the neighboring elements -- this is used as
c                       an additional criterion for forming blocking groups
c
      if ( parallel .and. (numprocs .gt. 1) ) then
         call trnlmetis (numprocs)
         if ( parallel_ebe ) then
            call trnlproc_access
         else
            do i = 1, numele
               elem_access(i) = 1
            end do
         end if
      else
         do i = 1, numele
            eleproc(i) = 1
            elem_access(i) = 1
         end do
      end if
c
c                       find groups of similar elements; all the elements
c                       in a group must be the same type and have the same
c                       configuration number. If we are constructing a model
c                       to run in parallel, then all elements in a group
c                       must also be on the same processor.
c
      call getgrp
c
c                       sort each group of similar elements into
c                       non-conflicting blocks and store these blocks
c                       into the appropriate global data structure.
c
c                       on a scalar machine, nonconflicting blocks are
c                       not necessary.  However, the blocks must still
c                       contain similar elements, so if there are
c                       multiple configuration numbers or element
c                       types, then we will still have to renumber
c                       the elements (warp requires all the elements
c                       in a block to have the same material properties,
c                       element types, etc.).  Output a warning
c                       if we have scalar blocking but must renumber
c                       anyway because we have multiple configuration
c                       numbers.
c
c
      nelblk = 0
      blk    = 0
c
      if ( scalar .and. (numgrp .gt. 1) ) write (termot,9100)
      if ( .not. scalar ) write (termot,9200)
c
c                       looping over groups
c
      if (debug) write (termot,*) 'looping over groups'
      do i = 1, numgrp
c
c                       ==   setting group properties
c
         nnode = elennd(grphead(i))
         nel   = grpnum(i)
         if (debug)
     &        write (termot,'("grp:",i4," nnode:",i2," elms:",i5)')
     &        i,nnode, nel
c
c                       ==  figure out how many blocks you need.
c                              if scalar -- just divide number of
c                                  elements by blocksize.
c                              if vector -- use the red-black
c                                  algorithm to assign the blocks
c
c
         lstblk = nelblk
         if ( scalar ) then
            nelblk = nelblk + (nel/blksz)
            if ( mod(nel,blksz) .ne. 0 ) nelblk = nelblk + 1
         else
            call gtblrb( nel, nnode, blksz, i, iblock )
         end if
         if (debug) write (termot,'(" need ",i4," blocks.")') nelblk
c
c                       ==  now build elblks structure -- an entry for
c                       ==  each block
c
c                                zero out array
c
         do j = lstblk+1, nelblk
            elblks(0,j) = 0
         end do
c
c                                if scalar, then fill the blocks in
c                                current group in same order as group
c
         if ( scalar ) then
c
c                                   fill next block
c
            blk = blk + 1
            start_blk = .true.
c
c                                   loop over elements in group. If we
c                                   fill current block, fill the next block.
c                                   if we are constructing an input file for
c                                   parallel computers, store the processor
c                                   number which owns the elements in the block
c                                   (in blkproc).
c
	    nxtele = grphead(i)
            do j = 1,nel
               elem = nxtele
	       if ( elem .eq. 0) then
		  write (termot,*) '>>>> ptr error in trnler'
		  stop
               end if
	       nxtele = grplst(elem)
               if ( parallel .and. start_blk ) then
                  blkproc ( blk ) = eleproc(elem) - 1
                  start_blk = .false.
               endif
               elblks(0,blk) = elblks(0,blk)+1
               elblks(elblks(0,blk),blk) = elem
               if ( (elblks(0,blk).eq.blksz) .and. (j.ne.nel) ) then
                  blk = blk + 1
                  start_blk = .true.
                  if ( blk .gt. mxnmbl ) then
                     write(termot,9000)
                     stop
                  end if
               end if
            end do
c
c                                if vectorized, then fill in blocks as
c                                decided in the red-block algorithm,
c                                using iblock.
c
         else
	    nxtele = grphead(i)
            do j = 1, nel
               elem = nxtele
	       if ( elem .eq. 0) then
		  write (termot,*) '>> ptr error in trnler'
		  stop
               end if
	       nxtele = grplst(elem)
               blk = iblock(elem)
               elblks(0,blk) = elblks(0,blk)+1
               elblks(elblks(0,blk),blk) = elem
            end do
         end if
c
c                        ==  end of looping over groups
c
      end do
c
c                        if debug, then print the eleblk structure.
c
      if (debug) then
         write (termot,*) 'the eleblk structure:'
         do blk = 1, nelblk
            write (termot,'("   in block ",i4," are elements:")')blk
            do elems = 1, elblks(0,blk)
               write (termot,'(6x,i5)')elblks(elems,blk)
            end do
         end do
      end if
c
c                        now create the reordered element incidences.
c                        if the blocking is scalar and there is
c                        only one element type and one configuration
c                        number, then the new numbering will be the
c                        same as the old numbering.  Otherwise, the
c                        new and old orderings will be different.
c
c                        elreno, the vector that holds the reordering,
c                        is structured as:
c                             elreno(new number) = old element number
c
c
      nel = 0
      do j = 1, nelblk
         span = elblks(0,j)
         do l = 1, span
            elem  = elblks(l,j)
            nnode = elennd(elem)
            nel   = nel + 1
            elreno (nel) = elem
         end do
      end do
c
c                        if blocking is vectorized, ask if the
c                        user wants the old element to new element
c                        ordering printed.
c
      if ( scalar .and. (numgrp.eq.1) ) go to 1000

c
c                        sixth question asked (q_num = 6)
c                        'print the new=>old element listing...?'
c                        ynanswers(6) = 1, => yes
c                        ynanswers(6) = 0, => no
c
c
      q_num = 6
      if ( ynanswers(q_num) .eq. 0 ) go to 1000
c
      print_to_file = .false.
c
c
c
c                         the file to print to was obtained in the
c                         'questions' subroutine
c                         now strip file name and open file(if specified)
c
 300  call stripf(elem_print_file)
      if (elem_print_file(1:4) .ne. ' ') then
         print_to_file = .true.
         if (elem_print_file(1:2) .eq. 'dd')
     &        elem_print_file(1:18)='old_to_new_ele.txt'
       open(unit=print_file_num,file=elem_print_file,status='unknown',
     &               err=400)
      else
         print_file_num = termot
      endif
      goto 500
c
c                         if can't open file, print to screen
c
 400  continue
      write(termot,*)'>> Error in opening file. Writing to terminal.'
      print_to_file = .false.
c
c                         now print the old to new element listing
c
 500  continue
      lines = numele/9
      if (lines.gt.0) then
         do i= 1, numele/9
            write(print_file_num,9050)(i-1)*9+1, (i-1)*9 + 9
            write(print_file_num,9060)(9*(i-1)+j,j=1,9)
            write(print_file_num,9070)(elreno(9*(i-1)+j),j=1,9)
         enddo
      endif
      left_over = mod(numele,9)
      if (left_over.gt.0) then
         write(print_file_num,9050) lines*9+1,lines*9+left_over
         write(print_file_num,9060)(9*lines+j,j=1,left_over)
         write(print_file_num,9070)(elreno(9*lines+j),j=1,left_over)
      endif
      if (print_to_file) close(print_file_num)
c
c
 1000 continue
c
c
      if (debug) write (termot,*) '<<<< leaving trnler'
      return
 9000 format (' >>>> exceeded maximum number of blocks: ',/,
     &        '       fatal error.  ending processing.')
 9040 format (9i8)
 9050 format (' >> new element listing: ',i8,' to ',i8)
 9060 format ('new:',9i8)
 9070 format ('old:',9i8)
 9080 format (/,7x,' >> using standard blocking -- element',
     &          ' reordering skipped')
 9100 format (/,7x,' >> Notes:',/,
     &        7x,'       >> Multiple element configurations',
     &                      ' caused a renumbering process.',/,
     &        7x,'       >> Elements of the same configuration ',
     &                      ' are now numbered sequentially.',/,
     &        7x,'       >> Patran element results files',
     &                      ' should be displayed using newly',/,
     &        7x,'              provided Patran neutral file for model.',
     &       / )
 9200 format (/,7x,' >> Notes: ',/,
     &          7x '       >> Elements are renumbered. ',/)
 9400 format (/,8x,'>> begin element reordering')
 9500 format (/,8x,'>> use WARP3D automatic blocking assignment')
      end
c
c ************************************************************************
c *                                                                      *
c *                      routine  trnlmetis                              *
c *                                                                      *
c *     this subroutine calls the metis partitioning library to assign   *
c *     elements to processors.  The ordering results are then used in   *
c *     getgrp to make groups of elements corresponding to processors.   *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlmetis (numprocs)
      use patwarp_data
      implicit integer (a-z)
      integer, allocatable, dimension (:) :: etype, incvec, dadjncy,
     &                                       elemwgt, options, dxadj,
     &                                       node_reorder
c
c
      allocate ( options(0:4),
     &           etype(numele),
     &           dxadj(numele+1),
     &           elemwgt(numele),
     &           incvec(numele*8),
     &           dadjncy(numele*8),
     &           node_reorder(numnod) )
c
      do i = 0, 4
         options(i) = 0
      enddo
c
c                          build incidence structure for use of Metis
c                          decomposition library.  Metis only accepts linear
c                          elems, and requires that the nodes be numbered
c                          contiguously.  Thus we have to create a nodal
c                          renumbering to remove the holes in the corner
c                          numbering.  This does not effect the actual
c                          node numbering, it is just a trick to get metis
c                          to work correctly.
c
c                              set node_reorder to 1 for a node if it is
c                              a corner node, 0 otherwise
c
      node_reorder(1:numnod) = 0
c
      do elem = 1, numele
         incptr = eleipt(elem) - 1
         if ((eletyp(elem).eq.11).or.(eletyp(elem).eq.12)) then
            cnt = 4
            etype(elem) = 2
         else
            cnt = 8
            etype(elem) = 3
         end if
         do node = 1, cnt
            node_reorder(eleinc(incptr+node)) = 1
         end do
      end do
c
c                              renumber all nodes which have a 1 in
c                              node_reorder
c
      node = 0
      do i = 1, numnod
         if ( node_reorder(i) .eq. 1) then
            node = node + 1
            node_reorder(i) = node
         end if
      end do
      tot_node = node
c
c                              use reordering in incidence list
c                              for Metis
c
      ptr = 0
      do elem = 1, numele
         incptr = eleipt (elem) - 1
         if (etype(elem) .eq. 2) then
            cnt = 4
         else
            cnt = 8
         end if
         do node = 1, cnt
            incvec(ptr+node) = node_reorder( eleinc(incptr+node) )
      if (ptr+node.gt.numele*8)  write(*,*) ptr+node
         end do
         ptr = ptr + cnt
      end do
c
c                          We must weight the graphs so that each
c                          processor will take similar amounts of time.
c
c                          call METIS library to create graph from mesh
c
      write (termot,1200)
      write (termot,1000)  tot_node
      call METIS_MeshToDual ( numele, tot_node, incvec, etype, 1,
     &     dxadj, dadjncy )
      write(termot,1210)
c
c                          loop thru element types to assign weights to
c                          the elements on the graphs.  Also make all
c                          edges of weight 1.
c
      do elem = 1, numele
         if ( eletyp(elem) .eq. 1 ) then
            elemwgt(elem) = 10
         else if ( eletyp(elem) .eq. 2) then
            elemwgt(elem) = 23
         else if ( eletyp(elem) .eq. 3) then
            elemwgt(elem) = 35
         else if ( eletyp(elem) .eq. 4) then
            elemwgt(elem) = 63
         else if ( eletyp(elem) .eq. 5) then
            elemwgt(elem) = 13
         else if ( eletyp(elem) .eq. 11) then
            elemwgt(elem) = 10
         else if ( eletyp(elem) .eq. 12) then
            elemwgt(elem) = 63
         end if
      end do
c
c                          call METIS library to partition our weighted
c                          graph.
c
      write(termot,1220)
      if ( numprocs .ge. 8 ) then
         call METIS_PartGraphKway (numele, dxadj, dadjncy, elemwgt,
     &        dum, 2, 1, numprocs, options, dum2, eleproc)
      else
         call METIS_PartGraphRecursive (numele, dxadj, dadjncy,
     &        elemwgt, dum, 2, 1, numprocs, options, dum2, eleproc)
      endif
c
      deallocate( etype,incvec,dadjncy,elemwgt,options,node_reorder,
     &            dxadj )
      write(termot,1230)
c
      return
 1000 format(8x,'>> total number of corner nodes:',i8)
 1200 format(8x,'>> using metis to build graph of mesh' )
 1210 format(8x,'>> completed graph of mesh')
 1220 format(8x,'>> using metis to patition graph')
 1230 format(8x,'>> completed graph partitioning')
c
      end
c
c ************************************************************************
c *                                                                      *
c *                      routine trnlproc_access                         *
c *                                                                      *
c *     this subroutine identifies elements which are on the border      *
c *     between processors. The border elements are then placed in       *
c *     groups corresponding to the processors which they access (have   *
c *     as neighbors). This permits efficient partitioning for the ebe   *
c *     preconditioner in parallel.                                      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlproc_access
      use patwarp_data
      implicit integer (a-z)
      dimension invinc(0:maxcon,numnod), proc_access(maxprocs),
     &     access_grps(0:maxcon, maxaccgrps)
      logical new, match_access, grp_match, local_debug
      data local_debug / .false. /
c
c           first form the inverse incidence table -- given a node, what
c           elements are attatched to it
c
      do i = 0, maxcon
         do j = 1, numnod
            invinc(i,j) = 0
         enddo
      enddo
c
      do elem = 1, numele
         incptr = eleipt ( elem ) - 1
         do node_loop = 1, elennd(elem)
            node = eleinc(incptr+node_loop)
            invinc(0,node) = invinc(0,node) + 1
            invinc(invinc(0,node),node) = elem
         enddo
      enddo
c
c           now loop over all elements and determine what processors own
c           its neighbors
c
c              loop over elements
c
      num_access_grps = 0
      do elem = 1, numele
         if (local_debug) write (termot,*) '>>>>>> elem, proc:',
     &        elem, eleproc(elem)
         num_access = 1
         proc_access(1) = eleproc(elem)

c
c                 loop over nodes of element
c
         incptr = eleipt ( elem ) - 1
         do node_loop = 1, elennd(elem)
            node = eleinc(incptr+node_loop)
c
c                    loop over elements connected to node -- skip current
c                    element
c
            do inv_loop = 1, invinc(0,node)
               neighbor_elem = invinc(inv_loop,node)
               if (elem .eq. neighbor_elem) goto 100
c
c                       find proc which owns element.  Add to proc_access
c                       list (don't duplicate processors in proc_access)
c
               proc = eleproc(neighbor_elem)
               new = .true.
               do i=1, num_access
                  if ( proc .eq. proc_access(i)) new = .false.
               enddo
               if (new) then
                  num_access = num_access + 1
                  proc_access(num_access) = proc
               endif
c
 100           continue
            end do
         end do
         if (local_debug) write(termot,9000) elem, eleproc(elem),
     &        (proc_access(j),j=1, num_access)
c
c                 we now have list of processors which access the neighboring
c                 elements to the current element.  identify if this set
c                 of neighboring processors is the same as previous elements
c                 or if it is new.
c
c                    If all neighbors were owned by same processor which
c                    owns the element, then set grp = 0 (internal element).
c
         if ( num_access .eq. 1 ) then
            elem_access(elem) = 0
            goto 200
         endif
c
c                    loop over access groups -- see if current proc accesses
c                    match with defined group
c
         matching_grp = 0
         do i=1, num_access_grps
            if ( num_access .ne. access_grps(0,i)) goto 150
c
            grp_match = .true.
            do j = 1, num_access
c
               match_access = .false.
               do k = 1, num_access
                  if (proc_access(j) .eq. access_grps(k,i)) then
                     match_access = .true.
                  endif
               enddo
c
               if ( .not. match_access ) then
                  grp_match = .false.
               endif
c
            enddo
c
            if ( grp_match) then
               matching_grp = i
            endif
c
  150       continue
         end do
c
c                    if matching grp found, set elem_access to the grp
c                    number. otherwise, greate a new access grp.
c
         if ( matching_grp .gt. 0 ) then
c
            elem_access (elem) = matching_grp
c
         else
c
            num_access_grps = num_access_grps + 1
            elem_access (elem) = num_access_grps
c
            if ( num_access_grps .gt. maxaccgrps) then
               write (termot,9100)
               stop
            endif
c
            access_grps(0,num_access_grps) = num_access
            do i=1, num_access
               access_grps(i, num_access_grps) = proc_access(i)
            end do
c
         endif
c
 200     continue
         if (local_debug) write (termot,'(9x,"group:",i6)')
     &        elem_access(elem)
      enddo
c
      if (debug) then
         write (termot,*) '>>> access groups:'
         do i=1, num_access_grps
            write (termot,*) '  -> grp:',i,'  accesses procs:'
            write (termot,'(10x,20i3)') (access_grps(j,i),j=1,
     &           access_grps(0,i))
         end do
      endif
c
 9000 format ("  elem:",i6," owner proc:",i6," access procs:",20i6)
 9100 format ('>>> FATAL ERROR: exceeded number of access groups.',/,
     &        '>>>      recompile patwarp with maxaccgrps set larger.',
     &        /)
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine getgrp                       *
c     *                                                              *
c     *     this subroutine partitions elements into groupings of    *
c     *     similar elements.                                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine getgrp
      use patwarp_data
      implicit integer (a-z)
      dimension grprm(4,mxnmgp)
      logical newgrp
c
c
c                       the element properties defining similarity are:
c
c                                1) type of element
c                                2) configuration number
c
c
      if (debug) write (termot,*) '>>> in getgrp'
c
c
      numgrp = 0
      do i = 1, mxnmgp
       grpnum(i)  = 0
       grphead(i) = 0
       grptail(i) = 0
      end do
      grplst(1:numele) = 0
c
c                    check over all elements
c
      do elem = 1, numele
         if (debug) write (termot,'("  working on elem ",i6)') elem
c
c                           check if element fits in a group
c
         if ( elem .eq. 1 ) then
            newgrp = .true.
         else
            do grp = 1, numgrp
               if ( (grprm(1,grp).eq.eletyp(elem)) .and.
     &              (grprm(2,grp).eq.elecfg(elem)) .and.
     &              (grprm(3,grp).eq.eleproc(elem)) .and.
     &              (grprm(4,grp).eq.elem_access(elem))) then
                  if (debug) write (termot,*) ' elem matches'
                  grpnum(grp)= grpnum(grp)+1
                  grplst(grptail(grp)) = elem
                  grptail(grp) = elem
                  newgrp = .false.
                  go to 10
               else
                  if (debug) write (termot,*) ' elem doesnt match'
                  newgrp = .true.
               end if
            end do
  10     end if
c
c                          if needed, create a new group
c
         if ( newgrp ) then
            numgrp = numgrp+1
            if (debug) write (termot,'("  making group ",i4)') numgrp
            if ( numgrp .gt. mxnmgp ) then
               write (termot,9000) mxnmgp
               stop
            end if
            grpnum(numgrp)  = 1
            grphead(numgrp) = elem
            grptail(numgrp) = elem
            grprm(1,numgrp) = eletyp(elem)
            grprm(2,numgrp) = elecfg(elem)
            grprm(3,numgrp) = eleproc(elem)
            grprm(4,numgrp) = elem_access(elem)
         end if
c
      end do
c
      if ( debug ) then
         do i=1, numgrp
            write (termot,*) '>>> for group ',i
            write (termot,'("   num_elems:",i4," type:",i3," cfg:",i3,
     &           " proc:",i3," access grp:",i3)') grpnum(i),
     &           grprm(1,i), grprm(2,i), grprm(3,i), grprm(4,i)
            write (termot,*) '  head, tail of group:',grphead(i),
     &           grptail(i)
         enddo
         write (termot,*) '>>> here is the pointer structure:'
         do i=1, numele
            write (termot,*) '     elem: ',i,' ptr:', grplst(i)
         enddo
      end if
c
c
 9000 format ('>>>>>> ERROR: this version of patwarp only supports ',
     &        i6,' groups.',
     &      /,'              Recompile with a larger value for mxnmgp.')

      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine gtblrb                       *
c     *                                                              *
c     *     this subroutine sets a vector of block numbers for each  *
c     *     element in a grouping of similar elements. elements in a *
c     *     given block are non-conflicting by red-black ordering.   *
c     *                                                              *
c     ****************************************************************
c
      subroutine gtblrb(nel,nnode,maxbsz,grp,iblock)
      use patwarp_data
      implicit integer (a-z)
      dimension iblock(1),inode(numnod)
      logical restart
c
c
      if (debug) write (termot,*) '>>>> in gtblrb'
      do i= 1,numele
         iblock(i)= 0
      enddo
      restart = .true.
c
 10   continue
c
c                find the first element in the group not
c                previously assigned to a block.
c
      if ( restart ) then
         elem = grphead(grp)
         restart = .false.
      else
         elem = grplst(elem)
c
c                    if elem=0, then all elements in group have been
c                     assigend to a block; return.
c
         if ( elem .eq. 0 ) goto 9999
      endif
c
c                    if elem already assigned, skip to next element in group
c
      if (iblock(elem).ne.0) goto 10
c
c                    start a new block.
c
      nelblk= nelblk+1
      iblock(elem)= nelblk
      ibsize= 1
c
      do i = 1,numnod
         inode(i) = 0
      enddo
c
      do i = 1,nnode
         inode(eleinc(eleipt(elem)+i-1))= 1
      endd o
c
c                 loop through remaining elements and test un-
c                 assigned elements for inclusion in the
c                 current block.
c
 20   continue
      elem = grplst(elem)
c
c                    if elem=0, then no more elements to add to block.
c                    check for a new block
c
      if ( elem .eq. 0 ) then
         restart = .true.
         goto 10
      endif
c
c                    if previously assigned, skip to next elem in group
c
      if (iblock(elem).ne.0) goto 20
c
c                    the current element is not blocked. test
c                    for common nodes with any element previously
c                    assigned to the current block.
c
      isum = 0
      do i = 1,nnode
c
c                       If transition element are being processed then
c                       zero incidences must not be counted as shared node
c
         if (eleinc(eleipt(elem)+i-1) .ne. 0) then
            isum = isum+inode(eleinc(eleipt(elem)+i-1))
         end if
      end do
c
c                    if isum is not zero, then skip to next elem in group
c
      if (isum.ne.0) goto 20
c
c                    the current element does not conflict. add
c                    to block.
c
      iblock(elem)= nelblk
      ibsize= ibsize+1
      do i = 1,nnode
         inode(eleinc(eleipt(elem)+i-1))= 1
      end do
c
c                    check for overflow of the maximum block
c                    size.
c
      if (ibsize.eq.maxbsz) then
         restart = .true.
         go to 10
      endif
c
c                    keep processing elements for this block
c
      goto 20
c
c
 9999 continue
c
      if (debug) then
         write (termot,*) ' writing iblock:'
         do i = 1, numele
            write (termot,'("iblock(",i5,") = ",i4)')i,iblock(i)
         enddo
         write (termot,*) '<<<< leaving gtblrb'
      endif
c
      return
      end

c ************************************************************************
c *                                                                      *
c *   routine  trnlin -- initialize variable for the warp3d trans.       *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlin
      use patwarp_data
      implicit integer (a-z)
c
c            initialize warp3d element type lists.
c
      nxtels = 1
      do i = 1, maxetp
          etypes(i) = ' '
          econfg(i) = 0
          ettypt(i) = 0
          etty_tail_pt(i) = 0
      enddo
      do i = 1, numele
          etypls(i,1) = 0
          etypls(i,2) = 0
      end do
c
      return
      end

c ************************************************************************
c *                                                                      *
c *  routine trnlnd -- write nodal coordinate data in warp3d format      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlnd
      use patwarp_data
      implicit integer (a-z)
      double precision zero
      data zero /0.0/
c
c            write nodal coordinate data in warp3d format.
c            if all z coordinates are zero (i.e, a 2-d)
c            problem, write only x and y values.
c
c            if writting coordinates to separate file, change
c            device number
c
      if( sep_file_control )then
         write(ofile,9200)coord_file(1:cdf_len)
         ofile = cd_file
      end if
c
c
      do i = 1, numnod
        if ( coord(3,i) .ne. zero ) go to 300
      end do
      write(ofile,9000)
      do i = 1, numnod
        write(ofile,9010) i, coord(1,i), coord(2,i)
      end do
      write(ofile,9100)
      ofile = ofile_store
      if( sep_file_control ) close( unit=cd_file )
      return
c
 300  continue
      write(ofile,9000)
      do i = 1, numnod
        write(ofile,9010) i, coord(1,i), coord(2,i), coord(3,i)
      end do
      write(ofile,9100)
      ofile = ofile_store
      if( sep_file_control ) close( unit=cd_file )
      return
c
 9000 format('c',/,'coordinates',/,'*echo off')
 9010 format(i7,3(1x,e16.9))
 9100 format('*echo on',/,'c ')
 9200 format('c ',/,'*input from ','''',a,'''',/,'c ')
c
      end

c ************************************************************************
c *                                                                      *
c *  routine trnlec -- check elements and drive the building of          *
c *                     the element list (listed by type)                *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlec
      use patwarp_data
      implicit integer (a-z)
c
c            go through each element -- get type and config #. then
c            call trnlae to build element list.
c            "elenum" means new element number. elreno converts new
c            element number to the old element number.
c
      do  elenum = 1, numele
         etype = eletyp(elreno(elenum))
         confg = elecfg(elreno(elenum))
         nnode = elennd(elreno(elenum))
c
c            Elements
c
         select case ( etype )
         case ( 1 )
             call trnlae( 'l3disop', elenum, confg)
         case ( 2 )
             call trnlae( 'ts12isop', elenum, confg)
         case ( 3 )
             call trnlae( 'ts15isop', elenum, confg)
         case ( 4 )
             call trnlae( 'q3disop', elenum, confg)
         case ( 5 )
             call trnlae( 'ts9isop', elenum, confg)
         case ( 11 )
            call trnlae( 'tet4', elenum, confg)
         case ( 12 )
            call trnlae( 'tet10', elenum, confg)
         case ( 21 )
            call trnlae( 'wedge6', elenum, confg)
         case ( 22 )
            call trnlae( 'wedge15', elenum, confg)
         case default
            write(termot,9030) elenum
         end select
c
      end do
c
      return
c
 9030 format(1x,'>>>> invalid element type. element ',i5,' skipped. ' )
c
      end
c ************************************************************************
c *                                                                      *
c *  function valid_etype -- determine if patwarp element type           *
c *                          is currently supported                      *
c *                                                                      *
c ************************************************************************
c
c
      logical function  valid_etype ( etype, hex, tet, wedge )
      implicit integer (a-z)
      logical  hex, tet, wedge
c
      hex   =  etype .ge. 1  .and. etype .le. 5
      tet   =  etype .ge. 11 .and. etype .le. 12
      wedge =  etype .ge. 21 .and. etype .le. 22
      valid_etype = hex .or. tet .or. wedge
c
      return
      end
c ************************************************************************
c *                                                                      *
c *  routine trnled -- write element incidences in warp3d format         *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnled
      use patwarp_data
      implicit integer (a-z)
      dimension  nod(50), tnod(50),  wedge15(15)
      logical    dupl, hex, tet, wedge, valid_etype
      data wedge15 / 1,2,3,4,5,6,7,8,9,13,14,15,10,11,12 /
c
c            write element incidence data in warp3d format.
c            the list of element types is also built for later printing.
c            "elenum" means new element number. elreno converts new
c            element number to the old element number.
c
c            if writting incidences to separate file, change
c            device number
c
c            nod() is actually written to file.
c            tnod() used for swapping patran->warp3d node numbers
c                   and for sorting to check for duplicate nodes
c
      if( sep_file_control ) then
         write(ofile,9200) incid_block_file(1:inbl_len)
         ofile = inbl_file
      end if
c
      write(ofile,9000)
      do elenum = 1, numele
         etype  = eletyp(elreno(elenum))
         if ( .not. valid_etype ( etype, hex, tet, wedge ) ) then
              write(termot,9040) elenum
              cycle
         end if
         nnode  = elennd(elreno(elenum))
         incptr = eleipt(elreno(elenum)) - 1
         do i = 1, nnode
           tnod(i)  = eleinc(incptr+i)
           nod(i)   = tnod(i)
         end do
c
c               if we have created hex transition elements this time, then
c               make sure the incidence list is correct for the
c               transition elements
c
         if ( hex .and. made_transition ) then
            if ( etype .eq. 4 ) call order1215( nod, nnode, 2 )
            etran = eletran(elreno(elenum))
            if ( etran .eq. 1 ) then
               etrptr = eletrpoin(elreno(elenum))
               if ( etype .eq. 5) then
                  nnode = 9
               elseif ( etype .eq. 2 ) then
                  nnode = 12
               elseif ( etype .eq. 3 ) then
                  nnode = 15
               end if
               iptr = eletrpoin(elreno(elenum))
               do i = 1, nnode
                  nod(i)  = eletrinc(iptr + i)
                  tnod(i) = nod(i)
               end do
            end if
         else
c
c              transition elements were not created this time.  If the
c              elements is not an 8 node brick, call the re-ordering
c              algorithm so that we output the proper incidences for
c              the transition elements.
c
            if (  hex .and. etype .ne. 1 )
     &             call order1215( nod, nnode, 2 )
c
         end if
c
c            convert tet and wedge ordering from patran to
c            warp3d ordering
c
         if ( tet ) continue
c
         if ( wedge ) then
            if ( etype .eq. 22 ) then
              do i = 1, 15
                tnod(wedge15(i)) = nod(i)
              end do
              do i = 1, 15
                nod(i) = tnod(i)
              end do
            end if
         end if
c
c            check for duplicate node numbers (except 0) in the incidence
c            list.  sort the node list, then look for duplicates.
c
         call sort( nnode, tnod )
         dupl = .false.
         do i = 2, nnode
            if ( tnod(i) .ne. tnod(i - 1) ) cycle
            if ( tnod(i) .eq. 0 ) cycle
            dupl = .true.
         end do
         if ( dupl ) write(termot,9070) elenum
c
c            write the element incidences
c
         if ( nnode .eq. 8 ) then
            write(ofile,9020) elenum, (nod(i), i=1,nnode)
          else
            write(ofile,9010) elenum, (nod(i), i=1,nnode)
         end if
c
      end do
      write(ofile,9100)
c
      ofile = ofile_store
c
      return
c
 9000 format('c',/,'incidences ',/,'*echo off' )
 9010 format(9i8,',')
 9020 format(9i8)
 9030 format(1x,'>>>> invalid hex element type  must have 8, 12, 15,'
     &   /,     '     20, 9 nodes:element ',i7,' skipped. ' )
 9040 format(1x,'>>>> invalid element type. Element ',i7,' skipped.' )
 9070 format(1x,'>>>> warning.  element ',i7,' has duplicate nodes in',
     &   /,     '               its incidences.')
 9080 format(1x,'>>>> warning.  element ',i7,' type ',a,' has a',
     &   /,     '               zero node in its incidence list')
 9100 format('*echo on',/,'c ')
 9200 format('c ',/,'*input from ','''',a,'''',/,'c ')
c
      end
c ************************************************************************
c *                                                                      *
c *  routine sort integer list in increasing order                       *
c *                                                                      *
c ************************************************************************
c
      subroutine sort( n, a )
      implicit integer (a-z)
      dimension  a(n)
      inc = 1
1     inc = 3*inc+1
      if( inc .le. n ) goto 1
2     continue
        inc = inc/3
        do 11 i = inc+1, n
          v = a(i)
          j = i
3         if( a(j-inc) .gt. v ) then
            a(j) = a(j-inc)
            j = j-inc
            if( j .le. inc ) go to 4
          go to 3
          end if
4         a(j) = v
11      continue
      if( inc .gt. 1 ) go to 2
      return
      end
c
c ************************************************************************
c *                                                                      *
c *   routine  trnlae --  add an element to the list of element types    *
c *                       for the warp3d translator                      *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlae( ename, elenum, config )
      use patwarp_data
      implicit integer (a-z)
      character * (*) ename
      character * 16   typin, slttyp
c
c
c            add element number 'elenum' of element type 'ename'
c            to the list of element numbers for the same type (note
c            that elenum is the new element number). the type
c            is a string of up to 16 characters of the element name
c            as used in warp3d.  element numbers of the same type are
c            stored in a linked list in the vector etypls.  the list
c            heads are stored in vector ettypt.  the forward pointer is
c            stored in the right half of the word and the element number
c            in the left half.  the next available space in the linked
c            list is kept in nxtels.  zero is used as an end list mark.
c
c
c            search the list of element types already processed.  if this
c            type not found, set up the linked list in the next available
c            row.
c
      avslot = 0
      typin  = ename
      do slot = 1, maxetp
         slttyp = etypes(slot)
         sltcfg = econfg(slot)
         if ( (typin.eq.slttyp) .and. (config.eq.sltcfg) ) go to 20
         if ( slttyp .eq. ' ' .and. avslot .eq. 0 ) avslot = slot
      end do
c
c            no elements of this type yet defined.  start this type's
c            linked list in row avslot.
c
      etypes(avslot) = typin
      econfg(avslot) = config
      slot           = avslot
c
c            elements of this type are in the one-way linked list
c            with head in row ettypt(slot) and tail in
c            etty_tail_pt(slot) slot. add elenum at end of the
c            list.
c
 20   continue
      nowptr = ettypt(slot)
      if ( debug ) write(termot,9000) slot, s1, s2, nowptr
      if ( nowptr .gt. 0 ) go to 30
c
c            no elements in the list yet.  add this element at head and
c            make it the tail as well.
c            note: element is stored as new element number
c
      etypls(nxtels,2)   = elenum
      ettypt(slot)       = nxtels
      etty_tail_pt(slot) = nxtels
      nxtels             = nxtels + 1
      return
c
c            elements are in the list. get current tail and
c            insert new element at tail.
c
 30   continue
      nowptr = etty_tail_pt(slot)
c
c            end of list at row nowptr.  add new element to list in next
c            available slot.  update previous element forward pointer.
c            update next available location. update tail of list.
c
      if ( etypls(nowptr,1) .ne. 0 ) then
       write(termot,*) '>> Fatal Error in trnlae @ 1'
       write(termot,*) '   program terminated...'
       stop
      end if
      etypls(nowptr,1)   = nxtels
      etypls(nxtels,2)   = elenum
      etty_tail_pt(slot) = nxtels
      nxtels = nxtels + 1
      return
c
 9000 format(1x,'>> routine trnlae -- slot, s1   = ',i3,1x,a16,
     &     /,1x,'                     s2, nowptr = ',a16,i4)
      end

c ************************************************************************
c *                                                                      *
c *   routine  trnlpe --  print element types in warp3d format           *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlpe
      use patwarp_data
      implicit integer (a-z)
      logical    comprs
c
c
c            build warp3d input data lines specifying the element
c            types.  the list of element types is scanned and output
c            data lines built and printed.
c
c            lists of consecutively numbered elements are compressed
c            into a finite integer list format to avoid a list overflow
c            in grammar store.  the following decision table applies:
c
c |comprs|nowptr|ilast+1|ifirst|
c |      |.ne.0 |=elenum|=ilast|   action
c |------|------|-------|------|----------------------------------------
c |  t   |   t  |   t   |  -   | add elenum to current list.  go to 10.
c |------|------|-------|------|----------------------------------------
c |  t   |   t  |   f   |  t   | write current 1-entry list.  start a
c |      |      |       |      | new list.  go to 10.
c |------|------|-------|------|----------------------------------------
c |  t   |   t  |   f   |  f   | write current multi-entry list. start a
c |      |      |       |      | new list.  go to 10.
c |------|------|-------|------|----------------------------------------
c |  t   |   f  |   t   |  -   | add elenum to current list & write it
c |      |      |       |      | out. go to 1000.
c |------|------|-------|------|----------------------------------------
c |  t   |   f  |   f   |  t   | write out current 1-entry list. write
c |      |      |       |      | out elenum. go to 1000.
c |------|------|-------|------|----------------------------------------
c |  t   |   f  |   f   |  f   | write out current multi-entry list.
c |      |      |       |      | write out elenum. go to 1000.
c |------|------|-------|------|----------------------------------------
c |  f   |   t  |   -   |  -   | start a new list.  go to 10.
c |------|------|-------|------|----------------------------------------
c |  f   |   f  |   -   |  -   | write out elenum.  go to 1000.
c |------|------|-------|------|----------------------------------------
c
c
      write(ofile,9000)
      do 1000 nowtyp = 1, maxetp
        if (debug ) write(termot,'("nowtyp, type, config:",
     &       i3,2x,a16,2x,i3)')nowtyp, etypes(nowtyp), econfg(nowtyp)
        if (etypes(nowtyp) .eq. '  ') go to 1000
        write(ofile,9060) econfg(nowtyp)
        nowptr = ettypt(nowtyp)
c
c            traverse the linked list containing element numbers for this
c            element type.  build each line of output element numbers and
c            print it.  nowptr = 0 means end of the linked list. note that
c            "elenum" is the new element number.
c
      comprs = .false.
c
 10   continue
      elenum = etypls(nowptr,2)
      nowptr = etypls(nowptr,1)
c
      if ( comprs ) then
        if ( nowptr .ne. 0 ) then
          if ( (ilast+1) .eq. elenum ) then
            ilast = elenum
            go to 10
          else
            if ( ifirst .eq. ilast ) then
              write( ofile, 9030 ) ilast
            else
              write( ofile, 9040 ) ifirst, ilast
            end if
            ifirst = elenum
            ilast  = elenum
            go to 10
          end if
        else
          if ( (ilast+1) .eq. elenum ) then
            ilast = elenum
            write( ofile, 9020 ) ifirst, ilast
            go to 900
          else
            if ( ifirst .eq. ilast ) then
              write( ofile, 9030 ) ilast
            else
              write( ofile, 9040 ) ifirst, ilast
            end if
            write( ofile, 9010 ) elenum
            go to 900
          end if
        end if
      else
        if ( nowptr .ne. 0 ) then
          ifirst = elenum
          ilast  = elenum
          comprs = .true.
          go to 10
        else
          write( ofile, 9010 ) elenum
          go to 900
        end if
      end if
c
c      now add element type and material at end
c
 900  continue
      if (etypes(nowtyp) .eq. 'tet4') then
         write(ofile,9048) etypes(nowtyp)
      else if (etypes(nowtyp) .eq. 'tet10') then
         write(ofile,9049) etypes(nowtyp)
      else
         write(ofile,9050) etypes(nowtyp)
      end if

c
 1000 continue
c
      return
c
 9000 format('c',/,'elements')
 9010 format( 6x, i7 , $             )
 9020 format( 6x, i7, ' -' , i7, $   )
 9030 format( 6x, i7,            ',' )
 9040 format( 6x, i7, ' -' , i7, ',' )
 9048 format (' type ',a10,' linear material default,',/,
     &        23x,'order 1pt_rule gausspts short')
 9049 format (' type ',a10,' linear material default,',/,
     &        23x,'order 4pt_rule gausspts short')
 9050 format (' type ',a10,' linear material default,',/,
     &        23x,'order 2x2x2 bbar center_output short')
 9060 format ('c   for config number ',i3)
c
      end
c
c ************************************************************************
c *                                                                      *
c *   routine  trnlcn -- warp3d constraint generator                     *
c *                                                                      *
c *         note that warp3d only handles 3 constraints, because it      *
c *         does not handle moments. thus the variables are changed      *
c *         from 6 dof to 3 dof.                                         *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trnlcn
      use patwarp_data
      implicit integer (a-z)
      character*2  mpclab(3)
      character * 3 labels(3), label_print(3)
      character * 80 outline
      dimension  dspbuf(3)
      real       dspbuf, tol, zero, constr_print(3)
      double precision mlt
      logical    bitchk
      data   labels / 'u', 'v', 'w' /
      data   mpclab / ' u', ' v', ' w' /
      data   tol    / 1.0e-06 /
c
c
c            generate warp3d constraint input lines if there
c            are nodes with constraints.
c
c            if writting constraints to separate file, change
c            device number
c
      if( sep_file_control )then
         write(ofile,9060)const_file(1:ctf_len)
         ofile = ct_file
      end if
c
c
      if( debug ) write(termot,9040)
      if( nxtnls .eq. 1)then
         if( sep_file_control ) close( unit=ct_file )
         ofile = ofile_store
         return
      end if
      limit = nxtnls - 1
      write(ofile,9000)
c
c            loop to process each constrained node.  class a constraints
c            only at this time.  pull the node number, the dof bit mask.
c            blank out the labels buffer.  count the number of dof actually
c            constrained (i.e., bits on).  build local displacement
c            buffer.
c
      do slot = 1, limit
         node = ndcnls(slot,2)
         ptr  = ndcnls(slot,1)
         mask = ndcndt(ptr)
         if( debug ) write(termot,9030) node, ptr, mask
         ptr  = ptr + 1
         cdof = 0
         do dof = 1, 3
            dspbuf(dof) = 0.0
            if( bitchk( mask, dof ) ) then
               cdof = cdof + 1
               dspbuf(cdof) = rndcndt(ptr)
               ptr = ptr + 1
               label_print(cdof) = labels(dof)
               constr_print(cdof) = dspbuf(cdof)
            endif
         end do
         if( debug ) write(termot,9050) cdof, dspbuf
         outline(1:80) = ' '
         if( cdof .ne. 0 ) then
           write(outline,9020) node, (label_print(i),
     &        constr_print(i), i=1,cdof )
           call squash_zeroes( outline )
           write(ofile,fmt='(a80)') outline
         end if
      end do
c
c            if mpcs exist, loop to process
c            note that the constant is always 0.0
c
      if( nmpc .gt. 0 ) then
         zero = 0.0
         write(ofile,9080)
         do mpc = 1, nmpc
            nod = mpcnod(mpc,1)
            dof = mpcdof(mpc,1)
            mlt = mpcmlt(mpc,1)
            write(ofile,9081) '   ', nod, mlt, mpclab(dof),','
            do trm = 2, mpctrm(mpc)
               nod = mpcnod(mpc,trm)
               dof = mpcdof(mpc,trm)
               mlt = mpcmlt(mpc,trm)
               write(ofile,9081) ' + ', nod, mlt, mpclab(dof),','
            end do
            write(ofile,9082) ' = ', zero
         end do
      end if
c
      write(ofile,9070)
c
      ofile = ofile_store
      if( sep_file_control )close( unit=ct_file )
c
      return
c
 9000 format('c',/,'*echo off',/,'constraints')
 9020 format(i8,3(1x,a8,1x,e13.6))
 9030 format(1x,'>> node, ptr, mask = ',i8,i5,i10)
 9040 format(1x,'>> entered routine trnlcn')
 9050 format(1x,'         cdof = ',i5,
     &     /,1x,'         dspbuf = ',3e16.9,
     &     /,1x,'                  ',3e16.9 )
 9060 format('c ',/,'*input from ','''',a,'''',/,'c ')
 9070 format('c ',/,'*echo on')
 9080 format('multipoint')
 9081 format (8x,a,i6,es16.7,a2,a)
 9082 format (8x,a,f15.4)
      end
c ****************************************************************************
c *                                                                          *
c *     routine squash_zeroes                                                *
c *                                                                          *
c *         replace the E format zeros in line of constraints output line    *
c *         with 0.0 for readability                                         *
c *                                                                          *
c ****************************************************************************
c
c
      subroutine squash_zeroes( outline )
      character * 80 outline
      if( outline(19:31) .eq. ' 0.000000E+00' )
     &  outline(19:31) =  ' 0.0         '
      if( outline(42:54) .eq. ' 0.000000E+00' )
     &  outline(42:54) =  ' 0.0         '
      if( outline(65:77) .eq. ' 0.000000E+00' )
     &  outline(65:77) =  ' 0.0         '
      return
      end
c
c ****************************************************************************
c *                                                                          *
c *     routine trnlnl -- warp3d nodal load generator                        *
c *                                                                          *
c *         note that warp3d only handles 3 load directions, because it      *
c *         does not handle moments. thus the variables are changed          *
c *         from 6 dof to 3 dof.                                             *
c *                                                                          *
c ****************************************************************************
c
c
      subroutine trnlnl
      use patwarp_data
      implicit integer (a-z)
      character * 8 labbuf(3), labels(3), ldname
      dimension lodbuf(3)
      real      lodbuf, temp
      logical   bitchk, nl_header_done
      data labels / 'force_x', 'force_y', 'force_z'/
c
c
c            generate warp3d nodal load input lines.
c
      if ( debug ) write(termot,9040)
c     if ( nxtnvl .eq. 1 ) return
c
c            loop to process nodal loads and temperatures for
c            each load set.
c
      do nowset = 1, maxset
         lodset = nlodtb(nowset,1)
         if ( lodset .eq. 0 ) go to 1000
c
c            issue new loading condition statement if required.
c            build character variable for load name from the load set
c            number.
c
         ldname = ' '
         ldname(1:4) = 'set_'
         if ( lodset .le. 9 ) then
                    ldname(5:5) = '0'
                    write(ldname(6:6),fmt='(i1)') lodset
             else if ( lodset .le. 99 ) then
                    write(ldname(5:6),fmt='(i2)') lodset
             else if ( lodset .le. 999 ) then
                    write(ldname(5:7),fmt='(i3)') lodset
             else if ( lodset .le. 9999 ) then
                    write(ldname(5:8),fmt='(i4)') lodset
         end if
         write(ofile,9000) ldname
c
c            extract and print the nodal loads.
c
         nl_header_done = .false.
         ptr = nlodtb(nowset,2)
         if ( ptr .eq. 0 ) go to 500
         write(ofile,9005)
         nl_header_done = .true.
 100     continue
         node = nodval(ptr)
         fptr = nodptr(ptr)
         ptr  = ptr + 1
         mask = nodval(ptr)
         ptr  = ptr + 1
c
c            build local vector of loads and labels.
c
         ldof = 0
         do  dof = 1, 3
             labbuf(dof) = ' '
             lodbuf(dof) = 0.0
             if ( bitchk( mask, dof ) ) then
                   ldof = ldof + 1
                   lodbuf(ldof) = rnodval(ptr)
                   ptr = ptr + 1
                   labbuf(ldof) = labels(dof)
             end if
         enddo
         if ( debug ) write(termot,9050) ldof, labbuf, lodbuf
         write(ofile,9020) node, (labbuf(dof),lodbuf(dof),dof=1,ldof)
c
c            process next node if one in in the list.
c
         ptr = fptr
         if ( ptr .gt. 0 ) go to 100
c
c            loop through the same load set to extract and
c            print the nodal temperatures.
c
 500     continue
         ptr = nlodtb(nowset,4)
         if ( ptr .eq. 0 ) go to 1000
         if ( .not. nl_header_done ) write(ofile,9005)
  600    continue
         node = nodval(ptr)
         fptr = nodptr(ptr)
         ptr  = ptr + 1
         temp = rnodval(ptr)
         write(ofile,9007) node, temp
         ptr  = fptr
         if ( ptr .gt. 0 ) go to 600

c
c            all done with loads and temperatures this load set.
c
 1000    continue
         call trneload(nowset,lodset)
      end do
c
c
c            all done outputting nodal load data
c
      return
c
 9000 format('c',/,'loading ',a8 )
 9005 format('  nodal loads' )
 9007 format(i8,' temperature ',e12.5)
 9010 format(6x,i5,1x,6(a3,1x),' = ',e12.5 )
 9020 format(i8,3(1x,a8,1x,e13.6) )
 9030 format(1x,'>> node, ptr, mask, load set = ',i5,i5,2i8)
 9040 format(1x,'>> entered routine trnlnl')
 9050 format(1x,'         ldof, labbuf = ',i5,6a8,
     &     /,1x,'         lodbuf = ',3e16.9,
     &     /,1x,'                  ',3e16.9 )
c
      end
c
c ************************************************************************
c *                                                                      *
c *      routine trneload -- write distributed loads on elements         *
c *                                                                      *
c ************************************************************************
c
c
      subroutine trneload(loading_case,fload)
      use patwarp_data
      implicit integer (a-z)
      logical   bitchk, local_debug, hex, tet, wedge, valid_etype
      character * 8 ldname
      integer   pspc(17), pointer, face_map(6), elnewo(numele)
      integer   nod_face(20),nod(30),nod1(30),nod2(30),list_tr(20)
      real      rvals
      dimension ivals(5), rvals(5)
      equivalence ( rvals, ivals )
c
      data face_map / 6,   5,   3,   4,   1,   2 /
c
      local_debug = .false.
c
c            elreno(new) ---> old
c            elnewo(old) ---> new
c
             do ielem=1, numele
                jelem=elreno(ielem)
                elnewo(jelem)=ielem
             end do
c
c            loop to process element pressure for each load set.
c            patran supports only constant pressure applied to an
c            element face.
c
      lodtyp = 0
      do i = 1, nxteld-1
         lodset  = eloads(i,1)
         if ( lodset .eq. 0 ) go to 1000
         if ( lodset .eq. loading_case) then
            if ( lodset .ne. lodtyp ) then
               if ( fload .eq. 0 ) then
                  ldname = ' '
                  ldname(1:4) = 'set_'
                  if ( lodset .le. 9 ) then
                     ldname(5:5) = '0'
                     write(ldname(6:6),fmt='(i1)') lodset
                  else if ( lodset .le. 99 ) then
                     write(ldname(5:6),fmt='(i2)') lodset
                  else if ( lodset .le. 999 ) then
                     write(ldname(5:7),fmt='(i3)') lodset
                  else if ( lodset .le. 9999 ) then
                     write(ldname(5:8),fmt='(i4)') lodset
                  end if
                  write(ofile,9000) ldname
               end if
               write(ofile,9005)
            end if
         end if
c
c              set up info, set element type, check validity,
c              unpack the bit maps
c
         lodtyp  = lodset
         pointer = eloads(i, 2)
         element = eloads(i, 3)
         map     = deload(pointer, 1)
         nfe     = deload(pointer, 2)
         npv     = deload(pointer + 1, 1)
         etype   = eletyp(element)
         if ( .not. valid_etype ( etype, hex, tet, wedge ) ) then
              write(termot,9040) element
              cycle
         end if
         do ij = 1, 17
            pspc(ij) = 0
            if ( bitchk(map,ij) ) pspc(ij) = 1
         end do
         do j = 1, npv
            ivals(j) = deload(pointer+1+j,1)
         end do
         do j = 1, 6
            if ( pspc(j+3) .eq. 1 ) then
               nonz = j
            end if
         end do
         patran_face = nfe
c
c              process loads applied to faces of tetrahedral elements at
c              this point. this is quite simple. hex elements can
c              involve transition elements and are more complex.
c
         if ( tet ) then
            if( lodset .eq. loading_case) then
             call tet_loaded_face( tface, pspc(10) )
             write(ofile,9100) elnewo(element), tface, rvals(1)
            end if
            cycle
         end if
c
c            *** processing a hex type element ***
c
c            convert patran distributed forces to
c            warp pressures. even faces change sign,
c            odd faces keep sign.
c
         if ( mod(patran_face,2) .eq. 0 ) then
            rvals(1) = -rvals(1)
         end if
         if (local_debug)
     &       write(termot,*) 'Elem', element,' load= ',rvals(1)
c
c            corner nodes of element faces
c
         if (patran_face .eq. 1) then
            nod_face(1) = 1
            nod_face(2) = 2
            nod_face(3) = 5
            nod_face(4) = 6
         elseif (patran_face .eq. 2) then
            nod_face(1) = 3
            nod_face(2) = 4
            nod_face(3) = 7
            nod_face(4) = 8
         elseif (patran_face .eq. 3) then
            nod_face(1) = 2
            nod_face(2) = 3
            nod_face(3) = 6
            nod_face(4) = 7
         elseif (patran_face .eq. 4) then
            nod_face(1) = 1
            nod_face(2) = 4
            nod_face(3) = 5
            nod_face(4) = 8
         elseif (patran_face .eq. 5) then
            nod_face(1) = 5
            nod_face(2) = 6
            nod_face(3) = 7
            nod_face(4) = 8
         elseif (patran_face .eq. 6) then
            nod_face(1) = 1
            nod_face(2) = 2
            nod_face(3) = 3
            nod_face(4) = 4
         end if
c
c           incidences of element
c
         etype  = eletyp(element)
         nnode  = elennd(element)
         incptr = eleipt(element) - 1
         do enode = 1, nnode
            nod2(enode) = eleinc(incptr+enode)
         end do
         etran = eletran(element)
c
c           if element is a non-transition
c
         if( etran .eq. 0 ) go to 100
c
c           if element is a transition element we need to know what are the
c           quadratic faces.
c
         etrptr = eletrpoin(element)
         if (etype .eq. 5) then
             nnode = 9
         elseif (etype .eq. 2) then
             nnode = 12
         elseif (etype .eq. 3) then
             nnode = 15
         end if
         iptr = eletrpoin(element)
         do enode = 1, nnode
            nod(enode) = eletrinc(iptr+enode)
         end do
         do enode = 1, 4
            nod1(enode) = nod2(nod_face(enode))
         end do
         if ( local_debug) then
            write(termot,*) (nod(enode), enode = 1, nnode)
            write(termot,'(4I5)') (nod1(enode),enode = 1, 4)
         end if
         inodes = 0
         do j = 1, 4
            do k = 1, nnode
               if ( nod1(j). eq . nod(k)) then
                  inodes = inodes + 1
                  list_tr(inodes) = k
               end if
            end do
         end do
         if (local_debug) write(termot,*) (list_tr(j), j = 1, 4)
         do j = 1, 4
            do k = 1, 4
               if (list_tr(j) .lt. list_tr(k)) then
                  call swap(list_tr(k),list_tr(j))
               end if
            end do
         end do
         face=1000*list_tr(1)+100*list_tr(2)+10*list_tr(3)+list_tr(4)
         if (local_debug) write(termot,*) 'face= ',face
         if     (face .eq. 1234) then
               warp_face = 1
         elseif (face .eq. 5678) then
               warp_face = 2
         elseif (face .eq. 2367) then
               warp_face = 3
         elseif (face .eq. 1458) then
               warp_face = 4
         elseif (face .eq. 1256) then
               warp_face = 5
         elseif (face .eq. 3478) then
               warp_face = 6
         end if
         goto 200
  100    continue
c
c            Hex element mapping:
c            Patran faces  1   2   3   4   5   6   correspond to
c            Warp faces    6   5   3   4   1   2
c
         warp_face = face_map(patran_face)
  200    continue
c
         if (local_debug) write(termot,*) 'warp_face', warp_face
         if ( debug ) then
            write(termot,'(3(i2,x),6i1,x,8i1)') pspc(1),pspc(2),pspc(3),
     &                        (pspc(ij),ij=4,9),(pspc(ij),ij=10,17)
            write(termot,*) lodset,pointer,element,map,nfe,npv
         end if
c
c             write the element, face and pressure value
c
         if( lodset .eq. loading_case ) then
             write(ofile,9100) elnewo(element), warp_face, rvals(1)
         end if
c
c            all done outputting element load data
c
 1000 continue
      end do
c
 9000 format('c',/,'loading ',a8 )
 9005 format('  element loads' )
 9100 format(2x,i8,' face ',i1,' pressure ',f12.3)
 9040 format(1x,'>>>> invalid element type. Element ',i7,' skipped.' )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  routine tet_loaded_face                     *
c     *                                                              *
c     *     given local corner nodes on a tet face, return the       *
c     *     loaded face number                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine tet_loaded_face( tface, flags )
      implicit integer (a-z)
      dimension flags(*)
      logical n1, n2, n3, n4
c
c             flags(1:4) indicate which (local) corner nodes are
c             on the loaded face. sort this out and return the
c             face number:
c
c               tet face:      corner nodes:
c               ---------      -------------
c                  1             1, 2, 3
c                  2             1, 2, 4
c                  3             2, 3, 4
c                  4             1, 3, 4
c
c             Patran and Warp3d face numbers are identical
c             as are the node ordering for tet elements
c
      n1 = flags(1) .eq. 1
      n2 = flags(2) .eq. 1
      n3 = flags(3) .eq. 1
      n4 = flags(4) .eq. 1

      if ( n1 .and. n2 .and. n3 ) then
         tface = 1
         return
      end if
c
      if ( n1 .and. n2 .and. n4 ) then
         tface = 2
         return
      end if
c
      if ( n2 .and. n3 .and. n4 ) then
         tface = 3
         return
      end if
c
      if ( n1 .and. n3 .and. n4 ) then
         tface = 4
         return
      end if
c
      write(*,*) '>> FATAL ERROR: routine tet_loaded_face...'
      write(*,*) '                patwarp aborted'
c
      stop
      end
c     ****************************************************************
c     *                                                              *
c     *           routine trnlbc -- originally part of genout        *
c     *                                                              *
c     *     this subroutine writes out the blocking command          *
c     *                                                              *
c     ****************************************************************
c
      subroutine trnlbc
      use patwarp_data
      implicit integer (a-z)
c
c           writes out the block numbers, the number of elements in each
c              block, and the first element in each block.  If the file
c              being generated is for parallel computation, then add the
c              processor number which owns the block as the last number.
c
c
      if( block_method .eq. 1 ) then
         write(ofile,9210) blksz
         write(ofile,9200)
         return
      end if
c
      if( scalar ) then
         write(ofile,9000)
      else
         write(ofile,9010)
      endif
      write(ofile,9220)
c
      felem = 0
      span  = 1
c
      do blk = 1,nelblk
c
         oldspn = span
         span   = elblks(0,blk)
         felem  = felem+oldspn
         if( parallel ) then
            write(ofile,9100) blk, span, felem, blkproc(blk)
         else
            write(ofile,9100) blk, span, felem
         endif
c
      end do
c
      write(ofile,9200)
c
      return
 9000 format ('c',/,'c',/,'*echo off',/,'blocking',
     &    ' $ scalar - elems in blk share nodes')
 9010 format ('c',/,'c',/,'*echo off',/,'blocking',
     &    ' $ vectorized - no elems in blk share nodes')
 9100 format (4i8)
 9200 format('c ',/,'*echo on')
 9210 format('c',/,'blocking automatic size = ',i3 )
 9220 format('c     <blk>  <# elems> <1st elem> <domain>')
      end
c
c     ****************************************************************
c     *                                                              *
c     *         routine trnlnf -- originally subroutine nuneut       *
c     *                                                              *
c     *     this subroutine creates a patran neutral file containing *
c     *     the reordered elements.                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine trnlnf
      use patwarp_data
      implicit integer (a-z)
      integer dspmsk(6), pspc(6), lodmsk(6)
      logical bitchk
      real dspbuf(6), lodbuf(6)
c
      out = 16
c
      if (nptfle.eq.'none') goto 9999
 10   continue
c
c
c                          eight question asked (q_num = 8)
c                          'make an updated patran neutral file...'
c                          ynanswers(8) = 0, => no
c
c                          if creating new neutral file, process file
c                          name and open file
c

      q_num = 8
c
      if ( ynanswers(q_num) .eq. 0 )go to 9999
c
      call stripf(nptfle)
      if (nptfle(1:1).eq.' ') nptfle = 'newpat.out'
c
      open(out,file=nptfle,status='unknown',err=9998)
c
c             start writing packets
c
c               => output title -- packet 25
c
      write(out,9000) 25,0,0,1,0,0,0,0,0
      write(out,9010) ustitl
c
c               => output creation date & time and release number
c                                        -- packet 26
c
      write(out,9000) 26,0,0,1,numnod,numele,0,0,0
      write(out,9020) date,time,versn
c
c               => output nodal coordinates -- packet 01
c
      do node = 1, numnod
         write(out,9000) 1,node,0,2,0,0,0,0,0
         write(out,9030) coord(1,node), coord(2,node), coord(3,node)
         do i=1,6
            pspc(i) = 0
            if (bitchk(nodpcn(node),i)) pspc(i)= 1
         enddo
         write(out,9040) nodcon(node), nodtyp(node), nodndf(node),
     &        nodcfg(node), nodcid(node), (pspc(i), i=1,6)
      end do
c
c               => output element incidences, etc -- packet 02. for
c                  9, 12, 15 node hex transition elements, write as
c                  8 node elements.
c
      do elem = 1, numele
         oldele     = elreno(elem)
         etype      = eletyp(oldele)
         num_enodes = elennd(oldele)
c
c            Types of elements: in warp all elements are solids
c             iv = 5   tets (patwarp element types 11,12)
c             iv = 7   wedges (patwarp element types 21, 22)
c             iv = 8   hexes (patwarp element types 1-5)
c
         select case ( etype )
         case ( 1:5 )
            pat_etype = 8
            if ( etype .eq. 4 ) then
               num_enodes = 20
            else
               num_enodes = 8
            endif
         case ( 11,12 )
            pat_etype = 5
         case ( 21, 22 )
            pat_etype = 7
         case default
            write(termot,9200)  elem
            stop
         end select
c
         nlines = 1 + (num_enodes-1)/10 + 1
         write(out,9000) 2,elem,pat_etype,nlines,elenad(oldele),0,0,0,0
         write(out,'(4i8,3e16.9)') num_enodes,elecfg(oldele),
     &        elepid(oldele),elecid(oldele),(eleang(i,oldele),i=1,3)
         write(out,'(10i8)') (eleinc(eleipt(oldele)+i-1),
     &        i=1,num_enodes)
      end do
c
c               => output nodal loadings -- packet 07
c                       Note: this skips the temp loads;
c                         ptrs = nlodtb(load set, {4,5})
c
      do nowset = 1,maxset
         lodset = nlodtb(nowset,1)
         if ( lodset .eq. 0) goto 20
         ptr = nlodtb(nowset,2)
         if ( ptr .eq. 0 ) go to 20
 15      continue
         node = nodval(ptr)
         fptr = nodptr(ptr)
         ptr  = ptr + 1
         mask = nodval(ptr)
         ptr  = ptr + 1
c
c                      build local vector of loads and labels.
c                         Note: warp3d dooesn't handle
c                         moments, but we run all six dofs
c                         anyway so the mask is set correctly
c
         ldof = 0
         do dof = 1, 6
             lodmsk(dof) = 0
             lodbuf(dof) = 0.0
             if ( bitchk( mask, dof ) ) then
                   ldof = ldof + 1
                   lodbuf(ldof) = rnodval(ptr)
                   ptr = ptr + 1
                   lodmsk(dof) = 1
             end if
         end do
         write(out,9000) 7, node, nowset, 2, 0,0,0,0,0
         write(out,9070) 0,(lodmsk(dof), dof=1,6)
         write(out,9100) (lodbuf(dof), dof=1,ldof)
c
c                        process next node if one in in the list.
c
         ptr = fptr
         if ( ptr .gt. 0 ) go to 15
 20   enddo
c
c               => output nodal constraints  -- packet 08
c                         note that we run all 6 dofs despite the
c                         fact warp3d only handles 3; this is to
c                         set the constraint mask correctly
c
c
      if ( nxtnls .eq. 1 ) goto 1000
      limit = nxtnls - 1
      do slot = 1, limit
         node = ndcnls(slot,2)
         ptr  = ndcnls(slot,1)
         mask = ndcndt(ptr)
         ptr = ptr + 1
         cdof = 0
         do dof = 1,6
            dspmsk(dof) = 0
            dspbuf(dof) = 0.0
            if (bitchk(mask,dof)) then
               cdof = cdof+1
               dspbuf(cdof) = rndcndt(ptr)
               ptr = ptr+1
               dspmsk(dof) = 1
            endif
         enddo
         if ( debug ) write(termot,9030) node, ptr, mask
         write(out,9000) 8,node,1,2,0,0,0,0,0
         write(out,9070) 0,(dspmsk(x),x=1,6)
         write(out,9100) (dspbuf(x),x=1,cdof)
      enddo
 1000 continue
c
c               => end the neutral file -- packet 99
c
      write(out,9000) 99,0,0,1,0,0,0,0,0
c
c
      close(out,status='keep')
      write(termot,8010) nptfle(1:50)
c
      go to 9999
c
 9998 write(termot,2000)nptfle
      write(termot,2100)
c
c
 9999 return
c
 2000 format (' >>unable to open file named ', a80)
 2100 format (/)
 8010 format(14x,'>> updated neutral file created..',
     & /,14x,    '   name: ',a50)
 9000 format (I2,8I8)
 9010 format (a80)
 9020 format (a12,a8,a12)
 9030 format (3e16.9)
 9040 format (i1,1a1,i8,i8,i8,2x,6i1)
 9050 format (i8,i8,i8,i8,3e16.9)
 9060 format (10i8)
 9070 format (i8,6i1)
 9080 format (1e16.9)
 9100 format (5e16.9)
 9200 format(1x,'>> element id: ',i8,' has unsupported type..' )
      END
c     ****************************************************************
c     *                                                              *
c     *         routine trnlpo                                       *
c     *                                                              *
c     *     this subroutine creates a patran output file showing     *
c     *     the element processor assignment and the blocking        *
c     *     of the model.                                            *
c     *                                                              *
c     ****************************************************************
c
      subroutine trnlpo
      use patwarp_data
      implicit integer (a-z)
c
      out = 16
c
 10   continue
c
c
c                    ninth question asked (q_num = 9)
c                    'make a patran-readable results file...'
c                    ynanswers(9) = 0, => no, otherwise => yes
c
c                    if yes, process title, open, create...
c
      q_num = 9
      if( ynanswers(q_num) .eq. 0 ) go to 9999
c
      call stripf(resfile)
c        rwh -- modification
      if (resfile(1:1).eq.' ') resfile = 'metis_pat_info.els'
c
c                write elem processor assignment file first
c
      open(out,file=resfile,status='unknown',err=9998)
c
      write (out, '(80a1)') (resfile(i:i),i=1,80)
      write (out, '(i5)') 2
      write (out, '(80a1)') (resfile(i:i),i=1,80)
      write (out, '(80a1)') (resfile(i:i),i=1,80)
c
      felem = 1
      do blk = 1, nelblk
         span = elblks(0,blk)
         do elem = felem, felem + span - 1
            write (out, '(2i8,/,2e13.7)') elem, 8,
     &           float(eleproc(elreno(elem))), float(blk)
         enddo
         felem = felem + span
      enddo
c
      close(out,status='keep')
c
      go to 9999
c
 9998 write(termot,*) ' unable to open file named ', resfile
      write(termot,*)
      goto 9999
c
c
 9999 return
c
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                    routine order1215                         *
c     *                                                              *
c     *     This subroutine performs two tasks:                      *
c     *     task = 1  Determines the number of non zero nodes only   *
c     *     task = 2  reorders the incidences for a 9, 12, 15        *
c     *               and 20 node element to match the warp          *
c     *               convention.                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine order1215 ( n, nnode, task )
      implicit integer (a-z)
      integer n(*), aux(20), face(6,4), k(6), mids(12),nnode,task,
     & list20(20), list(15), list1(12), list2(12), list3(12),
     & list4(12), list5(12), list6(12), list13(15), list14(15),
     & list15(15), list16(15), list23(15), list24(15), list25(15),
     & list26(15), list35(15), list36(15), list45(15), list46(15)
c
      data nfaces, nsides / 6, 4 /
      data list20 /1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16/
      data list1  /1,2,3,4,5,6,7,8,9,10,11,12/
      data list2  /6,5,8,7,2,1,4,3,13,16,15,14/
      data list3  /2,6,7,3,1,5,8,4,18,14,19,10/
      data list4  /5,1,4,8,6,2,3,7,17,12,20,16/
      data list5  /6,2,1,5,7,3,4,8,18,9,17,13/
      data list6  /3,7,8,4,2,6,5,1,19,15,20,11/
      data list13 /2,3,4,1,6,7,8,5,10,11,12, 9,14,18,19/
      data list14 /4,1,2,3,8,5,6,7,12, 9,10,11,16,20,17/
      data list15 /1,2,3,4,5,6,7,8, 9,10,11,12,13,17,18/
      data list16 /3,4,1,2,7,8,5,6,11,12, 9,10,15,19,20/
      data list23 /7,6,5,8,3,2,1,4,14,13,16,15,10,19,18/
      data list24 /5,8,7,6,1,4,3,2,16,15,14,13,12,17,20/
      data list25 /6,5,8,7,2,1,4,3,13,16,15,14, 9,18,17/
      data list26 /8,7,6,5,4,3,2,1,15,14,13,16,11,20,19/
      data list35 /2,6,7,3,1,5,8,4,18,14,19,10,17, 9,13/
      data list36 /7,3,2,6,8,4,1,5,19,10,18,14,20,15,11/
      data list45 /5,1,4,8,6,2,3,7,17,12,20,16,18,13, 9/
      data list46 /4,8,5,1,3,7,6,2,20,16,17,12,19,11,15/
c
c          How many nodes does the element have? (non-zeros)
c
      nnode = 20
      do i = 1, 20
         if ( n(i) .eq. 0 ) nnode = nnode - 1
      end do
      if ( task .eq. 1) return
c
c          Convert patran numbering to warp numbering. all elements
c          have 20 nodes at this point or we shouldn't be here.
c
      do i = 1, 20
         aux(i) = n(i)
      end do
      do i = 1, 20
         n(i) = aux(list20(i))
      end do
c
      if ( nnode .eq. 20 ) return
c
c          new ordering for 9-node element. this is the simplest
c          one to process.
c
      if ( nnode .ne. 9 ) go to 14
      do i = 1, 12
          mids(i) = n(i+8)
      end do
      do i = 1, 12
         if ( mids(i) .ne. 0 ) then
            quadratic_edge = i
         end if
      end do
      do i = 1, 20
         aux(i) = n(i)
      end do
      go to (1,2,3,4,5,6,7,8,9,10,11,12), quadratic_edge
c
   1        do i = 1, 9
               list(i) = list1(i)
            end do
            go to 13
c
   2        do i = 1, 9
               list(i) = list13(i)
            end do
            go to 13
c
   3        do i = 1, 9
               list(i) = list16(i)
            end do
            go to 13
c
   4        do i = 1, 9
               list(i) = list14(i)
            end do
            go to 13
c
   5        do i = 1, 9
               list(i) = list2(i)
            end do
            go to 13
c
   6        do i = 1, 9
               list(i) = list23(i)
            end do
            go to 13
c
   7        do i = 1, 9
               list(i) = list26(i)
            end do
            go to 13
c
   8        do i = 1, 9
               list(i) = list24(i)
            end do
            go to 13
c
   9        do i = 1, 9
               list(i) = list45(i)
            end do
            go to 13
c
   10       do i = 1, 9
               list(i) = list35(i)
            end do
            go to 13
c
   11       do i = 1, 9
               list(i) = list36(i)
            end do
            go to 13
c
   12       do i = 1, 9
               list(i) = list46(i)
            end do
c
   13 continue
      do i = 1, 9
         n(i) = aux(list(i))
      end do
      do i = 10, 20
        n(i) = 0
      end do
      return
c
c          process 12, 15 node elements here. we must find the quadratic
c          faces then re-order the element nodes so that the
c          quadratic faces are those required in warp.
c
c
 14   continue
c
c
c         face 1
c
      face(1,1) = n(9)
      face(1,2) = n(10)
      face(1,3) = n(11)
      face(1,4) = n(12)
c
c         face 2
c
      face(2,1) = n(13)
      face(2,2) = n(14)
      face(2,3) = n(15)
      face(2,4) = n(16)
c
c         face 3
c
      face(3,1) = n(10)
      face(3,2) = n(19)
      face(3,3) = n(14)
      face(3,4) = n(18)
c
c         face 4
c
      face(4,1) = n(12)
      face(4,2) = n(20)
      face(4,3) = n(16)
      face(4,4) = n(17)
c
c         face 5
c
      face(5,1) = n(9)
      face(5,2) = n(18)
      face(5,3) = n(13)
      face(5,4) = n(17)
c
c         face 6
c
      face(6,1) = n(11)
      face(6,2) = n(19)
      face(6,3) = n(15)
      face(6,4) = n(20)
c
c            which are the full faces?, i.e. find out which
c            faces are quadratic for the transition elements.
c
      face1 = 0
      face2 = 0
      do i = 1, nfaces
         k(i) = 0
         do j = 1, 4
            if ( face(i,j) .ne. 0 ) k(i) = k(i) + 1
         end do
         if ( k(i) .eq. nsides ) then
            if ( face1 .eq. 0 ) then
               face1 = i
            else
               face2 = i
            end if
         end if
      end do
      do i = 1, 20
         aux(i) = n(i)
      end do
      if ( face2 .ne. 0 ) then
         if ( face1 .gt. face2 ) then
            facex = face2
            face2 = face1
            face1 = facex
         end if
      end if
c
c            new order for 12-node element
c
      if ( face2 .eq. 0 ) then
         go to (15,25,35,45,55,65), face1
c
   15       do i = 1, 12
               list(i) = list1(i)
            end do
            go to 75
c
   25       do i = 1, 12
               list(i) = list2(i)
            end do
            go to 75
c
   35       do i = 1, 12
               list(i) = list3(i)
            end do
            go to 75
c
   45       do i = 1, 12
               list(i) = list4(i)
            end do
            go to 75
c
   55       do i = 1, 12
               list(i) = list5(i)
            end do
            go to 75
c
   65       do i = 1, 12
               list(i) = list6(i)
            end do
c
   75    continue
         do i = 1, 12
            n(i) = aux(list(i))
         end do
         do i = 13, 20
            n(i) = 0
         end do
         return
      end if
c
c                new order for 15-node element
c
      go to (100,200,300,400), face1
c
  100    go to (130,140,150,160), face2 - 2
c
  130         do i = 1, 15
                 list(i) = list13(i)
              end do
              go to 170
c
  140         do i = 1, 15
                 list(i) = list14(i)
              end do
              go to 170
c
  150         do i = 1, 15
                 list(i) = list15(i)
              end do
              go to 170
c
  160         do i = 1, 15
                 list(i) = list16(i)
              end do
c
  170    continue
         go to 500
c
  200    go to (230,240,250,260), face2 - 2
c
  230          do i = 1, 15
                  list(i) = list23(i)
               end do
               go to 270
c
  240          do i = 1, 15
                  list(i) = list24(i)
               end do
               go to 270
c
  250          do i = 1, 15
                  list(i) = list25(i)
               end do
               go to 270
c
  260          do i = 1, 15
                  list(i) = list26(i)
               end do
  270    continue
         go to 500
c
  300    go to (350,360), face2 - 4
c
  350         do i = 1, 15
                 list(i) = list35(i)
              end do
              go to 370
c
  360         do i = 1, 15
                 list(i) = list36(i)
              end do
  370    continue
         go to 500
c
  400    go to (450,460), face2 - 4
c
  450         do i = 1, 15
                 list(i) = list45(i)
              end do
              go to 470
c
  460         do i = 1, 15
                 list(i) = list46(i)
              end do
  470    continue
c
  500 continue
      do i = 1, 15
        n(i) = aux(list(i))
      end do
      n(16) = 0
      n(17) = 0
      n(18) = 0
      n(19) = 0
      n(20) = 0
      return
      end

