      PROGRAM MESH_PLOT
C     A SIMPLE POSTSCRIPT GENERATOR FOR HPPLOT TYPE FILES
C     ver 1.0     Date: 911112 (Modified by J.F. 960229)
C----67--1---------2---------3---------4---------5---------6---------7-!
C
C  INCLUDED SUBROUTINES & FUNCTIONS ARE
C
C  (1) CURVE1           (incl. in curve1.f  if makefile is used)
C  (2) NUM_FROM_STRING  (incl. in strings.f if makefile is used)
C  (3) CHAR_TO_INT      (incl. in strings.f if makefile is used)
C  (4) CHAR_TO_DP       (incl. in strings.f if makefile is used)
C  (5) LAST_NON_BLANK   (incl. in strings.f if makefile is used)
C  (6) POSX             (incl. in pos.f if makefile is used)
C  (7) POSY             (incl. in pos.f if makefile is used)
C  (8) ROTATE           (incl. in rotate.f if makefile is used)
C  (9) TEXT             (incl. in text.f if makefile is used)
C*************************************************************
      IMPLICIT NONE
      EXTERNAL LAST_NON_BLANK
      INTEGER LAST_NON_BLANK
      INTEGER I
      CHARACTER ROW*80, MESHPLOT_STUB*80, ANSW*3
      LOGICAL OK,L1,L2

C-------------------------------------------------------------
C NOTE!!! This part has to be changed at installation of the program
C
C DEFINE THE LOCAL DIRECTORY CONTAINING THE Mesh_plot.stub file here
C
       MESHPLOT_STUB='./Mesh_plot.stub'
C-------------------------------------------------------------

      WRITE(*,*) 'mesh-plot v1.0, 911112'
      OK=.FALSE.
      DO WHILE(.NOT.OK)
        WRITE(*,'(A,$)') 'Name of the plot-file without extension: '
        READ(*,'(A)') ROW
        I=LAST_NON_BLANK(ROW)
        INQUIRE(FILE=ROW(1:I)//'.inf',EXIST=L1)
        INQUIRE(FILE=ROW(1:I)//'.dat',EXIST=L2)
        OK=L1.AND.L2
      ENDDO
      INQUIRE(FILE=ROW(1:I)//'.ps',EXIST=L1)
      IF(L1)THEN
        WRITE(*,*) 'A FILE ',ROW(1:I)//'.ps',' already exists!'
        WRITE(*,'(T1,A,$)') ' Remove the old file (y/n): '
	READ(*,'(A)') ANSW
	IF ((ANSW.NE.'Y').AND.(ANSW.NE.'y')) STOP
      ENDIF
      OPEN(UNIT=11,STATUS='OLD',FILE=ROW(1:I)//'.inf')
      OPEN(UNIT=12,STATUS='OLD',FILE=ROW(1:I)//'.dat')
      OPEN(UNIT=21,STATUS='UNKNOWN',FILE=ROW(1:I)//'.ps')
      OPEN(UNIT=1,STATUS='OLD',FILE=MESHPLOT_STUB)
      DO WHILE (.TRUE.)      ! Move the .stub file into out-file
        READ(1,'(A)',END=10) ROW
        WRITE(21,'(T1,A)') ROW
      ENDDO
 10   CONTINUE
      CLOSE(1)

      WRITE(21,*) 'gsave'
      WRITE(21,*) ' 72 2.54 div dup scale   % cm scale'
      READ(11,'(A)') ROW
      DO WHILE (ROW(1:7).NE.'*FINISH')
        IF (ROW(1:1).EQ.'*') THEN
          IF(ROW(1:7).EQ.'*CURVE1') CALL CURVE1
          IF(ROW(1:7).EQ.'*ROTATE') CALL ROTATE
          IF(ROW(1:5).EQ.'*TEXT') CALL TEXT
        ENDIF
        READ(11,'(A)',END=999) ROW
      ENDDO
      WRITE(21,*) 'showpage'
      WRITE(21,*) 'grestoreall'
      STOP
 1001 FORMAT (G12.5,X,A)
 1002 FORMAT (G12.5,X,G12.5,X,A)
 1011 FORMAT (I6,X,A)
 1012 FORMAT (I6,X,I6,X,A)
  999 WRITE(*,*) 'END OF .inf FILE'
      END
C****************************************************************
C  S  U  B  R  O  U  T  I  N  E  S
C****************************************************************
      SUBROUTINE CURVE1
C----67--1---------2---------3---------4---------5---------6---------7-!
      IMPLICIT NONE
      DOUBLE PRECISION POSX,POSY
      DOUBLE PRECISION X0AT,Y0AT,X,Y,XMIN,XMAX,YMIN,YMAX,
     &                 X_LENGTH,Y_LENGTH
      DOUBLE PRECISION XLSTART,XLSP,DXTICK,YLSTART,YLSP,DYTICK,LSIZE
      INTEGER I,J,NCUR,NLIN,NUM
      DOUBLE PRECISION DPVEK(60)
      CHARACTER*80 ROW
      COMMON/SCALE/ XMIN,XMAX,YMIN,YMAX,X_LENGTH,Y_LENGTH

      READ(11,'(A)') ROW
      CALL CHAR_TO_DP(ROW,DPVEK,NUM)
      X0AT=DPVEK(1)
      Y0AT=DPVEK(2)
      XMIN=DPVEK(3)
      XMAX=DPVEK(4)
      YMIN=DPVEK(5)
      YMAX=DPVEK(6)
      X_LENGTH=DPVEK(7)
      Y_LENGTH=DPVEK(8)
      READ(11,'(A)') ROW
      CALL CHAR_TO_DP(ROW,DPVEK,NUM)
      XLSTART=DPVEK(1)
      XLSP=DPVEK(2)
      DXTICK=DPVEK(3)
      YLSTART=DPVEK(4)
      YLSP=DPVEK(5)
      DYTICK=DPVEK(6)
      LSIZE=DPVEK(7)

      WRITE(21,*) '% *CURVE1 ROUTINE BEGIN'
      WRITE(21,*) 'gsave'
      WRITE(21,1002) X0AT,Y0AT,'translate'
      WRITE(21,*) ' 0.02 setlinewidth'
      WRITE(21,*) ' 2 setlinecap  0 setlinejoin'
      WRITE(21,1002) 0.0, 0.0, 'moveto'
      WRITE(21,1002) POSX(XMAX),0.0,'lineto'
      WRITE(21,1002) POSX(XMAX),POSY(YMAX),'lineto'
      WRITE(21,1002) 0.0, POSY(YMAX), 'lineto closepath stroke'
      WRITE(21,*) '/Helvetica findfont ',LSIZE,' scalefont setfont'
      WRITE(21,*) '/fontsize',LSIZE,
     &            ' def    %Fontsize needed by subroutines'
      WRITE(21,*) ' 0.02 setlinewidth'
      WRITE(21,*) ' 2 setlinecap  0 setlinejoin'
      DO X = XLSTART,XMAX,XLSP
        WRITE(21,1002) POSX(X),POSY(YMIN)-0.25,'moveto'
        WRITE(21,*)    '0.0 0.3 rlineto'
        WRITE(21,1002) POSX(X),POSY(YMIN)-0.55
        WRITE(21,'(X,A,G12.5,A)')    '(',X,') centertop' 
c        WRITE(21,'(X,A,F5.1,A)')    '(',X,') centertop' 
      ENDDO
      DO Y = YLSTART,YMAX,YLSP
        WRITE(21,1002) POSX(XMIN)-0.25,POSY(Y),'moveto'
        WRITE(21,*)    '0.3 0.0 rlineto'
        WRITE(21,1002) POSX(XMIN)-0.55,POSY(Y)
        WRITE(21,'(X,A,G12.5,A)')    '(',Y,') rightcenter' 
c        WRITE(21,'(X,A,F5.1,A)')    '(',Y,') rightcenter' 
      ENDDO
      DO X = XLSTART,XMAX,XLSP/(1+DXTICK)
        WRITE(21,1002) POSX(X),POSY(YMIN)-0.075,'moveto'
        WRITE(21,*)    '0.0 0.15 rlineto'
      ENDDO
      DO Y = YLSTART,YMAX,YLSP/(1+DYTICK)
        WRITE(21,1002) POSX(XMIN)-0.075,POSY(Y),'moveto'
        WRITE(21,*)    '0.15 0.0 rlineto'
      ENDDO
      WRITE(21,*) 'stroke'

      READ(11,*) NCUR
      WRITE(21,*) ' 2 setlinecap  0 setlinejoin'
      WRITE(21,*) ' 0.01 setlinewidth gsave newpath',
     &            '   % prepare to set a clippath in the current graph'
      WRITE(21,1002) 0.0, 0.0, 'moveto'
      WRITE(21,1002) POSX(XMAX),0.0,'lineto'
      WRITE(21,1002) POSX(XMAX),POSY(YMAX),'lineto'
      WRITE(21,1002) 0.0,POSY(YMAX),'lineto closepath clip'
      DO I=1,NCUR
        READ(11,*) NLIN
        WRITE(21,*) '%  Processing curve #',i,' ,',NLIN-1,' lines.'
        READ(12,*) X,Y
        WRITE(21,1002) POSX(X),POSY(Y),'moveto'
        DO J=2,NLIN
          READ(12,*) X,Y
          WRITE(21,1002) POSX(X),POSY(Y),'lineto'
        ENDDO
        WRITE(21,*) 'stroke'
      ENDDO
      WRITE(21,*) 'grestore    % restore old clippath'
      WRITE(21,*) 'grestore'
      WRITE(21,*) '% *CURVE1 ROUTINE END'

      RETURN
 1001 FORMAT (G12.5,X,A)
 1002 FORMAT (G12.5,X,G12.5,X,A)
 1011 FORMAT (I6,X,A)
 1012 FORMAT (I6,X,I6,X,A)
      END
C****************************************************************
      DOUBLE PRECISION FUNCTION POSX(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,XMIN,XMAX,YMIN,YMAX,X_LENGTH,Y_LENGTH
      COMMON/SCALE/ XMIN,XMAX,YMIN,YMAX,X_LENGTH,Y_LENGTH
      POSX=(X-XMIN)/(XMAX-XMIN)*X_LENGTH
      RETURN
      END

      DOUBLE PRECISION FUNCTION POSY(Y)
      IMPLICIT NONE
      DOUBLE PRECISION Y,XMIN,XMAX,YMIN,YMAX,X_LENGTH,Y_LENGTH
      COMMON/SCALE/ XMIN,XMAX,YMIN,YMAX,X_LENGTH,Y_LENGTH
      POSY=(Y-YMIN)/(YMAX-YMIN)*Y_LENGTH
      RETURN
      END
C****************************************************************
       INTEGER FUNCTION NUM_FROM_STRING(ROW,I1,I2)
C--------------------------------------------------------------C
C        Extract the number from ROW beginning at I1 and       C
C        return index of character after number in I2.         C
C        NOTE ;  Sign is not handled !!
C--------------------------------------------------------------C
      IMPLICIT NONE
      INTEGER I1,I2,I,ICH
      CHARACTER ROW(80)
      
      I=0
      IF(I1.GT.80) THEN
         I2=0
      ELSE
         I2=I1
         ICH=ICHAR(ROW(I2))-48
         DO WHILE((I2.LE.80).AND.(ICH.GE.0).AND.(ICH.LE.9))
            I=I*10+ICH
            I2=I2+1
            ICH=ICHAR(ROW(I2))-48
         ENDDO 
      ENDIF 
      NUM_FROM_STRING=I
      RETURN
      END
      

      
C--------------------------------------------------------------C
      SUBROUTINE CHAR_TO_INT(ROW,IVEK,NUM)
C--------------------------------------------------------------C
C     --- The routine transform a character to an integer number !
      IMPLICIT NONE
      EXTERNAL IS_A_DIGIT
      LOGICAL IS_A_DIGIT
      INTEGER    IVEK(60),NUM,N,I,SIGN
      CHARACTER ROW(80)
      NUM=0
      N=1
      DO WHILE (N.LT.80)
         DO WHILE ((.NOT.IS_A_DIGIT(ROW(N))).AND.(N.LT.80) )
            N=N+1
         ENDDO 
         SIGN=1
         IF (N.GT.1) THEN
            IF (ROW(N-1).EQ.'-') THEN
               SIGN=-1
            ENDIF
         ENDIF
         IF (IS_A_DIGIT(ROW(N))) THEN
            NUM=NUM+1                 ! A DIGIT IS DETECTED
            I=SIGN*(ICHAR(ROW(N))-48)
            IF (N.LT.80) N=N+1
         ENDIF
C---  There is a number to convert
         DO WHILE (IS_A_DIGIT(ROW(N)).AND.(N.LT.80))
            I=10*I+SIGN*(ICHAR(ROW(N))-48)
            N=N+1
         ENDDO
         IF (IS_A_DIGIT(ROW(N))) THEN
            I=10*I+SIGN*(ICHAR(ROW(N))-48)
            N=N+1
         ENDIF
         IVEK(NUM)=I
      ENDDO
      END
C--------------------------------------------------------------C
      LOGICAL FUNCTION IS_A_DIGIT(CHR)
      CHARACTER CHR
      IS_A_DIGIT=((CHR.GE.'0').AND.(CHR.LE.'9'))
      RETURN
      END
C--------------------------------------------------------------C

C--------------------------------------------------------------C
      SUBROUTINE CHAR_TO_DP(ROW,DPVEK,NUM)
C--------------------------------------------------------------C
C     --- The routine transform a character to a DOUBLE PRECISION number !
      IMPLICIT  NONE                            ! Mod. 910721 Nisse J.
      EXTERNAL IS_A_DIGIT
      LOGICAL IS_A_DIGIT
      INTEGER    NUM,N,I,NDEC,SIGN
      DOUBLE PRECISION  MANTISSA,DPVEK(60)
      CHARACTER ROW(80)
      NUM=0
      N=1
      DO WHILE(N.LT.80)
         DO WHILE ( (.NOT.IS_A_DIGIT(ROW(N))).AND.(N.LT.80) )
            N=N+1    ! LOOK FOR A NUMBER
         ENDDO 
         IF (N.GE.80) GOTO 999
         NUM=NUM+1
         IF (N.GT.1) THEN
            IF (ROW(N-1).EQ.'.') THEN
               N=N-1
            ENDIF
         ENDIF 
         SIGN=1
         IF (N.GT.1) THEN
            IF (ROW(N-1).EQ.'-') THEN
               SIGN=-1
            ENDIF
         ENDIF
         IF (IS_A_DIGIT(ROW(N))) THEN
            I=SIGN*(ICHAR(ROW(N))-48)
            N=N+1
         ELSE
            I=0    ! NO MANTISSA ONLY A DECIMAL PERIOD
         ENDIF
C---  There is a mantissa number to convert
         DO WHILE (IS_A_DIGIT(ROW(N)).AND.(N.LT.80))
            I=10*I+SIGN*(ICHAR(ROW(N))-48)
            N=N+1
         ENDDO
         IF (IS_A_DIGIT(ROW(N))) THEN
            I=10*I+SIGN*(ICHAR(ROW(N))-48)
            N=N+1
         ENDIF
         MANTISSA=FLOAT(I)

C---  NOW LOOK FOR DECIMAL PART
         I=0
         IF (ROW(N).EQ.'.') THEN
            N=N+1
            NDEC=0
            DO WHILE (IS_A_DIGIT(ROW(N)).AND.(N.LT.80))
               I=10*I+SIGN*(ICHAR(ROW(N))-48)
               NDEC=NDEC+1
               N=N+1
            ENDDO
            IF (IS_A_DIGIT(ROW(N))) THEN
               I=10*I+SIGN*(ICHAR(ROW(N))-48)
               NDEC=NDEC+1
               N=N+1
            ENDIF
            MANTISSA=MANTISSA+DFLOAT(I)/(10**NDEC)
         ENDIF
         
C---  END OF MANTISSA

         I=0
         IF ((ROW(N).EQ.'E').OR.(ROW(N).EQ.'e').OR.
     $        (ROW(N).EQ.'D').OR.(ROW(N).EQ.'d')) THEN
            N=N+1
            IF (ROW(N).EQ.'-') THEN
               SIGN=-1
               N=N+1
            ELSEIF (ROW(N).EQ.'+') THEN
               SIGN=1
               N=N+1
            ELSE
               SIGN=1
            ENDIF
            DO WHILE (IS_A_DIGIT(ROW(N)).AND.(N.LT.80))
               I=10*I+SIGN*(ICHAR(ROW(N))-48)
               N=N+1
            ENDDO
            IF (IS_A_DIGIT(ROW(N))) THEN
               I=10*I+SIGN*(ICHAR(ROW(N))-48)
               N=N+1
            ENDIF
         ENDIF
         DPVEK(NUM)=MANTISSA*10.0**I
      ENDDO
 999  CONTINUE
      RETURN
      END
C	*******************************************************************
	integer function LAST_NON_BLANK(string)
	implicit none
 	character string*(*)
	integer   i,i2
C	*******************************************************************
C	***** Returns the position of the last non-blank character in *****
C       ***** the string. Returns zero if string is all blanks.       *****
C	*******************************************************************
	i=len(string)
	i2=0
	do while (i.gt.i2)
	  if(string(i:i).ne.' ') i2=i
	  i=i-1
	enddo
	last_non_blank=i2
	return
	end
C****************************************************************
      SUBROUTINE ROTATE
      IMPLICIT NONE
      WRITE(21,*) ' 90 rotate 0 -21 translate    % rotate paper'
      RETURN
      END
C****************************************************************
      SUBROUTINE TEXT
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,LSIZE,ANGLE
      INTEGER TYPE,LAST_NON_BLANK
      CHARACTER*80 ROW

      WRITE(21,*) '% *TEXT ROUTINE BEGIN'
      WRITE(21,*) 'gsave'
      READ(11,*) X,Y,LSIZE,ANGLE,TYPE
      WRITE(21,*) X,Y,' translate'
      WRITE(21,*) ANGLE,' rotate'
      WRITE(21,*) '/Helvetica findfont ',LSIZE,' scalefont setfont'
      WRITE(21,*) '/fontsize',LSIZE,
     &            ' def    %Fontsize needed by subroutines'
      READ(11,'(A)') ROW
      WRITE(21,*) 0.0,0.0
      WRITE(21,1001) ROW(1:LAST_NON_BLANK(ROW))
      IF(TYPE.EQ.1) THEN 
        WRITE(21,*) ' leftbottom'
      ELSEIF(TYPE.EQ.2) THEN 
        WRITE(21,*) ' leftcenter'
      ELSEIF(TYPE.EQ.3) THEN 
        WRITE(21,*) ' lefttop'
      ELSEIF(TYPE.EQ.4) THEN 
        WRITE(21,*) ' centerbottom'
      ELSEIF(TYPE.EQ.5) THEN 
        WRITE(21,*) ' centercenter'
      ELSEIF(TYPE.EQ.6) THEN 
        WRITE(21,*) ' centertop'
      ELSEIF(TYPE.EQ.7) THEN 
        WRITE(21,*) ' rightbottom'
      ELSEIF(TYPE.EQ.8) THEN 
        WRITE(21,*) ' rightcenter'
      ELSEIF(TYPE.EQ.9) THEN 
        WRITE(21,*) ' righttop'
      ELSE
        WRITE(*,*) ' Error in function text; Adjustment type ',
     &             TYPE,' is not implemented, using 1 instead.'
        WRITE(21,*) ' leftbottom'
      ENDIF
      WRITE(21,*) -ANGLE,' rotate'
      WRITE(21,*) -X,-Y,' translate'
      WRITE(21,*) 'grestore'
      WRITE(21,*) '% *TEXT ROUTINE END'
      RETURN
 1001 FORMAT(T2,'(',A,')')
      END
C****************************************************************
