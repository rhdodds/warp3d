c
c----67--1---------2---------3---------4---------5---------6---------7-2
c Convert all the characters in a string from lower case to upper case
c
      SUBROUTINE ToUpper( String )
      CHARACTER (LEN = *) String
      INTEGER I, Ismall, IBIG
      Ismall = ICHAR( 'a' )
      IBIG = ICHAR( 'A' )
c    
      DO I = 1, LEN( String )
        IF (String(I:I) >= 'a' .AND. String(I:I) <= 'z') THEN
         String(I:I) = CHAR( ICHAR( String(I:I) ) + IBIG - Ismall )
        END IF
      END DO
c
      END SUBROUTINE
