
      SUBROUTINE VPRIN ( N, V, TITLE, LENTIT, IOUNIT )
C
C ======================================================================
C     VPRIN ===>vprin    J.T. BETTS
C ======================================================================
C
C     ROUTINE TO PRINT ONE DOUBLE PRECISION VECTOR
C
      INTEGER IOUNIT, N, I, LENTIT
C
      CHARACTER(LEN=60) FMT
      CHARACTER(LEN=*) TITLE
C
      DOUBLE PRECISION V(*)
C
C
      FMT = '( /, 1X, A  ," IS BELOW- " )'
      WRITE( FMT(11:12) , '(I2.2)' ) LENTIT
C
      WRITE( IOUNIT, FMT ) TITLE
C
      IF ( N .LE. 99 ) THEN
        WRITE(IOUNIT,905) (  V(I), I, I=1, N )
905     FORMAT( 4( G14.6, '(', I2.2, ')' ) )
      ELSEIF( N .LE. 999 ) THEN
        WRITE(IOUNIT,915) ( V(I), I, I=1, N )
915     FORMAT( 4( G14.6, '(', I3.3, ')' ) )
      ELSE
        WRITE(IOUNIT,925) ( V(I), I, I=1, N )
925     FORMAT( 4( G14.6, '(', I4.4, ')' ) )
      ENDIF
C
      RETURN
      END
