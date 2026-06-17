      
      SUBROUTINE MATPRN( A, NR, IMAX, JMAX, TITLE, LENTIT, IOUNIT )
C
C ======================================================================
C     MATPRN===>matprn   J.T. BETTS
C ======================================================================
C
C     ROUTINE PRINTS DOUBLE PRECISION RECTANGULAR MATRICES
C
      CHARACTER(LEN=60) FMT
      CHARACTER(LEN=*) TITLE
      INTEGER IMAX, JMAX, I, J, NR, IOUNIT, LENTIT
      DOUBLE PRECISION A(NR,*)
C
C
      FMT = '( /, 1X, A  ," IS BELOW- " )'
      WRITE( FMT(11:12) , '(I2.2)' ) LENTIT
C
      WRITE( IOUNIT, FMT ) TITLE
C
      DO I=1, IMAX
        WRITE( IOUNIT, 1001 ) ( A(I,J), I, J, J=1, JMAX )
        WRITE( IOUNIT, * )
      ENDDO
C
1001  FORMAT( 5( 1X, G14.6, '(', I2.2, ',', I2.2, ')', 4X ) )
C
      RETURN
      END
