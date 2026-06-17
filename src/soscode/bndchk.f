      SUBROUTINE BNDCHK( NX, X, BL, BU, BNDTOL, BNDVIO )
C
C ======================================================================
C     BNDCHK===>bndchk   J.T. BETTS
C ======================================================================
C
C     ROUTINE TO DETERMINE IF BOUNDS HAVE BECOME
C     EXCESSIVELY VIOLATED. 
C
      DOUBLE PRECISION ONE
      PARAMETER ( ONE = 1.0D0 )
C
      LOGICAL BNDVIO
C
      INTEGER I, NX
C
      DOUBLE PRECISION BL(NX), BU(NX)
      DOUBLE PRECISION BNDTOL, X(NX)
C
      BNDVIO = .FALSE.
C
      DO I = 1, NX
        IF ( BL(I) - X(I) .GT. BNDTOL ) THEN
          BNDVIO = .TRUE.
        ENDIF
C
        IF ( X(I) - BU(I) .GT. BNDTOL ) THEN
          BNDVIO = .TRUE.
        ENDIF
C
      enddo
C
C
      RETURN
      END
