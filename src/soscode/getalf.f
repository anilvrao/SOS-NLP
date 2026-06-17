C
C
C
      SUBROUTINE GETALF( NX, ALFNOI, SMLSTP, PSTEP, X, BL, BU,
     .                   IUPTYP, IADD, ALPHA                     )
C
C ======================================================================
C     GETALF===>getalf   J.T. BETTS
C ======================================================================
C
C     ROUTINE TO DETERMINE WHICH VARIABLE BOUND TO
C     ADD (IF ANY) AND ALPHA SUCH THAT X + ALPHA * PSTEP
C     HITS THIS BOUND. (PSTEP IS THE STEP FROM X
C     CORRESPONDING TO THE KT SYSTEM SOLUTION.)
C
C     UPON LEAVING GETALF:
C                       { 1 IF NO BOUND IS HIT, I.E. ALPHA=1.0
C            IUPTYP =  {-1 IF VARIABLE IADD HITS IT LOWER BOUND
C                       {-2 IF VARIABLE IADD HITS IT UPPER BOUND
C
      DOUBLE PRECISION ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
C
      INTEGER IUPTYP, I, NX, IADD
C
      DOUBLE PRECISION PSTEP(NX), X(NX) 
      DOUBLE PRECISION BU(NX), BL(NX)
      DOUBLE PRECISION ALF, ALPHA, ALFNOI, SMLSTP
      DOUBLE PRECISION PALFMX
C
C
      IADD = 0
      IUPTYP = 1
      ALPHA = ONE
      PALFMX = ZERO
C
      DO I = 1, NX
C       ENSURE THAT PSTEP(I) IS "REAL", I.E. NOT NONZERO
C       DUE SOLELY TO INACCURACY IN THE SOLUTION PROCESS.
        IF ( ABS( PSTEP(I) ) .GT. SMLSTP ) THEN
C
          IF ( PSTEP(I) .GT. ZERO ) THEN
            ALF = MAX( ZERO, ( BU(I) - X(I) ) / PSTEP(I) )
C
            IF ( ALF .LT. ALFNOI ) THEN
C
C             WILL HIT I-TH CONSTRAINT WITH A, VIRTUALLY, ZERO STEP.
C             AS AN ATTEMPT TO INHIBIT ILL-CONDITIONING, AMONG THOSE
C             BOUNDS HIT WITH A NEARLY ZERO STEP, ADD THE ONE
C             WHICH HAS THE LARGEST MAGNITUDE PSTEP COMPONENT.

              IF ( ABS(PSTEP(I)) .GT. PALFMX) THEN
                PALFMX = ABS(PSTEP(I))
                ALPHA = ALF
                IADD = I
                IUPTYP = -2
              ENDIF
C               
            ELSEIF ( ALF .LT. ALPHA ) THEN
              ALPHA = ALF
              IADD = I
              IUPTYP = -2
            ENDIF
C
          ELSEIF ( PSTEP(I) .LT. ZERO ) THEN
            ALF = MAX( ZERO, ( X(I) - BL(I) ) / (-PSTEP(I)) )
C
            IF ( ALF .LT. ALFNOI ) THEN
C
C             WILL HIT I-TH CONSTRAINT WITH A, VIRTUALLY, ZERO STEP.
C             AS AN ATTEMPT TO INHIBIT ILL-CONDITIONING, AMONG THOSE
C             BOUNDS HIT WITH A NEARLY ZERO STEP, ADD THE ONE
C             WHICH HAS THE LARGEST MAGNITUDE PSTEP COMPONENT.

              IF ( ABS(PSTEP(I)) .GT. PALFMX) THEN
                PALFMX = ABS(PSTEP(I))
                ALPHA = ALF
                IADD = I
                IUPTYP = -1
              ENDIF
C               
            ELSEIF ( ALF .LT. ALPHA ) THEN
              ALPHA = ALF
              IADD = I
              IUPTYP = -1
            ENDIF
C
          ENDIF
C
        ENDIF
      enddo
C
C
      RETURN
      END
