
      SUBROUTINE SAMSET(NDIM,ISTATC,ISTATV,IOLDC,IOLDV,NSAME)
C
C ======================================================================
C     SAMSET===>SAMSET   J.T. BETTS 
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C         PURPOSE:  COMPARE THE ACTIVE SET FROM ITERATION TO ITERATION
C                   AND COMPUTE THE NUMBER OF STEPS WITH THE SAME
C                   ACTIVE SET
C
C             NDIM   NUMBER OF VARIABLES
C             ISTATC INTEGER CONSTRAINT STATUS (NCON)
C                    = 0  --- FREE (INACTIVE) INEQUALITY
C                    = 1  --- FIXED ON LOWER BOUND
C                    = 2  --- FIXED ON UPPER BOUND
C                    = 3  --- EQUALITY
C                    = 4  --- IGNORED CONSTRAINT
C             ISTATV INTEGER VARIABLE STATUS (NDIM)
C                    = 0  --- FREE VARIABLE 
C                    = 1  --- FIXED ON LOWER BOUND
C                    = 2  --- FIXED ON UPPER BOUND
C                    = 3  --- FIXED PERMANENTLY 
C                    = 4  --- IGNORED BOUND
C             IOLDC  OLD VALUE OF ISTATC
C             IOLDV  OLD VALUE OF ISTATV
C             NSAME  NUMBER OF STEPS WITH THE SAME ACTIVE SET
C                    MUST BE INITIALIZED TO ZERO.
C
      DIMENSION ISTATC(*),IOLDC(*)
      DIMENSION ISTATV(NDIM),IOLDV(NDIM)
C
      IDIFF = 0
      DO I = 1,MCON
        IDIFF = IDIFF + ISTATC(I) - IOLDC(I)
        IOLDC(I) = ISTATC(I)
      ENDDO
C
      DO I = 1,NDIM
        IDIFF = IDIFF + ISTATV(I) - IOLDV(I)
        IOLDV(I) = ISTATV(I)
      ENDDO
C
      IF(IDIFF.EQ.0) THEN
        NSAME = NSAME + 1
      ELSE
        NSAME = 0
      ENDIF
C
      RETURN 
      END
