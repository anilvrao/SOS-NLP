

      SUBROUTINE CONSAT(CLWR,CUPR,CVEC,XLWR,XUPR,XVEC,MCON,NDIM,
     $    CONTOL,CONTST)
C
C ======================================================================
C     CONSAT===>consat   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  TEST FOR CONSTRAINT SATISFACTION
C                   I.E. CHECK IF THE MCON CONSTRAINTS
C
C                      CLWR < CVEC < CUPR
C
C                   AND BOUNDS
C
C                      XLWR < XVEC < XUPR
C
C         INPUT:
C
C            CLWR   CONSTRAINT LOWER BOUNDS (MCON)
C            CUPR   CONSTRAINT UPPER BOUNDS (MCON)
C            CVEC   CONSTRAINTS AT XVEC (MCON)
C            XLWR   VARIABLE LOWER BOUNDS (NDIM)
C            XUPR   VARIABLE UPPER BOUNDS (NDIM)
C            XVEC   VARIABLES (NDIM)
C            MCON   NUMBER OF CONSTRAINTS
C            NDIM   NUMBER OF VARIABLES
C            CONTOL CONSTRAINT TOLERANCE
C
C         OUTPUT:
C
C            CONTST LOGICAL VARIABLE
C                   = .TRUE.   IF CONSTRAINTS AND BOUNDS ARE SATISFIED
C                   = .FALSE.  OTHERWISE
C
      DIMENSION CLWR(*),CUPR(*),CVEC(*),XLWR(NDIM),XUPR(NDIM),
     $    XVEC(NDIM)
      LOGICAL CONTST
C
      CONTST = .FALSE.
C
C             CHECK ABSOLUTE CONSTRAINT ERROR.
C
      DO I = 1,MCON
C
        IF(CLWR(I).EQ.CUPR(I)) THEN
C
C             EQUALITY CONSTRAINT
C
          IF(ABS(CLWR(I)-CVEC(I)).GT.CONTOL) RETURN
C
        ELSE
C
C             INEQUALITY CONSTRAINT
C
          IF(CVEC(I).LE.CLWR(I)-CONTOL) RETURN
          IF(CVEC(I).GE.CUPR(I)+CONTOL) RETURN
C
        ENDIF
C
      ENDDO
C
      DO I = 1,NDIM
C
        IF(XLWR(I).EQ.XUPR(I)) THEN
C
C             EQUALITY BOUND
C
          IF(ABS(XLWR(I)-XVEC(I)).GT.CONTOL) RETURN
C
        ELSE
C
C             INEQUALITY BOUND
C
          IF(XVEC(I).LE.XLWR(I)-CONTOL) RETURN
          IF(XVEC(I).GE.XUPR(I)+CONTOL) RETURN
C
        ENDIF
C
      ENDDO
C
      CONTST = .TRUE.
C
      RETURN
C
      END
