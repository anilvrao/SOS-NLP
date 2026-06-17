
      DOUBLE PRECISION FUNCTION COMPLM(MSUBB,PENMU,BVEC,VLAMDA)
C ======================================================================
C     COMPLM===>CMPLMT     J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BVEC(*)  ,VLAMDA(*)
      PARAMETER (ZERO = 0.D0)
C
C         COMPUTE THE BOUND CONSTRAINT CENTRAL PATH CONDITION ERROR
C         CMPLMT = ||BVEC(I)*VLAMDA(I) - PENMU||
C
      CMPLMT = ZERO
      DO I = 1,MSUBB
        BLAM = BVEC(I)*VLAMDA(I)
        CMPLMT = MAX(CMPLMT,ABS(BLAM - PENMU))
      enddo
      COMPLM = CMPLMT
C
      RETURN
      END
