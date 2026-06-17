
      SUBROUTINE CRSOVR( MAXCON ,MCON   ,AVEC   ,ALWR   ,AUPR   ,VECLAM    
     $          ,ISTATA ,AMAT   ,IROWA  ,JCOLA  ,NONZA  ,DELF   ,XVEC     
     $          ,XLWR   ,XUPR   ,VECNU  ,ISTATX ,NDIM )  
C
C
C ======================================================================
C     CRSOVR===>CRSOVR   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C *** PURPOSE
C
C        THIS ROUTINE CONSTRUCTS AN ESTIMATE OF THE ACTIVE SET FROM
C        THE LAGRANGE MULTIPLIER AND CONSTRAINT VALUES
C
C *** CALLING ARGUMENTS
C
C        MAXCON   I    MAXIMUM NUMBER OF CONSTRAINTS
C        MCON     I    NUMBER OF USER CONSTRAINTS
C        AVEC     I    CONSTRAINTS AT XVEC (MAXCON)
C        ALWR     I    CONSTRAINTS LOWER BOUNDS (MAXCON)
C        AUPR     I    CONSTRAINTS LOWER BOUNDS (MAXCON)
C        VECLAM   I    LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C        ISTATA   O    INTEGER CONSTRAINT STATUS (MAXCON)
C        AMAT     I    CONSTRAINT DERIVATIVES AT XVEC (NONZA)
C        IROWA    I    ROW INDICES OF JACOBIAN NONZEROS (NONZA)
C        JCOLA    I    COLUMN INDICES OF JACOBIAN NONZEROS (NONZA)
C        NONZA    I    NUMBER OF JACOBIAN NONZEROS
C        DELF     I    GRADIENT AT XVEC (NDIM)
C        XVEC     I    CURRENT POINT (NDIM)
C        XLWR     I    VARIABLE LOWER BOUNDS (NDIM)
C        XUPR     I    VARIABLE UPPER BOUNDS (NDIM)
C        VECNU    O    LAGRANGE MULTIPLIERS FOR BOUNDS (NDIM)
C        ISTATX   O    INTEGER VARIABLE STATUS (NDIM)
C        NDIM     I    NUMBER OF USER VARIABLES
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      DIMENSION AVEC(MAXCON)   ,ALWR(MAXCON)   ,AUPR(MAXCON)   
     &         ,VECLAM(MAXCON) ,ISTATA(MAXCON) ,AMAT(NONZA)    
     &         ,IROWA(NONZA)   ,JCOLA(NONZA)   ,DELF(NDIM)     
     &         ,XVEC(NDIM)     ,XLWR(NDIM)     ,XUPR(NDIM) 
     &         ,VECNU(NDIM)    ,ISTATX(NDIM)   
C
      PARAMETER (ZERO=0.0D0,HUNDRD=100.0D0)
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      vecnu(1:ndim) = zero
      IF(MCON.GT.0) THEN
C
C         COMPUTE (AMAT**T)(VECLAM) AND STORE IN VECNU
C
        CALL MVPSRC(11,NDIM,MCON,NONZA,AMAT,IROWA,JCOLA,VECLAM,VECNU)
C
      ENDIF
C
C         COMPLETE CALCULATION OF VECNU = DELF - (AMAT**T)(VECLAM)
C
      DO I = 1,NDIM
        VECNU(I) = DELF(I) - VECNU(I)
        IF(XLWR(I).EQ.XUPR(I)) THEN
C
C         EQUALITY CONSTRAINT
C
          ISTATX(I) = 3
C
        ELSEIF(ABS(XVEC(I)-XLWR(I)).LT.HUNDRD*CONTOL
     $         .AND.VECNU(I).GT.HUNDRD*CONTOL) THEN
C
C         LOWER BOUND INEQUALITY CONSTRAINT
C
          ISTATX(I) = 1
C
        ELSEIF(ABS(XVEC(I)-XUPR(I)).LT.HUNDRD*CONTOL
     $         .AND.VECNU(I).LT.-HUNDRD*CONTOL) THEN
C
C         UPPER BOUND INEQUALITY CONSTRAINT
C
          ISTATX(I) = 2
C
        ELSE
C
C         FREE VARIABLE
C
          ISTATX(I) = 0
C
        ENDIF
      enddo
C
      MLOOP: DO I = 1,MCON
C
        IF(ISTATA(I).EQ.4) CYCLE MLOOP
C
        IF(ALWR(I).EQ.AUPR(I)) THEN
C
C         EQUALITY CONSTRAINT
C
          ISTATA(I) = 3
C
        ELSEIF(ABS(AVEC(I)-ALWR(I)).LE.CONTOL.AND.
     $         VECLAM(I).GT.ZERO) THEN
C
C         LOWER BOUND INEQUALITY CONSTRAINT
C
          ISTATA(I) = 1
C
        ELSEIF(ABS(AVEC(I)-AUPR(I)).LE.CONTOL.AND.
     $         VECLAM(I).LT.ZERO) THEN
C
C         UPPER BOUND INEQUALITY CONSTRAINT
C
          ISTATA(I) = 2
C
        ELSE
C
C         FREE VARIABLE
C
          ISTATA(I) = 0
C
        ENDIF
C
      ENDDO MLOOP
C
      RETURN
C
      END
