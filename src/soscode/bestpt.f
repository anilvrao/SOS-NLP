
      SUBROUTINE BESTPT(XBAR,XLWR,XUPR,NDIM,CBAR,CLWR,CUPR,
     $    MCON,FBAR,XBEST,CBEST,FBEST) 
C
C ======================================================================
C     BESTPT===>bestpt   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  SAVE/UPDATE THE BEST FEASIBLE POINT FOR THE NLP
C                   COMPUTE THE VALUES OF THE NDIM VARIABLES X WHICH
C
C                             MINIMIZE F(X)
C
C                   SUBJECT TO THE CONSTRAINTS
C
C                             CLWR .LE. CBAR(X) .LE. CUPR  
C
C                   BOUNDS ON THE VARIABLES ARE OF THE FORM
C
C                            XLWR .LE. X .LE. XUPR
C
C         ARGUMENTS:
C
C            XBAR   INDEPENDENT VARIABLES (NDIM)
C            XLWR   VARIABLE LOWER BOUNDS (NDIM)
C            XUPR   VARIABLE UPPER BOUNDS (NDIM)
C            NDIM   NUMBER OF INDEPENDENT VARIABLES
C            CBAR   VECTOR OF CONSTRAINTS (MAXCON)
C            CLWR   CONSTRAINT LOWER BOUNDS (MAXCON)
C            CUPR   CONSTRAINT UPPER BOUNDS (MAXCON)
C            MCON   NUMBER OF CONSTRAINTS
C            FBAR   OBJECTIVE FUNCTION
C            XBEST  VARIABLES AT BEST POINT (NDIM)
C            CBEST  CONSTRAINTS AT BEST POINT (MAXCON)
C            FBEST  OBJECTIVE AT BEST POINT FUNCTION
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      DIMENSION XBAR(NDIM),XLWR(NDIM),XUPR(NDIM),XBEST(NDIM)
      DIMENSION CBAR(*),CLWR(*),CUPR(*),CBEST(*)
C
      LOGICAL FEEZ
C
C         CHECK IF THE POINT IS FEASIBLE 
C
      CALL CONSAT(CLWR,CUPR,CBAR,XLWR,XUPR,XBAR,MCON,NDIM,
     $    CONTOL,FEEZ)
C
      IF(FEEZ) THEN
C
C        THE POINT XBAR IS FEASIBLE
C
        IF(FBAR.GT.FBEST) RETURN
C
C         XBAR IS ALSO BETTER --- SAVE IT 
C         AND UPDATE CBEST, FBEST
C
        XBEST(1:NDIM) = XBAR(1:NDIM)
        IF(MCON.NE.0) CBEST(1:MCON) = CBAR(1:MCON)
        FBEST = FBAR
C
      ENDIF
C
      RETURN
      END
