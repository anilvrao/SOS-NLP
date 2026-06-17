
      SUBROUTINE SLAKER(BVEC,CVEC,YVEC,MAXBND,MAXCON,   
     $        MEQUAL,MINEQL,NALWR,NAUPR,NFREE,NVARMX,NXLWR,NXUPR)
C
C
C ======================================================================
C     SLAKER===>SLAKER   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C *** PURPOSE
C
C        THIS ROUTINE CONSTRUCTS THE SLACK VARIABLES AND MODIFIED
C        CONSTRAINTS AND BOUNDS FOR THE RELAXED FORMULATION GIVEN
C        INFORMATION FROM THE STANDARD FORMULATION
C
C *** CALLING ARGUMENTS
C
C        BVEC     IO   BOUND INEQUALITIES (MAXBND)
C        CVEC     IO   EQUALITY CONSTRAINTS AT YVEC (MAXCON)
C        YVEC     IO   CURRENT POINT (NVARMX)
C        MAXBND   I    MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C        MAXCON   I    MAXIMUM NUMBER OF CONSTRAINTS
C        MEQUAL   I    NUMBER OF USER EQUALITY CONSTRAINTS
C        MINEQL   I    NUMBER OF USER INEQUALITY CONSTRAINTS
C        NALWR    I    NUMBER OF CONSTRAINT LOWER BOUNDS
C        NAUPR    I    NUMBER OF CONSTRAINT UPPER BOUNDS
C        NFREE    I    NUMBER OF FREE USER VARIABLES
C        NVARMX   I    MAX(NVAR,NDIM)
C        NXLWR    I    NUMBER OF VARIABLE LOWER BOUNDS
C        NXUPR    I    NUMBER OF VARIABLE UPPER BOUNDS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      DIMENSION BVEC(MAXBND) ,CVEC(MAXCON) ,YVEC(NVARMX)       
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      LCSVEC = NFREE + 1
      LCTVEC = LCSVEC + MINEQL
      LCUVEC = LCTVEC + MEQUAL
      LCVVEC = LCUVEC + MEQUAL
      LCWVEC = LCVVEC + NALWR
C
C         INITIALIZE THE EQUALITY SLACKS (T AND U)
C
      FEASTL = FEATOL
      IF(NFREE.EQ.MEQUAL) FEASTL = ZEROOT
      DO I = 1,MEQUAL
        FEASTP = MAX(FEASTL,FEASTL*ABS(CVEC(I)))
        IF(CVEC(I).GT.ZERO) THEN
          YVEC(LCTVEC+I-1) = CVEC(I) + FEASTP
          YVEC(LCUVEC+I-1) = FEASTP
        ELSE
          YVEC(LCTVEC+I-1) = FEASTP
          YVEC(LCUVEC+I-1) = FEASTP - CVEC(I)
        ENDIF
      enddo
C
C         INITIALIZE THE INEQUALITY SLACKS (V AND W)
C
      DO I = 1,NALWR
        II = NXLWR + NXUPR + I
        FEASTP = MAX(FEASTL,FEASTL*ABS(BVEC(I)))
        YVEC(LCVVEC+I-1) = FEASTP
      enddo
C
      DO I = 1,NAUPR
        II = NXLWR + NXUPR + NALWR + I
        FEASTP = MAX(FEASTL,FEASTL*ABS(BVEC(I)))
        YVEC(LCWVEC+I-1) = FEASTP
      enddo
C
C         COMPLETE CONSTRUCTION OF CVEC 
C
      DO I = 1,MEQUAL
        CVEC(I) = CVEC(I) - YVEC(LCTVEC+I-1) 
     $                               + YVEC(LCUVEC+I-1)
      enddo
C
      LNSLAK = 2*MEQUAL + NALWR + NAUPR
      LC1 = 1
      LC2 = LC1 + NXLWR
      LC3 = LC2 + NXUPR
      LC4 = LC3 + NALWR
      LC5 = LC4 + NAUPR
C
      DO I = 1,NALWR
        II = LC3 + I - 1
        BVEC(II) = BVEC(II) + YVEC(LCVVEC+I-1)
      enddo
C
      DO I = 1,NAUPR
        II = LC4 + I - 1
        BVEC(II) = BVEC(II) + YVEC(LCWVEC+I-1)
      enddo
C
      BVEC(LC5:LC5+LNSLAK-1) = YVEC(LCTVEC:LCTVEC+LNSLAK-1)
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      RETURN
C
      END
