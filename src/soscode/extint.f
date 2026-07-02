
      SUBROUTINE EXTINT( MODE   ,RELAX  ,IREVRS ,IPRMC  ,IPRMG   
     $          ,IPRMX  ,IPRMH  ,MAXCON ,MCON   ,MEQUAL ,MIGNOR    
     $          ,AVEC   ,VECLAM ,AMAT   ,IROWA  ,JCOLA  ,NONZA    
     $          ,DELF   ,XVEC   ,NDIM   ,NFREE  ,HMAT   ,IROWH    
     $          ,JSTRH  ,NONZH  ,FBAR   ,CVEC   ,MSUBE  ,ALWRI       
     $          ,AUPRI  ,ETAVEC ,CMAT   ,IROWC  ,JCOLC  ,NONZC   
     $          ,GVEC   ,YVEC   ,NVAR   ,XLWRI  ,XUPRI  ,BVEC     
     $          ,MSUBB  ,NXUPR  ,NXLWR  ,NAUPR  ,NALWR  ,IBOUND  
     $          ,VLAMDA ,WMAT   ,IROWW  ,JSTRW  ,NONZW  ,FNLP     
     $          ,SCRTCH ,IFERR  ,IEREXI)
C
C
C ======================================================================
C     HDSnyy===>extint   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C *** PURPOSE
C
C        THIS ROUTINE TRANSFORMS THE USER DATA (EXTERNAL) FORMAT
C        TO THE BARRIER ALGORITHM (INTERNAL) FORMAT AND VICE VERSA.
C
C *** CALLING ARGUMENTS
C
C        MODE         MODE OF OPERATION
C                     = 1 CONVERT EXTERNAL TO INTERNAL FORMAT
C                     = 2 CONVERT INTERNAL TO EXTERNAL FORMAT
C        IREVRS       REVERSE COMMUNICATION CONTROL INDICATOR
C
C        PERMUTATION ARRAYS
C
C        IPRMC        INTEGER CONSTRAINT PERMUTATION ARRAY (MAXCON)
C        IPRMG        INTEGER JACOBIAN PERMUTATION ARRAY (NONZA)
C        IPRMX        INTEGER VARIABLE PERMUTATION ARRAY (NDIM)
C        IPRMH        INTEGER HESSIAN PERMUTATION ARRAY (NONZH)
C
C        EXTERNAL FORMAT
C
C        MAXCON       MAXIMUM NUMBER OF CONSTRAINTS
C        MCON         NUMBER OF USER CONSTRAINTS
C        MEQUAL       NUMBER OF USER EQUALITY CONSTRAINTS
C        MIGNOR       NUMBER OF USER IGNORED CONSTRAINTS
C        AVEC         CONSTRAINTS AT XVEC (MAXCON)
C        VECLAM       LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C        AMAT         CONSTRAINT DERIVATIVES AT XVEC (NONZA)
C        IROWA        ROW INDICES OF JACOBIAN NONZEROS (NONZA)
C        JCOLA        COLUMN INDICES OF JACOBIAN NONZEROS (NONZA)
C        NONZA        NUMBER OF JACOBIAN NONZEROS
C        DELF         GRADIENT AT XVEC (NDIM)
C        XVEC         CURRENT POINT (NDIM)
C        NDIM         NUMBER OF USER VARIABLES
C        NFREE        NUMBER OF FREE USER VARIABLES
C        HMAT         HESSIAN MATRIX AT XVEC (NONZH)
C        IROWH        ROW INDICES OF HESSIAN NONZEROS (NONZH)
C        JSTRH        COLUMN START INDICES OF HESSIAN (NDIM+1)
C        NONZH        NUMBER OF HESSIAN NONZEROS
C        FBAR         USER OBJECTIVE FUNCTION
C
C        INTERNAL FORMAT
C
C        CVEC         EQUALITY CONSTRAINTS AT YVEC (MSUBE)
C        MSUBE        NUMBER OF EQUALITY CONSTRAINTS MSUBE = MCON - MIGNOR
C        ALWRI        PERMUTED USER CONSTRAINT LOWER BOUNDS (MAXCON)
C        AUPRI        PERMUTED USER CONSTRAINT UPPER BOUNDS (MAXCON)
C        ETAVEC       LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MSUBE)
C        CMAT         CONSTRAINT DERIVATIVES AT YVEC (NONZC)
C        IROWC        ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C        JCOLC        COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C        NONZC        NUMBER OF JACOBIAN NONZEROS
C        GVEC         GRADIENT AT YVEC (NVAR)
C        YVEC         CURRENT POINT (NVAR)
C        NVAR         NUMBER OF INTERNAL VARIABLES;  NVAR = NFREE + NSLK
C        XLWRI        PERMUTED USER VARIABLE LOWER BOUNDS (NDIM)
C        XUPRI        PERMUTED USER VARIABLE UPPER BOUNDS (NDIM)
C        BVEC         BOUND INEQUALITIES (MSUBB)
C        MSUBB        NUMBER OF BOUNDS
C        NXUPR        NUMBER OF VARIABLE UPPER BOUNDS
C        NXLWR        NUMBER OF VARIABLE LOWER BOUNDS
C        NAUPR        NUMBER OF CONSTRAINT UPPER BOUNDS
C        NALWR        NUMBER OF CONSTRAINT LOWER BOUNDS
C        IBOUND       INTEGER BOUND DEFINITION (MSUBB)
C        VLAMDA       LAGRANGE MULTIPLIERS FOR BOUNDS (MSUBB)
C        WMAT         HESSIAN MATRIX AT YVEC (NONZW)
C        IROWW        ROW INDICES OF HESSIAN NONZEROS (NONZW)
C        JSTRW        COLUMN START INDICES OF HESSIAN (NVAR+1)
C        NONZW        NUMBER OF HESSIAN NONZEROS
C        FNLP         NLP OBJECTIVE FUNCTION
C
C        WORK ARRAYS
C
C        SCRTCH       SCRATCH ARRAY (MAXCON)
C
C        IFERR        FUNCTION ERROR FLAG
C        IEREXI       SUCCESS/ERROR CODE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      DIMENSION IREVRS(5)
      DIMENSION IPRMC(MAXCON)  ,IPRMG(NONZA)       ,IPRMX(NDIM)
     &         ,IPRMH(NONZH)
      DIMENSION AVEC(MAXCON)   ,VECLAM(MAXCON)     ,AMAT(NONZA)    
     &         ,IROWA(NONZA)   ,JCOLA(NONZA)       ,DELF(NDIM)     
     &         ,XVEC(NDIM)     ,HMAT(NONZH)        ,IROWH(NONZH)    
     &         ,JSTRH(NDIM+1)  ,CVEC(MSUBE)        ,ALWRI(MAXCON)
     &         ,AUPRI(MAXCON)  ,ETAVEC(MSUBE)      ,CMAT(NONZC)
     &         ,IROWC(NONZC)   ,JCOLC(NVAR+1)      ,GVEC(NVAR)
     &         ,YVEC(NVAR)     ,XLWRI(NDIM)        ,XUPRI(NDIM)
     &         ,BVEC(MSUBB)    ,IBOUND(MSUBB)      ,VLAMDA(MSUBB)
     &         ,WMAT(NONZW)    ,IROWW(NONZW)       ,JSTRW(NVAR+1)
     &         ,SCRTCH(MAXCON) 
C
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, POINT5=0.5D0)
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      LOGICAL FZMODE,RELAX
C
      IEREXI = 0
      FZMODE = MODE.GT.0
      LODE = ABS(MODE)
C
      IF(LODE.EQ.1) THEN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C --------------------EXTERNAL TO INTERNAL (RELAXATION)-----------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
        IF(IFERR.NE.0) RETURN
C
        MINEQL = MCON - MEQUAL - MIGNOR
        NSLK = NVAR - NFREE 
C
        LCSVEC = NFREE + 1
        LCTVEC = LCSVEC + MINEQL
        LCUVEC = LCTVEC + MEQUAL
        LCVVEC = LCUVEC + MEQUAL
        LCWVEC = LCVVEC + NALWR
        IF(FZMODE) THEN
          FNLP = ZERO
        ELSE
          FNLP = FBAR
        ENDIF
C
        IF(MCON.GT.0) THEN
C
C         COPY EXTERNAL QUANTITIES INTO FIRST MCON LOCATIONS
C
          CVEC(1:MCON) = AVEC(1:MCON)
C
C         RESTORE ALL QUANTITIES TO INTERNAL ORDER
C
          CALL HDPRMX(CVEC,MCON,IPRMC,IERP)
          IF(IERP.NE.0) THEN
            IEREXI = -147
            RETURN
          ENDIF
          IF(IREVRS(2).GE.1.OR.(IREVRS(3).GE.1.AND.IHESHN.NE.0)) THEN
            CALL HDPRMX(AMAT,NONZA,IPRMG,IERP)
            IF(IERP.NE.0) THEN
              IEREXI = -147
              RETURN
            ENDIF
            NZC = JCOLC(NFREE+1)-1
            CMAT(1:NZC) = AMAT(1:NZC)
            CALL HDPRMY(AMAT,NONZA,IPRMG,IERP)
            IF(IERP.NE.0) THEN
              IEREXI = -147
              RETURN
            ENDIF
          ENDIF
C
C         COMPLETE CONSTRUCTION OF CVEC 
C
          VIOLMX = ZERO
          DO I = 1,MEQUAL
            CVEC(I) = CVEC(I) - ALWRI(I) - YVEC(LCTVEC+I-1) 
     $                                   + YVEC(LCUVEC+I-1)
            VIOLMX = MAX(VIOLMX,ABS(CVEC(I)))
          enddo
C
          DO I = 1,MINEQL
            II = MEQUAL + I
            VIOL = MAX(ZERO,ALWRI(II)-CVEC(II),CVEC(II)-AUPRI(II))
            VIOLMX = MAX(VIOLMX,VIOL)
            CVEC(II) = CVEC(II) - YVEC(LCSVEC+I-1) 
          enddo
C
          IF(FZMODE) THEN
C
            M = MEQUAL + MINEQL
            FNLP = POINT5*DOT_PRODUCT(CVEC(1:M),CVEC(1:M))
C
C         CHECK FOR CONSTRAINT SATISFACTION 
C
            IF(VIOLMX.LT.CONTOL) IFERR = -101
C
          ENDIF
C
        ENDIF
C
        DO I = 1,NFREE
C
C         CONSTRUCT THE INTERNAL VARIABLES Y
C
          YVEC(I) = XVEC(IPRMX(I))
C
        enddo
C
C         COMPLETE CONSTRUCTION OF BVEC
C
        LC1 = 1
        LC2 = LC1 + NXLWR
        LC3 = LC2 + NXUPR
        LC4 = LC3 + NALWR
        LC5 = LC4 + NAUPR
        LC6 = LC5 + 2*MEQUAL + NALWR + NAUPR
C
        DO JJ = LC1,LC2-1
          II = IBOUND(JJ)
          BVEC(JJ) = YVEC(II) - XLWRI(II)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
        enddo
C
        DO JJ = LC2,LC3-1
          II = IBOUND(JJ)
          BVEC(JJ) = XUPRI(II) - YVEC(II)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
        enddo
C
        DO JJ = LC3,LC4-1
          KK = JJ - LC3 
          II = IBOUND(JJ)
          KCOL1 = LCSVEC + II - MEQUAL - 1
          KCOL2 = LCVVEC + KK 
          BVEC(JJ) = YVEC(KCOL1) + YVEC(KCOL2) - ALWRI(II)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
        enddo
C
        DO JJ = LC4,LC5-1
          KK = JJ - LC4 
          II = IBOUND(JJ)
          KCOL1 = LCSVEC + II - MEQUAL - 1
          KCOL2 = LCWVEC + KK 
          BVEC(JJ) = AUPRI(II) - YVEC(KCOL1) + YVEC(KCOL2)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
        enddo
C
        KCOL = LCTVEC - 1
        DO JJ = LC5,LC6-1
          KCOL = KCOL + 1
          BVEC(JJ) = YVEC(KCOL)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
        enddo
C
        IF(IREVRS(2).GE.1) THEN
C
          IF(FZMODE) THEN
C
C         COMPUTE GRADIENT = (CMAT**T)(CVEC) 
C
            IF(RELAX) THEN
              NVRGVC = NVAR
            ELSE
              NVRGVC = NFREE + MINEQL
            ENDIF
C
            CALL MVPSPR(11,NVRGVC,MSUBE,CMAT,IROWC,JCOLC,CVEC,GVEC)
C
          ELSE
C
            DO I = 1,NFREE
C
C         CONSTRUCT GRADIENT WITH RESPECT TO INTERNAL VARIABLES Y
C
              GVEC(I) = DELF(IPRMX(I))
C
            enddo
C
            IF(NSLK.GT.0) GVEC(NFREE+1:NFREE+NSLK) = ZERO
C
          ENDIF
C
          IF(RELAX.AND.NSLK.GT.0) THEN
            GVEC(NFREE+MINEQL+1:NFREE+NSLK) = ONE
          ENDIF
C
        ENDIF
C
        IF(IREVRS(3).GE.1) THEN
          WMAT(1:NONZH) = HMAT(1:NONZH)
          CALL HDPRMX(WMAT,NONZH,IPRMH,IERP)
          IF(IERP.NE.0) THEN
            IEREXI = -147
            RETURN
          ENDIF
C
C         ZERO OUT THE HESSIAN ELEMENTS CORRESPONDING TO SLACKS
C
          IF(NSLK.GT.0) WMAT(NONZW-NSLK+1:NONZW) = ZERO
C
        ENDIF
C
      ELSEIF(LODE.EQ.2) THEN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C --------------------INTERNAL TO EXTERNAL (RELAXATION)-----------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
        IF(MCON.GT.0) THEN
C
C         COPY INTERNAL QUANTITIES INTO FIRST MSUBE LOCATIONS
C
          IF(FZMODE) THEN
            IREVRS(5) = 4
            VECLAM(1:MSUBE) = -CVEC(1:MSUBE)
          ELSE
            VECLAM(1:MSUBE) = ETAVEC(1:MSUBE)
          ENDIF
          IF(MIGNOR.GT.0) VECLAM(MSUBE+1:MSUBE+MIGNOR) = ZERO
C
C         RESTORE ALL QUANTITIES TO EXTERNAL ORDER
C
          CALL HDPRMY(VECLAM,MCON,IPRMC,IERP)
          IF(IERP.NE.0) THEN
            IEREXI = -147
            RETURN
          ENDIF
C
        ENDIF
C
C         CONSTRUCT EXTERNAL FROM INTERNAL 
C
        DO I = 1,NFREE
C
C         CONSTRUCT THE EXTERNAL VARIABLES X
C
          XVEC(IPRMX(I)) = YVEC(I)

        enddo
C
      ELSE
        RETURN
      ENDIF
C
      RETURN
C
      END
