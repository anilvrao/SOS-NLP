
      SUBROUTINE SHUFLR( ISTATX ,NDIM   ,NVAR   ,NRES   ,IPRMX  
     $          ,RMAT   ,IROWR  ,JCOLR  ,NONZR  ,JCOLRF ,IPRMR  
     $          ,IWORK  ,RWORK  ,IERNLP )
C
C
C ======================================================================
C     SHUFLR===>shuflr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C *** PURPOSE
C
C        THIS ROUTINE APPLIES ONLY TO THE RESIDUAL JACOBIAN MATRIX.  IT
C 
C        (A) CONSTRUCTS THE TRANSFORMATION FROM THE USER DATA (EXTERNAL) 
C            FORMAT TO THE BARRIER ALGORITHM (INTERNAL) FORMAT
C        (B) APPLIES THE TRANSFORMATION AT THE INITIAL POINT
C
C        EXTERNAL FORMAT:
C
C        RMAT ---> NRES X NDIM, STORED AS TRIPLES 
C                  [RMAT(K),IROWR(K),JCOLR(K)];        K=1,...,NONZR
C
C        INTERNAL FORMAT:
C
C        RMAT ---> NRES X NVAR, STORED BY COLUMNS 
C                  [RMAT(K),IROWR(K),JCOLRF(JCOL)];   JCOL = 1,...,NVAR
C                  WHERE IMPLICITLY RMAT = 0 FOR JCOL > NFREE.
C
C *** CALLING ARGUMENTS
C
C        ISTATX   I    INTEGER VARIABLE STATUS (NDIM)
C        NDIM     I    NUMBER OF USER VARIABLES
C        NVAR     I    NUMBER OF VARIABLES (NFREE + NSLK)
C        NRES     I    NUMBER OF RESIDUALS
C        IPRMX    I    INTEGER VARIABLE WORK ARRAY (NDIM)
C        RMAT     IO   RESIDUAL JACOBIAN MATRIX AT XVEC (NONZR)
C        IROWR    IO   ROW INDICES OF RESIDUAL JACOBIAN NONZEROS (NONZR)
C        JCOLR    IO   COLUMN INDICES OF RESIDUAL JACOBIAN (NONZR)
C        NONZR    I    NUMBER OF RESIDUAL JACOBIAN NONZEROS
C        JCOLRF   O    COLUMN INDICES OF INTERNAL JACOBIAN NONZEROS (NVAR+1)
C        IPRMR    O    INTEGER RESIDUAL JACOBIAN PERMUTATION ARRAY (NONZR)
C        IWORK    IO   INTEGER WORK ARRAY (NVAR+1)
C        RWORK    IO   REAL WORK ARRAY (NONZR)
C        IERNLP   O    SUCCESS/ERROR CODE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      DIMENSION ISTATX(NDIM)   ,IPRMX(NDIM)        ,IPRMR(NONZR)
      DIMENSION RMAT(NONZR)    ,IROWR(NONZR)       ,JCOLR(NONZR)
     &         ,JCOLRF(NVAR+1) ,IWORK(NVAR+1)      ,RWORK(NONZR)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         COMPUTE THE DESIRED "NEW" LOCATION FOR THE VARIABLES AND 
C         TEMPORARILY SAVE IN IPRMX
C
      NFIX = 0
      DO I=1,NDIM
        IF (ISTATX(I).EQ.3)  NFIX = NFIX + 1
      ENDDO
      NFREE = NDIM - NFIX
C
      K1 = NFREE
      K2 = 0
      DO I = 1,NDIM
C
        IF(ISTATX(I).EQ.3) THEN
C
C         PUT FIXED VARIABLES AT END OF THE LIST
C
          K1 = K1 + 1
          IPRMX(I) = K1
C 
        ELSE
C
C         PUT FREE VARIABLES AT THE FRONT OF THE LIST
C
          K2 = K2 + 1
          IPRMX(I) = K2
C
        ENDIF
      enddo
C
C         LOOP OVER ALL RESIDUAL JACOBIAN NONZEROS AND COMPUTE NEW LOCATION NUMBER
C         ALSO COUNT THE NUMBER OF NONZEROS IN EACH COLUMN
C
      iwork(1:nvar+1) = 0
      XFREE = NFREE
      XNRES = NRES
      DO KK = 1,NONZR
        XROW = IROWR(KK)
        JCOL = IPRMX(JCOLR(KK))
        IWORK(JCOL) = IWORK(JCOL) + 1
        XCOL = JCOL
        RWORK(KK) = XROW + (XCOL-ONE)*XNRES
      enddo
C
C         ADD UP NONZEROS TO FORM JCOLRF
C
      jcolrf(1:nvar+1) = 0
      JCOLRF(1) = 1
      DO JCOL = 2,NFREE+1
        JCOLRF(JCOL) = JCOLRF(JCOL-1) + IWORK(JCOL-1)
      enddo
C
C         IMPLICITLY "ZERO" OUT THE COLUMNS FOR THE SLACK VARIABLES
C
      DO JCOL = NFREE+2,NVAR+1
        JCOLRF(JCOL) = JCOLRF(JCOL-1)
      enddo
C
C         SORT THE COLUMN INDEX ARRAY RWORK INTO ASCENDING ORDER
C         THEREBY CONSTRUCTING THE MASTER RESIDUAL JACOBIAN PERMUTATION ARRAY
C
      CALL HDSRTN(RWORK,NONZR,0,0,IPRMR,IERP)
C
      IF(IERP.NE.0) THEN
        IERNLP = -147
        RETURN
      ENDIF
C
C         APPLY THE PERMUTATION TO THE RESIDUAL JACOBIAN, 
C         ROW, AND COLUMN STORAGE ARRAYS CONVERTING THEM TO INTERNAL ORDER.
C
      CALL HDPRMX(RMAT,NONZR,IPRMR,IERP)
      IF(IERP.NE.0) THEN
        IERNLP = -147
        RETURN
      ENDIF
      CALL HJPRMX(IROWR,NONZR,IPRMR,IERP)
      IF(IERP.NE.0) THEN
        IERNLP = -147
        RETURN
      ENDIF
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
