 
      SUBROUTINE CONELM(GMAT,IROWG,JSTRG,NONZG,MCON,MEQUAL,NDIM,
     $    RWORK,LNRWRK,RHSWRK,WORK,NWORK,IWORK,LNIWRK,CLWR,CVEC,CBAR,
     $    CONTOL,DVEC,ISTATC,ISTATV,YVEC,YBAR,ITEMP,IFIXVR,IFREVR,
     $    SLOWTL,IOFLAG,IPU,IT,ITMAX,ITMX1,NPANIC,IFC,NCALL,IFERR,
     $    ICTERM)
C
C ======================================================================
C     CONELM===>conelm   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C
C         PURPOSE:  COMPUTE THE VALUES OF THE NDIM VARIABLES YVEC WHICH
C                   SATISFY THE MEQUAL CONSTRAINTS
C
C                               CVEC(YVEC) = CLWR,
C         
C                   WHERE MEQUAL.LT.NDIM AND THE TOTAL NUMBER OF CONSTRAINTS
C                   IS MCON.
C
C         INPUT:
C
C            GMAT   JACOBIAN MATRIX STORED AS AN ARRAY (NONZG)
C            IROWG  ROW INDEX OF NONZERO IN GMAT (NONZG)
C            JSTRG  COLUMN START INDEX OF NONZERO IN GMAT (NDIM+1)
C                   NONZG  = JSTRG(NDIM+1)-1 
C            MCON   NUMBER OF CONSTRAINTS
C            MEQUAL NUMBER OF EQUALITY CONSTRAINTS
C            NDIM   NUMBER OF VARIABLES
C            RWORK    WORK ARRAY (LNRWRK)
C            LNRWRK  LENGTH OF RWORK .GT. MAX(2*NDIM,MCON)
C            RHSWRK RIGHT HAND SIDE WORK ARRAY (4*NDIM)
C            WORK   WORK ARRAY (NWORK)
C            NWORK  LENGTH OF WORK ARRAY 
C            IWORK  INTEGER WORK ARRAY (LNIWRK)
C            LNIWRK LENGTH OF INTEGER WORK ARRAY .GT. 
C                   (2*(NONZG + NDIM) + 18*NEQNS + 30 + NEQNS)
C            CLWR   CONSTRAINT LOWER BOUND (MCON)
C            CVEC   CONSTRAINTS AT YVEC (MCON)
C            CBAR   CONSTRAINTS AT YBAR (MCON)
C            CONTOL CONSTRAINT TOLERANCE
C            DVEC   SEARCH VECTOR (NDIM)
C            ISTATC INTEGER CONSTRAINT STATUS ARRAY (MCON)
C            ISTATV INTEGER VARIABLE STATUS ARRAY (NDIM)
C            YVEC   ESTIMATE OF INDEPENDENT VARIABLES (NDIM)
C            YBAR   NEW ESTIMATE OF YVEC (NDIM)
C            ITEMP  INTEGER TEMPORARY ARRAY (MCON)
C            IFIXVR INTEGER FIXED VARIABLE ARRAY (NDIM)
C            IFREVR INTEGER FREE VARIABLE ARRAY (MCON)
C            SLOWTL SLOW STEP TOLERANCE
C            IOFLAG OUTPUT CONTROL FLAG
C                   = 0   NO OUTPUT
C                   = 10  NORMAL OUTPUT
C            IPU    OUTPUT UNIT NO.
C            IT     ITERATION NUMBER
C            ITMAX  MAXIMUM NUMBER OF ITERATIONS
C            ITMX1  MAXIMUM NUMBER OF ONE DIMENSIONAL SEARCH STEPS
C            NPANIC SLOW STEP TOLERANCE
C            IFC    =1 WHEN EVALUATING FUNCTIONS
C            NCALL  =1 WHEN THIS IS THE FIRST TIME THE ROUTINE IS CALLED
C                   =3 WHEN PREVIOUSLY FACTORED JACOBIAN IS TO BE USED
C            IFERR  =1 WHEN FUNCTION EVALUATION IS IMPOSSIBLE
C            ICTERM TERMINATION FLAG
C                   = 0   NORMAL INPUT VALUE
C
C         OUTPUT:
C
C            CVEC   THE CONSTRAINTS AT YVEC
C            ICTERM TERMINATION FLAG
C                   = 0       NO TERMINATION
C                   = 1       NORMAL TERMINATION
C                   = 2       SMALL STEPS
C                   = 3       MAXIMUM NO. OF ITERATIONS
C                   = 4       INVALID INPUT FOR NCALL 
C                   = 5       NOT ENOUGH STORAGE ALLOCATED
C                   = 6       MAX. NO. OF FUNCTION EVALS.
C            YVEC   THE SOLUTION POINT
C
C  ** SUBROUTINES REQUIRED: 
C
C        DSWAP
C        SPRMVP    SPARSE MATRIX-VECTOR PRODUCT.
C        MFRNLD    MULTIFRONTAL SPARSE LINEAR SYSTEM SOLVER
C
C
C        NOTE:   1.INPUT CONSTRAINTS CVEC EVALUATED AT YVEC MUST BE SUPPLIED
C                  ON FIRST ENTRY.  THEREAFTER CBAR EVALUATED AT YBAR
C                  IS SUPPLIED WHEN IFC = 1
C
C
C     ******************************************************************
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,
     $           ONEEP1=1.0D1,ONEEM2=1.0D-2,POINT2=2.0D-1,POINT5=5.0D-1)
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
C         INTERNAL: LOCAL VARIABLES NEEDED FROM ITERATION TO ITERATION
C                   ARE SAVED USING THE "SAVE" STATEMENT.
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
C
      DIMENSION  GMAT(NONZG),IROWG(NONZG),JSTRG(NDIM+1),RWORK(LNRWRK),
     $   RHSWRK(4*NDIM),WORK(NWORK),IWORK(LNIWRK),CLWR(MCON),
     $   CVEC(MCON),CBAR(MCON),DVEC(NDIM),ISTATC(MCON),YVEC(NDIM),
     $   YBAR(NDIM),ITEMP(MCON),IFIXVR(NDIM),IFREVR(MCON),ISTATV(NDIM)
C
C
C     ******************************************************************
C
      IF(IFC.EQ.1) GO TO 501
C
C ----------------------------------------------------------------------
C ----------------- INITIALIZATION -------------------------------------
C ----------------------------------------------------------------------
C
      IER = 0
      IT = 1
      NSLOW = 0
      SMALL = SQRT(ZEROMN)
      XMU = ONE
C
C         CONSTRAINT ERROR AT INITIAL POINT
C
C         COMPUTE PHI = ((CVEC-CLWR)**T)*(CVEC-CLWR)/CONTOL
C
      PHI = ZERO
      DO I = 1,MEQUAL
        CIA = CVEC(I)-CLWR(I)
        PHI = PHI + CIA*(CIA/CONTOL)
        CBAR(I) = CIA
      ENDDO
C
C         CHECK NCALL FLAG
C
      IF(NCALL.NE.1.AND.NCALL.NE.3) THEN
        ICTERM = 4
        GO TO 230
      ENDIF
C
C         COMPUTE LENGTH OF INTEGER WORK ARRAY 
C
      NEQNS = MEQUAL + NDIM
      IWRKCK = 2*(NONZG + NDIM) + 18*NEQNS + 30 + NEQNS
      IF(LNIWRK.LT.IWRKCK) THEN
        NEEDS = IWRKCK
        IF(IOFLAG.GE.10) WRITE(IPU,1002) LNIWRK,NEEDS
        ICTERM = 5
        GO TO 230
      ENDIF
C
C
C ----------------------------------------------------------------------
C ----------------- BEGIN ITERATION ------------------------------------
C ----------------------------------------------------------------------
C
 120  CONTINUE
C
C ----------------------------------------------------------------------
C ----------------- COMPUTE SEARCH DIRECTION AND MAGNITUDE--------------
C ----------------------------------------------------------------------
C
C             ITERATION PRINT
C
      IF(IOFLAG.GE.10)  WRITE(IPU,1001) IT,SQRT(PHI*CONTOL)
C
C         DEFINE CONSTRAINT POINTER ARRAY
C
      DO I = 1,MCON
        IF(ISTATC(I).EQ.3) THEN
          ITEMP(I) = -1
        ELSE
          ITEMP(I) = 0
        ENDIF
      ENDDO
C
C         DEFINE FIXED VARIABLE POINTER
C
      NFIXVR = 0
      DO I = 1,NDIM
        IF(ISTATV(I).EQ.3) THEN
          NFIXVR = NFIXVR + 1
          IFIXVR(I) = 1
        ELSE
          IFIXVR(I) = 0
        ENDIF
      ENDDO
C
C           LOAD ZERO IN GRADIENT VECTOR SLOT
C
      RWORK(1:NDIM) = ZERO
C
C         COMPUTE SEARCH VECTOR USING SPARSE UNDERDETERMINED SYSTEM
C         SOLVER 
C
C ----------------------------------------------------------------------
C
C         CALL SPARSE SYSTEM SOLVER TO COMPUTE 
C         SEARCH DIRECTION.  NOTE RWORK CONTAINS GRADIENT ON INPUT
C         AND LAGRANGE MULTIPLIER ON OUTPUT.
C
      CALL MFRNLD(MCON,NDIM,GMAT,IROWG,JSTRG,NONZG,RWORK
     $    ,CVEC,ITEMP,IFIXVR,DVEC,RWORK,RWORK(NDIM+1)
     $    ,RHSWRK,WORK,NWORK,IWORK
     $    ,LNIWRK,NEEDED,IFREVR
     $    ,NCALL,IPU,IOFLAG,CNDNUM,IER)
C
C
C             CHECK FOR SUCCESSFUL TERMINATION.
C
      IF (IER.EQ.0)
     $  THEN
C
C             SET NCALL TO 3, SO THAT SUBSEQUENT CALLS TO 
C             L.D.P. WILL USE THE SAME FACTORIZATION OF THE
C             COEFFICIENT MATRIX.  NOTE WORK ARRAY OF LENGTH 
C             NWORK MUST NOT BE CHANGED.
C
          NCALL = 3
C
C
          SNORM = ZERO
          DO I=1,NDIM
            SNORM = SNORM + DVEC(I)**2
          ENDDO
          SNORM = SQRT(SNORM)
C
C             GO TO END OF IF BLOCK AFTER STATEMENT 210 
C
        ELSE
C
C             IER .NE. 0 
C
          IF (IOFLAG.GE.10) THEN
              WRITE(IPU,1009) IER
              WRITE(IPU,1004)
          ENDIF
C
C             COMPUTE GRADIENT DIRECTION = 2*(GMAT**T)*(CVEC-CLWR)
C
          RWORK(MEQUAL+1:MCON) = ZERO
          DO I = 1,MEQUAL
            RWORK(I) = TWO*(CVEC(I)-CLWR(I))
          ENDDO
C
C             COMPUTE PRODUCT (GMAT**T)*RWORK = DVEC
C
          CALL MVPSPR(11,NDIM,MCON,GMAT,IROWG,JSTRG,RWORK,DVEC)
C
          GNORM = ZERO
          DO I=1,NDIM
            GNORM = GNORM + DVEC(I)**2
          ENDDO
          SNORM = SQRT(GNORM)
C
C             COMPUTE STEP LENGTH FROM SECOND DERIVATIVE
C             GMG = ((GMAT*DVEC)**T)(GMAT*DVEC).  FORM
C             (GMAT*DVEC) IN RWORK.
C
          CALL MVPSPR(1,MCON,NDIM,GMAT,IROWG,JSTRG,DVEC,RWORK)
C
C             COMPUTE SUM OF SQUARES OF EQUALITY ELEMENTS
C
          GMG = ZERO
          DO I = 1,MEQUAL
            GMG = GMG + RWORK(I)**2
          ENDDO
C
          IF (GMG.NE.ZERO)
     $      THEN
              STP = GNORM/GMG
              DVEC(1:NDIM) = STP*DVEC(1:NDIM)
              SNORM = STP*SNORM
          ENDIF
C
      ENDIF
C
C             END OF IER TEST.
C
      IF(IT.EQ.1) STEPL = SNORM
C
C ----------------------------------------------------------------------
C ----------------- CONVERGENCE TESTS ----------------------------------
C ----------------------------------------------------------------------
C
C             CHECK ABSOLUTE CONSTRAINT ERROR.
C
      DO I = 1,MEQUAL
        CTEST = ABS(CVEC(I)-CLWR(I))
        IF(CTEST.GT.CONTOL) GO TO 180
      ENDDO
C
C         SET NORMAL TERMINATION FLAG
C
      ICTERM = 1
      GO TO 230
C
 180  CONTINUE
C
C         CHECK FOR SLOW PROGRESS
C
      IF(NSLOW.GT.NPANIC.OR.ABS(XMU).LT.SMALL) THEN
        ICTERM = 2
        GO TO 230
      ENDIF
C
C             CHECK FOR MAXIMUM ITERATION COUNT.
C
      IF (IT.GT.ITMAX)
     $  THEN
          ICTERM = 3
          GO TO 230
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------- ONE DIMENSIONAL SEARCH -----------------------------
C ----------------------------------------------------------------------
C
C             INITIALIZATION.
C
C             ON FIRST STEP SET LENGTH TO AT MOST 1.
C
      XMU = MIN(ONE,ONEEP1*STEPL/MAX(ZEROMN,SNORM))
C
C         CHOOSE LYNIMP = 1 --- LIMITED TESTING SHOWED ACCURATE 
C         LINE SEARCH REDUCED RUN TIME (FEWER RHS SOLVES), BUT INCREASED
C         NO. OF F.E.
C
      LYNIMP = 1
      MAXFER = 5
      IFC = 0
C
 190  CONTINUE
C
      CALL LYNSRC(CLWR,CVEC,CBAR,DF0,D2F0,PHIBAR,
     $    PHI,CONTOL,XMU,LYNTRM,MEQUAL,LYNIMP,
     $    ITMX1,MAXFER,IFERR,IFC,IOFLAG,IPU)
C
C         IF THE RETURN IS NOT FOR AN EVALUATION TERMINATE
C
      IF(IFC.NE.1) GO TO 220
C
C
C ----------------------------------------------------------------------
C
C             FUNCTION EVALUATION.
C
C             COMPUTE NEW ESTIMATE OF VARIABLES.
C
      DO I = 1,NDIM
        YBAR(I) = YVEC(I) - XMU*DVEC(I)
      ENDDO
C
C             RETURN FOR FUNCTION EVALUATION.
C
      RETURN
 501  CONTINUE
C
C             THE FUNCTION EVALUATION IS COMPLETED.
C
      IF(IFERR.EQ.0) THEN
C
C         COMPUTE CONSTRAINT ERROR
C
      PHIBAR = ZERO
      DO I = 1,MEQUAL
        CIA = CBAR(I)-CLWR(I)
        PHIBAR = PHIBAR + CIA*(CIA/CONTOL)
      ENDDO
C
      ENDIF
C
      GO TO 190
C
C ----------------------------------------------------------------------
C ----------------- UPDATE INFORMATION ---------------------------------
C ----------------------------------------------------------------------
C
 220  CONTINUE
C
      IF(LYNTRM.EQ.-5) THEN
        ICTERM = 6
        GO TO 230
      ENDIF
C
C         COMPUTE ACTUAL REDUCTION RATIO
C
      ACTRED = SQRT(PHIBAR)/SQRT(PHI)
C
C         COUNT SUCCESSIVE SLOW ITERATIONS -- AN ITERATION IS CALLED 
C         'SLOW' IF THE CONSTRAINT ERROR DID NOT GO DOWN BY AT LEAST 
C         A FACTOR OF SLOWTL
C
      IF(ACTRED.GT.SLOWTL) THEN
        NSLOW = NSLOW + 1
      ELSE
        NSLOW = 0
      ENDIF
C
C             UPDATE POINT, YVEC.
C
      YVEC(1:NDIM) = YBAR(1:NDIM)
C
C         UPDATE CONSTRAINTS--STORE NEW VALUE IN CVEC
C
      CVEC(1:MEQUAL) = CBAR(1:MEQUAL)
C
C             UPDATE CONSTRAINT ERROR, PHI, STEP LENGTH, STEPL, AND
C             ITERATION COUNTER, IT.
C
      PHI = PHIBAR
      STEPL = ABS(XMU)*SNORM
      IT = IT + 1
C
C ----------------------------------------------------------------------
C ----------------- END OF ITERATION -----------------------------------
C ----------------------------------------------------------------------
C
      GO TO 120
C
C ----------------------------------------------------------------------
C ----------------- ALGORITHM TERMINATION PROCESSING -------------------
C ----------------------------------------------------------------------
C
 230  CONTINUE
C
C         GET ALL QUANTITIES IN SYNC AT THE FINAL POINT
C
      YBAR(1:NDIM) = YVEC(1:NDIM)
C
      CBAR(1:MEQUAL) = CVEC(1:MEQUAL)
C
      IF(ICTERM.EQ.1) THEN
        INSTAT(4) = MAX(INSTAT(4),IT)
        INSTAT(5) = INSTAT(5) + IT
        INSTAT(6) = INSTAT(6) + 1
      ELSEIF(ICTERM.GT.1) THEN
        INSTAT(7) = MAX(INSTAT(7),IT)
        INSTAT(8) = INSTAT(8) + IT
        INSTAT(9) = INSTAT(9) + 1
      ENDIF
C
      RETURN
C
 1001 FORMAT(T3,'*',T106,'*'/
     $    T3,'*',T16,'Constraint Iteration',I5,T50,'Error =',
     $    1PG16.8,T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'INTEGER WORK ARRAY FOR SPARSE S
     $EARCH DIMENSIONED ',I6,' NEED ',I6,T106,'*')
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'USE GRADIENT DIRECTION',
     $    T106,'*')
 1009 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'LINEAR SYSTEM ALGORITHM FAILED 
     $-- IER =',I5,T106,'*')
C
      END
