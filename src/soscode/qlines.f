

      SUBROUTINE QLINES(QZERO,QPRIME,QBAR,ALFA,ALFBND,PENMU,CBTCB,
     $    C0TC0,CNORM,IT,MAXIT,NWFLTR,IFUNER,IOFLAG,
     $    IPU,IFC,IERLIN,SAMSCL)
C
C
C ======================================================================
C     QLINES===>QLINES   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C
C         PURPOSE:  COMPUTE AN ESTIMATE ALFA, WHICH WILL EITHER
C
C                   (A) REDUCE THE LOG-BARRIER OBJECTIVE FUNCTION
C                       B(alfa) = f(alfa) - mu sum (ln b_i(alfa))
C
C                   (B) OR, REDUCE THE CONSTRAINT ERROR
C                       ||c(alfa)||
C
C                   THE SINGULARITY IN THE LOG BARRIER FUNCTION IS
C                   LOCATED AT ALFA = ALFBND
C
C
C         ARGUMENTS:
C
C            QZERO  LOG-BARRIER FUNCTION AT ALFA=0
C            QPRIME SLOPE OF LOG-BARRIER FUNCTION ALFA=0
C            QBAR   LOG-BARRIER FUNCTION AT ALFA
C            ALFA   AN ESTIMATE FOR THE STEP LENGTH
C            ALFBND AN ESTIMATE FOR THE NEAREST INEQUALITY BOUNDARY
C            PENMU  THE BARRIER PARAMETER, MU
C            CBTCB  THE DOT PRODUCT (cbar)^T(cbar)
C            C0TC0  THE DOT PRODUCT (czero)^T(czero)
C            CNORM  NORM OF CONSTRAINT ERROR
C            IT     ITERATION NUMBER. (INITIALIZED TO 1)
C            MAXIT  MAXIMUM NUMBER OF ITERATIONS
C            NWFLTR INTEGER FLAG DENOTING BLOCKING FILTER ENTRY
C                   (WHEN IFUNER=2)
C            IFUNER FUNCTION ERROR FLAG
C                   = 0    FUNCTION EVALUATED
C                   = 1    FUNCTION EVALUTION IMPOSSIBLE
C                   = 2    FUNCTION EVALUATED BUT VIOLATES FILTER
C                   .LT. 0 TERMINATE IMMEDIATELY
C            IOFLAG OUTPUT CONTROL FLAG
C                   =0   NO OUTPUT
C                   =10  NORMAL OUTPUT
C            IPU    OUTPUT UNIT NO.
C            IFC    INTEGER VARIABLE:  IFC = 1 WHEN RETURNING TO THE
C                   CALLING ROUTINE TO EVALUATE FBAR(ALFA)
C                   IFC = 0 WHEN ALGORITHM TERMINATES
C            IERLIN INTEGER ERROR FLAG (MUST BE CHECKED WHEN IFC=0)
C                   = 0     NORMAL RETURN, WHEN IFC=0
C                   ----WARNING ERRORS
C                   = 1     NEGATIVE CURVATURE AT ALFA=0
C                   = 2     SMALL STEPS (ALFA*SNORM .LE. ZEROOT)
C                   = 3     MAXIMUM NUMBER OF FUNCTION ERRORS (MAXFER
C                           SET IN PARAMETER STATEMENT BELOW)
C                   = 4     LINEAR OBJECTIVE FUNCTION
C                   = 5     MAXIMUM NUMBER OF ITERATIONS WITH NO IMPROVED 
C                           POINT
C                   = 7     MAX. NUMBER OF F.E.
C                   ----FATAL ERRORS
C                   = -1    QPRIME > 0
C                   = -2    ALFBND .LE. 0
C            SAMSCL LOGICAL FLAG; TRUE = SAME PRIMAL-DUAL STEP SCALING
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0,EIGHT=8.0D0,
     $    ONEEM1=1.0D-1,ONEEM5=1.0D-5,POINT5=5.0D-1)
      PARAMETER (MAXFER=5,REDTOL=ONEEM5)
C
      LOGICAL SAMSCL
C      
      CHARACTER(LEN=60) BLANK
C
      CHARACTER(LEN=18) MODEL
      CHARACTER(LEN=14) LYNE
      CHARACTER(LEN=32) TYTLE
      CHARACTER(LEN=40) FLTOUT
      DATA BLANK(1:60) / ' '/
C
C
C     ******************************************************************
C
      IF (IFC.EQ.1) THEN
        IF (IFESR.EQ.1)  GO TO 501
      ENDIF
C
C ----------------------------------------------------------------------
C
C         INITIALIZATION
C
      MODEL(1:18) = 'Quadratic         '
      LYNE(1:14) = 'Line Search   '
C
      IERLIN = 0
      IF(QPRIME.GT.ZERO.AND.C0TC0.LT.ZERO) THEN
        IERLIN = -1
        GO TO 130
      ENDIF
C
      IF(ALFBND.LE.ZERO) THEN
        IERLIN = -2
        GO TO 130
      ENDIF
C
C         COMPUTE FUNCTION ERROR REDUCTION FACTOR
C
      INTHAV = 0
      RMAXFE = DBLE(MAXFER)
      REDFAC = POINT5*(ONE+RMAXFE)*RMAXFE
      REDFAC = EXP(LOG(ALFA/ZEROMN)/REDFAC)
C
C ----------------------------------------------------------------------
C
C         BEGIN ITERATION
C
 110  CONTINUE
C
C         FUNCTION EVALUATION SEQUENCE
C
C
C         PRINT STEP LENGTH
C
      TYTLE(1:32) = MODEL(1:18)//LYNE(1:14)
      CALL HHCPRS(TYTLE,' ','.')
      IF(IOFLAG.GE.10) WRITE(IPU,1003) TYTLE(1:28),IT,ALFA
      IFC = 1
      IFESR = 1
C
C         RETURN FOR FUNCTION (AND GRADIENT) INFORMATION
C
      RETURN
 501  CONTINUE
C
C         THE FUNCTION HAS BEEN EVALUATED
C
      IFC = 0
C
C ----------------------------------------------------------------------
C
C         IF FUNCTION COULD NOT BE EVALUATED BYPASS FURTHER STEPS
C
      IF(IFUNER.EQ.1) THEN
C
C         >>> FUNCTION ERROR RESPONSE
C
        IF(IOFLAG.GE.10) WRITE(IPU,1001)
C
C
        INTHAV = INTHAV + 1
        MODEL(1:16) = 'Backtrack         '
        ALFA = ALFA/(REDFAC**INTHAV)
        IF(INTHAV.LT.MAXFER) THEN
          GO TO 110
        ELSE
          GO TO 120
        ENDIF
C 
C         >>> END OF FUNCTION ERROR RESPONSE
C
      ENDIF
C
C         ITERATION PRINT
C
 120  CONTINUE
C
      IT = IT + 1
C
      IF(IFUNER.NE.1) THEN
        FLTOUT(1:40) = '........................... (B,|c|) = ('
        IF(NWFLTR.EQ.-1) THEN
          FLTOUT(24:27) = '<NW>'
        ELSEIF(NWFLTR.EQ.-2) THEN
          FLTOUT(24:27) = '<SE>'
        ELSEIF(NWFLTR.GT.0) THEN
          FLTOUT(24:27) = '<  >'
          WRITE(FLTOUT(25:26),'(I2)') NWFLTR
        ENDIF
        IF(IOFLAG.GE.10) WRITE(IPU,1002) FLTOUT(1:40),QBAR,CNORM
      ENDIF
C
C         NORMAL TERMINATION
C
      IF(IFUNER.EQ.0) GO TO 130
C
C         TERMINATE WHEN FILTER REJECTS POINT WITH UNEQUAL PRIMAL-DUAL SCALING
C
      IF(IFUNER.EQ.2.AND..NOT.SAMSCL) THEN
        IERLIN = 8
        GO TO 130
      ENDIF
C
C         MAXIMUM ITERATION COUNT
C
      IF(IT.GE.MAXIT) THEN
        IERLIN = 5
        GO TO 130
      ENDIF
C
C         IF MAXIMUM NUMBER OF FUNCTION ERRORS, TERMINATE.
C
      IF(INTHAV.GE.MAXFER) THEN
        IERLIN = 3
        GO TO 130
      ENDIF
C
C ----------------------------------------------------------------------
C
C
C         CONSTRUCT QUADRATIC LOG-BARRIER MODEL
C
      ALFBIG = MIN(ALFBND,BIGBND)
      FACOEF = QZERO + PENMU*LOG(ALFBIG)
      FBCOEF = QPRIME - PENMU/ALFBIG
      IF(ZERO.LT.ALFA.AND.ALFA.LT.ALFBIG) THEN
        FCCOEF = QBAR - FACOEF - FBCOEF*ALFA + PENMU*LOG(ALFBIG-ALFA)
        FCCOEF = FCCOEF/ALFA**2
      ELSE
        FCCOEF = ZERO
      ENDIF
C
C         QUADRATIC LOG-BARRIER MINIMIZATION STEP ESTIMATE
C
      ALFALN = ZERO
      ALFAQD = ZERO
C
      IF(QPRIME.LT.ZERO) THEN
C
        IF(FCCOEF.EQ.ZERO) THEN
C
C         LINEAR LOG-BARRIER ESTIMATE
C
          ALFALN = ALFBIG + PENMU/FBCOEF
C
        ELSE
C
C         QUADRATIC LOG-BARRIER ESTIMATE
C
          ARG = (FBCOEF + TWO*FCCOEF*ALFBIG)**2 + EIGHT*PENMU*FCCOEF
          ALFAQD = (TWO*FCCOEF*ALFBIG - FBCOEF - SQRT(ARG))
     $              /(FOUR*FCCOEF)
C
        ENDIF
C
      ENDIF
C
C         COMPUTE LOG-BARRIER STEP
C
      ALFALB = MAX(ALFALN,ALFAQD)
C
C ----------------------------------------------------------------------
C
C         CONSTRAINT ERROR MODEL
C
      EACOEF = C0TC0/TWO
      EBCOEF = -C0TC0
      ECCOEF = (CBTCB/TWO - EACOEF - EBCOEF*ALFA)/(ALFA**2)
C
C         CONSTRAINT ERROR MINIMIZATION STEP ESTIMATE
C
      IF(ECCOEF.NE.ZERO) THEN
        ALFACN = -EBCOEF/(TWO*ECCOEF)
      ELSE
        ALFACN = ZERO
      ENDIF
C
C ----------------------------------------------------------------------
C
C
C         SELECT TRIAL STEP ESTIMATE
C
      IF(ZERO.LT.ALFACN.AND.ALFACN.LT..99D0*ALFA) THEN
        IF(ZERO.LT.ALFALB.AND.ALFALB.LT..99D0*ALFA) THEN
          ALFTRL = MAX(ALFALB,ALFACN)
        ELSE
          ALFTRL = ALFACN
        ENDIF
      ELSEIF(ZERO.LT.ALFALB.AND.ALFALB.LT..99D0*ALFA) THEN
        ALFTRL = ALFALB
      ELSE
        ALFTRL = ALFA/TWO
      ENDIF
C
C         LIMIT CONTRACTION OF TRIAL STEP
C
      ALFCON = ONEEM1*ALFA
      ALFA = MAX(ALFTRL,ALFCON)
C
      IF(ALFA.EQ.ALFAQD) THEN
        MODEL(1:18) = 'Quadratic Barrier '
      ELSEIF(ALFA.EQ.ALFALN) THEN
        MODEL(1:18) = 'Linear Barrier    '
      ELSEIF(ALFA.EQ.ALFACN) THEN
        MODEL(1:18) = 'Constraint Error  '
      ELSEIF(ALFA.EQ.ALFCON) THEN
        MODEL(1:18) = 'Contraction       '
      ELSEIF(ALFA.EQ.ALFTRL) THEN
        MODEL(1:18) = 'Filter            '
      ENDIF
C
      GO TO 110
C
C ----------------------------------------------------------------------
C
C         TERMINATION PROCEDURES
C
 130  CONTINUE
C
 1001 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.........................Functi
     $on Error',T106,'*')
 1002 FORMAT(T3,'*',T106,'*',/T3,'*',T15,A40,T55,SP,1PE14.6,',',
     $   SP,1PE14.6,' )',T106,'*')
 1003 FORMAT(T3,'*',T106,'*',/T3,'*',T15,A28,'Step:',I3,'.....ALFA =',
     $   1PG16.8,T106,'*')
      RETURN
      END
