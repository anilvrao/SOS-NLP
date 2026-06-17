
      SUBROUTINE LYNSRC(CLWR,CVEC,CBAR,DF0,D2F0,F,F0,CONTOL,XMU,
     $    ITERM,MEQUAL,MAXIMP,MAXIT,MAXFER,IFUNER,IFC,IOFLAG,IPU)
C
C ======================================================================
C     LYNSRC===>lynsrc   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C
C         PURPOSE:  COMPUTE AN ESTIMATE XMU, OF THE MINIMUM OF THE
C                   FUNCTION F(XMU) = (C(XMU)-CLWR)**T)*(C(XMU)-CLWR)
C
C         INPUT:
C
C            CLWR   CONSTRAINT LOWER BOUND (MEQUAL)
C            CVEC   CONSTRAINTS AT XMU = 0 (MEQUAL)
C            CBAR   CONSTRAINTS AT XMU (MEQUAL)
C            DF0    THE SLOPE OF F AT XMU=0, I.E. DF(XMU=0)/DXMU
C            D2F0   THE CURVATURE OF F AT XMU=0, I.E. D2F(XMU=0)/DXMU2
C            F      THE VALUE OF THE FUNCTION AT XMU.
C            F0     THE VALUE OF F AT THE INITIAL POINT XMU=0
C            CONTOL CONSTRAINT ERROR TOLERANCE
C            XMU    AN ESTIMATE (GT.0) OF THE OPTIMUM VALUE
C            ITERM  USER ABORT WHEN ITERM.LT.0
C            MEQUAL NUMBER OF EQUALITY CONSTRAINTS
C            MAXIMP MAXIMUM NUMBER OF IMPROVING STEPS BEFORE TERMINATION
C            MAXIT  MAXIMUM NUMBER OF ITERATIONS
C            MAXFER MAXIMUM NUMBER OF FUNCTION ERRORS
C            IFUNER FUNCTION ERROR FLAG
C                   = 0    FUNCTION EVALUATED
C                   = 1    FUNCTION EVALUTION IMPOSSIBLE
C            IOFLAG OUTPUT CONTROL FLAG
C                   =0   NO OUTPUT
C                   =10  NORMAL OUTPUT
C            IPU    OUTPUT UNIT NO.
C
C         OUTPUT:
C
C            F      THE VALUE OF F AT THE MINIMUM POINT
C            XMU    THE OPTIMUM POINT
C            IFC    INTEGER VARIABLE:  SET = 1 WHEN RETURNING TO THE
C                   CALLING ROUTINE TO EVALUATE F(XMU)
C
C         INTERNAL: QUANTITIES ARE SAVED FROM ITERATION TO ITERATION
C                   USING THE "SAVE" STATEMENT.  
C
      DIMENSION FCLOSE(3),FVEC(4),RCLOSE(3),XMUVEC(4)
C
      LOGICAL  INTRPL,KLUSTR
C
      CHARACTER(LEN=10)  TYTLE(3)
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,
     $    ONEEP1=1.0D1,ONEEP2=1.0D2,ONEEP3=1.0D3,ONEEP5=1.0D5,
     $    ONEEM1=1.0D-1,ONEEM2=1.0D-2,ONEEM3=1.0D-3,ONEEM5=1.0D-5,
     $    POINT2=2.0D-1,POINT5=5.0D-1,POINT8=8.0D-1,POINT9=9.0D-1,
     $    THREE=3.0D0)
C
      DIMENSION CLWR(MEQUAL),CVEC(MEQUAL),CBAR(MEQUAL)
C
      LOGICAL LTEST
C
C     ******************************************************************
C
      IF (IFC.EQ.1) THEN
        IF (IFESR.EQ.1) THEN
          GO TO 501
        ELSEIF (IFESR.EQ.2) THEN
          GO TO 502
        ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C         INITIALIZATION
C
      IT = 1
      ITERM = 0
      NIMPST = 0
      INTHAV = 0
      INTRPL = .TRUE.
      KLUSTR = .FALSE.
      ITYTL1 = 1
      TYTLE(1) =  'Bracket   '
      TYTLE(2) =  'Bisection '
      TYTLE(3) =  'Cluster   '
      FVEC(1:4) = ZERO
      XMUVEC(1:4) = ZERO
C
      FVEC(1) = F0
      FVEC(4) = FVEC(1)
      FOLD = FVEC(1)
      F = FVEC(1)
      FZERO = FVEC(1)
      FREDUC = ZERO
      XMULD = ZERO
      BNDMAX = BIGNUM
      BNDMIN = -BIGNUM
C
C         COMPUTE FUNCTION ERROR REDUCTION FACTOR
C
      RMAXFE = DBLE(MAXFER)
      REDFAC = POINT5*(ONE+RMAXFE)*RMAXFE
      REDFAC = EXP(LOG(XMU/ZEROMN)/REDFAC)
C
C ----------------------------------------------------------------------
C
C         BEGIN ITERATION
C
 110  CONTINUE
C
C         IF IT = 1 USE INPUT STEP XMU AND EVALUATE FUNCTION
C
      IF(IT.EQ.1) GO TO 130
C
C         ITERATION PRINT
C
      IF(IFUNER.EQ.0) THEN
        FREDUC = SQRT(F*CONTOL) - SQRT(F0*CONTOL) 
        FCHNG = F - FOLD
        IF(IOFLAG.GE.10) WRITE(IPU,1002) FREDUC
      ENDIF
C
C         CONVERGENCE TESTS
C
C         MAXIMUM ITERATION COUNT, OR USER ABORT
C
      IF(IT.GE.MAXIT) GO TO 230
      IF(ITERM.LT.0) GO TO 230
C
C         CONSTRAINT ERROR TOLERANCE TEST
C
      IF(IFUNER.EQ.0) THEN
        LTEST = F.LT.(FMIN + DBLE(MEQUAL)*CONTOL)
     $    .AND.D2F0.GT.ZEROMN.AND.IT.GT.2
        CBARMX = DAMAX(MEQUAL,CBAR,1)
        IF(LTEST.OR.CBARMX.LT.CONTOL) GO TO 230
      ENDIF
C
C         IMPROVED POINT
C
      IF (FCHNG.LE.ZEROMN)
     $  THEN
          FOLD = F
          IF(INTRPL) NIMPST = NIMPST + 1
          IF(NIMPST.GE.MAXIMP) GO TO 230
      ENDIF
C
C
C ----------------------------------------------------------------------
C
C         ESTIMATE STEP LENGTH XMU.
C
      XMULD = XMU
      XMU = ZERO
      IF(ABS(D2F0).GT.ZEROMN) XMU = -DF0/D2F0
C
C         LIMIT CONTRACTION OF STEP
C
      SMALST = ABS(XMULD)/ONEEP3
      IF (XMU.GT.ZERO)
     $  THEN
          XMU = MAX(SMALST,XMU)
        ELSEIF (XMU.LT.ZERO)  THEN
          XMU = MIN(-SMALST,XMU)
      ENDIF
C
C ----------------------------------------------------------------------
C
C         COMPUTE STEP EXPANSION FACTOR
C
      STPMXX = TWO*MAX(ABS(XMUVEC(1)),ABS(XMUVEC(2)),
     $                 ABS(XMUVEC(3)),ABS(XMULD))
C
C ----------------------------------------------------------------------
C
      STPMAX = MIN(STPMXX,BNDMAX) 
      STPMIN = MAX(-STPMXX,BNDMIN) 
      INTRPL = .TRUE.
      IF (D2F0.GT.ZEROMN)
     $  THEN
C
C             QUADRATIC IS POSITIVE DEFINITE.  CHECK FOR STEP LIMITS.
C
          IF (XMU.GT.STPMAX)
     $      THEN
              XMU = STPMAX
              INTRPL = .FALSE.
            ELSEIF (XMU.LT.STPMIN) THEN
              XMU = STPMIN
              INTRPL = .FALSE.
          ENDIF
C
        ELSEIF(D2F0.LT.-ZEROMN) THEN
C
C             QUADRATIC IS NEG. DEFINITE. PLACE POINT ON BOUND.
C
          INTRPL = .FALSE.
          AC = FVEC(1) - (DF0 + POINT5*D2F0*XMUVEC(1))*XMUVEC(1)
          FLEF = AC + (DF0 + POINT5*D2F0*STPMIN)*STPMIN
          FRT = AC + (DF0 + POINT5*D2F0*STPMAX)*STPMAX
          XMU = STPMAX
          IF(FRT.GT.FLEF) XMU = STPMIN
          IF(XMU.EQ.BNDMAX.OR.XMU.EQ.BNDMIN) INTRPL = .TRUE.
C
        ELSE
C
C            SECOND DERIVATIVE IS NEAR ZERO -- USE LINEAR MODEL
C
          IF(DF0.LT.ZERO) THEN
            XMU = STPMAX
            INTRPL = .FALSE.
          ELSE
            XMU = STPMIN
            INTRPL = .FALSE.
          ENDIF
C
      ENDIF
C
C         IF FUNCTION HAS ALREADY BEEN EVALUATED AT THE NEW
C         PREDICTION PLACE POINT IN MIDDLE OF BRACKET
C
      LTEST = (ABS(XMU-XMUVEC(1)).LT.ZEROMN).OR.
     $    (ABS(XMU-XMUVEC(2)).LT.ZEROMN).OR.
     $    ((ABS(XMU-XMUVEC(3)).LT.ZEROMN).AND.IT.GE.2)
      IF (LTEST)
     $  THEN
          XMU = POINT5*(BNDMIN + BNDMAX)
          ITYTL1 = 2
C
C             IF SLOPES
C             AT BOUNDS HAVE THE SAME SIGN TERMINATE.
C
          LTEST = (DF0+D2F0*STPMAX)*(DF0+D2F0*STPMIN).GT.ZERO
          LTEST = LTEST.AND.(STPMAX*STPMIN.GT.ZERO)
          IF(LTEST) GO TO 230
      ENDIF
      GO TO 130
 120  CONTINUE
C
C            FUNCTION VALUE AT XMU IS NOT COMPUTABLE
C
      IF(IOFLAG.GE.10) WRITE(IPU,1001)
C
      RMX = MAX(XMUVEC(1),XMUVEC(2),XMUVEC(3))
      RMN = MIN(XMUVEC(1),XMUVEC(2),XMUVEC(3))
C
      INTRPL = .TRUE.
      INTHAV = INTHAV + 1
C
C         IF MAXIMUM NUMBER OF FUNCTION ERRORS, TERMINATE.
C
      IF(INTHAV.GE.MAXFER) THEN
        ITERM = -4
        GO TO 230
      ENDIF
C
C         MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C
      IF(IFUNER.EQ.-100) THEN
        ITERM = -5
        GO TO 230
      ENDIF
C
      IF (XMU.GE.ZERO)
     $  THEN
          BNDMAX = RMX + (XMU-RMX)*(ONE-ONEEM2)
          XMU = RMX + (BNDMAX-RMX)/(REDFAC**INTHAV)
        ELSE
          BNDMIN = RMN + (XMU-RMN)*(ONE-ONEEM2)
          XMU = RMN + (BNDMIN-RMN)/(REDFAC**INTHAV)
      ENDIF
C
 130  CONTINUE
C
C         IF XMU IS CLOSE TO PREVIOUS POINTS, TERMINATE.
C
      LTEST = (ABS(XMU-XMUVEC(1)).LT.ZEROMN).OR.
     $    (ABS(XMU-XMUVEC(2)).LT.ZEROMN).OR.
     $    (ABS(XMU-XMUVEC(3)).LT.ZEROMN).OR.
     $    (ABS(XMU-XMULD).LT.ZEROMN)
      IF(LTEST) GO TO 230
C
C ----------------------------------------------------------------------
C
C         FUNCTION EVALUATION SEQUENCE
C
C
C         PRINT STEP LENGTH
C
      IF (IOFLAG.GE.10) WRITE(IPU,1003) IT,XMU
C
      IFC = 1
      IFESR = 1
C
C         RETURN FOR FUNCTION INFORMATION
C
      RETURN
 501  CONTINUE
C
C         THE FUNCTION HAS BEEN EVALUATED
C
      IFC = 0
C
      IF(IFUNER.NE.0) THEN
C
C         IF FUNCTION COULD NOT BE EVALUATED BYPASS FURTHER STEPS
C
        GO TO 120
C
      ENDIF
C
C         IF FUNCTION INCREASED TOO MUCH RESPOND LIKE A FUNCTION ERROR
C
      IF((F-FOLD)*ZEROMN.GT.ONE) GO TO 120
C
C
C ----------------------------------------------------------------------
C
C             UPDATE INFORMATION
C
      IF (IT.EQ.1)
     $  THEN
C
C             SAVE STEP LENGTH AND FUNCTION VALUE
C
          FVEC(2) = F
          XMUVEC(2) = XMU
C
C             SAVE SMALLEST F
C
          FVEC(4) = FVEC(1)
          XMUVEC(4) = XMUVEC(1)
          IF (FVEC(4).GE.FVEC(2))
     $      THEN
              FVEC(4) = FVEC(2)
              XMUVEC(4) = XMUVEC(2)
          ENDIF
          GO TO 210
      ENDIF 
C
      IF (IT.EQ.2)
     $  THEN
C
C             SECOND ITERATION.  SAVE FUNCTION AND STEP
C
          IF (XMU.LE.XMUVEC(1))
     $      THEN
              XMUVEC(3) = XMUVEC(2)
              XMUVEC(2) = XMUVEC(1)
              XMUVEC(1) = XMU
              FVEC(3) = FVEC(2)
              FVEC(2) = FVEC(1)
              FVEC(1) = F
            ELSEIF (XMU.LE.XMUVEC(2)) THEN
              XMUVEC(3) = XMUVEC(2)
              XMUVEC(2) = XMU
              FVEC(3) = FVEC(2)
              FVEC(2) = F
            ELSE 
              XMUVEC(3) = XMU
              FVEC(3) = F
          ENDIF
C
C             SAVE SMALLEST F
C
          FVEC(4) = FVEC(1)
          DO I = 1,3
            IF (FVEC(4).GE.FVEC(I))
     $        THEN
                FVEC(4) = FVEC(I)
                XMUVEC(4) = XMUVEC(I)
            ENDIF
          ENDDO
          GO TO 210
      ENDIF
C
C             ITERATIONS GREATER THAN 2.
C
C
C         STORE NEW POINT IN CORRECT POSITION
C
      DO I = 1,3
        IF (XMU.LT.XMUVEC(I))
     $    THEN
            IP1 = I + 1
            DO L = 4,IP1,-1
              LM1 = L - 1
              XMUVEC(L) = XMUVEC(LM1)
              FVEC(L) = FVEC(LM1)
            ENDDO
            GO TO 170
        ENDIF
      ENDDO
      I = 4
 170  CONTINUE
      XMUVEC(I) = XMU
      FVEC(I) = F
C
C             COMPARE PREDICTED BEHAVIOR WITH ACTUAL
C             TO DETERMINE MODEL FOR NEXT STEP
C
      IF (IT.GT.3)
     $  THEN
          FPRED = FZERO + (DF0 + POINT5*D2F0*XMU)*XMU
          FPREDC = FZEROC + (DF0C + POINT5*D2F0C*XMU)*XMU
          IF(ABS(F-FPREDC).LT.ABS(F-FPRED)) KLUSTR = .NOT.KLUSTR
          IF (KLUSTR)
     $      THEN
              ITYTL1 = 3
            ELSE
              ITYTL1 = 1
          ENDIF
      ENDIF
C
C             COMPUTE SMALLEST F AND SAVE CORRESPONDING XMU.
C
      FSML = FVEC(1)
      XMUSML = XMUVEC(1)
      DO I = 2,4
        LTEST = FSML.GE.FVEC(I) 
        IF (LTEST)
     $    THEN
            FSML = FVEC(I)
            XMUSML = XMUVEC(I)
        ENDIF
      ENDDO
C
C         SAVE POINTS CLUSTERED AROUND POINT WITH SMALLEST F
C
C         BEGIN BY DISCARDING POINT FARTHEST FROM XMUSML
C
      IP = 1
      IF(ABS(XMUVEC(1)-XMUSML).LT.ABS(XMUVEC(4)-XMUSML)) IP = 0
C
C         ELIMINATE POINT ON LEFT OF INTERVAL IF IP = 1
C         ELIMINATE POINT ON RIGHT OF INTERVAL IF IP = 0
C
      DO I = 1,3
        FCLOSE(I) = FVEC(I+IP)
        RCLOSE(I) = XMUVEC(I+IP)
      ENDDO
C
C         IF FEASIBLE POINT WITH SMALLEST F IS ON RIGHT END OF 
C         INTERVAL, ELIMINATE POINT ON LEFT (SMALLEST XMU).  
C         OTHERWISE ELIMINATE POINT ON RIGHT (LARGEST XMU).
C
      IF (XMUSML.GE.XMUVEC(3))
     $  THEN
C
C         ELIMINATE LOWEST (LEFT) POINT
C
          DO I = 1,3
            XMUVEC(I) = XMUVEC(I+1)
            FVEC(I) = FVEC(I+1)
          ENDDO
      ENDIF
C
C         SAVE SMALLEST VALUES IN FVEC(4),XMUVEC(4)
C
      FVEC(4) = FSML
      XMUVEC(4) = XMUSML
C
 210  CONTINUE
C
C         DEFINE BRACKET SIZE
C
      IF (FVEC(2).GT.FVEC(1))
     $  THEN
          BNDMAX = XMUVEC(2)
        ELSE
          BNDMIN = XMUVEC(1)
          IF (IT.GT.1)
     $      THEN
              IF (FVEC(3).GT.FVEC(2))
     $          THEN
                  BNDMAX = XMUVEC(3)
                ELSE
                  BNDMIN=  XMUVEC(2)
              ENDIF
          ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C         COMPUTE QUADRATIC FIT COEFFICIENTS
C
      IF (IT.EQ.1)
     $  THEN
C
C         DEGENERATE SITUATION (ONLY ONE POINT)
C
          DF0 = ZERO
          DO J = 1,MEQUAL
            DCDMU = (CBAR(J) - CVEC(J))/XMUVEC(2)
            DF0 = DF0 + DCDMU*(CVEC(J)-CLWR(J))
          ENDDO
          DF0 = TWO*DF0/CONTOL
C
          D2F0 = TWO*(FVEC(2) - F0 - DF0*XMUVEC(2))/XMUVEC(2)**2
          FZERO = FVEC(1) - XMUVEC(1)*(DF0 + POINT5*D2F0*XMUVEC(1))
        ELSE
          DIF31 = (FVEC(3)-FVEC(1))/(XMUVEC(3)-XMUVEC(1))
          DIF21 = (FVEC(2)-FVEC(1))/(XMUVEC(2)-XMUVEC(1))
          D2F0 = TWO*(DIF31-DIF21)/(XMUVEC(3)-XMUVEC(2))
C
          DF0 = DIF21 - POINT5*(XMUVEC(2)+XMUVEC(1))*D2F0
          FZERO = FVEC(1) - XMUVEC(1)*(DF0 + POINT5*D2F0*XMUVEC(1))
      ENDIF
C
C         COMPUTE PREDICTED MINIMUM VALUE
C
      IF(ABS(D2F0).GT.ZEROMN)
     $    FMIN = FZERO - POINT5*DF0**2/D2F0
      IF (IT.GT.2)
     $  THEN
C
C         REPEAT QUADRATIC CALCULATIONS FOR CLUSTERED POINTS
C
          DIF31 = (FCLOSE(3)-FCLOSE(1))/(RCLOSE(3)-RCLOSE(1))
          DIF21 = (FCLOSE(2)-FCLOSE(1))/(RCLOSE(2)-RCLOSE(1))
          D2F0C = TWO*(DIF31-DIF21)/(RCLOSE(3)-RCLOSE(2))
C
          DF0C = DIF21 - POINT5*(RCLOSE(2)+RCLOSE(1))*D2F0C
          FZEROC = FCLOSE(1) - RCLOSE(1)*(DF0C + POINT5*D2F0C*RCLOSE(1))
C
C         COMPUTE PREDICTED MINIMUM VALUE
C
          IF (ABS(D2F0).GT.ZEROMN) FMINC = FZEROC - POINT5*DF0C**2/D2F0C
C
C         IF CLUSTERING POINTS IS PREFERED MODEL SWAP COEFFICIENTS
C         OTHERWISE CONTINUE.
C
          IF (KLUSTR)
     $      THEN
              CALL DSWAP(1,D2F0,1,D2F0C,1)
              CALL DSWAP(1,DF0,1,DF0C,1)
              CALL DSWAP(1,FZERO,1,FZEROC,1)
              CALL DSWAP(1,FMIN,1,FMINC,1)
          ENDIF
      ENDIF
C
      IT = IT + 1
C
C ----------------------------------------------------------------------
C
C         END OF ITERATION
C
      GO TO 110
C
 230  CONTINUE
C
C         SET ITERM FOR LINEAR FUNCTIONS PROVIDED
C         IT HAS NOT ALREADY BEEN SET
C
      IF (ABS(D2F0).LT.ONEEP3*ZEROMN*ABS(DF0)
     $  .AND.ITERM.NE.-2) THEN
          ITERM = -1
          IF (IOFLAG.GE.10) WRITE(IPU,1004)
      ENDIF
C 
      IF (IFUNER.EQ.0) THEN
        LTEST = (F.EQ.FVEC(4)).AND.INTHAV.LT.MAXFER 
        IF (LTEST) GO TO 240
      ENDIF
C
C         REEVALUATE F AT SMALLEST VALUE.
C
      XMU = XMUVEC(4)
C
C         PRINT STEP LENGTH
C
      IF (IOFLAG.GE.10) THEN
        IF (IOFLAG.GE.20) WRITE(IPU,1005)
        WRITE(IPU,1003) IT,XMU
      ENDIF
      IFC = 1
      IFESR = 2
      RETURN
C
 502  CONTINUE
      IFC = 0
      FREDUC = SQRT(F*CONTOL) - SQRT(F0*CONTOL)
      IF(IOFLAG.GE.10) WRITE(IPU,1002) FREDUC
 240  CONTINUE
C
 1001 FORMAT(T3,'*',T16,'.................................Function Error
     $',T106,'*')
 1002 FORMAT(T3,'*',T16,'.................................Change =',
     $    1PG16.8,T106,'*')
 1003 FORMAT(T3,'*',T16,'Error Reduction Step',I3,'.....XMU =',1PG16.8,
     $    T106,'*')
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T16,'CONSTRAINT ERROR IS LINEAR',
     $    T106,'*')
 1005 FORMAT(T3,'*',T106,'*'/T3,'*',T16,'REEVALUATE CONSTRAINT ERROR AT
     $ BEST POINT',T106,'*')
      RETURN
      END
