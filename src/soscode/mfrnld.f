
      SUBROUTINE MFRNLD(NROW,NCOL,AMAT,IROW,JCST,NONZA,GVEC,
     $    CVEC,IVSTAT,IFIXVR,SVEC,PVEC,VECNU,RHSWRK,RWORK,
     $    LNRWRK,IWORK,NIWORK,NEEDED,IAROW,NCALL,IPU,IPC,CNDNUM,
     $    IER)
C
C ======================================================================
C     MFRNLD===>mfrnld   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C          PURPOSE:  LEAST DISTANCE PROBLEM -- MINIMIZE || X || 
C                    SUBJECT TO THE LINEAR EQUALITIES AX = B,
C                    AND THE SPECIAL VARIABLE FIXED CONSTRAINTS
C                    X[IFIXVR] = B[IFIXVR].
C                    THIS PROBLEM IS POSED AS THE SOLUTION TO
C                        
C                        | I   A**T ||X | = |0|
C                        | A    0   ||-P| = |B|
C  
C                    A SQUARE SYSTEM (NROW=NCOL) CAN BE SOLVED WITH
C 
C                        | 0   A**T ||X | = |0|
C                        | A    0   ||-P| = |B|
C                    
C                    THE PROJECTED GRADIENT X CAN BE COMPUTED FROM 
C                        
C                        | I   A**T ||X | = |G|
C                        | A    0   ||-P| = |C|
C  
C                    THE MATRIX A IS CONSTRUCTED FROM A 
C                    SUBSET OF THE ROWS OF THE INPUT MATRIX AMAT.
C                    THE ROWS ARE DEFINED BY THE VECTOR IAROW 
C                    ELEMENTS OF IVSTAT SPECIFY THE FREE
C                    AND FIXED ROWS (SEE INPUT). VARIABLES ARE 
C                    FIXED BY AUGMENTING AMAT WITH "ONE" IN THE 
C                    APPROPRIATE COLUMN.
C
C         INPUT:    
C
C             NROW   NUMBER OF ROWS IN AMAT
C             NCOL   NUMBER OF COLUMNS IN AMAT
C             AMAT   NONZERO ELEMENTS OF LINEAR SYSTEM (NONZA)
C             IROW   INTEGER ROW INDEX VECTOR (NONZA)
C             JCST   INTEGER COLUMN START VECTOR (NCOL+1)
C                    NONZA = JCST(NCOL+1)-1
C             GVEC   RIGHT HAND SIDE (GRADIENT) VECTOR (NCOL)
C             CVEC   RIGHT HAND SIDE (CONSTRAINT) VECTOR (NROW)
C             IVSTAT INTEGER VARIABLE STATUS INDICATOR (NROW)
C                    THE NROW ELEMENTS CORRESPOND TO CONSTRAINT 
C                    SLACKS.
C                    = -1 --- EQUALITY CONSTRAINT
C                    =  0 --- FREE VARIABLE
C                    =  1 --- FIXED VARIABLE AT LOWER BOUND
C                    =  2 --- FIXED VARIABLE AT UPPER BOUND
C             IFIXVR INTEGER FIXED VARIABLE INDICATOR (NCOL)
C                    =  1 --- FIXED VARIABLE
C                    =  0 --- FREE VARIABLE
C             RHSWRK RIGHT HAND SIDE WORK ARRAY (4*NCOL)
C             RWORK  WORK ARRAY (LNRWRK)
C             LNRWRK LENGTH OF RWORK ARRAY 
C             IWORK  INTEGER WORK ARRAY (NIWORK)
C             NIWORK LENGTH OF IWORK ( > NEQNS )
C             IAROW  INTEGER WORK ARRAY (NROW)
C             NCALL  NUMBER OF THE CALL 
C                    = 1 FIRST CALL,  DO THE FOLLOWING STEPS:
C                        (1) INPUT MATRIX STRUCTURE
C                        (2) ORDER THE MATRIX
C                        (3) SYMBOLIC FACTORIZATION
C                        (4) INPUT MATRIX VALUES
C                        (5) NUMERICAL FACTORIZATION
C                        (6) SOLVE
C                    = 2 SECOND CALL  DO STEPS 4-6
C                    = 3 THIRD CALL   DO STEP  6
C             IPU    OUTPUT UNIT NO.
C             IPC    OUTPUT CONTROL FLAG
C
C
C         OUTPUT:
C
C             SVEC   SEARCH DIRECTION PART OF SOLUTION TO THE 
C                    LINEAR SYSTEM (NCOL)
C             PVEC   LAGRANGE MULTIPLIER -- PART OF SOLUTION TO 
C                    THE LINEAR SYSTEM -- NOTE SIGN CONVENTION, 
C                    AND SYMBOL PVEC TO DISTIGUISH SECOND ORDER 
C                    ESTIMATE FROM THE FIRST ORDER MULTIPLIER 
C                    XLAM (NROW)
C             VECNU  LAGRANGE MULTIPLIERS FOR FIXED VARIABLES (NCOL)
C             RWORK  WORK ARRAY CONTAINING MULTIFRONTAL FACTORIZATION
C                    INFORMATION FOR SUBSEQUENT USE STORED AS FOLLOWS
C                    RWORK(1) -- SYMBOLIC FACTORIZATION DATA (LNSYMB)
C                    RWORK(LNSYMB+1) -- NUMERIC FACTORIZATION DATA (LNRWRK)
C             CNDNUM CONDITION NUMBER OF MATRIX
C             IER    ERROR RETURN CODE
C             NEEDED STORAGE REQUIRED WHEN LNRWRK OR NIWORK IS TOO SMALL
C
      INCLUDE '../commons/NLPSPR.CMN'
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      DIMENSION AMAT(NONZA),IROW(NONZA),JCST(NCOL+1),GVEC(NCOL),
     $    CVEC(NROW),SVEC(NCOL),PVEC(NROW),RHSWRK(4*NCOL),RWORK(LNRWRK),
     $    IWORK(NIWORK),IAROW(NROW),IVSTAT(NROW),IFIXVR(NCOL),
     $    VECNU(NCOL)
      DIMENSION INRTIA(3),TYME(6),OPCNTS(2)
      LOGICAL SQUARE
      INTEGER LNZLTA, XDSLNI
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,ONEEM2=1.0D-2,ONEEM3=1.0D-3)
C
C         CHECK CALL NUMBER
C
      IER = 0
      NEEDED = 0
C
C         BEGIN FACTORIZATION TIMING CLOCK
C
      CALL CLKBEG(6)
      CALL CLKBEG(13)
C
      IF(NCALL.LT.1.OR.NCALL.GT.3) THEN
        IER = -61
        GO TO 210
      ENDIF
C
C         CHECK ROW DIMENSION
C
      IF(NROW.LT.0) THEN
        IER = -62
        GO TO 210
      ENDIF
C
      IF(NCOL.LE.0) THEN
        IER = -63
        GO TO 210
      ENDIF
C
C         COMPUTE NUMBER OF ACTIVE CONSTRAINTS 
C
      MACTIV = 0
      ITOTL  = 0
      DO I=1,NROW
        IF (IVSTAT(I).NE.0)  MACTIV = MACTIV + 1
        IF (IVSTAT(I).LT.-1 .OR. IVSTAT(I).GT.2)  ITOTL = ITOTL + 1
      ENDDO
C
C         COMPUTE NUMBER OF FIXED VARIABLES
C
      NFIXVR = 0
      DO I=1,NCOL
        IF (IFIXVR(I).NE.0)  NFIXVR = NFIXVR + 1
        IF (IFIXVR(I).LT.0 .OR. IFIXVR(I).GT.1)  ITOTL = ITOTL + 1
      ENDDO
C
C         COMPUTE TOTAL NUMBER OF CONSTRAINTS
C
      MTOTAL = MACTIV + NFIXVR
C
C         CHECK NUMBER OF CONSTRAINTS
C
      IF(MTOTAL.GT.NCOL) THEN
        IER = -635
        GO TO 210
      ENDIF
C
C         CHECK FOR INPUT ERROR IN IVSTAT AND IFIXVR
C
      IF(ITOTL.GT.0) THEN
        IER = -64
        GO TO 210
      ENDIF
C
C         DEFINE NUMBER OF EQUATIONS IN AUGMENTED SYSTEM
C
      NEQNS = MTOTAL + NCOL 
C       
C     ------------------------------------------------------------
C        PREDICTED MEMORY REQUIREMENTS FOR SYMBOLIC PHASES OF 
C        BCSLIB-EXT
C
C        REQUIREMENT BELOW SHOULD ALLOW SYMBOLIC PHASES TO PROCEED
C        "IN-CORE". 
C     ------------------------------------------------------------
     
      LNZLTA = XDSLNI (NONZA + NCOL)
C
C              MINIMUM REAL WORK REQUIRED
C
      NWRKCK = 200 + 2*LNZLTA + 2*NEQNS
      NWRKCK = MAX(NWRKCK, 200+10*(NEQNS+1))
C
C         CHECK SIZE OF RWORK
C
      IF(NWRKCK.GT.LNRWRK) THEN
        NEEDED = NWRKCK
        IER = -65
        GO TO 210
      ENDIF
      IF(NEQNS.GT.NIWORK) THEN
        NEEDED = NEQNS
        IER = -66
        GO TO 210
      ENDIF
C
C         DEFINE POINTERS FOR THE INTEGER WORK ARRAY
C
C         IWORK(LCTEMP) --- MFRNLD INTEGER TEMPORARY, (NEQNS)
C
      LCTEMP = 1
C
C         SET SQUARE SYSTEM FLAG
C
      IF(NCOL.EQ.MTOTAL) THEN
        SQUARE = .TRUE.
      ELSE
        SQUARE = .FALSE.
      ENDIF
C
C         CONSTRUCT ACTIVE ROW NUMBER IN IAROW 
C
      IAR = 0
      DO I = 1,NROW
        IF(IVSTAT(I).NE.0) THEN
          IAR = IAR + 1
          IAROW(I) = IAR
        ELSE
          IAROW(I) = -1
        ENDIF
      ENDDO
C
      IF (NCALL.EQ.2) THEN
        GO TO 160
      ELSEIF (NCALL.EQ.3) THEN
        GO TO 190
      ENDIF
C
 120  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 1, 2, AND 3-----------------------------------------
C ------------------------------------------------------------------
C
C
C         INPUT THE MATRIX STRUCTURE TO THE MULTIFRONTAL CODE
C
      IERMFR = 0
C
C         -- SET IN-CORE OPTION FLAG: O (IN-CORE); 1 (OUT-OF-CORE)
C
      MFRSTR = 0
C
C         -- SET MULTIFRONTAL I/O UNITS
C
      CALL XDSLIN ( NEQNS, 'SYMMETRIC', IOFMFR, IPU, IPUMF6, IPUMF7, 
     $              IPUMF3, IPUMF4, 0, IPUMF5, RWORK, LNRWRK,IERMFR)
C
      IF(IERMFR.NE.0) GO TO 300
C
C          SET PANEL SIZE TO NEQNS UNTIL A BETTER SCHEME IS FOUND
C
      CALL XDSLSP ( 'panel size', NEQNS, DUMMY, RWORK,IERMFR)
      IF(IERMFR.NE.0) GO TO 300
C
      CALL XDSLSP ( 'pivot tolerance', 0, TOLPVT, RWORK,IERMFR)
C 
      IF (IERMFR.NE. 0 ) GO TO 300
C
      CALL XDSLSP ( 'limit fill', 0, TOLFIL, RWORK,IERMFR)
C 
      IF (IERMFR.NE. 0 ) GO TO 300
C
      CALL XDSLSP ( 'save original matrix', 0, ZERO, RWORK,IERMFR)
C 
      IF (IERMFR.NE. 0 ) GO TO 300
C
C         -- FIRST NCOL COLUMNS
C
      NFIX = 0
C
      DO K = 1,NCOL
C
        NZCOL = 0
C
C         ------ HESSIAN BLOCK COLUMNS
C
        IF(.NOT.SQUARE) THEN
C
          NZCOL = 1
          IWORK(LCTEMP) = K
C
        ENDIF
C
C         ------ JACOBIAN BLOCK COLUMNS
C
C         LOAD COLUMN K OF AMAT
C
        DO I =  JCST(K),JCST(K+1)-1
C
          IR = IROW(I) 
C
C         IF CONSTRAINT IR IS ACTIVE LOAD ROW 
C
          IF(IAROW(IR).GT.0) THEN
            NZCOL = NZCOL + 1
            IWORK(LCTEMP+NZCOL-1) = NCOL + IAROW(IR)
          ENDIF
C
        ENDDO
C
C         ------ FIXED VARIABLE BLOCK
C
        IF(IFIXVR(K).NE.0) THEN
          NZCOL = NZCOL + 1
          NFIX = NFIX + 1
          IWORK(LCTEMP+NZCOL-1) = NCOL + MACTIV + NFIX
        ENDIF
C
        CALL XDSLIC ('A', K, NZCOL, IWORK(LCTEMP), RWORK, LNRWRK,IERMFR)
C
        IF(IERMFR.NE.0) GO TO 300
C
      ENDDO
C
C         -- LAST MTOTAL COLUMNS
C
      DO K = 1,MTOTAL
C
        NZCOL = 0
C
C         DUMMY ROW INDEX FOR NONZERO IN COLUMN JCOL
C
        IWORK(LCTEMP) = NCOL + K
C
        CALL XDSLIC ( 'A', NCOL+K, NZCOL, IWORK(LCTEMP), RWORK, 
     $                 LNRWRK,IERMFR)
C
        IF(IERMFR.NE.0) GO TO 300
C
      ENDDO
C
      CALL XDSLIF ( RWORK,  LNRWRK, NEEDS,IERMFR)
      IF (IERMFR.NE. 0 ) GO TO 300
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(NEEDS.GT.LNRWRK) THEN
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
        NEEDED = NEEDS
        IER = -661
        GO TO 210
      ENDIF
C
C         ORDER THE MATRIX
C
      CALL XDSLOR ( RWORK, LNRWRK, NEEDS,IERMFR)
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(NEEDS.GT.LNRWRK) THEN
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
        NEEDED = NEEDS
        IER = -67
        GO TO 210
      ENDIF
C
      IF(IERMFR.NE.0) GO TO 300
C
C         PERFORM SYMBOLIC FACTORIZATION
C
      CALL XDSLSF ( RWORK, LNRWRK, NEEDS, NEEDMN,IERMFR)
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(MFRSTR.EQ.0) THEN
        NEEDST = NEEDS
      ELSE
        NEEDST = NEEDMN
      ENDIF
C
      IF(NEEDST.GT.LNRWRK) THEN
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDST,LNRWRK
        NEEDED = NEEDST
        IER = -68
        GO TO 210
      ENDIF
C
C         GET MULTIFRONTAL STATISTICS -- SAVE LENGTH NEEDED
C         FOR SYMBOLIC FACTORIZATION
C
      CALL XDSLSR ( RWORK, LNSYMB, MXUSED, TYME, OPCNTS )
C
      IF(IERMFR.NE.0) GO TO 300
C
 160  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 4 AND 5---------------------------------------------
C ------------------------------------------------------------------
C
C         INPUT MATRIX VALUES TO MULTIFRONTAL SOLVER
C
C         -- FIRST NCOL COLUMNS
C
      NFIX = 0
C
      DO K = 1,NCOL
C
        NZCOL = 0
C
C         ------ HESSIAN BLOCK COLUMNS
C
        IF(.NOT.SQUARE) THEN
C
          NZCOL = 1
          IWORK(LCTEMP) = K 
          RHSWRK(NZCOL) = ONE
C
        ENDIF
C
C         ------ JACOBIAN BLOCK COLUMNS
C
C         LOAD COLUMN K OF AMAT
C
        DO I =  JCST(K),JCST(K+1)-1
C
          IR = IROW(I) 
C
C         IF CONSTRAINT IR IS ACTIVE LOAD ROW 
C
          IF(IAROW(IR).GT.0) THEN
            NZCOL = NZCOL + 1
            IWORK(LCTEMP+NZCOL-1) = NCOL + IAROW(IR)
            RHSWRK(NZCOL) = AMAT(I)
          ENDIF
C
        ENDDO
C
C         ------ FIXED VARIABLE BLOCK
C
        IF(IFIXVR(K).NE.0) THEN
          NZCOL = NZCOL + 1
          NFIX = NFIX + 1
          IWORK(LCTEMP+NZCOL-1) = NCOL + MACTIV + NFIX
          RHSWRK(NZCOL) = ONE
        ENDIF
C
        CALL XDSLVC ( 'A', K, NZCOL, IWORK(LCTEMP), RHSWRK, 
     $                RWORK, LNRWRK,IERMFR)
C
        IF(IERMFR.NE.0) GO TO 300
C
      ENDDO
C
C         -- LAST MTOTAL COLUMNS ARE ZERO BLOCK; NO INPUT REQUIRED
C
      CALL XDSLVF ( RWORK, LNRWRK,IERMFR)
      IF (IERMFR.NE. 0 ) GO TO 300
C
C         FACTOR THE MATRIX
C
      CALL XDSLFA ( RWORK, LNRWRK, CNDNUM, INRTIA, NEEDS, IERMFR)
C
      IF(IERMFR.EQ.-501.OR.IERMFR.EQ.-507.OR.IERMFR.EQ.-508) THEN
        IER = -65
        GO TO 210
      ENDIF
C
      IF(IERMFR.NE.0) GO TO 300
C
C         CHECK FOR ILL-CONDITIONED KT MATRIX
C
      IF(CNDNUM.GT.BIGCND) THEN
        IF(IPC.GT.0) WRITE(IPU,1004) CNDNUM
        IER = 5000
      ENDIF
C
      IF(IOFMFR.GT.0) CALL XDSLPS(RWORK)
C
C         GET MULTIFRONTAL STATISTICS (fix this for new ext release)
C
      LNIWRK = 0
      MXIUSD = RWORK(6)
      MXUSED = RWORK(14)
      IF(MXUSED.LE.LNRWRK) INSTAT(10) = MIN(INSTAT(10),LNRWRK-MXUSED)
      IF(MXIUSD.LE.LNIWRK) INSTAT(30) = MIN(INSTAT(30),LNIWRK-MXIUSD)
C
      IF(IER.NE.0) GO TO 300
C
 190  CONTINUE
C
C ------------------------------------------------------------------
C --------STEP 6----------------------------------------------------
C ------------------------------------------------------------------
C
C         DEFINE RIGHT HAND SIDE
C
      RHSWRK(1:NCOL) = GVEC(1:NCOL)
      DO I = 1,NROW
C
        IF(IAROW(I).GT.0) THEN
          RHSWRK(NCOL+IAROW(I)) = CVEC(I)
        ENDIF
C
      ENDDO
C
C         LOAD ZERO RHS FOR FIXED VARIABLES
C
      RHSWRK(NCOL+MACTIV+1:NCOL+MACTIV+NFIXVR) = ZERO
C
      LCSOLN = 2*NCOL + 1
C
C         SOLVE LINEAR SYSTEM
C
      CALL REFITR(RHSWRK,RHSWRK(LCSOLN),NEQNS,
     $    RWORK,LNRWRK,NEEDS,IERMFR)
C
      IF(IERMFR.NE.0) GO TO 300
C
C ------------------------------------------------------------------
C
C         --- STORE THE SOLUTION INTO SVEC, PVEC, AND VECNU
C
      SVEC(1:NCOL) = RHSWRK(LCSOLN:LCSOLN+NCOL-1)
      PVEC(1:MACTIV) = RHSWRK(LCSOLN+NCOL:LCSOLN+NCOL+MACTIV-1)
      IOFFIX = LCSOLN + NCOL + MACTIV
      VECNU(1:NFIXVR) = RHSWRK(IOFFIX:IOFFIX+NFIXVR-1)
C
 210  CONTINUE
C
      IF(IER.NE.0) THEN
        IF(IPC.GT.0) WRITE(IPU,1003) IER
C
C         SET WORK ARRAY STORAGE REQUIRED PLUS A SMALL PAD FOR SLOP
C
        NEEDED = INT((ONE + TOLPVT)*DBLE(NEEDS))
C
      ENDIF
C
C         END FACTORIZATION TIMING CLOCK
C
      CALL CLKSUM(6)
      CALL CLKSUM(13)
C      
      RETURN 
C     -----------------------------------------------------------
C     ... TERMINATE ON UNEXPECTED ERROR CONDITION FROM BCSLIB-EXT
C     -----------------------------------------------------------
 300  CONTINUE
      IF(IPC.GT.0.AND.IERMFR.NE.0) WRITE(IPU,1002) IERMFR
      IER = -999
      GO TO 210

 1001 FORMAT(T3,'*',T106,'*'/T3,'*',T11,
     $  'STORAGE ERROR IN LDP;  NEEDS =',I6,'  NEEDED =',I6,T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,
     $  'MULTIFRONTAL ERROR;  IERMFR =',I6,T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,
     $  'MFRNLD ERROR;  IER =',I6,T106,'*')
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T11,
     $  'ILL-CONDITIONED KT MATRIX; CNDNUM =',G16.6,T106,'*')

      END
