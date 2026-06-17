      SUBROUTINE QPSKT (NROW,NCOL,AMAT,IROW,JCOLST,NONZA,IDHESS
     $    ,HMAT,IROWH,JSTRH,NONZH,SFEZ,GFEZ,GVEC,CVEC,IVSTAT
     $    ,PVEC,PIVEC,RTHSKT,SOLNKT,RWORK,LNRWRK,NEEDED
     $    ,IFREE,ITEMP,NCALL,IPU,IPC,CNDNUM,IERSKT,LNSYMB,MEQUAL
     $    ,MINEQL,NSLK,NFEZ,NVAR,NEQNS)
C
C ======================================================================
C     QPSKT ===>qpskt    J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C         PURPOSE:  SOLVE THE SET OF EQUATIONS RESULTING FROM 
C                   LINEARIZATION OF THE KUHN-TUCKER CONDITIONS.
C                   THE MULTIFRONTAL SPARSE SYMMETRIC INDEFINITE
C                   SOLVER IS APPLIED TO A SYSTEM CONSTRUCTED 
C                   FROM COLUMNS (AND ROWS) OF THE SYMMETRIC SYSTEM
C
C                       | 0                   || -ps | = | 0 |
C                       | 0    1              || -pf | = | gf|
C                       | 0    0    H         || -pr | = | gr|
C                       | 0    se   Ae   0    ||(pi)e| = | ce|
C                       | I    si   Ai   0   0||(pi)i| = | ci| , 
C
C                   WHERE THE "SUBSCRIPTS" s, f AND r REFER TO SLACK,
C                   FEASIBILITY AND REAL VARIABLES, RESPECTIVELY. NOTE
C                   THAT gs = 0 (SLACKS DON'T APPEAR IN THE OBJECTIVE)
C                   AND THAT gf = rho + xi SINCE THE FEASBILITY TERM
C                   IN THE OBJECTIVE IS rho*xi + 0.5*(xi**2). (Note:
C                   pf and xi both denote the feasibility variable.)
C
C                   NOTE ALSO THAT THE ROWS AND COLS OF THE KT MATRIX
C                   CORRESPONDING TO THE FEASIBILTY STUFF DO NOT APPEAR
C                   (AND NFEZ = 0 ) WHEN xi IS FIXED ON A BOUND.
C
C                   THE SYMMETRIC MATRIX H IS (NCOL X NCOL).
C                   THE RECTANGULAR MATRIX A IS (NROW X NCOL) AND IS
C                   PARTITIONED WITH 
C                       A = |Ae|
C                           |Ai|
C                   WHERE THE FIRST MEQUAL ROWS CORRESPOND TO EQUALITIES
C                   AND THE NEXT MINEQL ROWS CORRESPOND TO INEQUALITIES,
C                   NROW = (MEQUAL + MINEQL). THE IDENTITY MATRIX I IS
C                   (MINEQL X MINEQL), AND THE SCALED FEASIBILITY RESIDUAL
C                   VECTOR s = (se,si)**T HAS NROW ELEMENTS.
C
C                   THE MATRICES Af AND Hf ARE CONSTRUCTED FROM A 
C                   SUBSET OF THE COLUMNS OF THE INPUT MATRICES AMAT
C                   AND HMAT.  THE COLUMN SUBSETS ARE DEFINED BY THE 
C                   VECTOR IFREE = (is,if,ir) WHERE is IS A VECTOR OF
C                   LENGTH MINEQL, if IS A SCALAR, AND ir IS A VECTOR OF 
C                   LENGTH NCOL.  ELEMENTS OF IVSTAT SPECIFY THE FREE
C                   AND FIXED VARIABLES (SEE INPUT).
C
C                   A MINIMUM NORM UNDERDETERMINED SYSTEM CAN BE 
C                   SOLVED WHEN NSLK = 0, g = 0, c = b, AND H = I, 
C
C                       | I   A**T || -p | = |0|
C                       | A    0   || pi | = |b|
C
C                   A SQUARE SYSTEM (NROW=NCOL) CAN BE SOLVED WITH
C
C                       | 0   A**T || -p | = |0|
C                       | A    0   || pi | = |b|
C
C                   BY SETTING IDHESS = 2
C
C         INPUT:    
C
C             NROW   NUMBER OF ROWS IN AMAT
C             NCOL   NUMBER OF COLUMNS IN AMAT
C             AMAT   NONZERO ELEMENTS OF LINEAR SYSTEM (NONZA)
C             IROW   INTEGER ROW INDEX VECTOR (NONZA)
C             JCOLST INTEGER COLUMN START VECTOR (NCOL+1)
C             NONZA  DIMENSION OF AMAT (GE. JCOLST(NCOL+1) - 1)
C             IDHESS INTEGER FLAG SPECIFYING HESSIAN INPUT AS FOLLOWS:
C                    = 0  HESSIAN IS INPUT
C                    = 1  HESSIAN IS IDENTITY MATRIX
C                    = 2  HESSIAN IS ZERO MATRIX (SQUARE SYSTEM SOLVE)
C             HMAT   NONZERO ELEMENTS OF LOWER TRIANGLE OF HMAT (NONZH)
C             IROWH  INTEGER ROW INDEX VECTOR (NONZH)
C             JSTRH  INTEGER COLUMN START VECTOR (NCOL+1)
C             NONZH  DIMENSION OF HMAT (.GE. JSTRH(NCOL+1)-1 )
C             SFEZ   SCALED FEASIBILITY RESIDUAL VECTOR (NROW)
C             GFEZ   SCALAR FEASIBILITY GRADIENT.
C             GVEC   RIGHT HAND SIDE (GRADIENT) VECTOR (NCOL)
C             CVEC   RIGHT HAND SIDE (CONSTRAINT) VECTOR (NROW)
C             IVSTAT INTEGER VARIABLE STATUS INDICATOR (NROW+1+NCOL).
C                    IVSTAT PARTITION =
C                          ( equalities, inequalities, xi, real vars )
C                    THE FIRST NROW ELEMENTS CORRESPOND TO CONSTRAINT 
C                    SLACKS, FOLLOWED BY A SCALAR FOR THE FEASIBLITY SLACK
C                    AND THE REMAINING NCOL ELEMENTS TO REAL VARIABLES.
C                    = -1 --- EQUALITY CONSTRAINT
C                    =  0 --- FREE VARIABLE
C                    =  1 --- FIXED VARIABLE AT LOWER BOUND
C                    =  2 --- FIXED VARIABLE AT UPPER BOUND
C                    =  3 --- REAL VARIABLE PERMANENTLY FIXED
C             RWORK  REAL WORK ARRAY (LNRWRK)
C             LNRWRK LENGTH OF RWORK ARRAY (SEE BELOW)
C
C                    LOWER BOUNDS FOR THE WORK ARRAYS
C                    REAL:
C                      LNRWRK > NKT + 4*NEQN + 200
C                    WHERE
C                      NKT  > MINEQL + 1 + MCON + NONZH + NONZG
C                      NEQN > MINEQL + 1 + NDIM + MCON
C
C             NEEDED STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C             IFREE  INTEGER WORK ARRAY (MINEQL+1+NCOL)
C                    IFREE PARTITION = ( slacks, feas, real )
C                    IFREE(k) = { -1  if k-th variable is fixed
C                               { j   if k-th variable is the j-th
C                                     free var among all free vars.
C             ITEMP  INTEGER WORK ARRAY (NEQNS)
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
C             PVEC   SEARCH DIRECTION PART OF SOLUTION TO THE 
C                    LINEAR SYSTEM (NSLK+NFEZ+NVAR)
C             PIVEC  LAGRANGE MULTIPLIER -- PART OF SOLUTION TO 
C                    THE LINEAR SYSTEM -- NOTE SIGN CONVENTION, 
C                    AND SYMBOL PIVEC TO DISTINGUISH SECOND ORDER 
C                    ESTIMATE FROM THE FIRST ORDER MULTIPLIER 
C                    XLAM (NROW)
C             RTHSKT RIGHT HAND SIDE OF THE KT SYSTEM (NEQNS)
C             SOLNKT SOLUTION OF THE KT SYSTEM (NEQNS)
C             RWORK  WORK ARRAY CONTAINING MULTIFRONTAL FACTORIZATION
C                    INFORMATION FOR SUBSEQUENT USE STORED AS FOLLOWS
C                    RWORK(1) -- SYMBOLIC FACTORIZATION DATA (LNSYMB)
C                    RWORK(LNSYMB+1) -- NUMERIC FACTORIZATION DATA (LNRWRK)
C             CNDNUM CONDITION NUMBER OF MATRIX
C             IERSKT ERROR RETURN CODE
C
C                    0    NORMAL RETURN
C
C             THE CONVENTION FOR ERROR RETURN NUMBERS IS AS FOLLOWS:
C             -600 > IERSKT > -700       I/O ERRORS (INSUFFICIENT DISK SPACE)
C             -700 > IERSKT > -800       ERRORS REQUIRES SCHUR-QP RESPONSE
C             -800 > IERSKT > -900   --- INPUT ERRORS
C             -900 > IERSKT          --- STORAGE ERRORS REQUIRING
C                                        SCHUR-QP RESPONSE
C             THE UNITS DIGIT IS NUMBERED CHRONOLOGICALLY WITHIN QPSKT
C
C                    -801  NCALL.LT.1 OR NCALL.GT.3
C                    -802  NROW.LT.0
C                    -803  NCOL.LE.0
C                    -804  IVSTAT IS WRONG
C                    -906  REAL WORK ARRAY TOO SMALL
C                    -907  REAL WORK ARRAY TOO SMALL (XDSLIN, IER=-101)
C                    -808  NEQNS .LT. 1 (XDSLIN, IER=-102)
C                    -809  INCORRECT PATH (XDSLIC, IER=-100)
C                    -910  REAL WORK ARRAY TOO SMALL (XDSLIC, IER=-101)
C                    -811  ILLEGAL VALUE FOR NZCOL (XDSLIC, IER=-103)
C                    -812  ILLEGAL VALUE FOR JCOL (XDSLIC, IER=-104)
C                    -813  ILLEGAL VALUE FOR JROWIN (XDSLIC, IER=-106)
C                    -613  I/O ERROR ON WAFIL1 (XDSLIC, IER=-110)
C                    -814  INCORRECT PATH (XDSLIf, IER=-100)
C                    -915  REAL WORK ARRAY TOO SMALL (NEEDS.GT.LNRWRK)
C                    -617  I/O ERROR ON WAFIL1 (XDSLIf, IER=-110)
C                    -818  INCORRECT PATH (XDSLOR, IER=-200)
C                    -919  REAL WORK ARRAY TOO SMALL (XDSLOR, IER=-201)
C                    -620  I/O ERROR ON SQFILE (XDSLOR, IER=-202)
C                    -921  INTEGER WORK ARRAY TOO SMALL (IER=-203)
C                    -822  INCORRECT PATH (XDSLSF, IER=-300)
C                    -923  INTEGER WORK ARRAY TOO SMALL (XDSLSF, IER=-301)
C                    -924  real WORK ARRAY TOO SMALL (XDSLSF, IER=-301)
C                    -925  INTERNAL ERROR DURING SYMBOLIC FACT. (XDSLSF, IER=-302)
C                    -626  I/O ERROR ON SQFIL2 (XDSLSF, IER=-303)
C                    -927  REAL WORK ARRAY TOO SMALL (NEEDST.GT.LNRWRK)
C                    -928  INTEGER WORK ARRAY TOO SMALL
C                    -829  INCORRECT PATH (XDSLVC, IER=-400)
C                    -930  REAL WORK ARRAY TOO SMALL (XDSLVC, IER=-401)
C                    -831  JROWIN(I),JCOL PAIR NOT SPECIFIED (XDSLVC, IER=-402)
C                    -632  I/O ERROR ON waFILE (XDSLVC, IER=-410)
C                    -833  INCORRECT PATH (XDSLvf, IER=-400)
C                    -934  REAL WORK ARRAY TOO SMALL (XDSLVf, IER=-401)
C                    -635  I/O ERROR ON sqfile, waFILE (XDSLVf, IER=-410)
C                    -836  INCORRECT PATH (XDSLFA, IER=-500,-502)
C                    -937  REAL WORK ARRAY TOO SMALL (XDSLFA, IER=-501,-507,-508)
C                    -731  ZERO PIVOT; SINGULAR MATRIX (XDSLFA, IER=-503)
C                    -639  I/O ERROR ON SQFIL1 or sqfil2 (XDSLFA, IER=-504)
C                    -640  I/O ERROR ON WAFIL1 (XDSLFA, IER=-505,-511,-512)
C                    -641  I/O ERROR ON WAFIL2 (XDSLFA, IER=-506)
C                    -742  EXCESSIVE FILL DUE TO PIVOTING (XDSLFA, IER=-509)
C                    -744  WRONG INERTIA OF KT MATRIX
C                    -845  INCORRECT PATH (XDSLSL, IER=-600)
C                    -946  REAL WORK ARRAY TOO SMALL FOR SOLUTION (XDSLSL, IER=-601)
C                    -847  NRHS .LT. 1 (XDSLSL, IER=-602)
C                    -848  LDRHS .LE. NEQNS (XDSLSL, IER=-603)
C                    -649  I/O ERROR ON WAFIL1 (XDSLSL, IER=-604,605)
C                    -750  MAX. ITERATION DURING ITERATIVE REFINEMENT (REFITR, IER=-701)
C                    -852  INCORRECT INPUT WHEN REQUESTING KT FILE OUTPUT
C                    -853  TERMINATE AFTER SUCCESSFUL KT FILE OUTPUT
C                    -732  PIVOTING FAILURE -- NEED TO INCREASE PANEL SIZE 
C                          (XDSLFA, IER=-513)
C
C             LNSYMB LENGTH OF RWORK ARRAY NEEDED FOR SYMBOLIC FACTORIZATION
C             MEQUAL NUMBER OF EQUALITY CONSTRAINTS
C             MINEQL NUMBER OF INEQUALITY CONSTRAINTS
C             NSLK   NUMBER OF SLACK VARIABLES
C             NFEZ   NUMBER OF FEASIBILITY VARIABLES
C             NVAR   NUMBER OF FREE REAL VARIABLES 
C             NEQNS  NUMBER OF EQUATIONS IN AUGMENTED KT SYSTEM
C             NEEDED STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C
      INCLUDE '../commons/NLPSPR.CMN'
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      COMMON /INERVL/ INREQD
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
      DIMENSION AMAT(NONZA),IROW(NONZA),JCOLST(NCOL+1),HMAT(NONZH),
     $    IROWH(NONZH),JSTRH(NCOL+1),GVEC(NCOL),CVEC(*),PVEC(*),
     $    PIVEC(*),RTHSKT(NEQNS),SOLNKT(NEQNS),RWORK(LNRWRK),
     $    IFREE(*),ITEMP(NEQNS),IVSTAT(NROW+1+NCOL),
     $    SFEZ(*)
      DIMENSION INRTIA(3),TYME(6),OPCNTS(2)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,ONEEM2=1.0D-2,ONEEM3=1.0D-3)
      PARAMETER (NITREF=0)
      INTEGER LNZLTA, XDSLNI, MAX
C
C         BEGIN K-T FACTORIZATION TIMING
C
      CALL CLKBEG(6)
      CALL CLKBEG(13)
C
C         CHECK CALL NUMBER
C
      IERSKT = 0
C
      IF(NCALL.LT.1.OR.NCALL.GT.3) THEN
        IERSKT = -801
        GO TO 280
      ENDIF
C
C         CHECK ROW AND COLUMN DIMENSIONS
C
      IF(NROW.LT.0) THEN
        IERSKT = -802
        GO TO 280
      ENDIF
C
      IF(NCOL.LE.0) THEN
        IERSKT = -803
        GO TO 280
      ENDIF
C
C         CHECK FOR INPUT ERROR IN IVSTAT 
C
      NTOTL = NROW + 1 + NCOL
      ITOTL = 0
      DO I=1,NTOTL
        IF (IVSTAT(I).GT.3 .OR. IVSTAT(I).LT.-1)   ITOTL = ITOTL + 1
      ENDDO
      IF(ITOTL.GT.0) THEN
        IERSKT = -804
        GO TO 280
      ENDIF
C
C         ESTIMATE FOR NUMBER OF NONZEROS IN LOWER TRIANGULAR PART OF KT MATRIX
C
      NZLTA = MINEQL + 1 + NROW + NONZH + NONZA
C       
C     ------------------------------------------------------------
C        PREDICTED MEMORY REQUIREMENTS FOR SYMBOLIC PHASES OF 
C        BCSLIB-EXT
C
C        REQUIREMENT BELOW SHOULD ALLOW SYMBOLIC PHASES TO PROCEED
C        "IN-CORE". 
C     ------------------------------------------------------------
     
      LNZLTA = XDSLNI ( NZLTA)

      NWRKCK = 200 + 2*LNZLTA + 2*NEQNS
      NWRKCK = MAX(NWRKCK, 200+10*(NEQNS+1))
      NEEDED = 0
C
      IF(NWRKCK.GT.LNRWRK) THEN
        NEEDED = NWRKCK
        IERSKT = -906
        GO TO 280
      ENDIF
C
C         CONSTRUCT FREE ROW NUMBER IN IFREE
C
      KK = 0
      DO I = 1,MINEQL+1+NCOL
        IF(IVSTAT(MEQUAL+I).EQ.0) THEN
          KK = KK + 1
          IFREE(I) = KK
        ELSE
          IFREE(I) = -ABS(IVSTAT(MEQUAL+I))
        ENDIF
      ENDDO
C
      IF (NCALL.EQ.2) THEN
        GO TO 190
      ELSEIF (NCALL.EQ.3) THEN
        GO TO 250
      ENDIF
C
 120  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 1, 2, AND 3-----------------------------------------
C ------------------------------------------------------------------
C
C         INPUT THE MATRIX STRUCTURE TO THE MULTIFRONTAL CODE
C
C
C         -- SET IN-CORE OPTION FLAG: O (IN-CORE); 1 (OUT-OF-CORE)
C
      MFRSTR = 0
C
      CALL XDSLIN ( NEQNS, 'SYMMETRIC', IOFMFR, IPU, IPUMF1, IPUMF2, 
     $              IPUMF3, IPUMF4, 0, IPUMF5, RWORK, LNRWRK, IER )
C
      IF(IER.EQ.-101) THEN
        IERSKT = -907
      ELSEIF(IER.EQ.-102) THEN
        IERSKT = -808
      ENDIF
      IF(IERSKT.NE.0) GO TO 280
C
C          SET PANEL SIZE TO NEQNS UNTIL A BETTER SCHEME IS FOUND
C
      CALL XDSLSP ( 'panel size', NEQNS, DUMMY, RWORK, IER )
      IF(IER.NE.0) GO TO 280
C
      CALL XDSLSP ( 'pivot tolerance', 0, TOLPVT, RWORK, IER )
C 
      IF ( IER .NE. 0 ) GO TO 280
C
      CALL XDSLSP ( 'limit fill', 0, TOLFIL, RWORK, IER )
C 
      IF ( IER .NE. 0 ) GO TO 280
C
      CALL XDSLSP ( 'save original matrix', 0, ZERO, RWORK, IER )
C 
      IF ( IER .NE. 0 ) GO TO 280
C
C         -- FIRST NSLK COLUMNS
C
      JCOL = 0
      DO K = 1,MINEQL
C
C         IF SLACK VARIABLE K IS FREE LOAD COLUMN JCOL
C
        IF(IFREE(K).GT.0) THEN
C
          JCOL = JCOL + 1
          NZCOL = 1
C
C         ROW INDEX FOR NONZERO IN COLUMN JCOL
C
          ITEMP(1) = NSLK + NFEZ + NVAR + MEQUAL + K
C
          CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
          IF(IER.EQ.-100) THEN
            IERSKT = -809
          ELSEIF(IER.EQ.-101) THEN
            IERSKT = -910
          ELSEIF(IER.EQ.-103) THEN
            IERSKT = -811
          ELSEIF(IER.EQ.-104) THEN
            IERSKT = -812
          ELSEIF(IER.EQ.-106) THEN
            IERSKT = -813
          ELSEIF(IER.EQ.-110) THEN
            IERSKT = -613
          ENDIF
          IF(IERSKT.NE.0) GO TO 280
C
        ENDIF
C
      ENDDO
C
C         -- NEXT NFEZ COLUMN
C
      IF(NFEZ.EQ.1) THEN
C
        JCOL = JCOL + 1
        NZCOL = 1 + NROW
C
C         ------ HESSIAN BLOCK
C
        ITEMP(1) = NSLK + 1
C
C         ------ FEASIBILTY VECTOR BLOCK 
C
        DO I = 1,NROW
          ITEMP(I+1) = NSLK + 1 + NVAR + I
        ENDDO
C
          CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
          IF(IER.EQ.-100) THEN
            IERSKT = -809
          ELSEIF(IER.EQ.-101) THEN
            IERSKT = -910
          ELSEIF(IER.EQ.-103) THEN
            IERSKT = -811
          ELSEIF(IER.EQ.-104) THEN
            IERSKT = -812
          ELSEIF(IER.EQ.-106) THEN
            IERSKT = -813
          ELSEIF(IER.EQ.-110) THEN
            IERSKT = -613
          ENDIF
          IF(IERSKT.NE.0) GO TO 280
C
      ENDIF
C
C         -- NEXT NVAR COLUMNS
C
      DO 170 K = MINEQL+2,MINEQL+1+NCOL
C
C         IF REAL VARIABLE K IS FREE LOAD COLUMN JCOL
C
        KK = K - MINEQL - 1 
        IF(IFREE(K).GT.0) THEN
C
          JCOL = JCOL + 1
          NZCOL = 0
C
C         ------ HESSIAN BLOCK COLUMNS
C
          IF(IDHESS.EQ.2) THEN
C
            ITEMP(1) = NSLK + NFEZ + JCOL
C
          ELSEIF(IDHESS.EQ.1) THEN
C
            NZCOL = 1
            ITEMP(1) = NSLK + NFEZ + JCOL
C
          ELSE
C
C           LOAD COLUMN K OF HMAT
C
            DO I =  JSTRH(KK),JSTRH(KK+1)-1
C
              IF(NZCOL.GT.0.AND.HMAT(I).EQ.ZERO) CYCLE
C
              IR = IROWH(I) + MINEQL + 1
C
C             IF REAL VARIABLE IR IS FREE, LOAD ROW 
C
              IF(IFREE(IR).GT.0) THEN
                NZCOL = NZCOL + 1
                ITEMP(NZCOL) = IFREE(IR)
              ENDIF
C
            ENDDO
C
          ENDIF
C
C         ------ JACOBIAN BLOCK COLUMNS
C
C         LOAD COLUMN K OF AMAT
C
          IF(NROW.GT.0) THEN
            DO I =  JCOLST(KK),JCOLST(KK+1)-1
C
              IF(AMAT(I).EQ.ZERO) CYCLE
C
              NZCOL = NZCOL + 1
              ITEMP(NZCOL) = NSLK + NFEZ + NVAR + IROW(I) 
C
            ENDDO
          ENDIF
C
          CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
          IF(IER.EQ.-100) THEN
            IERSKT = -809
          ELSEIF(IER.EQ.-101) THEN
            IERSKT = -910
          ELSEIF(IER.EQ.-103) THEN
            IERSKT = -811
          ELSEIF(IER.EQ.-104) THEN
            IERSKT = -812
          ELSEIF(IER.EQ.-106) THEN
            IERSKT = -813
          ELSEIF(IER.EQ.-110) THEN
            IERSKT = -613
          ENDIF
          IF(IERSKT.NE.0) GO TO 280
C
        ENDIF
C
 170  CONTINUE
C
C         -- LAST NROW COLUMNS
C
      DO K = 1,NROW
C
          JCOL = JCOL + 1
          NZCOL = 0
C
C         DUMMY ROW INDEX FOR NONZERO IN COLUMN JCOL
C
          ITEMP(1) = NSLK + NVAR + NFEZ + K
C
          CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
          IF(IER.EQ.-100) THEN
            IERSKT = -809
          ELSEIF(IER.EQ.-101) THEN
            IERSKT = -910
          ELSEIF(IER.EQ.-103) THEN
            IERSKT = -811
          ELSEIF(IER.EQ.-104) THEN
            IERSKT = -812
          ELSEIF(IER.EQ.-106) THEN
            IERSKT = -813
          ELSEIF(IER.EQ.-110) THEN
            IERSKT = -613
          ENDIF
          IF(IERSKT.NE.0) GO TO 280
C
      ENDDO
C
      CALL XDSLIF ( RWORK, LNRWRK, NEEDS, IER )
C
      IF(IER.EQ.-100) THEN
        IERSKT = -814
      ELSEIF(IER.EQ.-101) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
          NEEDED = NEEDS
          IERSKT = -915
        ENDIF
      ELSEIF(IER.EQ.-110) THEN
        IERSKT = -617
      ENDIF
      IF ( IERSKT .NE. 0 ) GO TO 280
C
C         ORDER THE MATRIX
C
      CALL XDSLOR ( RWORK, LNRWRK, NEEDS, IER )
C
      IF(IER.EQ.-200) THEN
        IERSKT = -818
      ELSEIF(IER.EQ.-201) THEN
        IERSKT = -919
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
        NEEDED = NEEDS
      ELSEIF(IER.EQ.-202) THEN
        IERSKT = -620
      ENDIF
      IF(IERSKT.NE.0) GO TO 280
C
C         PERFORM SYMBOLIC FACTORIZATION
C
      CALL XDSLSF ( RWORK, LNRWRK, NEEDS, NEEDMN, IER )
C
      IF(IER.EQ.-300) THEN
        IERSKT = -822
      ELSEIF(IER.EQ.-301) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IERSKT = -924
          IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
          NEEDED = NEEDS
        ENDIF
      ELSEIF(IER.EQ.-302) THEN
        IERSKT = -925
      ELSEIF(IER.EQ.-303) THEN
        IERSKT = -626
      ENDIF
      IF(IERSKT.NE.0) GO TO 280
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(MFRSTR.EQ.0) THEN
        NEEDST = NEEDS
      ELSE
        NEEDST = NEEDMN
      ENDIF
C
      NEEDST = INT((ONE + TOLPVT)*DBLE(NEEDST)) 
      IF(NEEDST.GT.LNRWRK) THEN
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDST,LNRWRK
        NEEDED = NEEDST
        IERSKT = -927
        GO TO 280
      ENDIF
C
C         GET MULTIFRONTAL STATISTICS -- SAVE LENGTH NEEDED
C         FOR SYMBOLIC FACTORIZATION
C
      CALL XDSLSR(RWORK,LNSYMB,MXUSED,TYME,OPCNTS)
C
 190  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 4 AND 5---------------------------------------------
C ------------------------------------------------------------------
C
C         INPUT MATRIX VALUES TO MULTIFRONTAL SOLVER
C
C         -- FIRST NSLK COLUMNS
C
      JCOL = 0
      DO K = 1,MINEQL
C
C         IF SLACK VARIABLE K IS FREE LOAD COLUMN JCOL
C
        IF(IFREE(K).GT.0) THEN
C
          JCOL = JCOL + 1
          NZCOL = 1
C
C         ROW INDEX FOR NONZERO IN COLUMN JCOL
C
          ITEMP(1) = NSLK + NFEZ + NVAR + MEQUAL + K
          RTHSKT(NZCOL) = ONE
C
          CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHSKT, 
     $                  RWORK, LNRWRK, IER )
C
          IF(IER.EQ.-400) THEN
            IERSKT = -829
          ELSEIF(IER.EQ.-401) THEN
            IERSKT = -930
          ELSEIF(IER.EQ.-402) THEN
            IERSKT = -831
          ELSEIF(IER.EQ.-410) THEN
            IERSKT = -632
          ENDIF
          IF(IERSKT.NE.0) GO TO 280
C
        ENDIF
C
      ENDDO
C
C         -- NEXT NFEZ COLUMN
C
      IF(NFEZ.EQ.1) THEN
C
        JCOL = JCOL + 1
        NZCOL = 1 
C
C         ------ HESSIAN BLOCK
C
        ITEMP(1) = NSLK + 1
        RTHSKT(NZCOL) = ONE
C
C         ------ FEASIBILTY VECTOR BLOCK 
C
        DO I = 1,NROW
          ITEMP(I+1) = NSLK + 1 + NVAR + I
          NZCOL = NZCOL + 1
          RTHSKT(NZCOL) = SFEZ(I)
        ENDDO
C
        CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHSKT, 
     $                 RWORK, LNRWRK, IER )
C
        IF(IER.EQ.-400) THEN
          IERSKT = -829
        ELSEIF(IER.EQ.-401) THEN
          IERSKT = -930
        ELSEIF(IER.EQ.-402) THEN
          IERSKT = -831
        ELSEIF(IER.EQ.-410) THEN
          IERSKT = -632
        ENDIF
        IF(IERSKT.NE.0) GO TO 280
C
      ENDIF
C
C         -- NEXT NVAR COLUMNS
C
      DO 240 K = MINEQL+2,MINEQL+1+NCOL
C
C         IF REAL VARIABLE K IS FREE LOAD COLUMN JCOL
C
        KK = K - MINEQL - 1 
        IF(IFREE(K).GT.0) THEN
C
          JCOL = JCOL + 1
          NZCOL = 0
C
C         ------ HESSIAN BLOCK COLUMNS
C
          IF(IDHESS.EQ.2) THEN
C
            ITEMP(1) = NSLK + NFEZ + JCOL
C
          ELSEIF(IDHESS.EQ.1) THEN
C
            NZCOL = 1
            ITEMP(1) = NSLK + NFEZ + JCOL
            RTHSKT(NZCOL) = ONE
C
          ELSE
C
            DO I =  JSTRH(KK),JSTRH(KK+1)-1
C
              IF(NZCOL.GT.0.AND.HMAT(I).EQ.ZERO) CYCLE
C
              IR = IROWH(I) + MINEQL + 1
C
C             IF REAL VARIABLE IR IS FREE, LOAD ROW 
C
              IF(IFREE(IR).GT.0) THEN
                NZCOL = NZCOL + 1
                ITEMP(NZCOL) = IFREE(IR)
                RTHSKT(NZCOL) = HMAT(I)
              ENDIF
C
            ENDDO
C
          ENDIF
C
C         ------ JACOBIAN BLOCK COLUMNS
C
C         LOAD COLUMN K OF AMAT
C
          IF(NROW.GT.0) THEN
            DO I =  JCOLST(KK),JCOLST(KK+1)-1
C
              IF(AMAT(I).EQ.ZERO) CYCLE
C
              NZCOL = NZCOL + 1
              ITEMP(NZCOL) = NSLK + NFEZ + NVAR + IROW(I) 
              RTHSKT(NZCOL) = AMAT(I)
C
            ENDDO
          ENDIF
C
          CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHSKT, 
     $                  RWORK, LNRWRK, IER )
C
          IF(IER.EQ.-400) THEN
            IERSKT = -829
          ELSEIF(IER.EQ.-401) THEN
            IERSKT = -930
          ELSEIF(IER.EQ.-402) THEN
            IERSKT = -831
          ELSEIF(IER.EQ.-410) THEN
            IERSKT = -632
          ENDIF
          IF(IERSKT.NE.0) GO TO 280
C
        ENDIF
C
 240  CONTINUE
C
C         -- LAST NROW COLUMNS ARE ZERO BLOCK; NO INPUT REQUIRED
C
      CALL XDSLVF ( RWORK, LNRWRK, IER )
C
      IF(IER.EQ.-400) THEN
        IERSKT = -833
      ELSEIF(IER.EQ.-401) THEN
        IERSKT = -934
      ELSEIF(IER.EQ.-410) THEN
        IERSKT = -635
      ENDIF
      IF ( IERSKT .NE. 0 ) GO TO 280
C
C         FACTOR THE MATRIX
C
      CALL XDSLFA (RWORK, LNRWRK, CNDNUM, INRTIA, NEEDS, IER)
C
C
      IF(IER.EQ.-500) THEN
        IERSKT = -836
      ELSEIF(IER.EQ.-501.OR.IER.EQ.-507.OR.IER.EQ.-508) THEN
        NEEDED = NEEDS
        IERSKT = -937
      ELSEIF(IER.EQ.-502) THEN
        IERSKT = -836
      ELSEIF(IER.EQ.-503) THEN
        IERSKT = -731
      ELSEIF(IER.EQ.-513) THEN
        IERSKT = -732
      ELSEIF(IER.EQ.-504) THEN
        IERSKT = -639
      ELSEIF(IER.EQ.-505.OR.IER.EQ.-511.OR.IER.EQ.-512) THEN
        IERSKT = -640
      ELSEIF(IER.EQ.-506) THEN
        IERSKT = -641
      ELSEIF(IER.EQ.-509) THEN
        IERSKT = -742
      ENDIF
C
C         CHECK TO SEE IF OUT-OF-CORE MODE WAS USED
C
      IF(RWORK(55).NE.ZERO) INSTAT(28) = 1
C
      IF(IERSKT.NE.0) GO TO 280
C
C         CHECK FOR ILL-CONDITIONED KT MATRIX 
C
      IF(IPC.GT.0) WRITE(IPU,*) 'CNDNUM =',CNDNUM,
     $    '    INERTIA =',(INRTIA(I),I=1,3)
C
C         CHECK INERTIA OF KT MATRIX
C
      IF(INRTIA(1).NE.(NEQNS-INREQD-NROW).OR.INRTIA(2).NE.(NROW+INREQD)
     $    .OR.INRTIA(3).NE.0) THEN
        IERSKT = -744
      ENDIF
C
      IF(IOFMFR.GT.0) CALL XDSLPS(RWORK)
C
C         GET MULTIFRONTAL STATISTICS
C
      CALL XDSLSR(RWORK,INUSE,MXUSED,TYME,OPCNTS)
      INSTAT(10) = MAX(INSTAT(10),INUSE)
C
      IF(IERSKT.NE.0) GO TO 280
C
 250  CONTINUE
C
C ------------------------------------------------------------------
C --------STEP 6----------------------------------------------------
C ------------------------------------------------------------------
C
C         DEFINE RIGHT HAND SIDE
C
      RTHSKT(1:NSLK) = ZERO
C
      IF(NFEZ.EQ.1) RTHSKT(NSLK+1) = GFEZ
C
      DO I = 1,NCOL
        IFRIM = IFREE(I+MINEQL+1)
        IF(IFRIM.GT.0) THEN
          RTHSKT(IFRIM) = GVEC(I)
        ENDIF
      ENDDO
      DO KK = 0,NROW-1
        RTHSKT(1+NSLK+NFEZ+NVAR+KK) = CVEC(KK+1)
      ENDDO
C
C         SOLVE LINEAR SYSTEM
C
      CALL REFITR(RTHSKT,SOLNKT,NEQNS,
     $    RWORK,LNRWRK,NEEDS,IER)
C
      IF(IER.EQ.-600) THEN
        IERSKT = -845
      ELSEIF(IER.EQ.-601) THEN
        IERSKT = -946
        NEEDED = NEEDS
      ELSEIF(IER.EQ.-602) THEN
        IERSKT = -847
      ELSEIF(IER.EQ.-603) THEN
        IERSKT = -848
      ELSEIF(IER.EQ.-604.OR.IER.EQ.-605) THEN
        IERSKT = -649
      ELSEIF(IER.EQ.-701) THEN
      ENDIF
      IF(IERSKT.NE.0) GO TO 280
C
C ------------------------------------------------------------------
C
C         --- STORE THE SOLUTION INTO PVEC AND PIVEC, AND
C         LOAD WORK ARRAY WITH RIGHT HAND SIDE
C
      DO KK = 1,NSLK+NFEZ+NVAR
        PVEC(KK) = -SOLNKT(KK)
      ENDDO
      DO KK = 0,NROW-1
        PIVEC(KK+1) = SOLNKT(1+NSLK+NFEZ+NVAR+KK)
      ENDDO
C
 280  CONTINUE
C
      IF(IERSKT.NE.0.AND.IPC.GT.0) WRITE(IPU,1002) IERSKT
C
C         END K-T FACTORIZATION TIMING
C
      CALL CLKSUM(6)
      CALL CLKSUM(13)
C      
 1001 FORMAT(5X,'QPSKT STORAGE ERROR;  NEEDS =',I6,'  LNRWRK =',I6)
 1002 FORMAT(5X,'MULTIFRONTAL ERROR;  IERSKT =',I6)
      RETURN 
      END
