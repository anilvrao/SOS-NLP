      SUBROUTINE BARKKT( NVAR, NEQNS, NRES, MSUBE, MSUBB, NSLK, IPU,
     $                  IPC, NCALL, NONZR, IROWR, JCOLR, RMAT,
     $                  NONZB, IROWB, JCOLB, BMAT, 
     $                  NONZC, IROWC, JCOLC, CMAT, NONZW, IROWW,
     $                  JCOLW, WMAT, USEW, ADDHSS, MU, BVEC, CVEC, 
     $                  GRADL, ETA, LAMBDA,    
     $                  LNRWRK, RWORK,
     $                  ITEMP, NZROWB, ROWSCL, COLSCL, 
     $                  RTHBKT, SOLNKT, LNSYMB, CNDNUM, NEEDED,
     $                  DY, DZ, DETA, DLMBDA, IERBKT         ) 
C
C
C =====================================================================
C     BARKKT===>barkkt    J.T. Betts
C =====================================================================
C
      SAVE
C
C         PURPOSE:  SOLVE THE SET OF EQUATIONS RESULTING FROM 
C                   LINEARIZATION OF THE KUHN-TUCKER CONDITIONS
C                   WITH RELAXED COMPLEMENTARITY.
C                   THE MULTIFRONTAL SPARSE SYMMETRIC INDEFINITE
C                   SOLVER IS APPLIED TO A SYSTEM CONSTRUCTED FROM
C                   COLUMNS (AND ROWS) OF THE (SCALED) SYMMETRIC SYSTEM
C
C         NOTE: THIS CODE ACCOMODATES LEAST SQUARES AS WELL AS
C               "STANDARD" FORMULATIONS. FOR STANDARD FORMULATIONS
C               NRES=0, AND THERE ARE NO 2ND ROW AND 2ND COLUMN BLOCKS
C               IN THE PICTURE BELOW.
C
C                                                     
C |E1*W*V1                       | | V1(INV)*DY    | |-E1*GRADL       |
C |                              | |               | |                |
C |E2*R*V1     -I                | | V2(INV)*DZ    | |-E2*0           |
C |                              | |               | |                |
C |E3*C*V1      0    0           | |-V3(INV)*DETA  |=|-E3*CVEC        |
C |                              | |               | |                |
C |E4*LAM*B*V1  0    0  -E3*DB*V4| |-V4(INV)*DLMBDA| |-E4*(LAMBDA-PIB)|
C
C                   IN THE ABOVE:
C                      GRADL = (NVAR VECTOR) GRADIENT OF THE LAGRANGIAN
C                           AT PT Y, I.E. G(Y) -C(TR)*ETA -B(TR)*LAMBDA
C                      CVEC = (MSUBE VECTOR) EQUALITY CONSTR VECTOR
C                      PIB = (MSUBB VECTOR) WHERE PIB(I) = MU/BVEC(I),
C                             AND BVEC IS THE BOUNDS CONSTR VECTOR AT
C                             CURRENT PT Y. WHERE, MU = (SCALAR)
C                             RELAXATION OF COMPLEMENTARITY
C                      ETA = (MSUBE VECTOR) OF MULIPLIERS FOR THE
C                            EQUALITY CONSTRS
C                      LAMBDA = (MSUBB VECTOR) OF MULTIPLIERS FOR
C                                BOUNDS CONSTRS
C                      LAM = (MSUBB X MSUBB) DIAGONAL MATRIX OF LAMBDAS
C                      DY = (NVAR VECTOR) DELTA Y
C                      DZ = (NRES VECTOR) DELTA Z, WHERE Z IS A VECTOR
C                           OF EXTRA VARIABLES, INTRODUCED TO PUT A
C                           LEAST-SQUARES PROBLEM INTO "AUGMENTED" FORM
C                      DETA = (MSUBE VECTOR) DELTA ETA
C                      DLMBDA = (MSUBE VECTOR) DELTA LAMBDA
C                      W == WMAT = (NVAR X NVAR) HESSIAN OF THE
C                                   LAGRANGIAN AT PT Y.
C                      R == RMAT = (NRES X NVAR ) JACOBIAN OF LEAST
C                                  SQUARES EQUATIONS
C                      C == CMAT = (MSUBE X NVAR) JACOBIAN OF THE
C                                   EQUALITY CONSTRS
C                      B == BMAT = (MSUBB X NVAR) JACOBIAN OF THE
C                                   BOUNDS CONSTRS
C                      DB = DIAGONAL MATRIX OF BOUNDS CONSTRS, I.E.
C                           DB(I,I) = BVEC(I), FOR I=1,...,MSUBB.
C                      E1, E2, E3 & E4 ARE ROW SCALING MATRICES
C                      V1, V2, V3 & V4 ARE COLUMN SCALING MATRICES.
C                      NOTE THAT SYMMETRY REQUIRES THAT E1=V1, E2=V2,
C                      E3=V3, AND LAM*E4=V4.
C
C
C         INPUT:    
C
C             NVAR   NUMBER OF VARIABLES (I.E. THE DiMENSION OF Y
C                    VECTOR)
C             NEQNS  NUMBER OF EQUATIONS IN KKT SYSTEM=NVAR+NRES+MSUBE
C                    +MSUBB
C             NRES   NUMBER OF LEAST-SQUARES EQUATIONS
C             MSUBE  NUMBER OF EQUALITY CONSTRAINTS (ALSO, THE NUMBER
C                    OF ROWS IN CMAT)
C             MSUBB  NUMBER OF BOUNDS CONSTRAINTS (ALSO, THE NUMBER OF
C                    ROWS IN BMAT)
C             NSLK   NUMBER OF SLACK VARS, I.E. NUMBER OF INEQUALITITES
C                    IN THE EXTERNAL PROBLEM FORMULATION
C
C             IPU    OUTPUT UNIT NO.
C             IPC    OUTPUT CONTROL FLAG
C
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
C
C
C             NONZR  DIMENSION OF RMAT (GE. JCOLR(NVAR+1) - 1)
C             IROWR  INTEGER ROW INDEX VECTOR (NONZR)
C             JCOLR  INTEGER COLUMN START VECTOR (NVAR+1)
C             RMAT   NONZERO ELEMENTS OF LINEAR SYSTEM (NONZR)
C
C             NONZB  DIMENSION OF BMAT (GE. JCOLB(NVAR+1) - 1)
C             IROWB  INTEGER ROW INDEX VECTOR (NONZB)
C             JCOLB  INTEGER COLUMN START VECTOR (NVAR+1)
C             BMAT   NONZERO ELEMENTS OF LINEAR SYSTEM (NONZB)
C
C             NONZC  DIMENSION OF CMAT (GE. JCOLC(NVAR+1) - 1)
C             IROWC  INTEGER ROW INDEX VECTOR (NONZC)
C             JCOLC  INTEGER COLUMN START VECTOR (NVAR+1)
C             CMAT   NONZERO ELEMENTS OF LINEAR SYSTEM (NONZC)
C
C             NONZW  DIMENSION OF WMAT (.GE. JCOLW(NVAR+1)-1 )
C             IROWW  INTEGER ROW INDEX VECTOR (NONZW)
C             JCOLW  INTEGER COLUMN START VECTOR (NVAR+1)
C             WMAT   NONZERO ELEMENTS OF LOWER TRIANGLE OF WMAT (NONZW)
C             USEW   LOGICAL FLAG: TRUE WHEN USING WMAT, FALSE OTHERWISE
C
C             ADDHSS THE AMOUNT TO ADD TO THE DIAGONAL OF WMAT, TO
C                    IMPLEMENT THE LEVENBERG-MARQUARDT METHOD.
C
C             MU     AMOUNT OF RELAXATION OF COMPLEMENTARITY CONDITIONS
C
C             BVEC   BOUNDS CONSTR VALUES AT POINT Y (MSUBB VECTOR)
C             CVEC   EQUALITY CONSTRAINT VALUES AT Y (MSUBE VECTOR)
C
C             GRADL  GRADIENT OF THE LAGRANGIAN AT PT Y, 
C                    I.E. G(Y) -C(TR)*ETA -B(TR)*LAMBDA (NVAR
C                    VECTOR) 
C             ETA    MULIPLIERS FOR THE EQUALITY CONSTRS (MSUBE
C                    VECTOR)
C             LAMBDA MULIPLIERS FOR THE BOUNDS CONSTRS (MSUBB VECTOR)
C
C             LNRWRK LENGTH OF RWORK ARRAY (SEE BELOW)
C             RWORK  REAL WORK ARRAY (LNRWRK)
C
C                    LOWER BOUNDS FOR THE WORK ARRAYS
C                    REAL:
C                      LNRWRK > NKT + 4*NEQNS + 170
C                    WHERE
C                      NKT  > NRES + MSUBB + MSUBE + 1 + NONZW
C                             + NONZC + NONZR
C                      NEQNS = NVAR + NRES + MSUBE + MSUBB
C
C             ITEMP  INTEGER WORK ARRAY (NEQNS)
C             NZROWB INTEGER UTILITY ARRAY TO STORE COUNTS OF NONZEROS
C                    IN ROWS OF BMAT (MSUBB VECTOR)
C
C
C         OUTPUT:
C
C             ROWSCL ROW SCALING VECTOR FOR KKT SYSTEM (NEQNS)
C             COLSCL COLUMN SCALING VECTOR FOR KKT SYSTEM (NEQNS)
C             RTHBKT RIGHT HAND SIDE OF THE KKT SYSTEM (NEQNS)
C             SOLNKT SOLUTION OF THE KKT SYSTEM (NEQNS)
C
C             RWORK  WORK ARRAY CONTAINING MULTIFRONTAL FACTORIZATION
C                    INFORMATION FOR SUBSEQUENT USE STORED AS FOLLOWS
C                    RWORK(1) -- SYMBOLIC FACTORIZATION DATA (LNSYMB)
C                    RWORK(LNSYMB+1) -- NUMERIC FACTORIZATION DATA 
C                                       (LNRWRK)
C
C             LNSYMB LENGTH OF RWORK ARRAY NEEDED FOR SYMBOLIC
C                    FACTORIZATION
C             CNDNUM CONDITION NUMBER OF MATRIX
C             NEEDED STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C
C             DY     STEP IN THE VARIABLES, Y (NVAR VECTOR)
C             DZ     STEP IN THE LEAST SQUARES VARIABLES (NRES VECTOR).
C                    NOTE THAT DZ SHOULD = R*DY
C             DETA   STEP IN THE EQUALITY CONSTR MULTIPLIERS, 
C                    ETA (MSUBE VECTOR)
C             DLMBDA STEP IN THE BOUNDS CONSTR MULTIPLIERS, 
C                    LAMBDA (MSUBB VECTOR)
C
C
C             IERBKT ERROR RETURN CODE
C
C                    0    NORMAL RETURN
C
C             THE CONVENTION FOR ERROR RETURN NUMBERS IS AS FOLLOWS:
C             -330 > IERBKT > -370       I/O ERRORS (INSUFFICIENT DISK
C                                                    SPACE)
C             -370 > IERBKT > -800       BAD CONDITIONING, WRONG
C                                        INERTIA, EXCESSIVE FILL, OR
C                                        POOR SOLUTION ERRORS
C             -800 > IERBKT > -900   --- ERRORS IN INPUT
C             -900 > IERBKT          --- STORAGE ERRORS REQUIRING
C                                        A RESPONSE HIGHER UP THE
C                                        CALLING CHAIN
C
C             IERBKT IN CHRONOLOGICAL ORDER IN BARKKT:
C
C                    -801  NCALL.LT.1 OR NCALL.GT.3
C                    -802  MSUBE.LT.0
C                    -803  MSUBB.LT.0
C                    -804  NVAR.LE.0
C                    -806  NRES.LT.0
C                    -906  REAL WORK ARRAY TOO SMALL
C                    -907  REAL WORK ARRAY TOO SMALL (XDSLIN, IER=-101)
C                    -808  NEQNS .LT. 1 (XDSLIN, IER=-102)
C                    -809  INCORRECT PATH (XDSLIC, IER=-100)
C                    -910  REAL WORK ARRAY TOO SMALL (XDSLIC, IER=-101)
C                    -811  ILLEGAL VALUE FOR NZCOL (XDSLIC, IER=-103)
C                    -812  ILLEGAL VALUE FOR JCOL (XDSLIC, IER=-104)
C                    -813  ILLEGAL VALUE FOR JROWIN (XDSLIC, IER=-105)
C                    -914  REAL WORK ARRAY TOO SMALL (NEEDS.GT.LNRWRK),
C                          AFTER CALL TO XDSLIC
C                    -815  INCORRECT PATH (XDSLOR, IER=-170)
C                    -916  REAL WORK ARRAY TOO SMALL (XDSLOR, IER=-201)
C                    -617  I/O ERROR ON SQFILE (XDSLOR, IER=-202)
C                    -918  INTEGER WORK ARRAY TOO SMALL
C                          (XDSLOR, IER=-203)
C                    -819  INCORRECT PATH (XDSLSF, IER=-190)
C                    -920  INTEGER WORK ARRAY TOO SMALL
C                          (XDSLSF, IER=-301)
C                    -921  INTERNAL ERROR DURING SYMBOLIC FACT.
C                          (XDSLSF, IER=-302)
C                    -622  I/O ERROR ON SQFIL2 (XDSLSF, IER=-303)
C                    -923  REAL WORK ARRAY TOO SMALL
C                          (NEEDST.GT.LNRWRK), AFTER CALL TO XDSLSF
C                    -924  INTEGER WORK ARRAY TOO SMALL, AFTER
C                          CALL TO XDSLSF
C                    -850  THE DIAGONAL IS NOT THE 1ST ENTRY IN SOME
C                          COLUMN OF THE HESSIAN OF THE LAGRANGIAN
C                    -825  INCORRECT PATH (XDSLVC, IER=-210)
C                    -926  REAL WORK ARRAY TOO SMALL (XDSLVC, IER=-401)
C                    -827  JROWIN(I),JCOL PAIR NOT SPECIFIED
C                          (XDSLVC, IER=-402)
C                    -628  I/O ERROR ON SQFILE (XDSLVC, IER=-403)        
C                    -829  INCORRECT PATH (XDSLFA, IER=-230)             
C                    -930  REAL WORK ARRAY TOO SMALL (XDSLFA, IER=-501)
C                    -830  INSUFFICIENT STORAGE CAUSED BY EXTRA FILL-IN  
C                          DUE TO PIVOTING (XDSLFA, IER=-502)            
C                    -731  ZERO PIVOT; SINGULAR MATRIX
C                          (XDSLFA, IER=-503)
C                    -632  I/O ERROR ON SQFILE (XDSLFA, IER=-504)
C                    -633  I/O ERROR ON WAFIL1 (XDSLFA, IER=-505)        
C                    -634  I/O ERROR ON WAFIL2 (XDSLFA, IER=-506)
C                    -736  EXCESSIVE FILL DUE TO PIVOTING                
C                          (XDSLFA, IER=-509)                            
C                    -737  WRONG INERTIA OF KT MATRIX, AFTER CALL TO
C                          XDSLFA
C                    -838  INCORRECT PATH (XDSLSL, IER=-330)
C                    -939  REAL WORK ARRAY TOO SMALL FOR SOLUTION
C                          (XDSLSL, IER=-601)
C                    -840  NRHS .LT. 1 (XDSLSL, IER=-602)
C                    -841  LDRHS .LE. NEQNS (XDSLSL, IER=-603)
C                    -642  I/O ERROR ON WAFIL1 (XDSLSL, IER=-604)
C                    -743  MAX. ITERATION DURING ITERATIVE REFINEMENT
C                          (REFITR, IER=-701)
C                    -944  INTEGER WORK ARRAY TOO SMALL
C                          (XDSLSL, IER=-605)
C                    -852  INCORRECT INPUT WHEN REQUESTING KT FILE OUTPUT
C                    -853  TERMINATE AFTER SUCCESSFUL KT FILE OUTPUT
C                    -732  PIVOTING FAILURE -- NEED TO INCREASE PANEL SIZE
C                          (XDSLFA, IER=-513)
C
C
C     ...COMMON BLOCKS
C
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C
      INTEGER INSTAT
      DOUBLE PRECISION RLSTAT
      COMMON /STATIS/ INSTAT(30), RLSTAT(20)
C
      DOUBLE PRECISION ZEROMN, ZEROOT, BIGNUM, BGROOT, BIGBND, BIGCND
      COMMON /KONSTN/  ZEROMN, ZEROOT, BIGNUM, BGROOT, BIGBND, BIGCND 
C
C     ...SUBROUTINE CALLING ARGUMENTS
C
      INTEGER NVAR, NEQNS, NRES, MSUBE, MSUBB, NSLK
      INTEGER IPU, IPC, NCALL
      INTEGER LNSYMB, IERBKT
      INTEGER LNRWRK, NEEDED 
      INTEGER NONZR, NONZB, NONZC, NONZW
      INTEGER IROWR(NONZR), JCOLR(NVAR+1)
      INTEGER IROWB(NONZB), JCOLB(NVAR+1)
      INTEGER IROWC(NONZC), JCOLC(NVAR+1)
      INTEGER IROWW(NONZW), JCOLW(NVAR+1)
      INTEGER ITEMP(NEQNS)
      INTEGER NZROWB(*)
C
      LOGICAL USEW
      DOUBLE PRECISION ADDLEV, WMATI
      DOUBLE PRECISION ADDHSS, MU, CNDNUM
      DOUBLE PRECISION BVEC(*), CVEC(*), GRADL(NVAR), ETA(*), LAMBDA(*)
      DOUBLE PRECISION DY(NVAR), DZ(*), DETA(*), DLMBDA(*) 
      DOUBLE PRECISION RMAT(NONZR), BMAT(NONZB)
      DOUBLE PRECISION CMAT(NONZC), WMAT(NONZW)
      DOUBLE PRECISION RTHBKT(NEQNS), SOLNKT(NEQNS)
      DOUBLE PRECISION ROWSCL(NEQNS), COLSCL(NEQNS)
      DOUBLE PRECISION RWORK(LNRWRK)
C    
      LOGICAL OPUN,DUMPKT
C
C
C     ...LOCAL VARIABLES
C
      DOUBLE PRECISION ZERO, ONE, ONEEM2, ONEEM3
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,ONEEM2=1.0D-2,ONEEM3=1.0D-3)
C
      INTEGER MXERRS
      PARAMETER ( MXERRS = 120 )
      INTEGER NITREF
      PARAMETER (NITREF = 0)
C
      INTEGER I, JVAR, NWRKCK, II
      INTEGER JCOL, NZCOL, NEEDS, NEEDMN
      INTEGER JHSNAA, NEEDST, NRHS
      INTEGER NZLTA, MFRSTR, IER 
      INTEGER KKTOUT, LNZLTA
      INTEGER MXUSED, IROWKT
      INTEGER INRTIA(3)
      INTEGER XDSLNI, MAX
C
      DOUBLE PRECISION PIB,  DUMMY
      DOUBLE PRECISION TYME(6), OPCNTS(2)
C
C
C     FILL-IN THE VECTORS THAT MAP MULTIFRONTAL AND SOLVER
C     ERRORS TO BARKKT ERRORS (IERBKT).
C
      INTEGER KTERV(MXERRS), MFERV(MXERRS)

      DATA (MFERV(I), I=1, 39) / -100, -101, -102, -103, -104, -105,
     $                           -170, -201, -202, -203, -190, -301,
     $                           -302, -303, -210, -401, -402, -403,
     $                           -410, -230, -501, -502, -503, -504,
     $                           -505, -506, -509, -280, -511, -512,
     $                           -513, -330, -601, -602, -603, -604,
     $                           -605, -507, -508 / 
C
      DATA (KTERV(I), I=1, 39) / -809, -910, -808, -811, -812, -813,
     $                           -815, -916, -617, -918, -819, -920,
     $                           -921, -622, -825, -926, -827, -628,
     $                           -628, -829, -930, -829, -731, -632,
     $                           -633, -634, -736, -736, -633, -633,
     $                           -732, -838, -939, -840, -841, -642,
     $                           -944, -930, -930 /
C
C
C
C     BEGIN KKT FACTORIZATION TIMING
C
      CALL CLKBEG(6)
C
C     CHECK CALL NUMBER
C
      IERBKT = 0
C
      IF(NCALL.LT.1.OR.NCALL.GT.3) THEN
        IERBKT = -801
        GO TO 410
      ENDIF
C
C     CHECK ROW AND COLUMN DIMENSIONS
C
      IF(MSUBE.LT.0) THEN
        IERBKT = -802
        GO TO 410
      ENDIF
C
      IF(MSUBB.LT.0) THEN
        IERBKT = -803
        GO TO 410
      ENDIF
C
      IF(NRES.LT.0) THEN
        IERBKT = -806
        GO TO 410
      ENDIF
C
      IF(NVAR.LE.0) THEN
        IERBKT = -804
        GO TO 410
      ENDIF
C
C
C     ESTIMATE FOR NUMBER OF NONZEROS IN LOWER TRIANGULAR PART
C     OF KKT MATRIX
C
      NZLTA = NONZR + NONZW + NONZC + NONZB + NRES + MSUBE + MSUBB

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
C
      IF(NWRKCK.GT.LNRWRK) THEN
        NEEDED = NWRKCK
        IERBKT = -906
        GO TO 410
      ENDIF
C
C         To dump the KKT system and right hand side, the flag iofmfr
C         must be set negative.   To dump on a specific iteration the
C         value of iofmfr can be set based on the values of irevrs.
C         see socxoc_optcst for an example.  The KKT system is not
C         solved when dumping the quantities.
C
C
      DUMPKT = IOFMFR.LT.0.AND.IOFMFR.GT.-99
C
      IF (NCALL.EQ.2) THEN
        GO TO 220
      ELSEIF (NCALL.EQ.3) THEN
        GO TO 320
      ENDIF
C
 110  CONTINUE
C
C     SET INITIAL COUNT OF NONZEROS IN EACH ROW OF B MATRIX TO ZERO.   
C
      DO I=1, MSUBB
        NZROWB(I) = 0
      enddo
C
C ------------------------------------------------------------------
C --------STEPS 1, 2, AND 3-----------------------------------------
C ------------------------------------------------------------------
C
C     ...INPUT THE MATRIX STRUCTURE TO THE MULTIFRONTAL CODE
C
C
C     -- SET IN-CORE OPTION FLAG: 0 (IN-CORE); 1 (OUT-OF-CORE)
C
      MFRSTR = 0

      CALL XDSLIN ( NEQNS, 'SYMMETRIC', IOFMFR, IPU, IPUMF1, IPUMF2, 
     $              IPUMF3, IPUMF4, 0, IPUMF5, RWORK, LNRWRK, IER )

C
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
      IF(IERBKT.NE.0) GO TO 410
C
C          SET PANEL SIZE TO NEQNS UNTIL A BETTER SCHEME IS FOUND
C
      CALL XDSLSP ( 'panel size', NEQNS, DUMMY, RWORK, IER )
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
      IF(IERBKT.NE.0) GO TO 410
C
      CALL XDSLSP ( 'pivot tolerance', 0, TOLPVT, RWORK, IER )
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
      IF(IERBKT.NE.0) GO TO 410
C
      CALL XDSLSP ( 'limit fill', 0, TOLFIL, RWORK, IER )
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
      IF(IERBKT.NE.0) GO TO 410
C
      CALL XDSLSP ( 'save original matrix', 0, ZERO, RWORK, IER )
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
      IF(IERBKT.NE.0) GO TO 410
C
C
C     -- FIRST NVAR COLUMNS
C
      DO JCOL =1, NVAR
C
C       LOAD COLUMN JCOL
C
        NZCOL = 0
C
C       ...LOAD WMAT BLOCK
C
        wmloop: DO I=JCOLW(JCOL), JCOLW(JCOL+1) - 1
          IF(NZCOL.NE.0..AND.WMAT(I).EQ.ZERO) cycle wmloop
          NZCOL = NZCOL + 1
          ITEMP(NZCOL) = IROWW(I)
C
        enddo wmloop
C
C       ...LOAD RMAT BLOCK
C
        IF(NRES.GT.0) THEN
          colrloop: DO I=JCOLR(JCOL), JCOLR(JCOL+1) - 1
            IF(RMAT(I).EQ.ZERO) cycle colrloop
            NZCOL = NZCOL + 1
            ITEMP(NZCOL) = NVAR + IROWR(I)
          enddo colrloop
        ENDIF
C
C       ...LOAD CMAT BLOCK
C
        IF(MSUBE.GT.0) THEN
          colcloop: DO I=JCOLC(JCOL), JCOLC(JCOL+1) - 1
            IF(CMAT(I).EQ.ZERO) cycle colcloop
            NZCOL = NZCOL + 1
            ITEMP(NZCOL) = NVAR + NRES + IROWC(I)
          enddo colcloop
        ENDIF
C
C       ...LOAD BMAT BLOCK
C
        IF(MSUBB.GT.0) THEN
          DO I=JCOLB(JCOL), JCOLB(JCOL+1) - 1
            NZCOL = NZCOL + 1
            ITEMP(NZCOL) = NVAR + NRES + MSUBE + IROWB(I)
            NZROWB(IROWB(I)) = NZROWB(IROWB(I)) + 1
          enddo
        ENDIF
C
C        
C
        CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )

C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
C
C     ERROR CHECK ON THE COUNTS OF NONZEROS IN THE ROW OF THE B MATRIX.
C
      DO I=1, MSUBB
        IF ( NZROWB(I) .GT. 2 ) THEN
          WRITE(IPU,*) ' BARKKT: NONZERO COUNT FOR ROW NUMBER ', I,
     $                 ' OF B MATRIX IS = ', NZROWB(I)
          WRITE(IPU,*) ' THE NONZERO COUNT MUST LE 2: RUN STOPPED!'
C
          RETURN
        ENDIF
      enddo
C
C 
C     -- NEXT NRES COLUMNS, CORRESPONDING TO A NEGATIVE IDENTITY
C        ON THE DIAGONAL WITH TWO ZERO BLOCKS BELOW IT.
C
      DO I=1, NRES
C       NOTE: SINCE -I IS A DIAGONAL MATRIX ON A DIAGONAL BLOCK
C             OF THE KKT MATRIX THE ROW AND COLUMN INDICES ARE
C             THE SAME.
C
        JCOL = NVAR + I
        NZCOL = 1
        ITEMP(NZCOL) = JCOL
C
        CALL XDSLIC ('A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )

C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
C
C     -- NEXT MSUBE COLUMNS, CORRESPONDING TO TWO BLOCKS OF ZEROS.
C        MUST INPUT STRUCTURE FOR ZERO DIAGONAL BLOCK
C
      DO I=1, MSUBE
C       NOTE: FOR A DIAGONAL MATRIX ON A DIAGONAL BLOCK
C             OF THE KKT MATRIX, THE ROW AND COLUMN INDICES ARE
C             THE SAME.
C
        JCOL = NVAR + NRES + I
        NZCOL = 1
        ITEMP(NZCOL) = JCOL
C
        CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER )

C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
C     -- NEXT MSUBB COLUMNS, CORRESPONDING TO -R(D_b)V
C
      DO I=1, MSUBB
C       NOTE: SINCE D_b IS A DIAGONAL MATRIX ON A DIAGONAL BLOCK
C             OF THE KKT MATRIX, THE ROW AND COLUMN INDICES ARE
C             THE SAME.
C
        JCOL = NVAR + NRES + MSUBE + I
        NZCOL = 1
        ITEMP(NZCOL) = JCOL
C
        CALL XDSLIC ( 'A', JCOL, NZCOL, ITEMP, RWORK, LNRWRK, IER)

C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
      CALL XDSLIF (  RWORK, LNRWRK, NEEDS, IER )
C
      IF(IER.EQ.-100) THEN
        IERBKT = -809
      ELSEIF(IER.EQ.-101) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IF(IPC.GT.0) WRITE(IPU,1201) NEEDS,LNRWRK
          NEEDED = NEEDS
          IERBKT = -915
        ENDIF
      ELSEIF(IER.EQ.-110) THEN
        IERBKT = -617
      ENDIF
      IF ( IERBKT .NE. 0 ) GO TO 410
C
C         ORDER THE MATRIX
C
      CALL XDSLOR ( RWORK, LNRWRK, NEEDS, IER )
C
      IF(IER.EQ.-201) THEN
        IF(IPC.GT.0) WRITE(IPU,1201) NEEDS,LNRWRK
        NEEDED = NEEDS
      ENDIF
C
C     MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C     CODE IERBKT.
C
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
      IF(IERBKT.NE.0) GO TO 410
C
C
C
C     ...PERFORM SYMBOLIC FACTORIZATION
C
      CALL XDSLSF ( RWORK, LNRWRK, NEEDS, NEEDMN, IER )
C
C
C     MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C     CODE IERBKT.
C
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
      IF(IER.EQ.-301) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IF(IPC.GT.0) WRITE(IPU,1201) NEEDS,LNRWRK
          NEEDED = NEEDS
        ENDIF
      ENDIF
C
      IF(IERBKT.NE.0) GO TO 410
C
C
C     CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(MFRSTR.EQ.0) THEN
        NEEDST = NEEDS
      ELSE
        NEEDST = NEEDMN
      ENDIF
C
      NEEDST = INT((ONE + TOLPVT)*DBLE(NEEDST)) 
      IF(NEEDST.GT.LNRWRK) THEN
        IF(IPC.GT.0) WRITE(IPU,1201) NEEDST,LNRWRK
        NEEDED = NEEDST
        IERBKT = -923
        GO TO 410
      ENDIF
C
C     GET MULTIFRONTAL STATISTICS -- SAVE LENGTH NEEDED
C     FOR SYMBOLIC FACTORIZATION
C
      CALL XDSLSR(RWORK,LNSYMB,MXUSED,TYME,OPCNTS)
C
 220  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 4 AND 5---------------------------------------------
C ------------------------------------------------------------------
C
C     ...INPUT MATRIX VALUES TO MULTIFRONTAL SOLVER
C
C     -- SET DEFAULT ROW AND COLUMN SCALING TO 1.0
C
      rowscl(1:neqns) = one
      colscl(1:neqns) = one
C
C     -- FIRST NVAR COLUMNS
C
      DO JCOL =1, NVAR
C
C       LOAD COLUMN JCOL
C
        NZCOL = 0
C
C       ...LOAD WMAT BLOCK
C
        IF(JCOL .GT. NVAR - NSLK ) THEN
          ADDLEV = ZERO
        ELSE
          ADDLEV = ADDHSS
        ENDIF
        colwloop: DO I=JCOLW(JCOL), JCOLW(JCOL+1) - 1
          IF(NZCOL.NE.0..AND.WMAT(I).EQ.ZERO) cycle colwloop
          NZCOL = NZCOL + 1
          ITEMP(NZCOL) = IROWW(I)
C
          IF(USEW) THEN
            WMATI = WMAT(I)
          ELSE
            WMATI = ZERO
          ENDIF
C
          IF ( I .EQ. JCOLW(JCOL) ) THEN
C           FIRST ENTRY IN COL JCOL OF HESSIAN BLOCK (ONLY LOWER
C           TRIANGLE IS STORED).
C
            IF ( IROWW(I) .NE. JCOL ) THEN
C             ERROR EXIT: DIAGONAL MUST BE 1ST ENTRY IN EACH COLUMN.
C
              IERBKT = -850
C
              IF(IPC.GT.0) WRITE(IPU,1203) JCOL
              GO TO 410
            ENDIF
C
C           ADD ADDHSS TO DIAGONAL ENTRY OF HESSIAN OF LAGRANGIAN.
C           THEN, SCALE IT AND SET IT UP TO BE PUT IN COL JCOL.
C
            RTHBKT(NZCOL) = ROWSCL(ITEMP(NZCOL))*(WMATI+ADDLEV)*
     $                      COLSCL(JCOL)
          ELSE
C
C           AN ENTRY IN COL JCOL THAT'S SUPPOSED TO BE BELOW THE DIAG
C           SCALE IT AND SET IT UP TO BE PUT IN COL JCOL.
C
            RTHBKT(NZCOL) = ROWSCL(ITEMP(NZCOL))*WMATI*COLSCL(JCOL)
          ENDIF
        enddo colwloop
C
C       ...LOAD RMAT BLOCK
C
        IF(NRES.GT.0) THEN
          clrloop: DO I=JCOLR(JCOL), JCOLR(JCOL+1) - 1
            IF(RMAT(I).EQ.ZERO) cycle clrloop
            NZCOL = NZCOL + 1
            ITEMP(NZCOL) = NVAR + IROWR(I)
            RTHBKT(NZCOL) = ROWSCL(ITEMP(NZCOL))*RMAT(I)*COLSCL(JCOL)
          enddo clrloop
        ENDIF
C
C       ...LOAD CMAT BLOCK
C
        IF(MSUBE.GT.0) THEN
          clcloop: DO I=JCOLC(JCOL), JCOLC(JCOL+1) - 1
            IF(CMAT(I).EQ.ZERO) cycle clcloop
            NZCOL = NZCOL + 1
            ITEMP(NZCOL) = NVAR + NRES + IROWC(I)
            RTHBKT(NZCOL) = ROWSCL(ITEMP(NZCOL))*CMAT(I)*COLSCL(JCOL)
          enddo clcloop
        ENDIF
C
C       ...LOAD BMAT BLOCK
C
        IF(MSUBB.GT.0) THEN
          DO I=JCOLB(JCOL), JCOLB(JCOL+1) - 1
            NZCOL = NZCOL + 1
            ITEMP(NZCOL) = NVAR + NRES + MSUBE + IROWB(I)
            IROWKT = ITEMP(NZCOL)
C
C
C         COMPUTE ROW SCALING FOR ROW IROWKT OF KKT MATRIX.
C         DUE TO THE NEED FOR THE SCALING TO SYMMETRIZE THE NEWTON
C         SYSTEM, THIS WILL ALSO DETERMINE THE COLUMN SCALING FOR
C         COLUMN IROWKT OF THE KKT MATRIX.
C
            ROWSCL(IROWKT) = ONE/SQRT( LAMBDA(IROWB(I)) )
            COLSCL(IROWKT) = SQRT( LAMBDA(IROWB(I)) )
C
            RTHBKT(NZCOL) = ROWSCL(IROWKT)*LAMBDA(IROWB(I))*
     $                    BMAT(I)*COLSCL(JCOL)

          enddo
        ENDIF
C
        CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHBKT, 
     $                 RWORK, LNRWRK, IER )
C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
C 
C     -- NEXT NRES COLUMNS, CORRESPONDING TO A NEGATIVE IDENTITY AND
C        TWO BLOCKS OF ZEROS BELOW THE DIAGONAL.
C
      DO I=1, NRES
C       NOTE: FOR A DIAGONAL MATRIX ON A DIAGONAL BLOCK
C             OF THE KKT MATRIX, THE ROW AND COLUMN INDICES ARE
C             THE SAME.
C
        JCOL = NVAR + I
        NZCOL = 1
        ITEMP(NZCOL) = JCOL
        RTHBKT(NZCOL) = -1.0D0
C
        CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHBKT, 
     $                 RWORK, LNRWRK, IER )
C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
C 
C     -- NEXT MSUBE COLUMNS, CORRESPONDING TO TWO BLOCKS OF ZEROS
C        MUST INPUT THE DIAGONAL BLOCK OF ZEROS
C
      DO I=1, MSUBE
C       NOTE: FOR A DIAGONAL MATRIX ON A DIAGONAL BLOCK
C             OF THE KKT MATRIX, THE ROW AND COLUMN INDICES ARE
C             THE SAME.
C
        JCOL = NVAR + NRES + I
        NZCOL = 1
        ITEMP(NZCOL) = JCOL
        RTHBKT(NZCOL) = ZERO
C
        CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHBKT, 
     $                 RWORK, LNRWRK, IER )
C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
C
C     -- NEXT MSUBB COLUMNS, CORRESPONDING TO -R(D_b)V
C
      DO I=1, MSUBB
C       NOTE: SINCE D_b iS A DIAGONAL MATRIX ON A DIAGONAL BLOCK
C             OF THE KKT MATRIX, THE ROW AND COLUMN INDICES ARE
C             THE SAME.
C
        JCOL = NVAR + NRES + MSUBE + I
        NZCOL = 1
        ITEMP(NZCOL) = JCOL
        RTHBKT(NZCOL) = ROWSCL(ITEMP(NZCOL))*(-BVEC(I))*COLSCL(JCOL)
C
        CALL XDSLVC ( 'A', JCOL, NZCOL, ITEMP, RTHBKT, 
     $                 RWORK, LNRWRK, IER )
C
C       CHECK DIMENSION OF REAL WORK ARRAY
C
        IF(IER.EQ.-401) THEN
          IF(IPC.GT.0) WRITE(IPU,1201) NEEDS,LNRWRK
          NEEDED = NEEDS
        ENDIF
C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
        IF(IERBKT.NE.0) GO TO 410
C
      enddo
C
      CALL XDSLVF ( RWORK, LNRWRK, IER )
C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
      IF(IERBKT.NE.0) GO TO 410
C
C         BYPASS FURTHER PROCESSING WHEN DUMPING KKT SYSTEM
C
      IF(DUMPKT) GO TO 320
C
C     ...FACTOR THE MATRIX
C
      CALL XDSLFA  ( RWORK, LNRWRK, CNDNUM, INRTIA, NEEDS, IER )
C
C     CHECK TO SEE IF OUT-OF-CORE MODE WAS USED
C
      IF(RWORK(55).NE.ZERO) INSTAT(28) = 1
C
C
C     CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(IER.EQ.-501 .OR. IER.EQ.-502) THEN
        IF(IPC.GT.0) WRITE(IPU,1201) NEEDS,LNRWRK
        NEEDED = NEEDS
      ENDIF
C
C     MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C     CODE IERBKT.
C
      IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
C
      IF(IERBKT.NE.0) GO TO 410
C
      IF(IPC.GT.0) WRITE(IPU,*) 'CNDNUM =',CNDNUM,
     $    '    INERTIA =',(INRTIA(I),I=1,3)
C
C     ...CHECK INERTIA OF KKT MATRIX
C
      IF( INRTIA(1) .NE. (NEQNS-(NRES+MSUBE+MSUBB)) 
     $    .OR. INRTIA(2) .NE. (NRES+MSUBE+MSUBB)
     $    .OR. INRTIA(3) .NE. 0 ) THEN
        IERBKT = -737
      ENDIF
C
      IF(IOFMFR.GT.0) CALL XDSLPS(RWORK)
C
C     ...GET MULTIFRONTAL STATISTICS
C
      CALL XDSLSR(RWORK,INUSE,MXUSED,TYME,OPCNTS)
      INSTAT(10) = MAX(INSTAT(10),INUSE)
C
      IF(IERBKT.NE.0) GO TO 410
C
 320  CONTINUE
C
C ------------------------------------------------------------------
C --------STEP 6----------------------------------------------------
C ------------------------------------------------------------------
C
C     ...DEFINE RIGHT HAND SIDE
C
C     FIRST BLOCK ( NVAR ENTRIES )
C
      DO I=1, NVAR
        IROWKT = I
C                                                  T       T
C       GRADL IS GRADIENT OF THE LAGRANGIAN = g - C eta - B lambda
C
        RTHBKT(IROWKT) = -( ROWSCL(IROWKT)*GRADL(I) )
      enddo
C
C
C     SECOND BLOCK (NRES ENTRIES)
C
      DO I=1, NRES
        IROWKT = NVAR + I
C
        RTHBKT(IROWKT) = ZERO
      enddo
C
C
C     THIRD BLOCK (MSUBE ENTRIES)
C
      DO I=1, MSUBE
        IROWKT = NVAR + NRES + I
C
        RTHBKT(IROWKT) = -( ROWSCL(IROWKT)*CVEC(I) )
      enddo
C
C
C     FOURTH BLOCK (MSUBB ENTRIES)
C
      DO I=1, MSUBB
        IROWKT = NVAR + NRES + MSUBE + I
C
C       RECALL, PIB(I) = MU/BVEC(I)
C
        PIB = MU/BVEC(I)
        RTHBKT(IROWKT) = -BVEC(I)*( ROWSCL(IROWKT)*(LAMBDA(I) - PIB) )
      enddo
C
C
      IF  ( DUMPKT ) THEN
C
C         CHECK IF |IOFMFR| IS A VALID UNIT NUMBER FOR KT OUTPUT
C
        KKTOUT = ABS(IOFMFR)
        OPUN = .FALSE.
        INQUIRE(KKTOUT,OPENED=OPUN)
        IF(OPUN) THEN
          IERBKT = -852
        ELSE
          OPEN(KKTOUT,FILE='KKTMATRIX.FIL',STATUS='UNKNOWN')
C
          CALL XDSLWM ( 'ASCII','KKT MATRIX FROM SOCX', KKTOUT, 
     $                        RWORK, LNRWRK, IER )
          IF(IER.NE.0) THEN
            IERBKT = -852
          ELSE
            WRITE(KKTOUT,'(D25.17)') (RTHBKT(II),II=1,NEQNS)
            IERBKT = -853
          ENDIF
        ENDIF
C
        CLOSE(KKTOUT)
        GO TO 410
C
      ENDIF
C
C     ...SOLVE SYSTEM
C
      CALL REFITR(RTHBKT,SOLNKT,NEQNS,RWORK,LNRWRK,NEEDS,IER)
C
C
      IF ( IER .EQ. -601 .OR. IER .EQ. -605 ) THEN
        NEEDED = NEEDS
      ENDIF
C
      IF ( IER .NE. -701 ) THEN
C
C       MAP MULTI-FRONTAL ERROR CODE IER TO BARRIER KKT ERROR
C       CODE IERBKT.
C
        IERBKT =  JHSNAA( IPU, IER, MFERV, KTERV, MXERRS )
c
      ENDIF
C
      IF(IERBKT.NE.0) GO TO 410
C
C
C ------------------------------------------------------------------
C
C     ...MAP SOLNKT INTO STEP VECTORS DY, DETA, AND DLMBDA.
C        (SCALING AND SIGN CHANGES MUST BE ACCOUNTED FOR.)
C
C
      DO I=1, NVAR
        JVAR = I
        DY(I) = SOLNKT(JVAR)*COLSCL(JVAR)
      enddo
C
      DO I=1, NRES
        JVAR = NVAR + I
        DZ(I) = SOLNKT(JVAR)*COLSCL(JVAR)
      enddo
C
      DO I=1, MSUBE
        JVAR = NVAR + NRES + I
        DETA(I) = -SOLNKT(JVAR)*COLSCL(JVAR)
      enddo
C
      DO I=1, MSUBB
        JVAR = NVAR + NRES + MSUBE + I
        DLMBDA(I) = -SOLNKT(JVAR)*COLSCL(JVAR)
      enddo
C
C
C     -------------
C     ... TERMINATE
C     -------------

 410  CONTINUE
      IF(IERBKT.NE.0.AND.IPC.GT.0) WRITE(IPU,1202) IERBKT
C
C     END KKT FACTORIZATION TIMING
C
      CALL CLKSUM(6)
C      
 1201 FORMAT(5X,'BARKKT STORAGE ERROR;  NEEDS =',I6,'  LNRWRK =',I6)    
 1202 FORMAT(5X,'MULTIFRONTAL ERROR;  IERBKT =',I6)
 1203 FORMAT(5X,'ERROR: HESSIAN MISSING DIAGONAL IN COL NO. = ', I6)
C
      RETURN 
      END
