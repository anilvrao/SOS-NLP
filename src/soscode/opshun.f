

      SUBROUTINE OPSHUN(LNOPTN)
C
C ======================================================================
C     OPSHUN===>opshun   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  DISPLAY THE NLPSPR.CMN COMMON OPTIONS
C
C-------------------------------------------------------------
C
      PARAMETER (NRSYM=11,NISYM=35,NCSYM=3)
      PARAMETER (NSYMBL=NRSYM+NISYM+NCSYM)
      PARAMETER (NLINES=87)
C
C======================================================================
C Array equivalences into include file NLPSPR.CMN
C
      INCLUDE '../commons/NLPSPR.CMN'
C
C     COMMON /NPSPRR/ ALFLWR, ...      double precision
C     COMMON /NPSPRI/ INNPER, ...      integer
C     COMMON /NPSPRC/ ALGOPT, ...      character(len=6)
C     COMMON /NPSALG/ ALGNAM           character(len=6) (no equivalence)
C
      DOUBLE PRECISION RNPSPR(NRSYM)
      INTEGER          INPSPR(NISYM+1)
      CHARACTER(LEN=6) CNPSPR(NCSYM)
C
      EQUIVALENCE (RNPSPR(1),ALFLWR)
      EQUIVALENCE (INPSPR(1),INNPER)
      EQUIVALENCE (CNPSPR(1),ALGOPT)
C
      INCLUDE '../commons/BARNLP.CMN'
C
C     COMMON /BNPSPR/ BIGCON, ...      double precision
C     COMMON /BNPSPI/ IMAXMU, ...      integer
C
C End of array equivalences
C======================================================================
C
C-------------------------------------------------------------
C
C         DEFAULT VALUE COMMON
C
      COMMON /NPSDF1/ RDFLT(11),IDFLT(35)
      COMMON /NBSDF1/ RBDFLT(5),IBDFLT(3)
      COMMON /NPSDF2/ CDFLT(3),DFAWLT
      CHARACTER(LEN=6) CDFLT
      CHARACTER(LEN=17) DFAWLT
C
C-------------------------------------------------------------
C
      CHARACTER(LEN=50) DESCRP(NLINES)
      CHARACTER(LEN=10) OUT1,OUT2
      CHARACTER(LEN=1)  STAR
      LOGICAL SUMRY
      DIMENSION IPTSYM(NSYMBL+1)
      CHARACTER(LEN=6)  SYMBUL(NSYMBL)
      DIMENSION IBRIEF(NSYMBL)
C
C-------------------------------------------------------------
C
      DATA (SYMBUL(I),I=1,NSYMBL) /
     $    'ALFLWR', 'ALFUPR', 'CONTOL', 'EPSRLF', 'OBJTOL', 
     $    'PGDTOL', 'SLPTOL', 'SFZTOL', 'TOLFIL', 'TOLKTC', 
     $    'TOLPVT', 'IHESHN', 'IOFLAG', 'IOFLIN', 'IOFMFR',      
     $    'IOFPAT', 'IOFSHR', 'IOFSRC', 'IPUDRF', 'IPUFZF',   
     $    'IPUMF1', 'IPUMF2', 'IPUMF3', 'IPUMF4', 'IPUMF5',    
     $    'IPUMF6', 'IPUMF7', 'IPUNLP', 'IPUSTF', 'IRELAX',       
     $    'ITDRQP', 'ITFZQP', 'IT1MAX', 'JACPRM', 'LYNFNC',    
     $    'LYNOUT', 'LYNPLT', 'LYNPNT', 'LYNVAR', 'MAXLYN',    
     $    'MAXNFE', 'MNSAME', 'NEWTON', 'NITMAX', 'NITMIN', 
     $    'NORMAL', 'ALGOPT', 'KTOPTN', 'QPOPTN' /
      DATA (IBRIEF(I),I=1,NSYMBL) /
     $  2*0,1,0,3*1,5*0,1,27*0,1,0,0,1,0,0,1,0,0 /
C
C
      DATA (DESCRP(I),I=1,10) /
     A  'Lower Bound on ALFA (Line Search Diagnostic Plot) ',
     B  'Upper Bound on ALFA (Line Search Diagnostic Plot) ',
     C  'Constraint Tolerance                              ',
     D  'Relative Perturbation Size Parameter              ',
     E  'Objective Function Tolerance                      ',
     F  'Projected Gradient Tolerance                      ',
     G  'Convergence Requires:                             ',
     H  '(A) Max Absolute Error in Active Constraints and  ',
     I  '    Bounds .LT. CONTOL                            ',
     J  '(B) Max Absolute Error in KT Conditions .LT.      '/
      DATA (DESCRP(I),I=11,20) /
     A  '    PGDTOL*MAX(1,|DELF|)                          ',
     B  '(C) |F(X) - FMIN| .LT. OBJTOL                     ',
     C  '(D) |Steplength| = ALFA*|P| .LT.                  ',
     D  '    SQRT( OBJTOL/(1 + |F(X)|) )*(1 + |X|)         ',
     E  '(E) Correct Sign For All Lagrange Multipliers     ',
     F  'Slope Tolerance For SQP Line Search               ',
     G  'Line Search Requires Sufficient Reduction         ',
     H  '    F(ALFA) - F(0) .LT. 1.D-4*ALFA*F''(0)         ',
     I  'and Slope Reduction                               ',
     J  '    |F''(ALFA)| .LT. -SLPTOL*F''(0)               '/
      DATA (DESCRP(I),I=21,30) /
     A  'Slope Tolerance For Feasiblity Line Search        ',
     B  'Multifrontal Fill Tolerance                       ',
     C  'KT Condition Number Tolerance                     ',
     D  'Multifrontal Pivot Tolerance                      ',
     E  'Hessian Matrix Evaluation Option                  ',
     F  ' 0  = IHESHN      Finite Difference or Analytic   ',
     G  '|1| = IHESHN      SR1 (Symmetric Rank One)        ',
     H  '|2| = IHESHN      BFGS (Symmetric Positive Def.)  ',
     I  '|3| = IHESHN      SSQN (Self-Scaling Quasi-Newton)',
     J  'IHESHN < 0 Finite Difference Initialization       '/
      DATA (DESCRP(I),I=31,40) /
     A  'Output Level                                      ',
     B  '0  = IOFLAG       No Output                       ',
     C  '0  < IOFLAG < 10  Terse Output                    ',
     D  '9  < IOFLAG < 20  Standard Output                 ',
     E  '19 < IOFLAG < 30  Interpretive Output             ',
     F  '     IOFLAG = 30  Diagnostic Output               ',
     G  'Line Search Output Level                          ',
     H  '(0,10,20) Overrides Default (IOFLIN=IOFLAG)       ',
     I  'Multifrontal Output Level (0,1,2,3,4)             ',
     J  'Output Sparsity Pattern                           '/
      DATA (DESCRP(I),I=41,50) /
     A  'Check Pattern (.lt.0); Override Default (.ge.10)  ',
     B  'QP Output Level                                   ',
     C  'Feasibility Search Output Level                   ',
     D  '(0,10,20,30) Overrides Default (IOFSRC=IOFLAG)    ',
     E  'Output Unit for QP Dump From Optimization         ',
     F  'Output Unit for QP Dump From Feasibility          ',
     G  'Multifrontal I/O Unit                             ',
     H  'Multifrontal I/O Unit                             ',
     I  'Multifrontal I/O Unit                             ',
     J  'Multifrontal I/O Unit                             '/
      DATA (DESCRP(I),I=51,60) /
     A  'Multifrontal I/O Unit                             ',
     B  'Multifrontal I/O Unit                             ',
     C  'Multifrontal I/O Unit                             ',
     D  'Output Unit Number                                ',
     E  'Output Unit No. for QP Dump From QP START         ',
     F  'Constraint Relaxation Strategy Option             ',
     G  '=  0       No Relaxation                          ',
     H  '=  1       Relaxation                             ',
     I  'Dump QP on SQP Iteration No.                      ',
     J  'Dump QP on Feasible Search Iteration No.          '/
      DATA (DESCRP(I),I=61,70) /
     A  'Maximum Number of Steps for SQP Line Search       ',
     B  'Jacobian Permutation Option (External Order = 1)  ',
     C  'Function Number (Line Search Diagnostic Plot)     ',
     D  'Output Unit Number (Line Search Diagnostic Plot)  ',
     E  'Iteration Number For Line Search Diagnostic Plot  ',
     F  'No. of Plot Points (Line Search Diagnostic Plot)  ',
     G  'Variable Number (Line Search Diagnostic Plot)     ',
     H  'Maximum Line Limit for Array Output               ',
     I  'Maximum Number of Function Evaluations            ',
     J  'Switch to Equality QP after MNSAME Steps          '/
      DATA (DESCRP(I),I=71,80) /
     A  'Newton Option (0,1,2) = (Default,Newton,Gauss)    ',
     B  'Maximum Number of Iterations                      ',
     C  'Minimum Number of Iterations                      ',
     D  'Least Squares Normal Matrix Option                ',
     E  'Algorithm Control Option                          ',
     F  '=  FM      Find Feasible Point Then Minimize      ',
     G  '=  FME     Find Feasible Point Then Minimize      ',
     H  '           With Equalities Binding                ',
     I  '=  M       Minimize From The Initial Point        ',
     J  '=  F       Find Feasible Point Only               '/
      DATA (DESCRP(I),I=81,NLINES) /
     A  '=  LLSQ    Linear Least Squares                   ',
     B  'KT Matrix Factorization Option                    ',
     C  '=  SMALL   Active Constraints Only                ',
     D  '=  LARGE   All Equality and Inequality Constraints',
     E  'Quadratic Programming Algorithm Option            ', 
     F  '=  SPARSE  Schur-Complement QP                    ', 
     G  '=  DENSE   Nullspace Indefinite QP                '/ 
C
C        NOTE IPTSYM(NSYMBL+1) = NLINES+1
C
      DATA (IPTSYM(I),I=1,NSYMBL+1) /
     $  1,2,3,4,5,6,16,21,22,23,24,25,31,37,39,40,42,43,45,
     $  46,47,48,49,50,51,52,53,54,55,56,59,60,61,62,63,64,
     $  65,66,67,68,69,70,71,72,73,74,75,82,85,88 /
C
C-------------------------------------------------------------
C
      SUMRY = .FALSE.
C
C         NOTE: IPUNLP = INPSPR(18) = IPUOPS
C
      IPUOPS = INPSPR(18)
C
      LNOPIN = LNOPTN
      IF(LNOPIN.LE.0) THEN
C
        SUMRY = .TRUE.
        LNOPIN = 1
C
      ELSEIF(LNOPIN.EQ.1) THEN
C
C         SHORT OPTIONS TITLE
C
        WRITE(IPUOPS,1002) ALGNAM
        WRITE(IPUOPS,1004)
        WRITE(IPUOPS,1005)
C
      ELSEIF(LNOPIN.LE.16) THEN
C
C         FULL OPTIONS TITLE
C
        WRITE(IPUOPS,1003) ALGNAM
        WRITE(IPUOPS,1004)
        WRITE(IPUOPS,1005)
C
      ELSE
        PRINT *,'LNOPIN TOO BIG'
        STOP
      ENDIF
C
      IF(SUMRY) THEN
         WRITE(IPUOPS,1007) DFAWLT
      ELSE
         WRITE(IPUOPS,1004)
      ENDIF
C
C         REAL VARIABLES
C
      DO I = 1,NRSYM
        WRITE(OUT1(1:10),'(1PG10.3)') RDFLT(I)
        WRITE(OUT2(1:10),'(1PG10.3)') RNPSPR(I)
        CALL HHADJF(OUT1,' ',' ','R',NS,IERSH)
        CALL HHADJF(OUT2,' ',' ','R',NS,IERSH)
        KLWR = IPTSYM(I)
        KUPR = KLWR + MIN(LNOPIN,IPTSYM(I+1)-KLWR) - 1
        IF(OUT1.EQ.OUT2) THEN
          STAR = ' '
        ELSE
          STAR = '*'
        ENDIF
        IF(SUMRY) THEN
          IF(STAR.EQ.'*'.OR.IBRIEF(I).EQ.1) THEN
            WRITE(IPUOPS,1006) STAR,SYMBUL(I),OUT2,
     $        (DESCRP(K),K=KLWR,KUPR)
          ENDIF
        ELSE
          WRITE(IPUOPS,1001) STAR,SYMBUL(I),OUT1,OUT2,
     $    (DESCRP(K),K=KLWR,KUPR)
        ENDIF
C
        IF(I.EQ.2) THEN
          WRITE(OUT1(1:10),'(1PG10.3)') RBDFLT(1)
          WRITE(OUT2(1:10),'(1PG10.3)') BIGCON
          CALL HHADJF(OUT1,' ',' ','R',NS,IERSH)
          CALL HHADJF(OUT2,' ',' ','R',NS,IERSH)
          IF(OUT1.EQ.OUT2) THEN
            STAR = ' '
          ELSE
            STAR = '*'
          ENDIF
          IF(SUMRY.AND.STAR.EQ.'*') THEN
            WRITE(IPUOPS,1006) STAR,'BIGCON',OUT2,
     $    'Upper Bound on Total Constraint Violation         '
          ELSEIF(.NOT.SUMRY) THEN
            WRITE(IPUOPS,1001) STAR,'BIGCON',OUT1,OUT2,
     $    'Upper Bound on Total Constraint Violation         ',
     $    'BIGCON < 0    Merit Function Globalization        ',
     $    'BIGCON > 0    Filter Globalization                '
          ENDIF
        ENDIF
C
      enddo
C
      IF(.NOT.SUMRY) WRITE(IPUOPS,1004)
C
C         INTEGER VARIABLES
C
      DO I = 1,NISYM
        II = I + NRSYM
        WRITE(OUT1(1:10),'(I10)') IDFLT(I)
        WRITE(OUT2(1:10),'(I10)') INPSPR(I+1)
        KLWR = IPTSYM(II)
        KUPR = KLWR + MIN(LNOPIN,IPTSYM(II+1)-KLWR) - 1
        IF(I.NE.3) THEN
          IF(OUT1.EQ.OUT2) THEN
            STAR = ' '
          ELSE
            STAR = '*'
          ENDIF
        ELSEIF(I.EQ.3) THEN
          IF(INPSPR(I+1).GT.0.AND.INPSPR(I+1).NE.INPSPR(I)) THEN
            STAR = '*'
          ELSE
            STAR = ' '
          ENDIF
        ENDIF
        IF(SUMRY) THEN
          IF(STAR.EQ.'*'.OR.IBRIEF(II).EQ.1) THEN
            WRITE(IPUOPS,1006) STAR,SYMBUL(II),OUT2,
     $        (DESCRP(K),K=KLWR,KUPR)
          ENDIF
        ELSE
          WRITE(IPUOPS,1001) STAR,SYMBUL(II),OUT1,OUT2,
     $    (DESCRP(K),K=KLWR,KUPR)
        ENDIF
      enddo
C
      IF(.NOT.SUMRY) WRITE(IPUOPS,1004)
C
C         CHARACTER VARIABLES
C
      DO I = 1,NCSYM
        II = I + NRSYM + NISYM
        OUT1(1:10) = '          '
        OUT2(1:10) = OUT1(1:10)
        OUT1(5:10) = CDFLT(I)
        OUT2(5:10) = CNPSPR(I)
        CALL HHADJF(OUT1,' ',' ','R',NS,IERSH)
        CALL HHADJF(OUT2,' ',' ','R',NS,IERSH)
        KLWR = IPTSYM(II)
        KUPR = KLWR + MIN(LNOPIN,IPTSYM(II+1)-KLWR) - 1
        IF(OUT1.EQ.OUT2) THEN
          STAR = ' '
        ELSE
          STAR = '*'
        ENDIF
        IF(SUMRY) THEN
          IF(STAR.EQ.'*'.OR.IBRIEF(II).EQ.1) THEN
            WRITE(IPUOPS,1006) STAR,SYMBUL(II),OUT2,
     $        (DESCRP(K),K=KLWR,KUPR)
          ENDIF
        ELSE
          WRITE(IPUOPS,1001) STAR,SYMBUL(II),OUT1,OUT2,
     $    (DESCRP(K),K=KLWR,KUPR)
        ENDIF
      enddo
C
 1001 FORMAT(T3,'*',T8,A1,A6,T20,A10,T35,A10,
     $    15(T50,A50,T106,'*',:,/T3,'*'))
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T33,A6,  ' OPTIONS  (ALSO SEE ''FULL
     $'' OPTIONS)',T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T33,A6,   ' FULL OPTIONS  (ALSO SEE 
     $''OPTIONS'')',T106,'*')
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T9,91('-'),T106,'*'/T3,'*',T106,'*')
 1005 FORMAT(T3,'*',T9,'SYMBOL',T23,'DEFAULT',T40,
     $  'VALUE',T50,'DESCRIPTION',T106,'*')
 1006 FORMAT(T3,'*',T14,A1,A6,T25,A10,
     $    15(T40,A50,T106,'*',:,/T3,'*'))
 1007 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'ALGORITHM CONTROL PARAMETERS',
     * 2X,A17,T106,'*'/T3,'*',T106,'*')
      RETURN
      END
