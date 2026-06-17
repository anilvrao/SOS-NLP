


 
      SUBROUTINE BARSQP(      FBAR     ,FOBJ     ,YBAR     ,YVEC           
     $   ,DELYVC   ,GVEC     ,GRADL    ,GRDLNM   ,NVAR     ,NFREE       
     $   ,MAXRES   ,NRES     ,DELZVC   ,RMAT     ,IROWR    ,JCOLR       
     $   ,NONZR    ,CBAR     ,CVEC     ,ETABAR   ,ETAVEC   ,DELETA      
     $   ,FHFLTR   ,MXFLTR   ,MSUBE    ,MINEQL   ,MAXCON   ,CMAT            
     $   ,IROWC    ,JCOLC    ,NONZC    ,BBAR     ,BVEC     ,VLAMBR                  
     $   ,VLAMDA   ,DELLAM   ,MSUBB    ,MAXBND   ,BMAT     ,IROWB                     
     $   ,JCOLB    ,NONZB    ,WMAT     ,IROWW    ,JCOLW    ,NONZW                
     $   ,PENMU    ,PENRHO   ,SMTUVW   ,RWORK          
     $   ,LNRWRK   ,SCRTCH   ,LNSCRT   ,ISCRTC   ,LNISCR   ,NEEDED              
     $   ,IREVRS   ,IRVCOM   ,IFERR    ,FZMODE   ,RELAX    ,ITERM    
     $   ,IQPTRM)
C
C
C ======================================================================
C     BARSQP===>barsqp   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C        PURPOSE:     COMPUTE THE VALUES OF THE NVAR VARIABLES Y WHICH
C                     MINIMIZE THE LOGARITHMIC BARRIER FUNCTION
C
C                                           m_b
C                       B(y,mu) = f(y) - mu sum ln[b_{i}(y)] 
C                                           i=1
C
C                                 + rho sum [ |t| + |u| + |v| + |w| ]
C
C                     SUBJECT TO THE MSUBE EQUALITY CONSTRAINTS
C
C                             c(y) = 0
C
C                     THE MSUBB BOUND INEQUALITY CONSTRAINTS
C
C                             b(y) .ge. 0 
C
C                     ARE TREATED USING THE LOGARITHMIC BARRIER TERM.
C
C        ARGUMENTS:
C
C        FBAR         OBJECTIVE FUNCTION AT YBAR
C        FOBJ         OBJECTIVE FUNCTION AT YVEC 
C        YBAR         CURRENT POINT (NVAR)
C        YVEC         OLD POINT (NVAR)
C        DELYVC       YVEC SEARCH DIRECTION (NVAR)
C        GVEC         GRADIENT AT YVEC (NVAR)
C        GRADL        GRADIENT OF LAGRANGIAN (NVAR)
C        GRDLNM       NORM OF GRADIENT OF LAGRANGIAN 
C        NVAR         NUMBER OF INTERNAL VARIABLES; NVAR 
C        NFREE        NUMBER OF FREE REAL VARIABLES
C        MAXRES       MAXIMUM NUMBER OF RESIDUALS MAX(NRES,1)
C        NRES         NUMBER OF RESIDUALS
C        DELZVC       RESIDUAL STEP DIRECTION (MAXRES)
C        RMAT         RESIDUAL JACOBIAN DERIVATIVES AT PVEC (NONZR)
C        IROWR        ROW INDICES OF RESIDUAL JACOBIAN NONZEROS (NONZR)
C        JCOLR        COLUMN START INDICES OF RESIDUAL JACOBIAN NONZEROS (NVAR+1)
C        NONZR        NUMBER OF RESIDUAL JACOBIAN NONZEROS; NONZR = JCOLR(NVAR+1)-1
C        CBAR         EQUALITY CONSTRAINTS AT YBAR (MAXCON)
C        CVEC         EQUALITY CONSTRAINTS AT YVEC (MAXCON)
C        ETABAR       CONSTRAINT MULTIPLIERS AT YBAR (MAXCON)
C        ETAVEC       CONSTRAINT MULTIPLIERS AT YVEC (MAXCON)
C        DELETA       ETAVEC SEARCH DIRECTION (MAXCON)
C        FHFLTR       FILTER ARRAY (MXFLTR*5)
C        MXFLTR       MAXIMUM LENGTH OF FILTER ARRAY 
C        MSUBE        NUMBER OF EQUALITY CONSTRAINTS MSUBE 
C        MAXCON       MAXIMUM NUMBER OF CONSTRAINTS MAX(MSUBE,1)
C        CMAT         CONSTRAINT DERIVATIVES AT YVEC (NONZC)
C        IROWC        ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C        JCOLC        COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C        NONZC        NUMBER OF JACOBIAN NONZEROS; NONZC = JCOLC(NVAR+1)-1
C        BBAR         BOUND INEQUALITIES AT YBAR (MAXBND)
C        BVEC         BOUND INEQUALITIES AT YVEC (MAXBND)
C        VLAMBR       BOUND MULTIPLIERS AT YBAR (MAXBND)
C        VLAMDA       BOUND MULTIPLIERS AT YVEC (MAXBND)
C        DELLAM       VLAMDA SEARCH DIRECTION (MAXBND)
C        MSUBB        NUMBER OF BOUNDS
C        MAXBND       MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C        BMAT         BOUND DERIVATIVES AT YVEC (NONZB)
C        IROWB        ROW INDICES OF BOUND JACOBIAN NONZEROS (NONZB)
C        JCOLB        COLUMN INDICES OF BOUND JACOBIAN NONZEROS (NVAR+1)
C        NONZB        NUMBER OF BOUND JACOBIAN NONZEROS
C        WMAT         HESSIAN OF THE LAGRANGIAN (NONZW).
C                     NOTE THAT 
C                               wmat = | hmat  0 |
C                                      | 0     0 |
C                     WHERE HMAT IS THE HESSIAN WITH RESPECT TO THE NFREE
C                     REAL VARIABLES
C        IROWW        INTEGER ROW INDEX VECTOR FOR LOWER TRIANGLE OF WMAT (NONZW)
C        JCOLW        INTEGER COLUMN START VECTOR (NVAR+1)
C        NONZW        NUMBER OF NONZEROS IN THE HESSIAN (NONZW = JCOLW(NVAR+1)-1)
C        PENMU        PENALTY WEIGHT MU 
C        PENRHO       PENALTY WEIGHT RHO (AND MAXIMUM MULTIPLIER LIMIT) 
C        SMTUVW       THE TERM sum [ |t| + |u| + |v| + |w| ]:  MUST BE SET
C                     ZERO FOR PRIMARY MODE
C        RWORK        REAL WORK ARRAY (LNRWRK)
C        LNRWRK       LENGTH OF WORK 
C        SCRTCH       REAL SCRATCH ARRAY (LNSCRT)
C        LNSCRT       LENGTH OF SCRTCH (4*NEQNS + 2 + 3*NVAR + 2*MAXCON + 2*MAXBND + 5*MXQPFL) 
C        ISCRTC       INTEGER SCRATCH ARRAY (LNISCR)
C        LNISCR       LENGTH OF ISCRTC MAX(MXFLTR,NEQNS + MAXBND + 2) 
C        NEEDED       STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C        IREVRS       REVERSE COMMUNICATION FLAG
C                     (1)  =1 WHEN EVALUATING FUNCTIONS
C                          =0 OTHERWISE.
C                     (2)  =1 , OR 2 WHEN EVALUATING GRADIENTS
C                          =0 OTHERWISE.
C                     (3)   HESSIAN CALL FLAG--
C                          =2 WHEN HESSIAN IS REQUESTED
C                          =1 WHEN HESSIAN DIAGONAL IS REQUESTED
C                          =0 OTHERWISE.
C                     (4)  =1 FOR SYSTEM PRINT
C                     (5)  AUXILLIARY INFORMATION
C        IRVCOM       = 0 WHEN NONLINEAR PROGRAMMING ALGORITHM HAS
C                         FINISHED PROCESSING
C                     = 1 WHEN RETURNING FOR AN EVALUATION 
C        IFERR        = 0 WHEN FUNCTION IS EVALUATED
C                     = 1 WHEN FUNCTION EVALUATION IS IMPOSSIBLE
C                     = -100 WHEN MAXIMUM NUMBER OF EVALUATIONS EXCEEDED
C        FZMODE       FEASIBILITY MODE FLAG
C        RELAX        RELAXATION MODE FLAG
C        ITERM        TERMINATION FLAG
C
C                     ---------------------------------------------------
C                     = -13  BARKKT ERROR OTHER THAN BELOW
C                     = -12  INTERNAL ERROR DURING MULTIFRONTAL SYMBOLIC FAC.
C                     = -11  INPUT ERROR TO MULTIFRONTAL OR LINE PLOT
C                     = -10  EXCESSIVE CONDITION NUMBER OR BARKKT FAILED AND
C                            HESSIAN HAS ALREADY BEEN RESET TO IDENTITY
C                            (CONSTRAINTS MAY BE INCONSISTENT)
C                     = -8   REAL HOLD ARRAY TOO SMALL
C                     = -6   I/O ERROR DURING MULTIFRONTAL
C                     = -4   MAXIMUM NUMBER OF FUNCTION ERRORS
C                     ---------------------------------------------------
C                     = 0    NO TERMINATION
C                     = 1    NORMAL TERMINATION
C                     = 2    NORMAL TERMINATION ON CENTRAL PATH
C                     ---------------------------------------------------
C                     = 3    SMALL STEPS
C                     = 4    MAXIMUM NUMBER OF ITERATIONS
C                     = 5    FUNCTION ERROR ON GRADIENT EVALUATION
C                     = 6    MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C                     = 7    UPHILL DIRECTION DERIVATIVE OR MAX NO. OF STEPS
C                     = 8    INFEASIBLE CONSTRAINTS
C                     = 9    MULTIPLIER CALCULATION FAILED 
C                     ---------------------------------------------------
C
C        IQPTRM       QP TERMINATION FLAG
C
C                     ---------------------------------------------------
C                     = -13  BARKKT ERROR OTHER THAN BELOW
C                     = -12  INTERNAL ERROR DURING MULTIFRONTAL SYMBOLIC FAC.
C                     = -11  MULTIFRONTAL INPUT ERROR
C                     = -10  EXCESSIVE CONDITION NUMBER OR BARKKT FAILED AND
C                            HESSIAN REST.  CONSTRAINTS MAY BE INCONSISTENT
C                     = -8   REAL HOLD ARRAY TOO SMALL
C                     = -6   I/O ERROR DURING MULTIFRONTAL
C                     = -4   MAXIMUM NUMBER OF FUNCTION ERRORS
C                     ---------------------------------------------------
C                     = 0    NO TERMINATION
C                     = 1    NORMAL TERMINATION
C                     = 2    NORMAL TERMINATION ON QP CENTRAL PATH
C                     ---------------------------------------------------
C                     = 3    SMALL STEPS
C                     = 4    MAXIMUM NUMBER OF ITERATIONS
C                     = 6    MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C                     = 7    UPHILL DIRECTION DERIVATIVE OR MAX NO. OF STEPS
C                     = 8    INFEASIBLE CONSTRAINTS
C                     = 9    RANK DEFICIENT CONSTRAINTS
C                     ---------------------------------------------------
C
      COMMON /ITEREF/ MAXREF,MAXRFN,IREFIN
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      LOGICAL LINEAR,BADSTP,GOODST,FZMODE,RELAX,REEDOO,SAMSCL,USEW
      LOGICAL GAUSNW,OLDSTP,GOODGN,OLDNWT,OLDGUS
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0,
     $    TEN=1.0D1,ONEEP2=1.0D2,ONEEM1=1.0D-1,ONEEM4=1.0D-4,
     $    POINT5=5.0D-1)
      PARAMETER (ONQRTR = ONE/FOUR, THQRTR = THREE/FOUR)
C
      DIMENSION YBAR(NVAR),YVEC(NVAR),DELYVC(NVAR),GVEC(NVAR),
     $    GRADL(NVAR),DELZVC(MAXRES),RMAT(NONZR),IROWR(NONZR),
     $    JCOLR(NVAR+1),CBAR(MAXCON),CVEC(MAXCON),ETABAR(MAXCON),
     $    ETAVEC(MAXCON),DELETA(MAXCON),FHFLTR(MXFLTR*5),
     $    CMAT(NONZC),IROWC(NONZC),JCOLC(NVAR+1),BBAR(MAXBND),
     $    BVEC(MAXBND),VLAMBR(MAXBND),VLAMDA(MAXBND),DELLAM(MAXBND),
     $    BMAT(NONZB),IROWB(NONZB),JCOLB(NVAR+1),
     $    WMAT(NONZW),IROWW(NONZW),JCOLW(NVAR+1),
     $    RWORK(LNRWRK),SCRTCH(LNSCRT),ISCRTC(LNISCR)
C
      DIMENSION IREVRS(5)
C
      LOGICAL HSNSTR,REFPNT
      LOGICAL HESEVL,GRDEVL,RESET,GOODGR
      LOGICAL OLDSET
      CHARACTER(LEN=6)  CHITER
C
      PARAMETER (MAXLVN=15)
      COMMON /LEVNP1/ LEVNIT
      COMMON /LEVNP2/ LEVNDY
      CHARACTER(LEN=60) LEVNDY(MAXLVN)
      CHARACTER(LEN=60) BLANK,CHWORK
      COMMON /CONDLB/ CNDLBL
      CHARACTER(LEN=17) CNDLBL
      CHARACTER(LEN=16) LFTVAL,RITVAL
C
C         COMMONS FOR ITERATION LOG QUANTITIES
C
C     ------------------------------------------------------------------
      INCLUDE  '../commons/NLPOUT.CMN'
C     ------------------------------------------------------------------
C
      LOGICAL CLOSE,NEWMLT
C
      DATA BLANK(1:60) / ' '/
C
C     ******************************************************************
C
      IF (IMAX(4,IREVRS,1).GE.1) THEN
        IF (IFESR.EQ.1) THEN
          GO TO 501
        ELSEIF (IFESR.EQ.2) THEN
          GO TO 502
        ELSEIF (IFESR.EQ.3) THEN
          GO TO 503
        ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C         ALGORITHM INITIALIZATION
C
C         LOGICAL QUANTITIES
C
      GAUSNW = NRES.GT.0.AND.(NEWTON.EQ.0.OR.NEWTON.EQ.2)
      USEW = .NOT.GAUSNW
      OLDSTP = GAUSNW
      OLDNWT = .NOT.GAUSNW
      OLDGUS = GAUSNW
      SAMSCL = .FALSE.
      NEWMLT = .FALSE.
      LINEAR = .FALSE.
      REFPNT = .FALSE.
      RESET = .FALSE.
      OLDSET = RESET
C
C         SET GRADIENT AND HESSIAN EVALUATION FLAGS
C
      GRDEVL = .TRUE.
      HESEVL = .TRUE.
C
C         SET HESSIAN INITIALIZATION FLAG 
C
      IF(IHESHN.LT.0) THEN
        HSNSTR = .TRUE.
        IHESHN = -IHESHN
      ELSE
        HSNSTR = .FALSE.
      ENDIF
C
C         SET GOOD GRADIENT FLAG BASED ON HESSIAN OPTION
C
      IF(IHESHN.EQ.0) THEN
        GOODGR = .TRUE.
      ELSE
        GOODGR = .FALSE.
      ENDIF
C
C         INTEGER QUANTITIES
C
      IT = 1
      IT1 = 1
      ITERM = 0
      IT1OLD = 1
      ITSAME = 0
      ITLEVN = 1
      LNFLTR = 0
      NEQNS = NVAR + NRES + MSUBE + MSUBB
      NSLK = NVAR - NFREE
      MEQUAL = MSUBE - MINEQL
C
C         REAL QUANTITIES
C
      SMALLN = ONEEP2*ZEROMN
      FBAR = FOBJ
      FMIN = FOBJ
      IF(FZMODE) FMIN = ZERO
      ACTRED = BIGNUM
      ALFA = ONE
      CNDNUM = ONE
      DIAGNL = ZERO
      EIGMAX = ZERO
      EIGMIN = ZERO
      PGRATE = POINT5
      ETAOLD = -ONE
      ETARAT = ZERO
      CVECMX = DAMAX(MSUBE,CVEC,1)
      CMPLMT = COMPLM(MSUBB,PENMU,BVEC,VLAMDA)
      BIGFLT = MAX(CVECMX,CMPLMT)
      CMAX = MAX(BIGCON,BIGFLT+CONTOL)
      REDFAC = .9D0
C
C         INITIALIZE QP BARRIER PARAMETER TO NLP BARRIER PARAMETER
C
      QPENMU = PENMU
C
C             ALLOCATE STORAGE FROM THE REAL SCRATCH ARRAY
C     ---ALLOCATE THE ARRAYS (I.E. CONSTRUCT THE POINTERS)
C
      MXQPFL = IMAXMU*IT1MAX
      LCAZER = 1
      LCDLPV = LCAZER + NVAR
      LCPBAR = LCDLPV + NVAR
      LCCZER = LCPBAR + NVAR
      LCETAZ = LCCZER + MAXCON
      LCQPFL = LCETAZ + MAXCON
      LCBZER = LCQPFL + MXQPFL*5
      LCVLMZ = LCBZER + MAXBND
      LCSKTR = LCVLMZ + MAXBND
      LNSKTR = 4*NEQNS + 2
      LCRSCR = LCSKTR + LNSKTR
      NCRSCR = LNSCRT - LCRSCR + 1
C
      IF(NCRSCR.LT.0) THEN
        PRINT *,'REAL SCRATCH ARRAY TOO SMALL; NEED =', LCRSCR - 1
        STOP
      ENDIF
C
C             ALLOCATE STORAGE FROM THE INTEGER SCRATCH ARRAY
C     ---ALLOCATE THE ARRAYS (I.E. CONSTRUCT THE POINTERS)
C
      LCITMP = 1
      LCNZRB = LCITMP + NEQNS
      LCISCR = LCNZRB + MAXBND
      NCISCR = LNISCR - LCISCR +1
      IF(NCISCR.LT.0) THEN
        PRINT *,'INTEGER SCRATCH ARRAY TOO SMALL; NEED =', LCISCR-1
        STOP
      ENDIF
C
C             INITIALIZE TERSE OUTPUT IF NECESSARY
C
      IF(FZMODE) THEN
        NPRFLG = -1
      ELSE
        NPRFLG = -2
      ENDIF
C
      CNDLBL = 'Cond(K)..........'
C
C ----------------------------------------------------------------------
C
C         BEGIN ITERATION
C
 110  CONTINUE
C
C ----------------------------------------------------------------------
C
C         COMPUTE GRADIENT INFORMATION
C
      IREVRS(1) = 0
C
C         SET GRADIENT CALL FLAG BASED ON THE FLAG GOODGR.  IF GOODGR 
C         IS TRUE REQUEST GOOD (I.E CENTRAL DIFFERENCE) GRADIENTS.  IF
C         GOODGR IS FALSE REQUEST CHEAP (I.E. FORWARD DIFFERENCE) GRADIENTS.
C         SET HESSIAN CALL FLAG TO COMPUTE THE HESSIAN DIAGONAL (OR GET
C         READY FOR A SUBSEQUENT FULL HESSIAN CALL).
C
      IF(.NOT.GRDEVL) GO TO 501
C
      IF(GOODGR) THEN
        IREVRS(2) = 2
        IREVRS(5) = 3
      ELSEIF(HSNSTR.AND.IT.EQ.1) THEN
        IREVRS(2) = 2
        IREVRS(5) = 5
      ELSE
        IREVRS(2) = 1
      ENDIF
C
      IF(HESEVL) THEN
        IREVRS(3) = 1
      ELSE
        IREVRS(3) = 0
      ENDIF
C
      IF(RESET) IREVRS(5) = 2
C
      IF (IOFLAG.GE.10) THEN
        IREVRS(4) = 1
      ELSE
        IREVRS(4) = 0
      ENDIF
      IFESR = 1
C
C         RETURN FOR FUNCTION AND GRADIENT EVALUATION
C
      RETURN
C
C         FUNCTION EVALUATION RETURN POINT
C
 501  CONTINUE
C
      irevrs(1:4) = 0
C
      IFUNER = IFERR
C
      IF(IFERR.EQ.-100) THEN
C
C         MAX. NO. OF FUNCTION EVALS.
C
        ITERM = 6
        GO TO 200
C
      ENDIF
C
C ----------------------------------------------------------------------
C
C         CONVERGENCE TESTS
C
C         ---FIRST COMPUTE PRIMAL-DUAL NECESSARY CONDITION QUANTITIES
C
      CALL MPDNEC(IFUNER,CVEC,MSUBE,MAXCON,ETAVEC,CMAT,IROWC,JCOLC,
     $    NONZC,BVEC,MSUBB,MAXBND,VLAMDA,BMAT,IROWB,JCOLB,NONZB,
     $    GVEC,NVAR,PENMU,SCRTCH,LNSCRT,GRADL,GRDLNM,CLAMNM,ERPTHB,
     $    DELFNM,ERREQL,ETANRM)
C
      IF(IFUNER.EQ.0) THEN
C
        ERRKTC = MAX(GRDLNM/MAX(ONE,DELFNM,CLAMNM),ERREQL,ERPTHB)
C
        IF(FZMODE) THEN
C
C         COMPUTE INFINITY NORM OF (CMAT**T)
C         FIRST COMPUTE ABSOLUTE COLUMN SUMS OF CMAT STORED IN SCRTCH
C
          CALL MRNSPR(2,NRES,NVAR,RMAT,IROWR,JCOLR,SCRTCH)
C
C         COMPUTE MAX ABSOLUTE SUMS OVER ALL NRES ROWS
C
          GRADNM = DAMAX(NVAR,GRADL,1)
          CTNORM = DAMAX(NRES,SCRTCH,1)
          IF(MSUBB.EQ.0) THEN
            EPSG = MAX(ZEROMN,CTNORM*CONTOL)
            EPSOBJ = POINT5*DBLE(NRES)*CONTOL**2
            EPSB = CONTOL
            EPSC = CONTOL
          ELSEIF(IT.EQ.1) THEN
            EPSG = MAX(ZEROMN,CTNORM*CONTOL)
            EPSOBJ = POINT5*DBLE(NRES)*CONTOL**2
            REDTOL = ONEEM1
            EPSG = MAX(EPSG,MIN(ONE,REDTOL*GRADNM))
            EPSOBJ = MAX(EPSOBJ,REDTOL*FOBJ)
            EPSB = MAX(CONTOL,MIN(ONE,REDTOL*ERRKTC))
            EPSC = MAX(CONTOL,MIN(ONE,REDTOL*ERRKTC))
          ENDIF
C
        ELSE
C
          EPSG = PGDTOL
          EPSOBJ = OBJTOL
          GRADNM = GRDLNM
          EPSB = CONTOL
          EPSC = CONTOL
C
        ENDIF
C
      ENDIF
C
C         ---PERFORM TESTS
C
      CALL BCNVRG(GRADNM,EPSG,CLAMNM,ERREQL,EPSC,ERPTHB,EPSB,
     $    PENMU,CNDNUM,ETARAT,DELYVC,YVEC,DELETA,ETAVEC,DELLAM,VLAMDA,
     $    ALFA,DELFNM,FOBJ,FMIN,EPSOBJ,NVAR,MSUBE,MSUBB,MAXCON,
     $    MAXBND,IT,NITMAX,NITMIN,IFUNER,IOFLAG,IPUNLP,FZMODE,RELAX,
     $    ITERM)
C
C         IF GRADIENTS ARE SMALL AND POSSIBLY INACCURATE REDO THE
C         GRADIENT EVALUATION WITH ACCURATE GRADIENTS
C
      IF(GOODGR.AND.IT.LE.9999) THEN
        CHITER = '(    )'
        WRITE(CHITER(2:5),'(I4)') IT
      ELSEIF(HSNSTR.AND.IT.EQ.1.AND.IT.LE.9999) THEN
        CHITER = '(    )'
        WRITE(CHITER(2:5),'(I4)') IT
      ELSE
        CHITER = '      '
        WRITE(CHITER(1:6),'(I6)') IT
        IF(GRDLNM.LT.PGDTOL) THEN
          GOODGR = .TRUE.
          GO TO 110
        ENDIF
      ENDIF
      CALL HHSDEL(CHITER,' ',' ',6,NDEL,IERDEL)
C
C         IF A FUNCTION ERROR OCCURED DURING GRADIENT EVALUATION TERMINATE
C
      IF(ITERM.EQ.5) GO TO 200
C
C         TRANSFER TO TERMINATION PROCEDURE BLOCK IF CONVERGENCE
C         TESTS ARE SATISFIED
C
      IF((ITERM.NE.0.AND.ITERM.NE.2) .OR. (ITERM.EQ.2.AND.FZMODE)) THEN
C
C         ITERATION PRINT AT FINAL POINT
C
        IF(IOFLAG.GE.10) THEN
          CALL HHCPRS(CHITER,' ','-')
          WRITE(IPUNLP,1001) CHITER
          WRITE(LFTVAL,'(SP,1PG14.6)') FOBJ + PENRHO*SMTUVW
          WRITE(RITVAL,'(SP,1PG14.6)') GRDLNM
          CALL DBLSTR('Objective Function',LFTVAL,
     $               'Gradient Of Lagrangian',RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') FMRT
          WRITE(RITVAL,'(SP,1PG14.6)') ERREQL
          CALL DBLSTR('Log-Barrier Function',LFTVAL,
     $               'Equality Error',RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') PENMU
          WRITE(RITVAL,'(SP,1PG14.6)') ERPTHB
          CALL DBLSTR('Barrier Parameter',LFTVAL,
     $               'Complementarity',RITVAL,IPUNLP)
          IF(IOFLAG.GE.20) THEN
            WRITE(LFTVAL,'(SP,1PG14.6)') DIAGNL
            WRITE(RITVAL,'(SP,1PG14.6)') CNDNUM
            CALL DBLSTR('Levenberg Parameter',LFTVAL,
     $               CNDLBL,RITVAL,IPUNLP)
            WRITE(LFTVAL,'(SP,1PG14.6)') EIGMIN
            WRITE(RITVAL,'(SP,1PG14.6)') EIGMAX
            CALL DBLSTR('Min. Eigenvalue',LFTVAL,
     $               'Max. Eigenvalue',RITVAL,IPUNLP)
          ENDIF
          CALL FLTRPR('   ',PENMU,PENRHO,FHFLTR,LNFLTR,MXFLTR)
        ENDIF
        IF(IOFLAG.GT.0.AND.IOFLAG.LT.10) THEN
          CALL BTURSE(NPRFLG,CHITER,LNFLTR,CNDNUM,ALFA,DIAGNL,ERPTHB,
     $    PENMU,ERREQL,GRDLNM,FOBJ,FMRT)
        ENDIF
C
C       ----------------------------------------------------------------
C
C         SAVE ITERATION LOG QUANTITIES
C
        ITB    = IT    
C
C         ALGORITHM KEY INDICATOR
C         KEYMB = 0    ----FEASIBILITY PHASE
C         KEYMB = 1    ----OPTIMIZATION PHASE
C
        KEYMB  = 1
        BFOBJ  = FOBJ + PENRHO*SMTUVW 
        BGRDLN = GRDLNM
        BFMRT  = FMRT  
        BERREQ = ERREQL
        BPENMU = PENMU
        BERPTH = ERPTHB
        BDIAGN = DIAGNL
        BCNDNU = CNDNUM
        BEIGMI = EIGMIN
        BEIGMA = EIGMAX
C
C       ----------------------------------------------------------------
C
        GO TO 200
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         COMPUTE FINITE DIFFERENCE HESSIAN MATRIX
C
      IREVRS(1) = 0
      IREVRS(2) = 0
      IREVRS(4) = 0
C
C        SET HESSIAN CALL FLAG
C
      IREVRS(3) = 2
C
      IFESR = 2
C
C         IF HESSIAN IS NOT NEEDED GO ON
C
      IF(.NOT.HESEVL) GO TO 502
C
C         RETURN FOR HESSIAN EVALUATION
C
      RETURN
C
C         FUNCTION EVALUATION RETURN POINT
C
 502  CONTINUE
C
      IREVRS(1:4) = 0
      IREVRS(5) = 1
C
      IF(IFERR.EQ.-100) THEN
C
C         MAX. NO. OF FUNCTION EVALS.
C
        ITERM = 6
        GO TO 200
C
      ENDIF
C
      IFUNER = IFERR
      IF(IFUNER.NE.0.AND.IFUNER.NE.-101) THEN
        ITERM = 5
        GO TO 200
      ENDIF
C
C ----------------------------------------------------------------------
C
      LEVNIT = 0
      IF(IT.EQ.ITLEVN) THEN
C
        PGNOLD = GRDLNM
C
      ELSE
C
C         ---TRUST REGION UPDATE FOR LEVENBERG CORRECTION 
C
        ITLEVN = ITLEVN + 1
C
C         COMPARE PREDICTED REDUCTION TO ACTUAL REDUCTION
C
        BADSTP = ABS(ACTRED-PRERED).GE.THQRTR*ABS(PRERED)
        BADSTP = BADSTP.AND.(IT1OLD.GT.2)
        IF(ACTRED.GT.ZERO) THEN
          GOODST = ACTRED.GE.THQRTR*PRERED
        ELSE
          GOODST = ABS(ACTRED-PRERED).LE.ONQRTR*ABS(PRERED)
        ENDIF
C
C          ---NOTE:  LEVENBERG CORRECTION IS INCLUDED EVEN WHEN WMAT IS 
C          NOT USED, TO COMPENSATE FOR RANK DEFICIENCIES.
C
        IF(BADSTP) THEN
C
C         BAD STEP WITH CURRENT TRUST RADIUS -- REDUCE THE SIZE
C         OF THE TRUST RADIUS (INCREASE HESSIAN DIAGONAL) 
C         FOR SUBSEQUENT STEPS
C
          DIAGOL = DIAGNL
          DIAGNL = TWO*MAX(DIAGNL,ZEROOT)
          DIAGNL = MIN(DIAGNL,ONE)
C
          IF(IOFLAG.GE.20.AND.LEVNIT.LT.MAXLVN.AND.DIAGOL.LT.ONE) THEN
            LEVNIT = LEVNIT + 1
            LEVNDY(LEVNIT) = BLANK
            WRITE(LEVNDY(LEVNIT)(1:16),'(SP,1PG16.6)') DIAGNL
            LEVNDY(LEVNIT)(20:49) = 'Bad step...reduce trust radius'
          ENDIF
C 
        ELSEIF(GOODST) THEN
C
C         GOOD STEP WITH CURRENT TRUST RADIUS -- INCREASE THE
C         SIZE (REDUCE HESSIAN DIAGONAL) FOR SUBSEQUENT STEPS
C
          REDUCT = MIN(POINT5,GRDLNM/(PGNOLD+ZEROOT))
          DIAGOL = DIAGNL
          DIAGNW = DIAGNL*REDUCT
          IF(DIAGNW.GE.ZEROOT) THEN
            DIAGNL = DIAGNW
          ELSE
            DIAGNL = ZERO
          ENDIF
C
          IF(IOFLAG.GE.20.AND.LEVNIT.LT.MAXLVN.AND.DIAGOL.GT.ZERO) THEN
            LEVNIT = LEVNIT + 1
            LEVNDY(LEVNIT) = BLANK
            WRITE(LEVNDY(LEVNIT)(1:16),'(SP,1PG16.6)') DIAGNL
            LEVNDY(LEVNIT)(20:52) = 'Good step...increase trust radius'
          ENDIF
C
        ENDIF
C
        PGRATE = GRDLNM/(PGNOLD+ZEROOT)
        PGNOLD = GRDLNM
C
C           COMPARE PREDICTED AND ACTUAL REDUCTION FOR GAUSS-NEWTON MODEL
C
        IF(ACTRED.GT.ZERO) THEN
          GOODGN = ACTRED.GE.THQRTR*GNPRED
        ELSE
          GOODGN = ABS(ACTRED-GNPRED).LE.ONQRTR*ABS(GNPRED)
        ENDIF
C
        IF(OLDSTP) THEN
C
C           CURRENT STEP IS GAUSS-NEWTON
C
          IF(OLDNWT.AND.GOODST.AND.(.NOT.OLDGUS.OR..NOT.GOODGN)) THEN
C
C           TWO SUCCESSIVE GOOD NEWTON STEPS, AND AT LEAST ONE BAD 
C           GAUSS-NEWTON ---> SWITCH TO NEWTON 
C
            IF(NEWTON.NE.2) GAUSNW = .FALSE.
            IF(.NOT.GAUSNW.AND.IOFLAG.GE.20) WRITE(IPUNLP,1008)
          ELSEIF(DIAGNL.NE.ZERO) THEN
            IF(NEWTON.NE.2) GAUSNW = .FALSE.
            IF(.NOT.GAUSNW.AND.IOFLAG.GE.20) WRITE(IPUNLP,1008)
          ELSE
            GAUSNW = .TRUE.
          ENDIF
C
        ELSE
C
C           CURRENT STEP IS NEWTON
C
          IF(OLDGUS.AND.GOODGN.AND.(.NOT.OLDNWT.OR..NOT.GOODST)) THEN
C
C           TWO SUCCESSIVE GOOD GAUSS-NEWTON STEPS, AND AT LEAST ONE BAD 
C           NEWTON ---> SWITCH TO GAUSS-NEWTON 
C
            GAUSNW = NRES.GT.0.AND.(NEWTON.EQ.0.OR.NEWTON.EQ.2)
            IF(GAUSNW.AND.IOFLAG.GE.20) WRITE(IPUNLP,1007)
          ELSE
            GAUSNW = .FALSE.
          ENDIF
C
        ENDIF
C
        OLDNWT = GOODST
        OLDGUS = GOODGN
        OLDSTP = GAUSNW
C
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         ---BARRIER PARAMETER UPDATE
C
      IF(IT.GT.1) THEN
C
        GRDTST = PGDTOL*MAX(ONE,DELFNM)*(ONE + ABS(FOBJ))
C
        IF(MSUBB.GT.0) THEN
C
          TOLPTH = MIN(TEN*PENMU,PTHTOL)
          PENOLD = PENMU
C
C         FIRST CHECK TO SEE THAT LEVENBERG MODIFICATION IS SMALL AND/OR
C         UNCONSTRAINED MINIMIZATION WITH CURRENT BARRIER PARAMETER 
C         IS NEAR COMPLETION (I.E. GRDLNM IS SMALL)
C
          ERRCON = MAX(ERREQL,ERPTHB)
          CLOSE = ALFA.GE.ONEEM1.AND.(ERRKTC.LT.TOLPTH.OR.
     $    ERRCON.LT.SMALLN.OR.GRDLNM.LT.GRDTST.OR.DIAGNL.LE.TEN*ZEROOT)
C
          IF(CLOSE.OR.ITSAME.EQ.IMAXMU) THEN
C
C           RESET THE TERMINATION FLAG IF NECESSARY
C
            IF(ITERM.EQ.2) ITERM = 0
C
            IF(ERRKTC.LT.TOLPTH.OR.ERRCON.LT.SMALLN) THEN
C
C               CURRENT POINT IS "NEAR" THE CENTRAL PATH
C               ---BARRIER PARAMETER REDUCTION
C
              IF(PENMU.LT.ONEEM4) THEN
                HATMU = TEN*PENMU**2
              ELSE
                HATMU = PENMU/TEN
              ENDIF
C
C           COMPUTE DIVERGENCE RATE TEST QUANTITIES
C
              IF(ETAOLD.GT.ZEROMN) ETARAT = (ETANRM - ETAOLD)/ETAOLD
              ETAOLD = ETANRM
C
C           UPDATE THE BARRIER PARAMETER
C
              IF(PENMU.GT.ZEROOT) THEN
                PENRAT = HATMU/PENMU
              ELSE
                PENRAT = ZERO
              ENDIF
              PENMU = HATMU
C
C           REFRESH THE FILTER AND ADJUST MAX VIOLATION BOUND
C
              LNFLTR = 0
C
              BIGFLT = MAX( DAMAX(MSUBE,CVEC,1),
     $                      COMPLM(MSUBB,PENMU,BVEC,VLAMDA) )
              CMAX = MAX(CMAX,BIGFLT+CONTOL)
C
C           RESET THE LEVENBERG PARAMETER 
C
              ITLEVN = 1
              DIAGOL = DIAGNL
              DIAGNW = DIAGNL*PENRAT
              IF(DIAGNW.GE.ZEROOT) THEN
                DIAGNL = DIAGNW
              ELSE
                DIAGNL = ZERO
              ENDIF
              IF(IOFLAG.GE.20.AND.LEVNIT.LT.MAXLVN
     $          .AND.DIAGOL.NE.DIAGNL) THEN
                LEVNIT = LEVNIT + 1
                LEVNDY(LEVNIT) = BLANK
                WRITE(LEVNDY(LEVNIT)(1:16),'(SP,1PG16.6)') DIAGNL
                LEVNDY(LEVNIT)(20:44) = 'Update Barrier Parameter'
              ENDIF
C
            ELSEIF(ITSAME.EQ.IMAXMU) THEN
C
              IF(PENMU.GT.ONE) THEN
                HATMU = ONEEM1*PENMU
              ELSE
                HATMU = REDFAC*PENMU
              ENDIF
C
C             COMPUTE DIVERGENCE RATE TEST QUANTITIES
C
              ETARAT = ZERO
              ETAOLD = ETANRM
C
C             UPDATE THE BARRIER PARAMETER
C
              PENMU = HATMU
C
C             REFRESH THE FILTER AND ADJUST MAX VIOLATION BOUND
C
              LNFLTR = 0
C
              BIGFLT = MAX( DAMAX(MSUBE,CVEC,1),
     $                      COMPLM(MSUBB,PENMU,BVEC,VLAMDA) )
              CMAX = MAX(CMAX,BIGFLT+CONTOL)
C
C             RESET THE LEVENBERG PARAMETER
C
              ITLEVN = 1
              DIAGOL = DIAGNL
              DIAGNW = DIAGNL*REDFAC
              IF(DIAGNW.GE.ZEROOT.AND.PENMU.GT.ZEROOT) THEN
                DIAGNL = DIAGNW
              ELSE
                DIAGNL = ZERO
              ENDIF
              IF(IOFLAG.GE.20.AND.LEVNIT.LT.MAXLVN
     $          .AND.DIAGOL.NE.DIAGNL) THEN
                LEVNIT = LEVNIT + 1
                LEVNDY(LEVNIT) = BLANK
                WRITE(LEVNDY(LEVNIT)(1:16),'(SP,1PG16.6)') DIAGNL
                LEVNDY(LEVNIT)(20:44) = 'Update Barrier Parameter'
              ENDIF
C
            ENDIF
C
          ENDIF
C
          IF(PENMU.EQ.PENOLD) THEN
C
C           ITERATION COUNT WITH UNCHANGED PENMU
C
            ITSAME = ITSAME + 1
C
          ELSE
C
C             RESET MU ITERATION COUNTER
C
            ITSAME = 0
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         IF CURRENT MULTIPLIER ESTIMATES ARE EXCESSIVELY LARGE
C         TRY TO RECOMPUTE NEW ESTIMATES.
C
      IF(ETAMAG.GT.PENRHO) NEWMLT = .TRUE.
C
C         FOR BARRIER SQP MODE, RECOMPUTE MULTIPLIERS AND BARRIER PARAMETER
C         AFTER EACH QP SUBPROBLEM (UNTIL NEAR THE SOLUTION)  
C
      IF(IT.GT.1.AND.MXQPIT.GT.1.AND.PENMU.GT.ONEEM4) NEWMLT = .TRUE.
C
      IF(NEWMLT) THEN
C
C             COMPUTE NEW FIRST ORDER MULTIPLIERS
C
        NEWMLT = .FALSE.
        REEDOO = .TRUE.
C
        BSTMU = PENMU
        ETABAR(1:MSUBE) = ETAVEC(1:MSUBE)
        VLAMBR(1:MSUBB) = VLAMDA(1:MSUBB)
        MUOPTN = 2
C
C            COMPUTE OPTIMUM MU 
C
 120    CONTINUE
C
        NCALL = 0
        CALL LSQMLT(MUOPTN,NCALL,CVEC,MSUBE,MAXCON,ETAVEC,CMAT,
     $           IROWC,JCOLC,NONZC,BVEC,MSUBB,MAXBND,VLAMDA,BMAT,IROWB,
     $           JCOLB,NONZB,GVEC,NVAR,PENMU,SCRTCH,RWORK,     
     $           LNRWRK,IPUNLP,IOFLAG,
     $           CNDNUM,GRDLNM,IERLSM,NEEDED)
C      
        IF(IERLSM.NE.0.AND.IERLSM.NE.-999) THEN
          PRINT *,'IERLSM =',IERLSM
          ITERM = 9
          GO TO 200
        ENDIF
C
C           COMPUTE ERROR IN KT CONDITIONS
C
        PMUCVC = ERREQL/(TEN*1.1D0)
C
C         CHECK THE MINIMUM VALUES FOR THE MULTIPLIERS AND BARRIER
C         PARAMETER
C
        IF(PENMU.LT.BSTMU.AND.MSUBB.GT.0.AND.REEDOO) THEN
          REEDOO = .FALSE.
          IF(PENMU.GT.PMUCVC) THEN
            PENMU = MAX(ZEROOT,PENMU)
          ELSE
            PENMU = MAX(ZEROOT,ONEEM1*BSTMU,MIN(PMUCVC,BSTMU))
          ENDIF
          MUOPTN = 1
          GO TO 120
        ELSEIF(PENMU.GE.BSTMU) THEN
          PENMU = BSTMU
        ENDIF
C
        ETAMAG = ZERO
        DO I=1,MSUBE
          ETABAR(I) = ETAVEC(I)
          ETAMAG = MAX(ETAMAG,ABS(ETAVEC(I)))
        ENDDO
C
        VLAMAG = ZERO
        DO I=1,MSUBB
          VLAMBR(I) = VLAMDA(I)
          VLAMAG = MAX(VLAMAG,ABS(VLAMDA(I)))
        ENDDO
C
      ENDIF
C
      IF(MIN(ETAMAG,VLAMAG,DELFNM).GT.PENRHO) THEN
        ITERM = 9
        GO TO 200
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C     
C         SOLVE WITHOUT ITERATIVE REFINEMENT UNLESS NEAR THE ANSWER
C
      IF(GRDLNM.LT.GRDTST.AND.ERREQL.LT.CONTOL) THEN
        IREFIN = 2
      ELSE
        IREFIN = 0
      ENDIF
C
C         INITIALIZE QP BARRIER PARAMETER TO NLP BARRIER PARAMETER
C
      QPENMU = PENMU
C
      IF(MXQPIT.EQ.1) THEN
C
        LCRSCL = 1
        LCCSCL = LCRSCL + NEQNS
        LCRTKT = LCCSCL + NEQNS
        LCSNKT = LCRTKT + NEQNS
        USEW = .NOT.GAUSNW
        LINEAR = GAUSNW
C
C         COMPUTE SEARCH DIRECTION
C
        CALL BRSRCH( NVAR, NEQNS, NRES, MSUBE, MSUBB, NSLK, 
     $       DIAGNL, LINEAR, USEW, IPUNLP, IOFLAG , NONZR, IROWR, 
     $       JCOLR, RMAT, NONZB, IROWB, JCOLB, BMAT, NONZC,
     $       IROWC, JCOLC, CMAT, NONZW, IROWW, JCOLW, WMAT,  QPENMU,
     $       BVEC, MAXBND, CVEC, MAXCON, GRADL, GVEC, ETAVEC, VLAMDA,  
     $       LNRWRK, RWORK,ISCRTC(LCITMP),ISCRTC(LCNZRB),
     $       SCRTCH(LCRSCL), SCRTCH(LCCSCL),SCRTCH(LCRTKT), 
     $       SCRTCH(LCSNKT), LNSYMB, CNDNUM, NEEDED, DELYVC,
     $       DELZVC, MAXRES, DELETA, DELLAM, EIGMIN,EIGMAX,  
     $       ADDHSS, RESET,IQPTRM  ) 
C
      ELSE
C
C         COMPUTE SEARCH DIRECTION
C
        CALL BRLPQP(            FBAR     ,FOBJ     ,SCRTCH(LCPBAR)         
     $     ,DELYVC   ,SCRTCH(LCDLPV)     ,GVEC     ,GRADL          
     $     ,SCRTCH(LCAZER)     ,NVAR     ,NFREE    ,MAXRES   ,NRES
     $     ,DELZVC   ,RMAT     ,IROWR    ,JCOLR    ,NONZR    ,CBAR          
     $     ,CVEC     ,SCRTCH(LCCZER)     ,ETABAR   ,ETAVEC     
     $     ,SCRTCH(LCETAZ)     ,DELETA   ,SCRTCH(LCQPFL)     ,MXQPFL         
     $     ,MSUBE    ,MINEQL   ,MAXCON   ,CMAT     ,IROWC    ,JCOLC                    
     $     ,NONZC    ,BBAR     ,BVEC     ,SCRTCH(LCBZER)     ,VLAMBR    
     $     ,VLAMDA   ,SCRTCH(LCVLMZ)     ,DELLAM   ,MSUBB    ,MAXBND             
     $     ,BMAT     ,IROWB    ,JCOLB    ,NONZB    ,WMAT     ,IROWW     
     $     ,JCOLW    ,NONZW    ,QPENMU   ,PENRHO   ,SMTUVW           
     $     ,RWORK    ,LNRWRK   ,SCRTCH(LCSKTR)     ,LNSKTR      
     $     ,ISCRTC   ,LNISCR   ,NEEDED   ,EIGMIN   ,EIGMAX   ,ADDHSS      
     $     ,DIAGNL   ,RESET    ,IQPTRM)      
C
      ENDIF
C
      IF(IQPTRM.LT.0) THEN
        ITERM = IQPTRM
        GO TO 200
      ENDIF
C
C         TRANSFER TO TERMINATION PROCEDURE BLOCK IF SEARCH
C         DIRECTION CALCULATION FAILED
C
      RESNOR = DAMAX(NVAR,DELYVC,1) + DAMAX(MSUBE,DELETA,1) 
     $         + DAMAX(MSUBB,DELLAM,1)
C
C         CHECK FOR SMALL STEP
C
      IF(RESNOR.LE.SMALLN) THEN
        ITERM = 3
        GO TO 200
      ENDIF
C
      IF(IQPTRM.EQ.9.OR.IQPTRM.EQ.8) THEN
        PRINT *,'BARRIER QP FAILED--TRY STEP'
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C
C         INITIALIZE LOG-BARRIER FUNCTION FOR UNIVARIATE SEARCH
C
      CALL BMERUT(CVEC,ETAVEC,BVEC,VLAMDA,FOBJ,SMTUVW,PENMU,
     $    PENRHO,MSUBE,MSUBB,SMLNBI,FMRT,CMRT,XLAGR0,C0TC0)
C
      CBTCB = C0TC0
C
C         RECOMPUTE COMPLEMENTARITY ERROR (WITH POSSIBLY NEW PENMU, VLAMDA)
C
      ERPTHB = ZERO
      DO I = 1,MSUBB
        ERPTHB = MAX(ERPTHB,ABS(BVEC(I)*VLAMDA(I)-PENMU))
      enddo
C
      ETAMAG = DAMAX(MSUBE,ETAVEC,1)
      VLAMAG = DAMAX(MSUBB,VLAMDA,1)
      IF(LNFLTR.EQ.0) THEN
        CALL FILTER(FOBJ,CMRT,CMAX,SMLNBI,SMTUVW,ETAMAG,PENMU,PENRHO,
     $    FHFLTR,LNFLTR,MXFLTR,NWFLTR,ISCRTC,LNISCR)
      ENDIF
C
C         ITERATION PRINT
C
      IF(IOFLAG.GE.10) THEN
        CALL HHCPRS(CHITER,' ','-')
        WRITE(IPUNLP,1001) CHITER
        WRITE(LFTVAL,'(SP,1PG14.6)') FOBJ + PENRHO*SMTUVW
        WRITE(RITVAL,'(SP,1PG14.6)') GRDLNM
        CALL DBLSTR('Objective Function',LFTVAL,
     $               'Gradient Of Lagrangian',RITVAL,IPUNLP)
        WRITE(LFTVAL,'(SP,1PG14.6)') FMRT
        WRITE(RITVAL,'(SP,1PG14.6)') ERREQL
        CALL DBLSTR('Log-Barrier Function',LFTVAL,
     $               'Equality Error',RITVAL,IPUNLP)
        WRITE(LFTVAL,'(SP,1PG14.6)') PENMU
        WRITE(RITVAL,'(SP,1PG14.6)') ERPTHB
        CALL DBLSTR('Barrier Parameter',LFTVAL,
     $               'Complementarity',RITVAL,IPUNLP)
        IF(IOFLAG.GE.20) THEN
          WRITE(LFTVAL,'(SP,1PG14.6)') DIAGNL
          WRITE(RITVAL,'(SP,1PG14.6)') CNDNUM
          CALL DBLSTR('Levenberg Parameter',LFTVAL,
     $               CNDLBL,RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') EIGMIN
          WRITE(RITVAL,'(SP,1PG14.6)') EIGMAX
          CALL DBLSTR('Min. Eigenvalue',LFTVAL,
     $               'Max. Eigenvalue',RITVAL,IPUNLP)
        ENDIF
      ENDIF
C
C       ----------------------------------------------------------------
C
C         SAVE ITERATION LOG QUANTITIES
C
      ITB    = IT    
C
C         ALGORITHM KEY INDICATOR
C         KEYMB = 0    ----FEASIBILITY PHASE
C         KEYMB = 1    ----OPTIMIZATION PHASE
C
      KEYMB  = 1
      BFOBJ  = FOBJ + PENRHO*SMTUVW 
      BGRDLN = GRDLNM
      BFMRT  = FMRT  
      BERREQ = ERREQL
      BPENMU = PENMU
      BERPTH = ERPTHB
      BDIAGN = DIAGNL
      BCNDNU = CNDNUM
      BEIGMI = EIGMIN
      BEIGMA = EIGMAX
C
C       ----------------------------------------------------------------
C
C
      IF(RESET) THEN
C
C         TIGHTEN THE SOUTHEAST CORNER ON THE FILTER
C
        SMLCON = SQRT(CONTOL)
        CMAX = MAX( SMLCON, DAMAX(MSUBE,CVEC,1),
     $                COMPLM(MSUBB,PENMU,BVEC,VLAMDA) )
        LNFLTR = 0
C        
      ENDIF
C
      IF(OLDSET.AND.RESET) ITERM = -10
      OLDSET = RESET
      IF(IOFLAG.GE.20.AND.LEVNIT.LT.MAXLVN.AND.LEVNIT.GT.0) THEN
        WRITE(IPUNLP,1004)
        CHWORK = BLANK
        CHWORK(4:16) = 'New Levenberg'
        CHWORK(20:25) = 'Reason'
        WRITE(IPUNLP,1005) CHWORK
        WRITE(IPUNLP,1005) (LEVNDY(II),II = 1,LEVNIT)
        WRITE(IPUNLP,1004)
      ENDIF
C
      IF(IQPTRM.LT.0) THEN
        ITERM = IQPTRM
        GO TO 200
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
 140  CONTINUE
C
C         LINE SEARCH.  
C
C
C         ---------------------------------------------------
C         COMPUTE PREDICTED MINIMUM, PREDICTED REDUCTION, AND
C         ESTIMATES OF THE FIRST STEPLENGTH
C         ---------------------------------------------------
C
      LCCDLY = 1
      LCBDLY = LCCDLY + MAXCON
      LCWDLY = LCBDLY + MAXBND
      CALL MSLOPE(CVEC,MSUBE,MAXCON,ETAVEC,CMAT,IROWC,JCOLC,
     $    NONZC,BVEC,MSUBB,MAXBND,VLAMDA,BMAT,IROWB,JCOLB,NONZB,     
     $    WMAT,IROWW,JCOLW,NONZW,GVEC,NVAR,NFREE,NRES,MAXRES,ADDHSS,
     $    PENMU,SCRTCH,SCRTCH(LCBDLY),SCRTCH(LCWDLY),DELYVC,DELZVC,
     $    DELETA,DELLAM,FMRT,SAMSCL,FMIN,SLOPEL,DELTAH,DELTAR,
     $    DELTAX,SLOPEB,ALFEST,ALFLIM,GAMEST)
C
C         PREDICTED REDUCTION 
C
      PRERED = -SLOPEL - POINT5*(DELTAH + DELTAR + ADDHSS*DELTAX)
C
C         GAUSS-NEWTON MODEL PREDICTION
C
      GNPRED = -SLOPEL - POINT5*(DELTAR + ADDHSS*DELTAX)
C
C         DEFINE INITIAL STEP LENGTH
C
      ALFA = MIN(ONE,ALFEST)
      GAMA = MIN(ONE,GAMEST)
      ALFA0 = ALFA
C
C         SET PRIMAL-DUAL SCALING FLAG
C
      SAMSCL = ALFA.EQ.GAMA
C
C         ITERATION PRINT
C
      IF(IOFLAG.GE.10.AND.IERLIN.NE.8) THEN
        CALL FLTRPR('   ',PENMU,PENRHO,FHFLTR,LNFLTR,MXFLTR)
        IF(RESET) WRITE(IPUNLP,1003)
      ENDIF
C
      IF(IOFLAG.GT.0.AND.IOFLAG.LT.10) THEN
        CALL BTURSE(NPRFLG,CHITER,LNFLTR,CNDNUM,ALFA,DIAGNL,ERPTHB,
     $    PENMU,ERREQL,GRDLNM,FOBJ,FMRT)
      ENDIF
C
      IFUNER = 0
C
 150  CONTINUE
C
      IF(IT.EQ.LYNPLT) THEN
C
        CALL PLTLYN(NVAR,MCON,ALFA,FBRMRT,FBAR,DELYVC,YVEC,YBAR,CBAR,
     $    VLAMBR,VECSBR,IREVRS(1),IERPLT)
C
        IF(IREVRS(1).EQ.0) THEN
          ITERM = -11
          GO TO 200
        ENDIF
C
      ELSE
C
        IF(MSUBB.EQ.0) THEN
          ALFLMI = TWO
          PENMUI = ZERO
        ELSE
          ALFLMI = ALFLIM 
          PENMUI = PENMU
        ENDIF
C
        CALL BLINES(FMRT,SLOPEB,FBRMRT,ALFA,ALFLMI,PENMUI,CBTCB,
     $    C0TC0,CBRMRT,IT1,IT1MAX,NWFLTR,IFUNER,IOFLIN,IPUNLP,
     $    IREVRS(1),IERLIN,SAMSCL)
C
      ENDIF
C
      IF(IREVRS(1).EQ.0) GO TO 190
C
C         COMPUTE LOCATION OF POINT 
C         YBAR = YVEC + ALFA*DELYVC 
C
      DO I = 1,NVAR
        YBAR(I) = YVEC(I) + ALFA*DELYVC(I)
      enddo
C
      RATIO = GAMA*(ALFA/ALFA0)
      DO I = 1,MSUBE
        ETABAR(I) = ETAVEC(I) + RATIO*DELETA(I)
      enddo
C
      DO I = 1,MSUBB
        VLAMBR(I) = VLAMDA(I) + RATIO*DELLAM(I)
        VLAMBR(I) = MAX(ZEROMN,VLAMBR(I))
      enddo
C
C         RETURN FOR FUNCTION EVALUATION
C
      IREVRS(4) = 0
      IF(IOFLAG.GE.10.AND.LYNPLT.NE.IT) IREVRS(4) = 1
      IFESR = 3
      RETURN
C
C         FUNCTION EVALUATION RETURN POINT
C
 503  CONTINUE
C
      IFUNER = IFERR
      IF(IFERR.EQ.0.OR.IFERR.EQ.-101) THEN
C
C         CONSTRUCT LOG-BARRIER FUNCTION FOR UNIVARIATE SEARCH
C
        CALL BMERUT(CBAR,ETABAR,BBAR,VLAMBR,FBAR,SMTUVW,PENMU,
     $    PENRHO,MSUBE,MSUBB,SMLNBI,FBRMRT,CBRMRT,XLAGRN,CBTCB)
C
C         CHECK CURRENT POINT AGAINST FILTER
C
        ETAMAG = DAMAX(MSUBE,ETAVEC,1)
        CALL FILTER(FBAR,CBRMRT,CMAX,SMLNBI,SMTUVW,ETAMAG,PENMU,PENRHO,
     $    FHFLTR,LNFLTR,MXFLTR,NWFLTR,ISCRTC,LNISCR)
        IF(NWFLTR.NE.0) IFUNER = 2
C
      ELSE
C
        FBRMRT = HDMCON(1)
C
      ENDIF
C
C         FUNCTION EVALUATED. REENTER ONE DIMENSIONAL SEARCH
C
      GO TO 150
C
C ----------------------------------------------------------------------
C
C         UPDATE INFORMATION
C
 190  CONTINUE
C
      IF(IERLIN.EQ.8) THEN
C
C         FILTER REJECTS POINT AND PRIMAL-DUAL STEP SCALING IS UNEQUAL
C         FORCE EQUAL STEP SCALING AND REPEAT 
C
        SAMSCL = .TRUE.
        IF(IOFLAG.GE.20) WRITE(IPUNLP,1006)
        GO TO 140
C
      ENDIF
C
C         SET ITERM FLAG BASED ON LINE SEARCH ERROR INDICATOR IERLIN
C
      IF(IERLIN.EQ.3) THEN
        ITERM = -4
      ELSEIF(IFERR.EQ.-1) THEN
        ITERM = 0
      ELSEIF(IFERR.EQ.-100) THEN
        ITERM = 6
      ELSEIF(IFERR.EQ.-101) THEN
        ITERM = 0
      ELSEIF(IERLIN.EQ.5) THEN
        ITERM = 7
      ELSEIF(IERLIN.EQ.2) THEN
        ITERM = 3
      ELSEIF(IERLIN.EQ.-1) THEN
        ITERM = 7
      ENDIF
C
C         COMPUTE ACTUAL REDUCTION
C
      ACTRED = XLAGR0 - XLAGRN
      IT1OLD = IT1
C
C         SET GRADIENT EVALUATION FLAG FOR NEXT ITERATION
C
      GRDEVL = .TRUE.
C
C         SET HESSIAN EVALUATION FLAG FOR NEXT ITERATION
C         THE OLD HESSIAN WILL BE USED IF THE ITERATION 
C         APPEARS TO BE PROGRESSING WELL--OTHERWISE A NEW HESSIAN
C         WILL BE COMPUTED.  GOOD PROGRESS IS INDICATED BY 
C         A) FULL LENGTH STEPS (ALFA = 1)
C         B) POSITIVE DEFINITE HESSIAN WITHOUT MODIFICATION (DIAGNL = 0)
C         C) REASONABLE REDUCTION IN THE PROJECTED GRADIENT (PGRATE .LE. 1/2)
C         D) HESSIAN HAS NOT BEEN RESET TO THE IDENTITY (RESET = .FALSE.)
C         E) A NEWTON METHOD HAS NOT BEEN SPECIFIED BY INPUT (NEWTON = 1)
C
      HESEVL = ALFA.LT.ONE.OR.DIAGNL.NE.ZERO.OR.PGRATE.GT.POINT5
     $        .OR.RESET.OR.NEWTON.EQ.1
C
C         CHECK FOR DISCONTINUOUS BEHAVIOR
C
      IF(ABS(ALFA)*RESNOR.LE.SMALLN.AND.ABS(ACTRED).GT.SMALLN
     $    .AND.IOFLAG.GE.10) WRITE(IPUNLP,1002)
C      
C         SET THE GRADIENT ACCURACY FOR THE NEXT ITERATION
C
      IF(.NOT.GOODGR) THEN
        SMALCH = MAX(OBJTOL,ZEROOT)
        IF(ABS(ACTRED).LT.SMALCH) THEN
          ACTRED = BIGNUM
          IF(ITERM.EQ.3) ITERM = 0
          GOODGR = .TRUE.
        ENDIF
        DELYNM = ZERO
        DO I=1,NVAR
          DELYNM = DELYNM + DELYVC(I)**2
        ENDDO
        DELYNM = SQRT(DELYNM)
        IF(ALFA*DELYNM.LT.SMALLN) THEN
          ALFA = TWO
          IF(ITERM.EQ.3) ITERM = 0
          GOODGR = .TRUE.
        ENDIF
      ENDIF
C
      IT = IT + 1
      IT1 = 1
      FOBJ = FBAR
      FMRT = FBRMRT
      REFPNT = .TRUE.
C
      SAMSCL = .FALSE.
C
      YVEC(1:NVAR) = YBAR(1:NVAR)
      DO I=1,MSUBE
        ETAVEC(I) = ETABAR(I)
        CVEC(I) = CBAR(I)
      ENDDO
      DO I=1,MSUBB
        VLAMDA(I) = VLAMBR(I)
        BVEC(I) = BBAR(I)
      ENDDO
C
      IF(ITERM.LE.-1.OR.ITERM.EQ.6) GO TO 200
C
C ----------------------------------------------------------------------
C
C         END OF ITERATION
C
      GO TO 110
C
C ----------------------------------------------------------------------
C
C         TERMINATION PROCEDURES
C
 200  CONTINUE
C
      FBAR = FOBJ
      YBAR(1:NVAR) = YVEC(1:NVAR)
C
      IF((NSLK-MINEQL).GT.0.AND..NOT.FZMODE) THEN
        TUVMAX = MAXVAL(YBAR(NFREE+MINEQL+1:NFREE+NSLK))
        IF((IQPTRM.EQ.1.OR.IQPTRM.EQ.2).AND.TUVMAX.GT.CONTOL) ITERM = 8
      ENDIF
C
C
 1001 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'-------------------------------
     $------ Iteration ',A6,'----------------------------------',
     $   T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....FUNCTION MAY BE DISCONTINU
     $OUS',T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....RESET HESSIAN AND LAGRANGE
     $ MULTIPLIERS',T106,'*')
 1004 FORMAT(T3,'*',T11,88('-'),T106,'*')
 1005 FORMAT(T3,'*',T11,'|',T20,A60,T98,'|',T106,'*')
 1006 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....SWITCH TO EQUAL PRIMAL-DUA
     $L STEP SCALING',T106,'*')
 1007 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....SWITCH TO GAUSS-NEWTON MET
     $HOD',T106,'*')
 1008 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....SWITCH TO NEWTON METHOD',
     $    T106,'*')
      RETURN
      END
