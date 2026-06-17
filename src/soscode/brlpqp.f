      SUBROUTINE BRLPQP(      QBAR     ,QOBJ     ,PBAR     ,PVEC           
     $   ,DELPVC   ,GVEC     ,GRADL    ,AZERO    ,NVAR     ,NFREE             
     $   ,MAXRES   ,NRES     ,DELZVC   ,RMAT     ,IROWR    ,JCOLR    
     $   ,NONZR    ,CBAR     ,CVEC     ,CZERO    ,ETABAR   ,ETAVEC     
     $   ,ETAZER   ,DELETA   ,QPFLTR   ,MXQPFL   ,MSUBE    ,MINEQL         
     $   ,MAXCON   ,CMAT     ,IROWC    ,JCOLC    ,NONZC    ,BBAR                      
     $   ,BVEC     ,BZERO    ,VLAMBR   ,VLAMDA   ,VLAMZR   ,DELLAM                     
     $   ,MSUBB    ,MAXBND   ,BMAT     ,IROWB    ,JCOLB    ,NONZB                  
     $   ,WMAT     ,IROWW    ,JCOLW    ,NONZW    ,QPENMU   ,PENRHO           
     $   ,SMTUVW   ,RWORK    ,LNRWRK   ,SCRTCH              
     $   ,LNSCRT   ,ISCRTC   ,LNISCR   ,NEEDED   ,EIGMIN   ,EIGMAX     
     $   ,ADDHSS   ,DIAGNL   ,RESET    ,IQPTRM)
C
C
C ======================================================================
C     BRLPQP===>brlpqp   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C        PURPOSE:     COMPUTE THE VALUES OF THE NVAR VARIABLES P WHICH
C                     MINIMIZE THE LOGARITHMIC BARRIER FUNCTION
C
C                                           M_B
C                       B(p,mu) = q(p) - mu sum ln[b_{i}(p)] 
C                                           i=1
C
C                                 + rho sum [ |t| + |u| + |v| + |w| ]
C
C                     subject to the msube equality constraints
C
C                             c(p) = 0
C
C                     The msubb bound inequality constraints
C
C                             b(p) .ge. 0 
C
C                     are treated using the logarithmic barrier term.
C                     The objective function is a quadratic of the form
C
C                       q(p) = .5*(p^T)(W + tau*I)p + (a0^T)p + q0
C
C                     with linear constraints of the form
C
C                       c(p) = Cp + c0
C
C                     and bound inequalities of the form
C
C                       b(p) = Bp + b0
C
C
C        ARGUMENTS:
C
C        QBAR         OBJECTIVE FUNCTION AT PBAR
C        QOBJ         OBJECTIVE FUNCTION AT PVEC 
C        PBAR         CURRENT POINT (NVAR)
C        PVEC         OLD POINT (NVAR)
C        DELPVC       PVEC SEARCH DIRECTION (NVAR)
C        GVEC         GRADIENT AT PVEC (NVAR)
C        GRADL        GRADIENT OF LAGRANGIAN (NVAR)
C        AZERO        INITIAL VALUE OF GRADIENT (NVAR)
C        NVAR         NUMBER OF INTERNAL VARIABLES; NVAR 
C        NFREE        NUMBER OF FREE REAL VARIABLES
C        MAXRES       MAXIMUM NUMBER OF RESIDUALS MAX(NRES,1)
C        NRES         NUMBER OF RESIDUALS
C        DELZVC       RESIDUAL STEP DIRECTION (MAXRES)
C        RMAT         RESIDUAL JACOBIAN DERIVATIVES AT PVEC (NONZR)
C        IROWR        ROW INDICES OF RESIDUAL JACOBIAN NONZEROS (NONZR)
C        JCOLR        COLUMN START INDICES OF RESIDUAL JACOBIAN NONZEROS (NVAR+1)
C        NONZR        NUMBER OF RESIDUAL JACOBIAN NONZEROS; NONZR = JCOLR(NVAR+1)-1
C        CBAR         EQUALITY CONSTRAINTS AT PBAR (MAXCON)
C        CVEC         EQUALITY CONSTRAINTS AT PVEC (MAXCON)
C        CZERO        INITIAL VALUE OF EQUALITY CONSTRAINTS (MAXCON)
C        ETABAR       CONSTRAINT MULTIPLIERS AT PBAR (MAXCON)
C        ETAVEC       CONSTRAINT MULTIPLIERS AT PVEC (MAXCON)
C        ETAZER       INITIAL MULTIPLIERS (MAXCON)
C        DELETA       ETAVEC SEARCH DIRECTION (MAXCON)
C        QPFLTR       FILTER ARRAY (MXQPFL*5)
C        MXQPFL       MAXIMUM LENGTH OF FILTER ARRAY 
C        MSUBE        NUMBER OF EQUALITY CONSTRAINTS MSUBE 
C        MINEQL       NUMBER OF INEQUALITIES
C        MAXCON       MAXIMUM NUMBER OF CONSTRAINTS MAX(MSUBE,1)
C        CMAT         CONSTRAINT DERIVATIVES AT PVEC (NONZC)
C        IROWC        ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C        JCOLC        COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C        NONZC        NUMBER OF JACOBIAN NONZEROS; NONZC = JCOLC(NVAR+1)-1
C        BBAR         BOUND INEQUALITIES AT PBAR (MAXBND)
C        BVEC         BOUND INEQUALITIES AT PVEC (MAXBND)
C        BZERO        INITIAL BOUND INEQUALITIES (MAXBND)
C        VLAMBR       BOUND MULTIPLIERS AT PBAR (MAXBND)
C        VLAMDA       BOUND MULTIPLIERS AT PVEC (MAXBND)
C        VLAMZR       INITIAL BOUND MULTIPLIERS (MAXBND)
C        DELLAM       VLAMDA SEARCH DIRECTION (MAXBND)
C        MSUBB        NUMBER OF BOUNDS
C        MAXBND       MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C        BMAT         BOUND DERIVATIVES AT PVEC (NONZB)
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
C        QPENMU       PENALTY WEIGHT MU:  INPUT INITIAL GUESS; OUTPUT FINAL VALUE 
C        PENRHO       PENALTY WEIGHT RHO (AND MAXIMUM MULTIPLIER LIMIT) 
C        SMTUVW       THE TERM sum [ |t| + |u| + |v| + |w| ]:  MUST BE SET
C                     ZERO FOR PRIMARY MODE
C        RWORK        REAL WORK ARRAY (LNRWRK)
C        LNRWRK       LENGTH OF WORK 
C        SCRTCH       REAL SCRATCH ARRAY (LNSCRT)
C        LNSCRT       LENGTH OF SCRTCH (4*NEQNS + 2) 
C        ISCRTC       INTEGER SCRATCH ARRAY (LNISCR)
C        LNISCR       LENGTH OF ISCRTC MAX(MXQPFL,NEQNS + MAXBND + 2) 
C        NEEDED       STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C        EIGMIN       MINIMUM EIGENVALUE ESTIMATE
C        EIGMAX       MAXIMUM EIGENVALUE ESTIMATE
C        ADDHSS       THE AMOUNT ADDED TO THE UPPER BLOCK (CORRES TO 
C                       TO NON-SLACK VARS) OF THE HESSIAN MATRIX WMAT
C        DIAGNL       LEVENBERG PARAMETER 
C                     ADDHSS = DIAGNL*(ONE - MIN(EIGMIN,ZERO))
C        RESET        LOGICAL VARIABLE; TRUE WHEN HESSIAN IS RESET TO I
C        IQPTRM       TERMINATION FLAG
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
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,
     $    TEN=1.0D1,ONEEM1=1.0D-1,ONEEM4=1.0D-4,POINT5=5.0D-1)
C
      DIMENSION PBAR(NVAR),PVEC(NVAR),DELPVC(NVAR),GVEC(NVAR),
     $    GRADL(NVAR),AZERO(NVAR),DELZVC(MAXRES),RMAT(NONZR),
     $    IROWR(NONZR),JCOLR(NVAR+1),CBAR(MAXCON),CVEC(MAXCON),
     $    CZERO(MAXCON),ETABAR(MAXCON),ETAVEC(MAXCON),ETAZER(MAXCON),
     $    DELETA(MAXCON),QPFLTR(MXQPFL*5),CMAT(NONZC),IROWC(NONZC),
     $    JCOLC(NVAR+1),BBAR(MAXBND),BVEC(MAXBND),BZERO(MAXBND),
     $    VLAMBR(MAXBND),VLAMDA(MAXBND),VLAMZR(MAXBND),DELLAM(MAXBND),
     $    BMAT(NONZB),IROWB(NONZB),JCOLB(NVAR+1),WMAT(NONZW),
     $    IROWW(NONZW),JCOLW(NVAR+1),RWORK(LNRWRK),
     $    SCRTCH(LNSCRT),ISCRTC(LNISCR)
C
      LOGICAL RESET
      LOGICAL OLDSET
      CHARACTER(LEN=6)  CHITER
C
      PARAMETER (MAXLVN=15)
      COMMON /LEVNP1/ LEVNIT
      COMMON /LEVNP2/ LEVNDY
      CHARACTER(LEN=60) LEVNDY(MAXLVN)
      CHARACTER(LEN=60) BLANK
      COMMON /CONDLB/ CNDLBL
      CHARACTER(LEN=17)  CNDLBL
      CHARACTER(LEN=16) LFTVAL,RITVAL
C
      LOGICAL LINEAR,CLOSE,NEWMLT,SAMSCL,USEW
C
      DATA BLANK(1:60) / ' '/
C
C     ******************************************************************
C
C ----------------------------------------------------------------------
C
C         ALGORITHM INITIALIZATION
C
C         LOGICAL QUANTITIES
C
      SAMSCL = .FALSE.
      NEWMLT = .FALSE.
      LINEAR = .FALSE.
      RESET = .FALSE.
      OLDSET = RESET
C
C         INTEGER QUANTITIES
C
      IT = 1
      IT1 = 1
      IQPTRM = 0
      IT1OLD = 1
      ITSAME = 0
      ITLEVN = 1
      LNQPFL = 0
      NEQNS = NVAR + NRES + MSUBE + MSUBB
      NSLK = NVAR - NFREE
      MEQUAL = MSUBE - MINEQL
C
C         INITIALIZE OUTPUT CONTROL FOR QP ALGORITHM
C
      IF(IOFSHR.EQ.0) THEN
        IPCSQP = 0
        IF(IOFLAG.GE.30) IPCSQP = IOFLAG
      ELSE
        IPCSQP = IOFSHR
      ENDIF
C
C         REAL QUANTITIES
C
      QBAR = QOBJ
      QMIN = QOBJ
      ACTRED = BIGNUM
      ALFA = ONE
      CNDNUM = ONE
      EIGMAX = ZERO
      EIGMIN = ZERO
      PGRATE = POINT5
      ETAOLD = -ONE
      ETARAT = ZERO
      CMAG = DAMAX(MSUBE,CVEC,1)
      CMAX = MAX(BIGCON,CMAG+ZEROOT)
      CMPLMT = ZERO
      DO I = 1,MSUBB
        CMPLMT = MAX(CMPLMT,ABS(BVEC(I)*VLAMDA(I)-QPENMU))
      enddo
      CMAX = MAX(BIGCON,CMAG+ZEROOT,CMPLMT+ZEROOT)
      REDFAC = .9D0
C
C             ALLOCATE STORAGE FROM THE REAL SCRATCH ARRAY
C     ---ALLOCATE THE ARRAYS (I.E. CONSTRUCT THE POINTERS)
C
      LCRSCL = 1
      LCCSCL = LCRSCL + NEQNS
      LCRTKT = LCCSCL + NEQNS
      LCSNKT = LCRTKT + NEQNS
      LCRSCR = LCSNKT + NEQNS
      NCRSCR = LNSCRT - LCRSCR + 1
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
      NCISCR = LNISCR - LCISCR + 1
      IF(NCISCR.LT.0) THEN
        PRINT *,'INTEGER SCRATCH ARRAY TOO SMALL; NEED =', LCISCR - 1
        STOP
      ENDIF
C
C             ---SAVE THE INITIAL VALUES 
C
      FZERO = QOBJ
      AZERO(1:NVAR) = GVEC(1:NVAR)
      DO I=1,MSUBE
        CZERO(I) = CVEC(I)
        ETAZER(I) = ETAVEC(I)
      ENDDO
      DO I=1,MSUBB
        BZERO(I) = BVEC(I)
        VLAMZR(I) = VLAMDA(I)
      ENDDO
C
C         ---INITIALIZE THE STEP P 
C
      pvec(1:nvar) = zero
C
C             INITIALIZE TERSE OUTPUT IF NECESSARY
C
      NPRFLG = -1
C
      CNDLBL = 'Cond(K)..........'
C
C ----------------------------------------------------------------------
C
C         BEGIN ITERATION
C
 120  CONTINUE
C
C ----------------------------------------------------------------------
C
      IFUNER = 0
C
C ----------------------------------------------------------------------
C
C         CONVERGENCE TESTS
C
C         ---FIRST COMPUTE PRIMAL-DUAL NECESSARY CONDITION QUANTITIES
C
      CALL MPDNEC(IFUNER,CVEC,MSUBE,MAXCON,ETAVEC,CMAT,IROWC,JCOLC,
     $    NONZC,BVEC,MSUBB,MAXBND,VLAMDA,BMAT,IROWB,JCOLB,NONZB,
     $    GVEC,NVAR,QPENMU,SCRTCH,LNSCRT,GRADL,GRDLNM,CLAMNM,ERPTHB,
     $    DELFNM,ERREQL,ETANRM)
C
C         ---PERFORM TESTS
C
      CALL QCNVRG(GRDLNM,PGDTOL,ERREQL,CONTOL,ERPTHB,CONTOL,
     $    QPENMU,CNDNUM,ETARAT,DELPVC,PVEC,DELETA,ETAVEC,DELLAM,VLAMDA,
     $    ALFA,DELFNM,QOBJ,QMIN,OBJTOL,NVAR,MSUBE,MSUBB,MAXCON,
     $    MAXBND,IT,MXQPIT,IQPTRM)
C
      CHITER = '      '
      WRITE(CHITER(2:4),'(I3)') IT
      CALL HHSDEL(CHITER,' ',' ',6,NDEL,IERDEL)
C
C         INITIALIZE LOG-BARRIER FUNCTION FOR UNIVARIATE SEARCH
C
      CALL BMERUT(CVEC,ETAVEC,BVEC,VLAMDA,QOBJ,SMTUVW,QPENMU,
     $    PENRHO,MSUBE,MSUBB,SMLNBI,QMRT,CMRT,XLAGR0,C0TC0)
C
      CBTCB = C0TC0
C
      ETAMAG = DAMAX(MSUBE,ETAVEC,1)
      VLAMAG = DAMAX(MSUBB,VLAMDA,1)
      IF(LNQPFL.EQ.0) THEN
        CALL FILTER(QOBJ,CMRT,CMAX,SMLNBI,SMTUVW,ETAMAG,QPENMU,PENRHO,
     $    QPFLTR,LNQPFL,MXQPFL,NWFLTR,ISCRTC,LNISCR)
      ENDIF
C
C         ITERATION PRINT
C
      IF(IPCSQP.GE.10) THEN
        CALL HHCPRS(CHITER,' ','-')
        WRITE(IPUNLP,1001) CHITER
        WRITE(LFTVAL,'(SP,1PG14.6)') QOBJ + PENRHO*SMTUVW
        WRITE(RITVAL,'(SP,1PG14.6)') GRDLNM
        CALL QPPRNT('Objective Function',LFTVAL,
     $               'Grad L',RITVAL,IPUNLP)
        WRITE(LFTVAL,'(SP,1PG14.6)') QMRT
        WRITE(RITVAL,'(SP,1PG14.6)') ERREQL
        CALL QPPRNT('Log-Barrier Function',LFTVAL,
     $               'Equality Error',RITVAL,IPUNLP)
      ENDIF
C
C         TRANSFER TO TERMINATION PROCEDURE BLOCK IF CONVERGENCE
C         TESTS ARE SATISFIED
C
      IF(IQPTRM.NE.0) THEN
        IF(IPCSQP.GE.10) THEN
          WRITE(LFTVAL,'(SP,1PG14.6)') QPENMU
          WRITE(RITVAL,'(SP,1PG14.6)') ERPTHB
          CALL QPPRNT('Barrier Parameter',LFTVAL,
     $               'Complementarity',RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') DIAGNL
          WRITE(RITVAL,'(SP,1PG14.6)') CNDNUM
          CALL QPPRNT('Levenberg Parameter',LFTVAL,
     $               CNDLBL,RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') EIGMIN
          WRITE(RITVAL,'(SP,1PG14.6)') EIGMAX
          CALL QPPRNT('Min. Eigenvalue',LFTVAL,
     $               'Max. Eigenvalue',RITVAL,IPUNLP)
          CALL FLTRPR(' QP',QPENMU,PENRHO,QPFLTR,LNQPFL,MXQPFL)
        ENDIF
        IF(IPCSQP.GT.0.AND.IPCSQP.LT.10) THEN
          CALL BTURSE(NPRFLG,CHITER,LNQPFL,CNDNUM,ALFA,DIAGNL,ERPTHB,
     $    QPENMU,ERREQL,GRDLNM,QOBJ,QMRT)
        ENDIF
        GO TO 220
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      IFUNER = 0
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         ---BARRIER PARAMETER UPDATE
C
      IF(IT.GT.1) THEN
C
        GRDTST = PGDTOL*MAX(ONE,DELFNM)*(ONE + ABS(QOBJ))
C
        IF(MSUBB.GT.0) THEN
C
          ERRKTC = MAX(GRDLNM/MAX(ONE,DELFNM,CLAMNM),ERREQL,ERPTHB)
          TOLPTH = MIN(TEN*QPENMU,PTHTOL)
          PENOLD = QPENMU
C
C         FIRST CHECK TO SEE THAT LEVENBERG MODIFICATION IS SMALL AND/OR
C         UNCONSTRAINED MINIMIZATION WITH CURRENT BARRIER PARAMETER 
C         IS NEAR COMPLETION (I.E. GRDLNM IS SMALL)
C
          CLOSE = ALFA.GE.ONEEM1.AND.(ERRKTC.LT.TOLPTH.OR.
     $            GRDLNM.LT.GRDTST.OR.DIAGNL.LE.TEN*ZEROOT)   
C
          IF(CLOSE.OR.ITSAME.EQ.IMAXMU) THEN
C
            IF(ERRKTC.LT.TOLPTH) THEN
C
C               CURRENT POINT IS "NEAR" THE CENTRAL PATH
C               ---BARRIER PARAMETER REDUCTION
C
              IF(QPENMU.LT.ONEEM4) THEN
                HATMU = TEN*QPENMU**2
              ELSE
                HATMU = QPENMU/TEN
              ENDIF
C
C               COMPUTE ACCELERATED ESTIMATE
C
              FKAPPA = 1.1D0*TEN
              CALL BQPSTR(CVEC,MSUBE,MAXCON,ETAVEC,CMAT,
     $          IROWC,JCOLC,NONZC,BVEC,MSUBB,MAXBND,BMAT,IROWB,JCOLB,
     $          NONZB,GVEC,NVAR,FKAPPA,PFAST,SCRTCH,LNSCRT)
C
              HATMU = MIN(HATMU,PFAST)
C
C           COMPUTE DIVERGENCE RATE TEST QUANTITIES
C
              IF(ETAOLD.GT.ZEROMN) ETARAT = (ETANRM - ETAOLD)/ETAOLD
              ETAOLD = ETANRM
C
C           UPDATE THE BARRIER PARAMETER
C
              QPENMU = HATMU
C
C           RECOMPUTE FILTER QUANTITIES
C
              CALL BMERUT(CVEC,ETAVEC,BVEC,VLAMDA,QOBJ,SMTUVW,QPENMU,
     $        PENRHO,MSUBE,MSUBB,SMLNBI,QMRT,CMRT,XLAGR0,C0TC0)
C
              CBTCB = C0TC0
C
C           REFRESH THE FILTER
C
              LNQPFL = 0
              CALL FILTER(QOBJ,CMRT,CMAX,SMLNBI,SMTUVW,ETAMAG,QPENMU,
     $          PENRHO,QPFLTR,LNQPFL,MXQPFL,NWFLTR,ISCRTC,LNISCR)
C
            ELSEIF(ITSAME.EQ.IMAXMU) THEN
C
              HATMU = REDFAC*QPENMU
C
C             COMPUTE DIVERGENCE RATE TEST QUANTITIES
C
              ETARAT = ZERO
              ETAOLD = ETANRM
C
C             UPDATE THE BARRIER PARAMETER
C
              QPENMU = HATMU
C
C           RECOMPUTE FILTER QUANTITIES
C
C
              CALL BMERUT(CVEC,ETAVEC,BVEC,VLAMDA,QOBJ,SMTUVW,QPENMU,
     $          PENRHO,MSUBE,MSUBB,SMLNBI,QMRT,CMRT,XLAGR0,C0TC0)
C
              CBTCB = C0TC0
C
C             REFRESH THE FILTER
C
              LNQPFL = 0
              CALL FILTER(QOBJ,CMRT,CMAX,SMLNBI,SMTUVW,ETAMAG,QPENMU,
     $          PENRHO,QPFLTR,LNQPFL,MXQPFL,NWFLTR,ISCRTC,LNISCR)
C
            ENDIF
C
          ENDIF
C
          IF(QPENMU.EQ.PENOLD) THEN
C
C           ITERATION COUNT WITH UNCHANGED QPENMU
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
      IF(ETAMAG.GT.PENRHO) NEWMLT = .TRUE.
C
      IF(NEWMLT) THEN
C
C             COMPUTE NEW FIRST ORDER MULTIPLIERS
C
        NEWMLT = .FALSE.
C
        BSTMU = QPENMU
        ETABAR(1:MSUBE) = ETAVEC(1:MSUBE)
        VLAMBR(1:MSUBB) = VLAMDA(1:MSUBB)
        MUOPTN = 2
C
C            COMPUTE OPTIMUM MU 
C
 130    CONTINUE
C
        NCALL = 0
        CALL LSQMLT(MUOPTN,NCALL,CVEC,MSUBE,MAXCON,ETAVEC,CMAT,
     $           IROWC,JCOLC,NONZC,BVEC,MSUBB,MAXBND,VLAMDA,BMAT,IROWB,
     $           JCOLB,NONZB,GVEC,NVAR,QPENMU,SCRTCH,RWORK,     
     $           LNRWRK,IPUNLP,IPCSQP,
     $           CNDNUM,GRDLNM,IERLSM,NEEDED)
C      
        IF(IERLSM.NE.0) THEN
          PRINT *,'IERLSM =',IERLSM
          IQPTRM = 9
          GO TO 220
        ENDIF
C
C           COMPUTE ERROR IN KT CONDITIONS
C
        PMUCVC = ERREQL/(TEN*1.1D0)
C
C         CHECK THE MINIMUM VALUES FOR THE MULTIPLIERS AND BARRIER
C         PARAMETER
C
        IF(QPENMU.LT.BSTMU.AND.MSUBB.GT.0.AND.MUOPTN.EQ.2) THEN
          IF(QPENMU.GT.PMUCVC) THEN
            QPENMU = MAX(ZEROOT,QPENMU)
          ELSE
            QPENMU = MAX(ZEROOT,MIN(PMUCVC,BSTMU))
          ENDIF
          MUOPTN = 1
          GO TO 130
        ENDIF
C
        ETAMAG = ZERO
        DO I=1,MSUBE
          ETAVEC(I) = ETABAR(I)
          ETAMAG = MAX(ETAMAG,ABS(ETABAR(I)))
        ENDDO
C
        VLAMAG = ZERO
        DO I=1,MSUBB
          VLAMDA(I) = VLAMBR(I)
          VLAMAG = MAX(VLAMAG,ABS(VLAMBR(I)))
        ENDDO
C
      ENDIF
C
      IF(ETAMAG.GT.PENRHO) THEN
        IQPTRM = 9
        GO TO 220
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C     
C         SOLVE WITHOUT ITERATIVE REFINEMENT UNLESS NEAR THE ANSWER
C
      IF(GRDLNM.LT.GRDTST.AND.ERREQL.LT.CONTOL.AND.QPENMU.LT.ZEROOT) 
     $    THEN
        IREFIN = 2
      ELSE
        IREFIN = 0
      ENDIF
C
C         COMPUTE SEARCH DIRECTION
C
      USEW = .TRUE.
      CALL BRSRCH( NVAR, NEQNS, NRES, MSUBE, MSUBB, NSLK, 
     $       DIAGNL, LINEAR, USEW, IPUNLP, IPCSQP, NONZR, IROWR, 
     $       JCOLR, RMAT, NONZB, IROWB, JCOLB, BMAT, NONZC,
     $       IROWC, JCOLC, CMAT, NONZW, IROWW, JCOLW, WMAT, QPENMU,
     $       BVEC, MAXBND, CVEC, MAXCON, GRADL, GVEC, ETAVEC, VLAMDA,  
     $       LNRWRK, RWORK,ISCRTC(LCITMP),ISCRTC(LCNZRB),
     $       SCRTCH(LCRSCL), SCRTCH(LCCSCL),SCRTCH(LCRTKT), 
     $       SCRTCH(LCSNKT), LNSYMB, CNDNUM, NEEDED, DELPVC,
     $       DELZVC, MAXRES, DELETA, DELLAM, EIGMIN,EIGMAX,  
     $       ADDHSS, RESET,IQPTRM  ) 
C
C         TRANSFER TO TERMINATION PROCEDURE BLOCK IF SEARCH
C         DIRECTION CALCULATION FAILED
C
      IF(OLDSET.AND.RESET) IQPTRM = -10
      OLDSET = RESET
      IF(IQPTRM.NE.0.AND.IQPTRM.NE.2) THEN
        GO TO 220
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
     $    QPENMU,SCRTCH,SCRTCH(LCBDLY),SCRTCH(LCWDLY),DELPVC,DELZVC,
     $    DELETA,DELLAM,QMRT,SAMSCL,QMIN,SLOPEL,DELTAH,DELTAR,
     $    DELTAX,SLOPEB,ALFEST,ALFLIM,GAMEST)
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
      IF(IPCSQP.GE.10) THEN
          WRITE(LFTVAL,'(SP,1PG14.6)') QPENMU
          WRITE(RITVAL,'(SP,1PG14.6)') ERPTHB
          CALL QPPRNT('Barrier Parameter',LFTVAL,
     $               'Complementarity',RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') DIAGNL
          WRITE(RITVAL,'(SP,1PG14.6)') CNDNUM
          CALL QPPRNT('Levenberg Parameter',LFTVAL,
     $               CNDLBL,RITVAL,IPUNLP)
          WRITE(LFTVAL,'(SP,1PG14.6)') EIGMIN
          WRITE(RITVAL,'(SP,1PG14.6)') EIGMAX
          CALL QPPRNT('Min. Eigenvalue',LFTVAL,
     $               'Max. Eigenvalue',RITVAL,IPUNLP)
          CALL FLTRPR(' QP',QPENMU,PENRHO,QPFLTR,LNQPFL,MXQPFL)
        IF(RESET) WRITE(IPUNLP,1003)
      ENDIF
C
      IF(IPCSQP.GT.0.AND.IPCSQP.LT.10) THEN
        CALL BTURSE(NPRFLG,CHITER,LNQPFL,CNDNUM,ALFA,DIAGNL,ERPTHB,
     $    QPENMU,ERREQL,GRDLNM,QOBJ,QMRT)
      ENDIF
C
      IFUNER = 0
C
 150  CONTINUE
C
      IF(MSUBB.EQ.0) THEN
        ALFLMI = TWO
        QPNMUI = ZERO
      ELSE
        ALFLMI = ALFLIM 
        QPNMUI = QPENMU
      ENDIF
C
      CALL QLINES(QMRT,SLOPEB,QBRMRT,ALFA,ALFLMI,QPNMUI,CBTCB,
     $    C0TC0,CBRMRT,IT1,IT1MAX,NWFLTR,IFUNER,IOFLIN,IPUNLP,
     $    IFCALL,IERLIN,SAMSCL)
C
      IF(IFCALL.EQ.0) GO TO 210
C
C         COMPUTE LOCATION OF POINT 
C         PBAR = PVEC + ALFA*DELPVC 
C
      DO I = 1,NVAR
        PBAR(I) = PVEC(I) + ALFA*DELPVC(I)
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
C         ---COMPUTE Wp AND SAVE IN SCRTCH
C
      CALL SYMMVP(WMAT,IROWW,JCOLW,NVAR,PBAR,SCRTCH)
C
C         ---COMPUTE THE OBJECTIVE AND GRADIENT, 
C              qbar = .5(p^T)(Wp + tau p) + (azero^T)p + fzero
C              g = (Wp + tau p) + azero 
C
      QBAR = FZERO
      DO I = 1,NVAR
        GVEC(I) = SCRTCH(I) + ADDHSS*PBAR(I) + AZERO(I)
        WPAZP = (POINT5*SCRTCH(I) + POINT5*ADDHSS*PBAR(I) 
     $           + AZERO(I))*PBAR(I)
        QBAR = QBAR + WPAZP
      enddo
C
C         ---COMPUTE Cp AND SAVE IN SCRTCH
C
      CALL MVPSPR(1,MSUBE,NVAR,CMAT,IROWC,JCOLC,PBAR,SCRTCH)
C
C         ---COMPUTE EQUALITIES, c = Cp + czero 
C
      CBAR(1:MSUBE) = SCRTCH(1:MSUBE) + CZERO(1:MSUBE)
C
C         ---COMPUTE Bp AND SAVE IN SCRTCH
C
      CALL MVPSPR(1,MSUBB,NVAR,BMAT,IROWB,JCOLB,PBAR,SCRTCH)
C
C         ---COMPUTE BOUNDS, b = Bp + bzero
C
      DO I = 1,MSUBB
        BBAR(I) = SCRTCH(I) + BZERO(I)
        BBAR(I) = MAX(ZEROMN,BBAR(I))
      enddo
C
      IFUNER = 0
C
C         CONSTRUCT LOG-BARRIER FUNCTION FOR UNIVARIATE SEARCH
C
      CALL BMERUT(CBAR,ETABAR,BBAR,VLAMBR,QBAR,SMTUVW,QPENMU,
     $    PENRHO,MSUBE,MSUBB,SMLNBI,QBRMRT,CBRMRT,XLAGRN,CBTCB)
C
C         CHECK CURRENT POINT AGAINST FILTER
C
      ETAMAG = DAMAX(MSUBE,ETAVEC,1)
      CALL FILTER(QBAR,CBRMRT,CMAX,SMLNBI,SMTUVW,ETAMAG,QPENMU,PENRHO,
     $    QPFLTR,LNQPFL,MXQPFL,NWFLTR,ISCRTC,LNISCR)
      IF(NWFLTR.NE.0) IFUNER = 2
C
C         FUNCTION EVALUATED. REENTER ONE DIMENSIONAL SEARCH
C
      GO TO 150
C
C ----------------------------------------------------------------------
C
C         UPDATE INFORMATION
C
 210  CONTINUE
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
C         SET IQPTRM FLAG BASED ON LINE SEARCH ERROR INDICATOR IERLIN
C
      IF(IERLIN.EQ.3) THEN
        IQPTRM = -4
      ELSEIF(IERLIN.EQ.5) THEN
        IQPTRM = 7
      ELSEIF(IERLIN.EQ.2) THEN
        IQPTRM = 3
      ELSEIF(IERLIN.EQ.-1) THEN
        IQPTRM = 7
      ENDIF
C
C         COMPUTE ACTUAL REDUCTION
C
      ACTRED = XLAGR0 - XLAGRN
      IT1OLD = IT1
C      
      IT = IT + 1
      IT1 = 1
      QOBJ = QBAR
      QMRT = QBRMRT
C
      PVEC(1:NVAR) = PBAR(1:NVAR)
      DO I=1,MSUBE
        ETAVEC(I) = ETABAR(I)
        CVEC(I) = CBAR(I)
      ENDDO
      DO I=1,MSUBB
        VLAMDA(I) = VLAMBR(I)
        BVEC(I) = BBAR(I)
      ENDDO
C
C         IF IQPTRM =   9 RANK DEFICIENT CONSTRAINTS
C         IF IQPTRM =   8 INCONSISTENT CONSTRAINTS
C         IF IQPTRM =   6 MAX. NO. OF FUNCTION EVALS.
C         IF IQPTRM =  -1 (REDUCED OBJECTIVE FUNCTION IS LINEAR) OR;
C         IF IQPTRM =  -6 BARKKT FAILED WITH UNEXPECTED ERROR
C         IF IQPTRM =  -7 SINGULAR MATRIX PREVENTED SEARCH DIRECTION 
C                        FROM BEING CALCULATED;
C         IF IQPTRM =  -8 INSUFFICIENT REAL STORAGE FOR SRCHDR
C         IF IQPTRM =  -9 SEARCH DIRECTION FAILED; SLOPE CONDITION
C         IF IQPTRM = -10 SEARCH DIRECTION FAILED; MAX. DIAGONAL VALUE.
C         IF IQPTRM = -12 INSUFFICIENT INTEGER STORAGE FOR SRCHDR
C         IF IQPTRM = -13 I/O ERROR (INSUFFICIENT DISK SPACE)
C 
C         THEN, TERMINATE IMMEDIATELY
C
      IF(IQPTRM.LE.-1.OR.IQPTRM.EQ.6) GO TO 220
C
C ----------------------------------------------------------------------
C
C         END OF ITERATION
C
      GO TO 120
C
C ----------------------------------------------------------------------
C
C         TERMINATION PROCEDURES
C
 220  CONTINUE
C
      QOBJ = FZERO
      QBAR = QOBJ
      PBAR(1:NVAR) = PVEC(1:NVAR)
C
C         ---COMPUTE THE TOTAL QP STEPS 
C
      DELETA(1:MSUBE) = ETABAR(1:MSUBE) - ETAZER(1:MSUBE)
      DELLAM(1:MSUBB) = VLAMBR(1:MSUBB) - VLAMZR(1:MSUBB)
C
C         ---RESTORE THE INITIAL VALUES
C
      GVEC(1:NVAR) = AZERO(1:NVAR)
      DO I=1,MSUBE
        CVEC(I) = CZERO(I)
        ETAVEC(I) = ETAZER(I)
      ENDDO
      DO I=1,MSUBB
        BVEC(I) = BZERO(I)
        VLAMDA(I) = VLAMZR(I)
      ENDDO
C
      IF((NSLK-MINEQL).GT.0) THEN
        TUVMAX = MAXVAL(PBAR(NFREE+MINEQL+1:NFREE+NSLK))
        IF((IQPTRM.EQ.1.OR.IQPTRM.EQ.2).AND.TUVMAX.GT.CONTOL) IQPTRM = 8
      ENDIF
C
C
 1001 FORMAT(T3,'*',T106,'*'/T3,'*',T15,'-------------------------------
     $ QP Iteration ',A6,'-----------------------------',
     $   T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....RESET HESSIAN AND LAGRANGE
     $ MULTIPLIERS',T106,'*')
 1006 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....SWITCH TO EQUAL PRIMAL-DUA
     $L STEP SCALING',T106,'*')
      RETURN
      END
