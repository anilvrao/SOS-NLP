      SUBROUTINE BARYER(      FNLPBR   ,GVEC     ,YBAR     ,NVARN
     $   ,NVARR    ,NFIX     ,NFREE    ,NXLWR    ,NXUPR    ,IFERR                  
     $   ,MAXRES   ,NRES     ,RMAT     ,IROWR    ,JCOLR    ,NONZR           
     $   ,CBAR     ,MSUBE    ,MIGNOR   ,MEQUAL   ,MINEQL   ,NALWR           
     $   ,NAUPR    ,MAXCON   ,ETABAR   ,CMAT     ,IROWC    ,JCOLC       
     $   ,NONZCN   ,NONZCR   ,BBAR     ,MSUBBN   ,MSUBBR   ,MAXBND       
     $   ,VLAMBR   ,BMAT     ,IROWB    ,JCOLB    ,NONZBN   ,NONZBR        
     $   ,WMAT     ,IROWW    ,JCOLW    ,NONZWN   ,NONZWR   ,IWORK      
     $   ,LNIWRK   ,RWORK    ,LNRWRK   ,NEEDED   ,IREVRS   ,IRVCOM   
     $   ,RELAX    ,FZMODE   ,IERNLP)
C
C
C ======================================================================
C     BARYER===>baryer   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C        PURPOSE:     Compute the values of the nvar variables y which
C
C                             Minimize f(y)
C
C                     subject to the msube equality constraints
C
C                             c(y) = 0
C
C                     and the msubb bound inequality constraints
C
C                             b(y) .ge. 0 
C
C                     using an interior point (logarithmic barrier)
C                     method.
C
C        ARGUMENTS:
C
C        FNLPBR       OBJECTIVE FUNCTION AT YBAR 
C        GVEC         GRADIENT AT YBAR (NVAR)
C        YBAR         CURRENT POINT (NVAR)
C        NVARN        NUMBER OF VARIABLES (NORMAL); NVARN = NFREE + NSLK
C        NVARR        NUMBER OF VARIABLES (RELAXATION); NVARR = NFREE + NSLK
C        NFIX         NUMBER OF FIXED REAL VARIABLES
C        NFREE        NUMBER OF FREE REAL VARIABLES
C        MSUBX        NUMBER OF LOWER AND UPPER VARIABLE BOUNDS 
C        IFERR        = 0 WHEN FUNCTION IS EVALUATED
C                     = 1 WHEN FUNCTION EVALUATION IS IMPOSSIBLE
C                     = -100 WHEN MAXIMUM NUMBER OF EVALUATIONS EXCEEDED
C                     = -101 WHEN POINT IS FEASIBLE (IN FEASIBILITY MODE)
C        MAXRES       MAXIMUM NUMBER OF RESIDUALS MAX(NRES,1)
C        NRES         NUMBER OF RESIDUALS
C        RMAT         RESIDUAL JACOBIAN DERIVATIVES AT PVEC (NONZR)
C        IROWR        ROW INDICES OF RESIDUAL JACOBIAN NONZEROS (NONZR)
C        JCOLR        COLUMN START INDICES OF RESIDUAL JACOBIAN NONZEROS (NVAR+1)
C        NONZR        NUMBER OF RESIDUAL JACOBIAN NONZEROS; NONZR = JCOLR(NVAR+1)-1
C        CBAR         EQUALITY CONSTRAINTS AT YBAR (MAXCON)
C        MSUBE        NUMBER OF EQUALITY CONSTRAINTS MSUBE = MCON - MIGNOR
C        MIGNOR       NUMBER OF IGNORED CONSTRAINTS 
C        MEQUAL       NUMBER OF REAL EQUALITY CONSTRAINTS 
C        MINEQL       NUMBER OF REAL INEQUALITY CONSTRAINTS 
C        MSUBS        NUMBER OF LOWER AND UPPER CONSTRAINT BOUNDS 
C        MAXCON       MAXIMUM NUMBER OF CONSTRAINTS MAX(MSUBE,1)
C        ETABAR       LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C        CMAT         CONSTRAINT DERIVATIVES AT YBAR (NONZC)
C        IROWC        ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C        JCOLC        COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C        NONZCN       NUMBER OF JACOBIAN NONZEROS (NORMAL); NONZC = JCOLC(NVAR+1)-1
C        NONZCR       NUMBER OF JACOBIAN NONZEROS (RELAXATION); NONZC = JCOLC(NVAR+1)-1
C        BBAR         BOUND INEQUALITIES (MAXBND)
C        MSUBBN       NUMBER OF BOUNDS (NORMAL)
C        MSUBBR       NUMBER OF BOUNDS (RELAXATION)
C        MAXBND       MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C        VLAMBR       LAGRANGE MULTIPLIERS FOR BOUNDS (MAXBND)
C        BMAT         BOUND DERIVATIVES AT YVEC (NONZB)
C        IROWB        ROW INDICES OF BOUND JACOBIAN NONZEROS (NONZB)
C        JCOLB        COLUMN INDICES OF BOUND JACOBIAN NONZEROS (NVAR+1)
C        NONZBN       NUMBER OF BOUND JACOBIAN NONZEROS (NORMAL)
C        NONZBR       NUMBER OF BOUND JACOBIAN NONZEROS (RELAXATION)
C        WMAT         HESSIAN OF THE LAGRANGIAN WITH RESPECT TO THE NVAR VARIABLES (NONZW)
C        IROWW        INTEGER ROW INDEX VECTOR FOR LOWER TRIANGLE OF WMAT (NONZW)
C        JCOLW        INTEGER COLUMN START VECTOR (NVAR+1)
C        NONZWN       NUMBER OF NONZEROS IN THE HESSIAN (NORMAL) (NONZW = JCOLW(NVAR+1)-1)
C        NONZWR       NUMBER OF NONZEROS IN THE HESSIAN (RELAXATION) (NONZW = JCOLW(NVAR+1)-1)
C        IWORK        INTEGER WORK ARRAY (LNIWRK)
C        LNIWRK       LENGTH OF IWORK 
C        RWORK        REAL WORK ARRAY (LNRWRK)
C        LNRWRK       LENGTH OF WORK 
C        NEEDED       STORAGE REQUIRED WHEN LNRWRK OR LNIWRK IS TOO SMALL
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
C        RELAX        LOGICAL FLAG:  TRUE FOR RELAXATION MODE, FALSE OTHERWISE
C        IERNLP       = 0 WHEN YBAR IS THE SOLUTION TO THE STATED PROBLEM;
C                     IF IERNLP .NE. O AND IRVCOM = 0 THE ALGORITHM
C                     TERMINATED ABNORMALLY
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TEN=1.0D1,
     $    ONEEP6=1.0D6,POINT5=5.0D-1,ONEP12=1.D12)
      DIMENSION GVEC(NVARR),YBAR(NVARR),RMAT(NONZR),IROWR(NONZR),
     $    JCOLR(NVARN+1),CBAR(MAXCON),ETABAR(MAXCON),CMAT(NONZCR),
     $    IROWC(NONZCR),JCOLC(NVARR+1),BBAR(MAXBND),VLAMBR(MAXBND),
     $    BMAT(NONZBR),IROWB(NONZBR),JCOLB(NVARR+1),WMAT(NONZWR),
     $    IROWW(NONZWR),JCOLW(NVARR+1),IWORK(*),RWORK(*)
      DIMENSION IREVRS(5)
      PARAMETER (MXFLGB=10, LNISGB=10)
      DIMENSION FLTRGB(MXFLGB,5),ISCRGB(LNISGB)
C
      LOGICAL STOPIT,RELAX,REEDOO
      LOGICAL FZMODE,FEEZPT,OPMODE,BIGMLT
      LOGICAL OPTPRB,FRSTRY
C
C         REVERSE COMMUNICATION RETURN BRANCH POINT -- BY CONVENTION ALL 
C         RETURN POINTS ARE NUMBERED BEGINNING WITH 501
C
      IF (IMAX(4,IREVRS,1).GE.1) THEN
        SELECT CASE (IFESR)
          CASE(1)
            GO TO 501
          CASE(2)
            GO TO 502
          CASE(3)
            GO TO 503
          CASE(4)
            GO TO 504
        END SELECT
      ENDIF
C
C ----------------------------------------------------------------------
C
C         ALGORITHM INITIALIZATION
C
C
C         ---LOGICAL VARIABLES
C
C            .....ALGORITHM TERMINATION FLAG: TERMINATE WHEN TRUE
      STOPIT = .FALSE.
C     FZMODE .....FEASIBILITY MODE FLAG: FINDING A FEASIBLE POINT WHEN TRUE
C            .....FEASIBLE POINT FOUND FLAG: SET TRUE WHEN A FEASIBLE POINT EXISTS
      FEEZPT = .FALSE.
C            .....OPTIMIZATION PROBLEM FLAG: A REAL OBJECTIVE EXISTS
      OPTPRB = ALGOPT(1:1).EQ.'M'.OR.ALGOPT(2:2).EQ.'M'
C            .....OPTIMIZATION MODE FLAG: OPTIMIZING WHEN TRUE
      OPMODE = ALGOPT(1:1).EQ.'M'
      IF(NFREE.EQ.MEQUAL.AND.IRELAX.NE.2) THEN
        FZMODE = .TRUE.
        OPMODE = .NOT.FZMODE
      ENDIF
      IF(IRELAX.EQ.2) THEN
        FZMODE = .FALSE.
        OPMODE = .NOT.FZMODE
      ENDIF
C            .....FIRST TRY FLAG
      FRSTRY = .TRUE.
C
C         ---INTEGER VARIABLES
C
      ITERM = 0
      MXFLTR = IMAXMU*IT1MAX
      MUREDN = 0
      CALL FLIPER(RELAX,
     $     NVARN,NONZCN,MSUBBN,NONZBN,NONZWN,
     $     NVARR,NONZCR,MSUBBR,NONZBR,NONZWR,
     $     NVAR,NONZC,MSUBB,NONZB,NONZW)
C
      NDIM = NFIX + NFREE
      MCON = MIGNOR + MEQUAL + MINEQL
      NSLK = NVAR - NFREE
      LNFLGB = 0
C
C         ---REAL VARIABLES
C
      FMIN = BIGNUM
      SMLRHO = RHOLWR
C
      IF(IRELAX.EQ.0) THEN 
        PMUBND = PMULWR
      ELSEIF(IRELAX.EQ.1) THEN
        IF(RELAX) THEN
          PMUBND = PMULST
        ELSE
          PMUBND = PMULWR
        ENDIF
      ELSEIF(IRELAX.EQ.2) THEN
        PMUBND = PMULWR
      ELSE
        PRINT *,'IRELAX IS WRONG'
        RETURN
      ENDIF
C
C
C ----------------------------------------------------------------------
C
C             ALLOCATE STORAGE BASED ON RELAXATION PROBLEM DIMENSIONS
C             STORAGE ALLOCATION FOR THE REAL ARRAY
C     ---ALLOCATE THE ARRAYS (I.E. CONSTRUCT THE POINTERS)
C
      LCYVEC = 1
      LCDLYV = LCYVEC + NVARR
      LCGRDL = LCDLYV + NVARR
      LCDLZV = LCGRDL + NVARR
      LCCVEC = LCDLZV + MAX(MAXCON,MAXRES)
      LCETAV = LCCVEC + MAXCON
      LCDLET = LCETAV + MAXCON
      LCFLTR = LCDLET + MAXCON
      LCBVEC = LCFLTR + 5*MXFLTR
      LCVLAM = LCBVEC + MAXBND
      LCDLLM = LCVLAM + MAXBND
      LNSCRT = 4*(NVARR+MAXRES+MAXCON+MAXBND) + 2
      MXQPFL = IMAXMU*IT1MAX
      LNSCRT = LNSCRT + 3*NVARR + 2*MAXCON + 2*MAXBND + 5*MXQPFL
      LCSCRT = LCDLLM + MAXBND
      LCRHLD = LCSCRT + LNSCRT
      LNRHLD = LNRWRK - LCRHLD + 1
C
      IF(LNRHLD.LT.0) THEN
        NEEDED = LCRHLD + 1
        IERNLP = -131
        GO TO 200
      ENDIF
C
C     ---------------
C
C             STORAGE ALLOCATION FOR THE INTEGER ARRAY
C     ---ALLOCATE THE ARRAYS (I.E. CONSTRUCT THE POINTERS)
C
      LNISCR = MAX(MXFLTR,NVARR + NRES + MAXCON + 2*MAXBND + 2)
C
      LCISCR = 1
      LCIWRK = LCISCR + LNISCR
      NCIWRK = LNIWRK - LCIWRK + 1
C
      IF(NCIWRK.LT.0) THEN
        NEEDED = LCIWRK - 1
        IERNLP = -132
        GO TO 200
      ENDIF
C
C        SPECIAL TEST FOR FEASIBLE INITIAL POINT---CHECK INITIAL FNLPBR
C
      IF(FZMODE) THEN
C
C         SET FEASIBILITY MODE TERMINATION FLAGS
C
        FEEZPT = IFERR.EQ.-101
C
        IF(FEEZPT.AND.OPTPRB) THEN
C
C         BEGIN OPTIMIZATION MODE
C
          ITERM = 0
          FZMODE = .FALSE.
          OPMODE = .NOT.FZMODE
C
C         COMPUTE REAL OBJECTIVE AND GRADIENT, AND THEN REENTER OPTIMIZATION
C
          GO TO 190
C
        ELSEIF(FEEZPT.AND..NOT.OPTPRB) THEN
C
C         TRANSFER TO ALGORITHM TERMINATION TESTS
C
          ITERM = 1
          GO TO 130
C
        ENDIF
C
      ENDIF
C
 110  CONTINUE
C
C         PRINT HEADING
C
      IF(IOFLAG.GE.10) THEN
        IF(FZMODE.EQV.OPMODE) THEN
          PRINT *,'FZMODE=OPMODE'
          RETURN
        ELSEIF(FZMODE) THEN
          WRITE(IPUNLP,1055)
        ELSE
          WRITE(IPUNLP,1005)
        ENDIF
      ENDIF
C
C         INITIALIZE ALL THE "PREVIOUS" POINT QUANTITIES
C
      FNLP = FNLPBR
C
      RWORK(LCYVEC:LCYVEC+NVAR-1) = YBAR(1:NVAR)
      RWORK(LCCVEC:LCCVEC+MSUBE-1) = CBAR(1:MSUBE)
      RWORK(LCBVEC:LCBVEC+MSUBB-1) = BBAR(1:MSUBB)
C
      MUOPTN = ABS(MUCALC)
      REEDOO = MUCALC.GT.ZERO
      SMTOLD = BIGNUM
C
C         COMPUTE INITIAL ESTIMATE FOR MULTIPLIERS AND BARRIER PARAMETER
C
 120  CONTINUE
C
      IF(RELAX.AND.(NSLK-MINEQL).GT.0) THEN
        GVEC(NFREE+MINEQL+1:NFREE+NSLK) = ZERO
      ENDIF
C
      IF(MSUBB.GT.0.AND.NFREE.NE.MEQUAL) THEN
        PENMU = PMUBND
      ELSE
        MUOPTN = 1
        PENMU = ZEROOT
      ENDIF
C
      NCALL = 0
C
      IF(FZMODE) THEN
        MSUBEM = 0
      ELSE
        MSUBEM = MSUBE
      ENDIF
C
C         COMPUTE THE MULTIPLIER ESTIMATES UNLESS VALUES HAVE BEEN INPUT
C
      IF(MUCALC.LT.4) THEN
C
        CALL LSQMLT(MUOPTN,NCALL,CBAR,MSUBEM,MAXCON,ETABAR,CMAT,
     $    IROWC,JCOLC,NONZC,BBAR,MSUBB,MAXBND,VLAMBR,BMAT,IROWB,
     $    JCOLB,NONZB,GVEC,NVAR,PENMU,RWORK(LCSCRT),RWORK(LCRHLD),
     $    LNRHLD,IPUNLP,IOFLAG,CNDNUM,
     $    GRDLNM,IERLSM,NEEDED)
C
        IF(IERLSM.EQ.-67.OR.IERLSM.EQ.-68.OR.IERLSM.EQ.-906.OR.
     $   IERLSM.EQ.-924) THEN
          NEEDED = LCRHLD - 1 + NEEDED
          IERNLP = -131
          GO TO 200
        ENDIF
C
      ENDIF
C
C           COMPUTE ERROR IN KT CONDITIONS
C
      ERRCVC = DAMAX(MSUBEM,CBAR,1)
      IF(RELAX) ERRCVC = MAX(ERRCVC,GRDLNM)
      PMUCVC = ERRCVC/(TEN*1.1D0)
C
C         CHECK THE MINIMUM VALUES FOR THE MULTIPLIERS AND BARRIER
C         PARAMETER
C
      IF((PENMU.LT.PMUBND.AND.MSUBB.GT.0
     $   .OR.(IERLSM.EQ.-999)).AND.REEDOO) THEN
        REEDOO = .FALSE.
        IF(PENMU.GT.PMUCVC) THEN
          PMUBND = MAX(ZEROOT,PENMU)
        ELSE
          PMUBND = MAX(ZEROOT,PMUCVC,PMUBND)
        ENDIF
        IF(MUCALC.GT.0) THEN
          PMUBND = MAX(PMUBND,PMULWR)
          MUOPTN = 1
          GO TO 120
        ENDIF
      ENDIF
C
C         TRUNCATE THE BARRIER PARAMETER ESTIMATE
C
      IF(MSUBB.GT.0) PENMU = MAX(PMUBND,PENMU)
C
      BIGMLT = DAMAX(MSUBEM,ETABAR,1).GT.ONEP12
      IF(IERLSM.EQ.-999.OR.BIGMLT) THEN
C
C         INITIAL MULTIPLIER CALCULATION FAILED:
C
        IF(FRSTRY.AND..NOT.RELAX) THEN
C
          FRSTRY = .FALSE.
          ITERM = 0
          FZMODE = .TRUE.
          OPMODE = .NOT.FZMODE
C
          IF(IOFLAG.GE.10) THEN
            WRITE(IPUNLP,1016)
          ENDIF
C
          GO TO 190
C
        ELSEIF(.NOT.RELAX) THEN
C
C         IF INITIAL MULTIPLIER CALCULATION FAILED AND RELAXATION HAS
C         NOT BEEN TRIED, TRY IT
C
          IF(MSUBB.GT.0) THEN
            PMULST = PENMU
          ELSE
            PMULST = PMULWR
          ENDIF
C
          IF(IOFLAG.GE.10) THEN
            WRITE(IPUNLP,1015)
          ENDIF
C
          RELAX = .TRUE.
          CALL FLIPER(RELAX,
     $      NVARN,NONZCN,MSUBBN,NONZBN,NONZWN,
     $      NVARR,NONZCR,MSUBBR,NONZBR,NONZWR,
     $      NVAR,NONZC,MSUBB,NONZB,NONZW)
          NSLK = NVAR - NFREE
C
C         MODIFY VARIABLES AND CONSTRAINTS
C
          CALL SLAKER(BBAR,CBAR,YBAR,MAXBND,MAXCON,MEQUAL,MINEQL,   
     $       NALWR,NAUPR,NFREE,NVARR,NXLWR,NXUPR)   
C
C         MODIFY JACOBIAN AND HESSIAN
C
          JCOLC(NVAR+1) = NONZC + 1
          JCOLB(NVAR+1) = NONZB + 1
          JCOLW(NVAR+1) = NONZW + 1
C
          IF(IOFLAG.GE.10) THEN
C
            WRITE(IPUNLP,1006)
C
            MSUBX = NXLWR + NXUPR
            MSUBS = NALWR + NAUPR
            IF(RELAX) THEN
              NVARP = NVARR
              NSLKP = NSLK
              MSUBBP = MSUBBR
            ELSE
              NVARP = NVARN
              NSLKP = NSLK
              MSUBBP = MSUBBN
            ENDIF
            WRITE(IPUNLP,1007) NDIM,NVARP,
     $                     NFIX,NSLKP,
     $                     NFREE,NFREE,
     $                     MCON,MSUBE,
     $                     MEQUAL,MSUBBP,
     $                     MINEQL,MSUBX,
     $                     MIGNOR,MSUBS
            WRITE(IPUNLP,1006)
          ENDIF
C
          GO TO 110
C
        ELSE
C
C         MULTIPLIER CALCULATION FAILED AND RELAXATION IS ACTIVE
C
          ITERM = 9
C
          GO TO 130
C
        ENDIF
C
      ELSEIF(IERLSM.NE.0) THEN
        PRINT *,'IERLSM =',IERLSM
        RETURN
      ENDIF
C
      RWORK(LCETAV:LCETAV+MSUBEM-1) = ETABAR(1:MSUBEM)
      RWORK(LCVLAM:LCVLAM+MSUBB-1) = VLAMBR(1:MSUBB)
C
C         IF RELAXATION MODE IS ACTIVE, DEFINE EXACT PENALTY PARAMETER
C
      IF(RELAX) THEN
C
        RHOFAC = TEN
        PENRHO = MAX(SMLRHO,RHOFAC*DAMAX(MSUBEM,ETABAR,1))
C
C         RELAXATION OBJECTIVE FUNCTION
C         MODIFY THE INITIAL GRADIENT
C
        SMTUVW = 0.0D0
        DO I=NFREE+MINEQL+1,NFREE+NSLK
          SMTUVW = SMTUVW + YBAR(I)
          GVEC(I) = PENRHO
        ENDDO
C
        IF(IOFLAG.GE.10) WRITE(IPUNLP,1023) PENRHO
C
      ELSE
C
        PENRHO = ONEEP6
        SMTUVW = ZERO
C
      ENDIF
C
 130  CONTINUE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------  CONSTRAINED OPTIMIZATION BLOCK  --------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         ALGORITHM TERMINATION TESTS
C
      IF(ITERM.NE.0) THEN
C
        IF(ITERM.EQ.1) THEN
C
C         NORMAL TERMINATION FROM BARRIER MINIMIZATION,
C         TERMINATE THE ALGORITHM NORMALLY
C
          IERNLP = 0
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.2) THEN
C
C         WEAK SOLUTION -- OPTIMALITY CONDITIONS SATISFIED, BUT
C         MULTIPLIERS ARE NEAR ZERO
C
          IERNLP = +101
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.3) THEN
C
C          SMALL STEP TERMINATION
C
          IERNLP = +105
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.4) THEN
C
C         MAX NO. OF ITERATIONS IN SPRBAR
C
          IERNLP = +106
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.5) THEN
C
C         FUNCTION ERROR DURING GRADIENT EVALUATION
C
          IERNLP = -129
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.6) THEN
C
C         MAX NO. OF FUNCTION EVALUATIONS
C
          IERNLP = +104
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.7) THEN
C
C         MAXIMUM NUMBER OF LINESEARCH STEPS WITH NO ACCEPTABLE FILTER POINT
C
          IERNLP = +116
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.8) THEN
C
C         FEASIBLE POINT NOT FOUND
C
          IERNLP = +108
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.9) THEN
C
C         RANK DEFICIENT CONSTRAINTS
C
          IERNLP = -133
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.10) THEN
C
C         ALGORITHM CYCLING WITH NO REDUCTION IN KKT CONDITIONS
C
          IERNLP = -133
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-1) THEN
C
C         TEST FOR LINEAR REDUCED OBJECTIVE FUNCTION
C
          IERNLP = +117
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-4) THEN
C
C         TEST FOR MAXIMUM NUMBER OF INTERVAL HALVES IN LINE SEARCH
C
          IERNLP = +109
          IF(IOFLAG.GE.10) WRITE(IPUNLP,1014)
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-6) THEN
C
C         UNEXPECTED ERROR IN SRCHDR
C
          IERNLP = +114
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-8) THEN
C
C         INSUFFICIENT REAL STORAGE FOR SRCHDR 
C
          NEEDED = LCRHLD - 1 + NEEDED
          IERNLP = -131
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-10) THEN
C
C         STOP AFTER SEARCH DIRECTION CALCULATION FAILED
C
          IERNLP = +111
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-11) THEN
C
C         STOP AFTER DIAGNOSTIC LINE SEARCH
C
          IERNLP = +119
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-12) THEN
C
C         INSUFFICIENT INTEGER STORAGE FOR SRCHDR 
C
          NEEDED = LCISCR - 1 + NEEDED
          IERNLP = -132
          STOPIT = .TRUE.
C
        ELSEIF(ITERM.EQ.-13) THEN
C
C         I/O ERROR (INSUFFICIENT DISK SPACE)
C
          IERNLP = -153
          STOPIT = .TRUE.
C
        ELSE
          PRINT *,'SPRBAR ERROR RETURN--ITERM',ITERM
          RETURN
        ENDIF
C
      ENDIF
C
C         IF STOPIT IS TRUE TERMINATE ALGORITHM
C
      IF(STOPIT) GO TO 200
C
C ----------------------------------------------------------------------
C
C     
C         SET ALGORITHM CONTROL PARAMETERS
C
      IREVRS(5) = 1
C
 140  CONTINUE
C
C         CALL BARRIER SQP MINIMIZATION ALGORITHM
C
      IF(FZMODE) THEN
C
        IF(MEQUAL.EQ.NFREE) THEN
          MSUBBF = 0
        ELSE
          MSUBBF = MSUBB
        ENDIF
C
        CALL BARSQP( FNLPBR ,FNLP  ,YBAR   ,RWORK(LCYVEC)             
     $   ,RWORK(LCDLYV)   ,GVEC  ,RWORK(LCGRDL) ,GRDLNM ,NVAR ,NFREE        
     $   ,MAXCON   ,MSUBE ,RWORK(LCDLZV) ,CMAT  ,IROWC  ,JCOLC
     $   ,NONZC    ,CBAR  ,RWORK(LCCVEC) ,ETABAR        ,RWORK(LCETAV)  
     $   ,RWORK(LCDLET)   ,RWORK(LCFLTR) ,MXFLTR ,0     ,0        
     $   ,MAXCON ,CMAT    ,IROWC ,JCOLC ,NONZC  ,BBAR   ,RWORK(LCBVEC)  
     $   ,VLAMBR ,RWORK(LCVLAM)  ,RWORK(LCDLLM) ,MSUBBF ,MAXBND        
     $   ,BMAT    ,IROWB  ,JCOLB ,NONZB  ,WMAT   ,IROWW ,JCOLW     
     $   ,NONZW   ,PENMU  ,PENRHO,SMTUVW     
     $   ,RWORK(LCRHLD)   ,LNRHLD,RWORK(LCSCRT) ,LNSCRT ,IWORK(LCISCR) 
     $   ,LNISCR  ,NEEDED ,IREVRS,IRVCOM ,IFERR ,FZMODE  
     $   ,RELAX   ,ITERM  ,IQPTRM)
C
      ELSE
C
        CALL BARSQP( FNLPBR ,FNLP  ,YBAR   ,RWORK(LCYVEC)             
     $   ,RWORK(LCDLYV)   ,GVEC  ,RWORK(LCGRDL) ,GRDLNM ,NVAR ,NFREE        
     $   ,MAXRES   ,NRES  ,RWORK(LCDLZV) ,RMAT  ,IROWR  ,JCOLR
     $   ,NONZR    ,CBAR  ,RWORK(LCCVEC) ,ETABAR        ,RWORK(LCETAV)  
     $   ,RWORK(LCDLET)   ,RWORK(LCFLTR) ,MXFLTR ,MSUBE ,MINEQL        
     $   ,MAXCON ,CMAT    ,IROWC ,JCOLC ,NONZC  ,BBAR  ,RWORK(LCBVEC)  
     $   ,VLAMBR ,RWORK(LCVLAM)  ,RWORK(LCDLLM) ,MSUBB ,MAXBND        
     $   ,BMAT    ,IROWB  ,JCOLB ,NONZB  ,WMAT   ,IROWW ,JCOLW     
     $   ,NONZW   ,PENMU  ,PENRHO,SMTUVW     
     $   ,RWORK(LCRHLD)   ,LNRHLD,RWORK(LCSCRT) ,LNSCRT ,IWORK(LCISCR) 
     $   ,LNISCR  ,NEEDED ,IREVRS,IRVCOM ,IFERR ,FZMODE  
     $   ,RELAX   ,ITERM  ,IQPTRM)
C
      ENDIF
C
C         CHECK REVERSE COMMUNICATION FLAG
C
      IF(IREVRS(1).EQ.1) THEN
C
C         FUNCTION EVALUATION REQUESTED
C
        GO TO 170
C
      ELSEIF(IREVRS(2).GE.1) THEN
C
C         GRADIENT EVALUATION REQUESTED
C
        INSTAT(2) = INSTAT(2) + 1
        GO TO 160
C
      ELSEIF(IREVRS(3).GE.1) THEN
C
C         HESSIAN EVALUATION REQUESTED
C
        GO TO 150
C
      ELSE
C
C         BARRIER MINIMIZATION HAS TERMINATED
C         KEY:                              
C         ITERM =  9  RANK DEFICIENT CONSTRAINTS (PENMU SMALL, |ETA| LARGE)
C         ITERM =  8  CONSTRAINTS ARE INFEASIBLE (PENMU SMALL, |C| LARGE)
C         ITERM =  7  MAX. LINESEARCH STEPS WITH NO ACCEPTABLE FILTER POINT
C         ITERM =  6  MAX. NO. OF FUNC. EVALS.
C         ITERM =  5  FUNC. ERROR ON GRADIENT CALL (SPRBAR)
C         ITERM =  4  MAX. NO. ITERATIONS (SPRBAR)
C         ITERM =  3  SMALL STEP (SPRBAR)
C         ITERM =  2  WEAK SOLUTION OBTAINED
C         ITERM =  1  NORMAL TERMINATION FROM SPRBAR
C         ITERM =  0  CONTINUE ITERATIONS 
C         ITERM = -1  LINEAR REDUCED OBJECTIVE FUNCTION 
C         ITERM = -4  MAX. NO. INTERVAL HALVES IN LINE SEARCH 
C         ITERM = -6  SCHUR-QP FAILED WITH UNEXPECTED ERROR 
C         ITERM = -8  INSUFFICIENT REAL STORAGE FOR SRCHDR 
C         ITERM = -9  SEARCH DIR. FAILED--SLOPE CONDITION 
C         ITERM = -10 SEARCH DIR. FAILED--INCONS. CONSTR. 
C         ITERM = -11 TERMINATION AFTER DIAGNOSTIC (LINE SEARCH)
C         ITERM = -12 INSUFFICIENT INTEGER STORAGE FOR SRCHDR 
C         ITERM = -13 I/O ERROR (INSUFFICIENT DISK SPACE) 
C
        GO TO 180
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ------------------- HESSIAN EVALUATION SEQUENCE ----------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         THE FOLLOWING SEQUENCE OF OPERATIONS IS PERFORMED WHENEVER
C         IREVRS(3) .GE. 1 TO EVALUATE THE HESSIAN OF THE LAGRANGIAN 
C
 150  CONTINUE
C
C         TURN OFF PRINT EXCEPT WHEN IOFLAG = 30 
C
      IF(IOFLAG.LE.20) THEN
        IREVRS(4) = 0
      ELSE
        IREVRS(4) = 1
      ENDIF
C
C         SET HES RETURN POINT FLAG (IFESR = 1) AND RETURN
C
      IFESR = 1
      IRVCOM = 1
C
      GO TO 10000
C
C ******* FUNCTION EVALUATION RETURN POINT
C
 501  CONTINUE
C
      IF(RELAX) THEN
C
C         RELAXATION OBJECTIVE FUNCTION CONTRIBUTION
C
        SMTUVW = 0.0D0
        DO I=NFREE+MINEQL+1,NFREE+NSLK
          SMTUVW = SMTUVW + YBAR(I)
        ENDDO
C
      ELSE
C
        SMTUVW = ZERO
C
      ENDIF
C
C         FUNCTION EVALUATION COMPLETED.
C         RETURN TO APPROPRIATE PLACE IN THE PROGRAM
C
      GO TO 140
C
C ----------------------------------------------------------------------
C -------------------- END OF HESSIAN EVALUATION SEQUENCE --------------
C ----------------------------------------------------------------------
C
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ------------------- GRADIENT EVALUATION SEQUENCE ---------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         THE FOLLOWING SEQUENCE OF OPERATIONS IS PERFORMED WHENEVER
C         IREVRS(2) .GE. 1 TO EVALUATE THE GRADIENT OF THE OBJECTIVE 
C         AND CONSTRAINT FUNCTIONS
C
 160  CONTINUE
C
C         TURN OFF PRINT EXCEPT WHEN IOFLAG = 30 
C
      IF(IOFLAG.LE.20) THEN
        IREVRS(4) = 0
      ELSE
        IREVRS(4) = 1
      ENDIF
C
C         SET RETURN POINT FLAG (IFESR = 2) AND RETURN
C
      IFESR = 2
      IRVCOM = 1
      GO TO 10000
C
C ******* GRADIENT EVALUATION RETURN POINT
C
 502  CONTINUE
C 
C         GRADIENT EVALUATION COMPLETED.
C
      IF(RELAX.AND.(NSLK-MINEQL).GT.0) THEN
        GVEC(NFREE+MINEQL+1:NFREE+NSLK) = PENRHO
      ENDIF
C
C         RETURN TO APPROPRIATE PLACE IN THE PROGRAM
C
      GO TO 140
C
C ----------------------------------------------------------------------
C -------------------- END OF GRADIENT EVALUATION SEQUENCE -------------
C ----------------------------------------------------------------------
C
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ------------------- FUNCTION EVALUATION SEQUENCE ---------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
 170  CONTINUE
C
C
C         TURN OFF PRINT WHEN IOFLAG = 10
C
      IF(IOFLAG.LE.10) IREVRS(4) = 0
C
C         EVALUATE FUNCTIONS AT INITIAL POINT.  
C
      ICES = 1
      IREVRS(2) = 0
C
C         SET CEES RETURN POINT FLAG (IFESR = 3) AND RETURN
C
      IFESR = 3
      IRVCOM = 1
      GO TO 10000
C
C ******* FUNCTION EVALUATION RETURN POINT
C
 503  CONTINUE
C
      IF(RELAX) THEN
C
C         RELAXATION OBJECTIVE FUNCTION CONTRIBUTION
C
        SMTUVW = 0.0D0
        DO I=NFREE+MINEQL+1,NFREE+NSLK
          SMTUVW = SMTUVW + YBAR(I)
        ENDDO
C
      ELSE
C
        SMTUVW = ZERO
C
      ENDIF
C
C         RETURN TO APPROPRIATE PLACE IN THE PROGRAM
C
      GO TO 140
C
C ----------------------------------------------------------------------
C -------------------- END OF FUNCTION EVALUATION SEQUENCE -------------
C ----------------------------------------------------------------------
C
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ------------------- BARRIER SQP ALGORITHM TERMINATION SEQUENCE -------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
 180  CONTINUE
C
C         CHECK FOR ABNORMAL TERMINATION CONDITIONS
C
      IF(ITERM.EQ.-1.OR.ITERM.EQ.-4.OR.ITERM.EQ.-6
     $   .OR.ITERM.EQ.-8.OR.ITERM.EQ.-11.OR.ITERM.EQ.-12
     $   .OR.ITERM.EQ.-13) GO TO 130
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      IF(FZMODE) THEN
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         FEASIBILITY MODE TERMINATED
C         SET FEASIBILITY MODE TERMINATION FLAGS
C
        IF(ITERM.EQ.1.OR.(ITERM.EQ.2.AND.MEQUAL.NE.NFREE)) THEN
C
          IF(IOFLAG.GE.10) THEN
            WRITE(IPUNLP,1004)
            WRITE(IPUNLP,1001)
          ENDIF
C
C         NORMAL TERMINATION AT MINIMUM SUM OF VIOLATIONS
C
          FEEZPT = IFERR.EQ.-101
C
        ELSE
C
          IF(IOFLAG.GE.10) THEN
            WRITE(IPUNLP,1003)
            WRITE(IPUNLP,1001)
          ENDIF
C
        ENDIF
C
C         WHAT NEXT?
C
        IF(FEEZPT.OR.(ITERM.EQ.2.AND.MEQUAL.NE.NFREE)) THEN
C
C         ---- POINT IS FEASIBLE AND/OR SUFFICIENTLY GOOD TO PROCEED ----
C
          IF(OPTPRB.AND.NFREE.NE.MEQUAL) THEN
C
C           OPTIMIZATION IS REQUIRED --> BEGIN OPTIMIZATION MODE
C
            ITERM = 0
            FZMODE = .FALSE.
            OPMODE = .NOT.FZMODE
C
C           RESET RELAXATION FLAG IF POSSIBLE
C
            IF(IRELAX.NE.2.AND.RELAX) THEN
C
              RELAX = .FALSE.
C
              CALL FLIPER(RELAX,
     $          NVARN,NONZCN,MSUBBN,NONZBN,NONZWN,
     $          NVARR,NONZCR,MSUBBR,NONZBR,NONZWR,
     $          NVAR,NONZC,MSUBB,NONZB,NONZW)
              NSLK = NVAR - NFREE
C
C                 MODIFY JACOBIAN AND HESSIAN
C
              JCOLC(NVAR+1) = NONZC + 1
              JCOLB(NVAR+1) = NONZB + 1
              JCOLW(NVAR+1) = NONZW + 1
C
              IF(IOFLAG.GE.10) THEN
C
                WRITE(IPUNLP,1006)
C
                MSUBX = NXLWR + NXUPR
                MSUBS = NALWR + NAUPR
                IF(RELAX) THEN
                  NVARP = NVARR
                  NSLKP = NSLK
                  MSUBBP = MSUBBR
                ELSE
                  NVARP = NVARN
                  NSLKP = NSLK
                  MSUBBP = MSUBBN
                ENDIF
                WRITE(IPUNLP,1007) NDIM,NVARP,
     $                   NFIX,NSLKP,
     $                   NFREE,NFREE,
     $                   MCON,MSUBE,
     $                   MEQUAL,MSUBBP,
     $                   MINEQL,MSUBX,
     $                   MIGNOR,MSUBS
                WRITE(IPUNLP,1006)
C
              ENDIF
C
            ENDIF
C
C           COMPUTE REAL OBJECTIVE AND GRADIENT, AND THEN REENTER OPTIMIZATION
C
            GO TO 190
C
          ELSE
C
C           OPTIMIZATION IS NOT REQUIRED --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
            GO TO 130
C
          ENDIF
C
        ELSE
C
C         ---- POINT IS NOT FEASIBLE ----
C
          IF(RELAX) THEN
C
C           RELAXATION IS ACTIVE AND DID NOT YIELD A FEASIBLE POINT
C               --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
            IF(ITERM.NE.5) ITERM = 8
            GO TO 130
C
C
          ELSE
C
C           RELAXATION IS NOT CURRENTLY ACTIVE
C
            IF(IRELAX.GT.0.AND.MSUBE.NE.0.AND.MSUBB.NE.0) THEN
C
C             RELAXATION CAN BE ATTEMPTED
C
              RELAX = (ITERM.EQ.9.OR.ITERM.EQ.8.OR.ITERM.EQ.7.OR
     $           .ITERM.EQ.-9.OR.ITERM.EQ.-10 ).AND.NFREE.NE.MEQUAL
              IF(ITERM.EQ.2.AND.NFREE.EQ.MEQUAL) THEN
                RELAX = .TRUE.
                SMLRHO = MAX(RHOLWR,DAMAX(MSUBE,CBAR,1))
                OPMODE = .TRUE.
                FZMODE = .FALSE.
              ENDIF
C
              IF(RELAX) THEN
C
                IF(IOFLAG.GE.10) THEN
                  WRITE(IPUNLP,1021)
                ENDIF
C
                ITERM = 0
                IF(MSUBB.GT.0) THEN
                  PMULST = PENMU
                ELSE
                  PMULST = PMULWR
                ENDIF
C
                CALL FLIPER(RELAX,
     $            NVARN,NONZCN,MSUBBN,NONZBN,NONZWN,
     $            NVARR,NONZCR,MSUBBR,NONZBR,NONZWR,
     $            NVAR,NONZC,MSUBB,NONZB,NONZW)
                NSLK = NVAR - NFREE
C
C                 MODIFY VARIABLES AND CONSTRAINTS
C
                CALL SLAKER(BBAR,CBAR,YBAR,MAXBND,MAXCON,MEQUAL,MINEQL,   
     $            NALWR,NAUPR,NFREE,NVARR,NXLWR,NXUPR)   
C
C                 MODIFY JACOBIAN AND HESSIAN
C
                JCOLC(NVAR+1) = NONZC + 1
                JCOLB(NVAR+1) = NONZB + 1
                JCOLW(NVAR+1) = NONZW + 1
C
                IF(IOFLAG.GE.10) THEN
C
                  WRITE(IPUNLP,1006)
C
                  MSUBX = NXLWR + NXUPR
                  MSUBS = NALWR + NAUPR
                  IF(RELAX) THEN
                    NVARP = NVARR
                    NSLKP = NSLK
                    MSUBBP = MSUBBR
                  ELSE
                    NVARP = NVARN
                    NSLKP = NSLK
                    MSUBBP = MSUBBN
                  ENDIF
                  WRITE(IPUNLP,1007) NDIM,NVARP,
     $                     NFIX,NSLKP,
     $                     NFREE,NFREE,
     $                     MCON,MSUBE,
     $                     MEQUAL,MSUBBP,
     $                     MINEQL,MSUBX,
     $                     MIGNOR,MSUBS
                  WRITE(IPUNLP,1006)
                ENDIF
C
                GO TO 110
C
              ELSE
C
C               RELAXATION NOT WORTHWHILE
C               --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
                IF(ITERM.NE.5) ITERM = 8
                GO TO 130
C
              ENDIF
C
            ELSE
C
C             RELAXATION CANNOT BE ATTEMPTED
C             --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
              IF(ITERM.NE.5) ITERM = 8
              GO TO 130
C
            ENDIF
C
          ENDIF
C
        ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      ELSEIF(OPMODE) THEN
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         OPTIMIZATION MODE TERMINATED
C
        IF(IOFLAG.GE.10) THEN
          WRITE(IPUNLP,1002)
          WRITE(IPUNLP,1001)
        ENDIF
C
        IF(ITERM.EQ.1) THEN
C
          IF(RELAX) THEN
C
            IF(SMTUVW.GT.DBLE(NSLK)*CONTOL) THEN
C
C             NORMAL TERMINATION AT NONZERO SUM OF VIOLATIONS
C
              IF(.NOT.OPTPRB) THEN
                ITERM = 8
                GO TO 130
              ENDIF
C
            ELSE
C
C             CHECK FOR WEAK SOLUTION (REQUIRES RELAXATION)
C
              ITERM = 2
              GO TO 130
C
            ENDIF
C
          ELSE
C
C             NORMAL TERMINATION AT THE SOLUTION
C
            IF(OPTPRB) GO TO 130
C
          ENDIF
C
        ENDIF
C
C       ---------------------------------------------------
C       ABNORMAL TERMINATION IN OPTIMIZATION MODE
C       ---------------------------------------------------
C
C           CHECK KKT CONDITIONS
C
        CERR = ZERO
        DO I=1,MSUBE
          CERR = MAX(CERR,ABS(CBAR(I)))
        ENDDO
C
C           CHECK CONSTRAINT FEASIBILITY
C
        IF(CERR.GT.CONTOL.AND.GRDLNM.LT.CERR.AND..NOT.RELAX) THEN
C
C           IS IT WORTH FINDING A FEASIBLE POINT?
C
C           COMPUTE GLOBAL FILTER ENTRIES
C
          ETAMAX = DAMAX(MSUBE,ETABAR,1)
C
          CALL FILTER(GRDLNM,CERR,BIGCON,ZERO,ZERO,ETAMAX,ZERO,
     $        ZERO,FLTRGB,LNFLGB,MXFLGB,NWFLGB,ISCRGB,LNISGB)
C
C           CHECK POINT AGAINST GLOBAL FILTER.  IF IT IS NOT ACCEPTED STOP
C
          IF(NWFLGB.NE.0.OR.LNFLGB.EQ.10) THEN
            IF(IOFLAG.GE.20) THEN
              WRITE(IPUNLP,1008)
              CALL FLTRPR('ALG',ZERO,ZERO,FLTRGB,LNFLGB,MXFLGB)
            ENDIF
            ITERM = 10
            GO TO 130
          ENDIF
          IF(ITERM.EQ.4.OR.ITERM.EQ.5.OR.ITERM.EQ.6.OR. 
     $         (ITERM.LT.0.AND.ITERM.NE.-10) ) THEN
C
C               TRYING TO FIND A FEASIBLE POINT IS NOT WORTHWHILE
C               --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
            GO TO 130
C
          ELSE
C
C              TRY TO FIND A FEASIBLE POINT
C
            ITERM = 0
            FZMODE = .TRUE.
            OPMODE = .NOT.FZMODE
C
            IF(IOFLAG.GE.10) THEN
              WRITE(IPUNLP,1017)
            ENDIF
C
            GO TO 190
C
          ENDIF
C
        ENDIF
C
        IF(RELAX) THEN
C
C           RELAXATION IS ACTIVE
C
          IF(ITERM.EQ.3.OR.ITERM.EQ.7) THEN
C
C               INCREASING RELAXATION PENALTY NOT WORTHWHILE
C               --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
            GO TO 130
C
          ELSEIF(PENRHO.LT.ONEEP6.AND.SMTUVW.LT.POINT5*SMTOLD) THEN
C
C             INCREASE PENALTY PARAMETER AND RESET ITERM
C
            MUOPTN = 2
            REEDOO = .TRUE.
            SMTOLD = SMTUVW
            ITERM = 0
            SMLRHO = MIN(ONEEP6,TEN*PENRHO)
            IF(IOFLAG.GE.10) WRITE(IPUNLP,1022)
C
            GO TO 120
C
          ENDIF
C
        ELSE
C
C           RELAXATION IS NOT CURRENTLY ACTIVE
C           SWITCH TO RELAXATION MODE IF GRAD L > |C|
C
          IF(IRELAX.GT.0.AND.MSUBE.NE.0.AND.MSUBB.NE.0
     $       .AND.GRDLNM.GT.CERR) THEN
C
C             RELAXATION CAN BE ATTEMPTED
C
            RELAX = ITERM.EQ.9.OR.ITERM.EQ.8.OR.ITERM.EQ.7.OR
     $                 .ITERM.EQ.-9.OR.ITERM.EQ.-10
C
            IF(RELAX) THEN
C
              IF(IOFLAG.GE.10) THEN
                WRITE(IPUNLP,1021)
              ENDIF
C
              ITERM = 0
              IF(MSUBB.GT.0) THEN
                PMULST = PENMU
              ELSE
                PMULST = PMULWR
              ENDIF
C
              CALL FLIPER(RELAX,
     $          NVARN,NONZCN,MSUBBN,NONZBN,NONZWN,
     $          NVARR,NONZCR,MSUBBR,NONZBR,NONZWR,
     $          NVAR,NONZC,MSUBB,NONZB,NONZW)
              NSLK = NVAR - NFREE
C
C                 MODIFY VARIABLES AND CONSTRAINTS
C
              CALL SLAKER(BBAR,CBAR,YBAR,MAXBND,MAXCON,MEQUAL,MINEQL,   
     $          NALWR,NAUPR,NFREE,NVARR,NXLWR,NXUPR)   
C
C                 MODIFY JACOBIAN AND HESSIAN
C
              JCOLC(NVAR+1) = NONZC + 1
              JCOLB(NVAR+1) = NONZB + 1
              JCOLW(NVAR+1) = NONZW + 1
C
              IF(IOFLAG.GE.10) THEN
C
                WRITE(IPUNLP,1006)
C
                MSUBX = NXLWR + NXUPR
                MSUBS = NALWR + NAUPR
                IF(RELAX) THEN
                  NVARP = NVARR
                  NSLKP = NSLK
                  MSUBBP = MSUBBR
                ELSE
                  NVARP = NVARN
                  NSLKP = NSLK
                  MSUBBP = MSUBBN
                ENDIF
                WRITE(IPUNLP,1007) NDIM,NVARP,
     $                   NFIX,NSLKP,
     $                   NFREE,NFREE,
     $                   MCON,MSUBE,
     $                   MEQUAL,MSUBBP,
     $                   MINEQL,MSUBX,
     $                   MIGNOR,MSUBS
                WRITE(IPUNLP,1006)
              ENDIF
C
              GO TO 110
C
            ELSE
C
C               RELAXATION NOT WORTHWHILE
C               --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
              GO TO 130
C
            ENDIF
C
          ELSE
C
C             RELAXATION CANNOT BE ATTEMPTED
C             --> TRANSFER TO ALGORITHM TERMINATION TESTS
C
            GO TO 130
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
 190  CONTINUE
C
      IFESR = 4
      IRVCOM = 1
      IREVRS(1) = 1
      IREVRS(2) = 1
      IREVRS(4) = 0
      GO TO 10000
C
C ******* FUNCTION EVALUATION RETURN POINT
C
 504  CONTINUE
      IREVRS(1) = 0
      IREVRS(2) = 0
      IREVRS(4) = 0
C
      GO TO 110
C
C
C
C ----------------------------------------------------------------------
C
C         ALGORITHM TERMINATION PROCESSING
C
 200  CONTINUE
C
C
C ----------------------------------------------------------------------
C
C
      IRVCOM = 0
C
10000 CONTINUE
C
 1001 FORMAT(T3,'*',T106,'*',/2X,104('*'))
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'CONSTRAINED OPTIMIZATION COMPLE
     $TED',T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'NONLINEAR SEARCH TERMINATED ABN
     $ORMALLY',T106,'*')
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'CONSTRAINT SATISFACTION COMPLET
     $ED',T106,'*')
 1005 FORMAT(T3,'*',T106,'*',/T3,'*',T11,'CONSTRAINED OPTIMIZATION',
     $    T106,'*')
 1055 FORMAT(T3,'*',T106,'*',/T3,'*',T11,'CONSTRAINT SATISFACTION'
     $    ,T106,'*')
 1006 FORMAT(T3,'*',T106,'*',/2X,104('*'))
 1007 FORMAT(T3,'*',T106,'*',
     $  /T3,'*',T11,85('-'),T106,'*',
     $  /T3,'*',T11,'External Representation',T54,'|'
     $  ,T59,'Internal Representation',T106,'*',
     $  /T3,'*',T11,85('-'),T106,'*',
     $  /T3,'*',T11,'Variables',T41,'=',I6,T54,'|'
     $  ,T59,'Variables',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,'  Fixed',T41,'=',I6,T54,'|'
     $  ,T59,'  Slack',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,'  Free',T41,'=',I6,T54,'|'
     $  ,T59,'  Free',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,'Constraints',T41,'=',I6,T54,'|'
     $  ,T59,'Constraints',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,'  Equality',T41,'=',I6,T54,'|'
     $  ,T59,'Bounds',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,'  Inequality',T41,'=',I6,T54,'|'
     $  ,T59,'  Variable Bounds',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,'  Ignored',T41,'=',I6,T54,'|'
     $  ,T59,'  Slack Bounds',T89,'=',I6,T106,'*',
     $  /T3,'*',T11,85('-'),T106,'*')
 1008 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....PROGRESS REJECTED BY FILTE
     $R--TERMINATE ALGORITHM',T106,'*')
 1014 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....MAXIMUM NUMBER OF INTERVAL
     $ HALVES IN LINE SEARCH--TERMINATE ALGORITHM',T106,'*')
 1015 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....INITIAL MULTIPLIER CALCULA
     $TION FAILED;  SWITCH TO RELAXATION MODE',T106,'*')
 1016 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....INITIAL MULTIPLIER CALCULA
     $TION FAILED;  FIND FEASIBLE POINT',T106,'*')
 1017 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....OPTIMIZATION FAILED;  FIND
     $ FEASIBLE POINT',T106,'*')
 1021 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....SWITCH TO RELAXATION MODE'
     $   ,T106,'*')
 1022 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....RELAXATION FAILED;  INCREA
     $SE PENALTY PARAMETER',T106,'*')
 1023 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'Relaxation Penalty Parameter ='
     $  ,1PG16.8,T106,'*')
      RETURN
      END
