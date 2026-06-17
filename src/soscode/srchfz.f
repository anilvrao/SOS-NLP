
      SUBROUTINE SRCHFZ(GMAT,IROWG,JSTRG,NONZG,MCON,NDIM,
     $    WORK,NWORK,IWORK,NIWORK,NEEDED,ISTATC,CLWR,CUPR,
     $    CVEC,CBAR,COLD,SVEC,ISTATV,XLWR,XUPR,XVEC,XBAR,IPCSRC,
     $    IPU,IREVRS,IFERR,ISTART,ISTERM)
C
C ======================================================================
C     SRCHFZ===>srchfz   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C
C         PURPOSE:  COMPUTE THE VALUES OF THE NDIM VARIABLES XVEC WHICH
C                   SATISFY THE MCON CONSTRAINTS
C
C                      CLWR < CVEC < CUPR
C
C                   AND BOUNDS
C
C                      XLWR < XVEC < XUPR
C
C         INPUT:
C
C            GMAT   JACOBIAN MATRIX STORED AS AN ARRAY (NONZG)
C            IROWG  ROW INDEX OF NONZERO IN GMAT (NONZG)
C            JSTRG  COLUMN START INDEX (NDIM+1)
C            NONZG  NUMBER OF NONZEROS IN GMAT
C                   NOTE: NONZG  = JSTRG(NDIM+1)-1
C            MCON   NUMBER OF CONSTRAINTS  (GT.0)
C            NDIM   NUMBER OF VARIABLES
C            WORK   WORK ARRAY (NWORK)
C            NWORK  LENGTH OF WORK ARRAY (GT. 
C                   7*NDIM + 11*MCON + NONZG )
C            IWORK  INTEGER WORK ARRAY (NIWORK)
C            NIWORK LENGTH OF INTEGER WORK ARRAY (GT. 5*NDIM + 6*MCON 
C                   NONZG + 2 + MAX(NONZG+MCON,NDIM+MCON+1),MCON+8,9) 
C            ISTATC CONSTRAINT STATUS ARRAY (MCON)
C            CLWR   CONSTRAINT LOWER BOUND (MCON)
C            CUPR   CONSTRAINT UPPER BOUND (MCON)
C            CVEC   CONSTRAINTS AT XVEC (MCON)
C            CBAR   CONSTRAINTS AT XBAR (MCON)
C            COLD   CONSTRAINTS AT OLD POINT (MCON)
C            SVEC   SEARCH VECTOR (NDIM)
C            ISTATV VARIABLE STATUS ARRAY (NDIM)
C            XLWR   VARIABLE LOWER BOUND (NDIM)
C            XUPR   VARIABLE UPPER BOUND (NDIM)
C            XVEC   ESTIMATE OF INDEPENDENT VARIABLES (NDIM)
C            XBAR   NEW ESTIMATE OF XVEC (NDIM)
C            IPCSRC OUTPUT CONTROL FLAG
C                   = 0   NO OUTPUT
C                   = 10  NORMAL OUTPUT
C            IPU    OUTPUT UNIT NO.
C            IREVRS (1) =1 WHEN EVALUATING FUNCTIONS
C                   (2) .GE.1 WHEN EVALUATING GRADIENT INFORMATION
C                   (4) =1 WHEN PRINTING OUTPUT
C            IFERR  =1 WHEN FUNCTION EVALUATION IS IMPOSSIBLE
C            ISTART QP START OPTION
C            ISTERM TERMINATION FLAG
C                   = 0   NORMAL INPUT VALUE
C
C         OUTPUT:
C
C            CVEC   THE CONSTRAINTS AT XVEC
C            ISTART QP START OPTION
C            ISTERM TERMINATION FLAG
C                   = 0       NO TERMINATION
C                   = 1       NORMAL TERMINATION
C                   = 2       SMALL STEPS
C                   = 3       MAXIMUM NO. OF ITERATIONS
C                   = 4       FUNCTION ERROR ON GRADIENT 
C                   = 5       NOT ENOUGH REAL STORAGE ALLOCATED
C                   = 6       SINGULAR JACOBIAN ON SUCCESSIVE ITERS.
C                   = 7       NOT ENOUGH INTEGER STORAGE ALLOCATED
C                   = 8       MAX. NO. OF FUNC. EVALS.
C                   = 9       UPHILL DIRECTION IN LINE SEARCH
C                   = 10      FUNCTION ERROR AT INITIAL POINT
C                   = 11      I/O ERROR (INSUFFICIENT DISK SPACE)
C                   = 12      USER EXTERNAL KILL
C            NEEDED STORAGE REQUIRED WITH NWORK IS TOO SMALL (ISTERM=5)
C                   OR NIWORK IS TOO SMALL (ISTERM=7)
C            XVEC   THE SOLUTION POINT
C
C  ** SUBROUTINES REQUIRED: 
C
C        NOTE:   1.INPUT CONSTRAINTS CVEC EVALUATED AT XVEC MUST BE SUPPLIED
C                  ON FIRST ENTRY.  THEREAFTER CBAR EVALUATED AT XBAR
C                  IS SUPPLIED WHEN IREVRS(1) = 1
C
C
C     ******************************************************************
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,ONEEM4=1.0D-4,
     $    ONEEP2=1.0D2,ONEEP1=1.0D1,ONEEM1=1.0D-1,ONEEM3=1.0D-3,
     $    POINT2=2.0D-1,POINT5=5.0D-1)
C
      INCLUDE '../commons/NLPSPR.CMN'
C
      COMMON /ITEREF/ MAXREF,MAXRFN,IREFIN
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      COMMON /PERCOM/ ISQPER(20)
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
C         COMMONS FOR ITERATION LOG QUANTITIES
C
C     ------------------------------------------------------------------
      INCLUDE  '../commons/NLPOUT.CMN'
C     ------------------------------------------------------------------
C
      DIMENSION  GMAT(NONZG),IROWG(NONZG),JSTRG(NDIM+1),WORK(NWORK),
     $    IWORK(NIWORK),CVEC(MCON),CBAR(MCON),COLD(MCON),
     $    SVEC(NDIM),XVEC(NDIM),XBAR(NDIM),CLWR(MCON),CUPR(MCON),
     $    XLWR(NDIM),XUPR(NDIM),ISTATC(MCON),ISTATV(NDIM)
      DIMENSION IREVRS(5)
C
      LOGICAL CONTST,GOODGR
      CHARACTER(LEN=6) CHITER
      CHARACTER(LEN=22) CNDLBL
      CHARACTER(LEN=16) LFTVAL,RITVAL
C
C     ******************************************************************
C
      IF (IMAX(4,IREVRS,1).GE.1) THEN
        IF (IRTRN.EQ.1) THEN
          GO TO 501
        ELSEIF (IRTRN.EQ.2) THEN
          GO TO 502
        ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------- INITIALIZATION -------------------------------------
C ----------------------------------------------------------------------
C
      IER = 0
      IT = 1
      SMALL = ZEROOT
      OLDSTP = BIGNUM
      ALFA = ONE
      FREDUC = BIGNUM
      GOODGR = .FALSE.
C        LNRSTP < 0 ---> NO LINEAR CONSTRAINT STEP
C        LNRSTP = 0 ---> LINEAR CONSTRAINT STEP
      LNRSTP = 0
      IF(SFZTOL.LT.ZERO) LNRSTP = -1
      IF(QPOPTN.EQ.'SPARSE') THEN
        CNDLBL = 'Cond(K)...............'
      ELSE
        CNDLBL = 'Cond(G)...............'
      ENDIF
C
C         SEARCH MODE FLAG
C         = 1  STANDARD (LDP) SEARCH DIRECTION
C         = 2  RELAXED SUBPROBLEM
C
      ISRFLG = 1
C
C         SET SCHUR-QP OUTPUT CONTROL FLAG
C
      IF(IOFSHR.EQ.0) THEN
        IPCSHR = 0
        IF(IPCSRC.GE.30) IPCSHR = IPCSRC
      ELSE
        IPCSHR = IOFSHR
      ENDIF
C
C         SET SRCHFZ OUTPUT CONTROL FLAG
C
      IF(IOFSRC.EQ.0) THEN
        IPCFEZ = IPCSRC
      ELSE
        IPCFEZ = IOFSRC
      ENDIF
C
C             INITIALIZE TERSE OUTPUT IF NECESSARY
C
      NPRFLG = -1
C
C         NOTE THE FIRST NDIM+MCON ELEMENTS OF THE WORK
C         ARRAY ARE USED FOR BOTH THE QP OBJECTIVE FUNCTION
C         DATA AND WORKING STORAGE IN THE LINE SEARCH
C
      IF(IFERR.EQ.0) THEN
C
C         COMPUTE CONSTRAINT ERROR AT INITIAL POINT
C
        CALL FMERIT(CVEC,CLWR,CUPR,XVEC,XLWR,XUPR,MCON,NDIM,CONTOL,
     $    PHI,IFERR)
C
      ELSE
C
        ISTERM = 10
        GO TO 200
C
      ENDIF
C
C         INITIALIZE SEARCH DIRECTION FOR FIRST QP CALL -- SUBSEQUENT
C         CALLS TO SCHUR-QP WILL BE INITIALIZED WITH THE CURRENT VALUE OF
C         SVEC
C
      SVEC(1:NDIM) = ZERO
C
C ----------------------------------------------------------------------
C ----------------- BEGIN ITERATION ------------------------------------
C ----------------------------------------------------------------------
C
 110  CONTINUE
C
C             CHECK FOR CONSTRAINT SATISFACTION
C
      CALL CONSAT(CLWR,CUPR,CVEC,XLWR,XUPR,XVEC,MCON,NDIM,
     $    CONTOL,CONTST)
C
C             COMPUTE GRADIENT INFORMATION
C
      IRTRN = 1
      IF(IPCFEZ.GE.10.AND.IPCFEZ.LT.20.AND.IT.NE.1) THEN
        IREVRS(4) = 1
      ELSE
        IREVRS(4) = 0
      ENDIF
      IF(CONTST.OR.IT.EQ.1) THEN
        IREVRS(2) = 0
      ELSEIF(OLDSTP.LT.SMALL.OR.ALFA.LT.SMALL.OR.
     $    ABS(FREDUC).LT.SMALL) THEN
        IREVRS(2) = 2
        GOODGR = .TRUE.
      ELSE
        IREVRS(2) = 1
        GOODGR = .FALSE.
      ENDIF
      IF(IREVRS(2).EQ.0.AND.IREVRS(4).EQ.0) GO TO 501
      RETURN
C
 501  CONTINUE
      IREVRS(2) = 0
      IREVRS(4) = 0
C
      IF (IFERR.EQ.1) THEN
C  
C             GRADIENT CANNOT BE EVALUATED AT POINT WHERE FUNCTION
C             CAN, TERMINATE ALGORITHM.
C
        ISTERM = 4
        GO TO 200
C
      ELSEIF(IFERR.EQ.-100) THEN
C
C             MAX. NO. OF. FUNCTION EVALS.
C
        ISTERM = 8
        GO TO 200
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------- COMPUTE SEARCH DIRECTION AND MAGNITUDE--------------
C ----------------------------------------------------------------------
C
      IF(GOODGR.AND.IT.LE.9999) THEN
        CHITER = '(    )'
        WRITE(CHITER(2:5),'(I4)') IT
        IREFIN = 2
      ELSE
        CHITER = '      '
        WRITE(CHITER(1:6),'(I6)') IT
        IREFIN = 0
      ENDIF
      CALL HHSDEL(CHITER,' ',' ',6,NDEL,IERDEL)
C
      IF(CONTST.AND.IT.GT.NITMIN) THEN
C
C             ITERATION PRINT (AT FINAL POINT)
C
        SRPHI = SQRT(PHI*CONTOL)
        MACTIV = 0
        DO I=1,NDIM
          IF (ISTATV(I).GT.0)  MACTIV = MACTIV + 1
        ENDDO
        DO I=1,MCON
          IF (ISTATC(I).GE.1 .AND. ISTATC(I).LE.3)  MACTIV = MACTIV + 1
        ENDDO
        NDOF = NDIM - MACTIV
C
        IF(IPCFEZ.GE.10) THEN
          CALL HHCPRS(CHITER,' ','-')
          WRITE(IPU,1001) CHITER
          WRITE(LFTVAL,'(SP,1PG14.6)') SRPHI
          WRITE(RITVAL,'(SP,1PG14.6)') CNDNUM
          CALL DBLSTR('Constraint Error',LFTVAL,
     $               CNDLBL,RITVAL,IPU)
          WRITE(LFTVAL,'(I8)') MACTIV
          WRITE(RITVAL,'(I8)') NDOF
          CALL DBLSTR('Active Constraints',LFTVAL,
     $               'Degrees of Freedom',RITVAL,IPU)
          WRITE(LFTVAL,'(I8)') 0
          WRITE(RITVAL,'(I8)') 0
          CALL DBLSTR('QP Iterations',LFTVAL,
     $               'Matrix Factorizations',RITVAL,IPU)
        ELSEIF(IPCFEZ.GT.0) THEN
          CALL WRTURS(NPRFLG,CHITER,0,NDOF,ALFA,SNORM,SRPHI,0,XDUM,
     $    SRPHI,XDUM,XDUM,CNDNUM,IPU) 
        ENDIF
C
C       ----------------------------------------------------------------
C
C         SAVE ITERATION LOG QUANTITIES
C
        ITP    = IT  
        KEYMP  = 0 
        PPGNOR = ZERO
        PFNOM  = SRPHI  
        PERREQ = ZERO
        PERRIN = ZERO
        PFMRT  = SRPHI  
        PPENMA = ZERO
        PDIAGN = ZERO
        PCNDNU = CNDNUM
        PEIGMI = ZERO
        PEIGMA = ZERO
        MACTIP = MACTIV
        NDOFP  = NDOF  
        NQPITP = 0
        NUMKTP = 0
C
C       ----------------------------------------------------------------
C
C
C         SET NORMAL TERMINATION FLAG
C
        ISTERM = 1
        GO TO 200
C
      ENDIF
C
C         INITIALIZE K-T FACTORIZATION COUNT CLOCK
C
      IF(QPOPTN.EQ.'SPARSE') CALL CLKSET(13)
C
C         INITIALIZE NUMBER OF QP ITERATIONS 
C
      NQPITR = 1
C
 120  CONTINUE
C
      ISQPER(1:20) = 0
      INSTAT(12) = INSTAT(12) + 1
C
      IF(ISRFLG.EQ.1) THEN
C
C ----------------- LEAST DISTANCE PROGRAMMING METHOD ------------------
C
        CALL LSDPQP(GMAT,IROWG,JSTRG,NONZG,MCON,NDIM,
     $    WORK,NWORK,IWORK,NIWORK,NEEDED,CLWR,CUPR,CVEC,
     $    SVEC,XLWR,XUPR,XVEC,IPCSHR,IPU,ISTART,
     $    IT,ISTATC,ISTATV,CNDNUM,NUMKTF,IERLDP)
C
        IER = IERLDP
C
      ELSE    
C
C ----------------- RELAXATION STRATEGY --------------------------------
C
        CALL RELXQP(GMAT,IROWG,JSTRG,NONZG,MCON,NDIM,
     $    WORK,NWORK,IWORK,NIWORK,NEEDED,CLWR,CUPR,CVEC,
     $    SVEC,XLWR,XUPR,XVEC,IPCSHR,IPU,ISTART,IT,
     $    ISTATC,ISTATV,CNDNUM,NUMKTF,IERELX)
C
        IER = IERELX
C
      ENDIF
C
C         SCHUR-QP TIMING STATISTICS
C
      IF(ISQPER(4).GT.INSTAT(14)) THEN
        INSTAT(14) = ISQPER(4)
        INSTAT(15) = -IT
      ENDIF
      NQPITR = NQPITR + ISQPER(2)
      INSTAT(26) = INSTAT(26) + NQPITR
      INSTAT(27) = MAX(INSTAT(27),NQPITR)
      IF(ISQPER(2).GT.INSTAT(16)) THEN
        INSTAT(16) = ISQPER(2)
        INSTAT(17) = -IT
      ENDIF
      IF(ISQPER(5).GT.INSTAT(18)) THEN
        INSTAT(18) = ISQPER(5)
        INSTAT(19) = -IT
      ENDIF
      IF(ISQPER(10).GT.INSTAT(20)) THEN
        INSTAT(20) = ISQPER(10)
        INSTAT(21) = -IT
      ENDIF
C
      IF(IER.NE.0) INSTAT(11) = INSTAT(11) + 1
C
C             CHECK FOR SUCCESSFUL TERMINATION.
C
      IF (IER.EQ.0.OR.IER.EQ.1.OR.IER.EQ.2.OR.
     $         (IER.EQ.-1106.AND.ISRFLG.EQ.2)) THEN
C
C         SUCCESSFUL STEP -- USE WARM START FOR SUBSEQUENT CALLS
C
        ISTART = 0
C
      ELSEIF(IER.EQ.3.AND.QPOPTN.NE.'SPARSE') THEN
C
C         MINIMUM SUM OF INFEASIBILITIES BY DENSE QP
C
        ISTART = 1
C
      ELSEIF(IER.EQ.-1.OR.IER.EQ.-1014.OR.IER.EQ.-1015) THEN
C
C           INSUFFICIENT REAL STORAGE
C
        ISTERM = 5
        GO TO 200
C
      ELSEIF(IER.EQ.-2.OR.IER.EQ.-1019) THEN
C
C           INSUFFICIENT INTEGER STORAGE
C
        ISTERM = 7
        GO TO 200
C
      ELSEIF(IER.EQ.-1111) THEN
C
C           I/O ERROR (INSUFFICIENT DISK SPACE)
C
        ISTERM = 11
        GO TO 200
C
      ELSEIF(IER.EQ.-1999) THEN
C
C           USER EXTERNAL KILL
C
        ISTERM = 12
        GO TO 200
C
      ELSE
C
C           IER .NE. 0 -- QP SUBPROBLEM FAILED
C
        IF (IPCFEZ.GE.10) THEN
          IF(IER.EQ.-1103.OR.IER.EQ.-1104) THEN
            WRITE(IPU,1007)
          ELSEIF(IER.EQ.-1106.OR.IER.EQ.-1108) THEN
            WRITE(IPU,1008)
          ELSEIF(IER.EQ.3) THEN
            WRITE(IPU,1009)
          ELSE
            WRITE(IPU,1006) IER
          ENDIF
        ENDIF
C
        IF(ISRFLG.EQ.1.AND.ISTART.EQ.0) THEN
C
C         SEARCH DIRECTION CALCULATION FAILED -- TRY AGAIN WITH COLD START
C
          ISTART = 1
          IF(IPCFEZ.GE.10) WRITE(IPU,1002)
C
          GO TO 120
C
        ELSEIF(ISRFLG.EQ.1.AND.IRELAX.GT.0) THEN
C
C         SEARCH DIRECTION CALCULATION FAILED -- USE COLD START FOR 
C         SUBSEQUENT CALLS AND SWITCH TO BACKUP STRATEGY
C
          ISTART = 1
          ISRFLG = 2
          IF(IPCFEZ.GE.10) WRITE(IPU,1005)
C
          GO TO 120
C
        ELSE
C
C           RELAXATION PROCEDURE FAILED
C
          IF(IPCFEZ.GE.10.AND.IRELAX.GT.0) WRITE(IPU,1003)
          ISTERM = 6
          GO TO 200
C
        ENDIF
C
      ENDIF
C
C           COMPUTE GRADIENT VECTOR = 2*(GMAT**T)*CHAT
C           WHERE CHAT = C FOR ACTIVE CONSTRAINTS AND 
C           VIOLATED INEQUALITIES.  FORM CHAT IN WORK.
C           AND GRADIENT IN WORK(MCON+1)
C
      DO I = 1,MCON
C
        IF(CLWR(I).EQ.CUPR(I)) THEN
C
C         EQUALITY CONSTRAINT
C
          WORK(I) = CLWR(I)-CVEC(I)
C
        ELSE
C
C         INEQUALITY CONSTRAINT
C
          IF(CVEC(I).LT.CLWR(I)-CONTOL) THEN
            WORK(I) = CLWR(I)-CVEC(I)
          ELSEIF(CVEC(I).GT.CUPR(I)+CONTOL) THEN
            WORK(I) = CUPR(I)-CVEC(I)
          ELSE
            WORK(I) = ZERO
          ENDIF
C
        ENDIF
C
        WORK(I) = TWO*WORK(I)
C
      ENDDO
C
C         COMPUTE PRODUCT OF GMAT**T WITH CHAT.
C
      CALL MVPSPR(11,NDIM,MCON,GMAT,IROWG,JSTRG,WORK,WORK(MCON+1))
C
C         AUGMENT GRADIENT WITH BOUND CONTRIBUTION
C
      DO I = 1,NDIM
C
        IF(XLWR(I).EQ.XUPR(I)) THEN
C
C         EQUALITY BOUND
C
          TERM = XLWR(I)-XVEC(I)
C
        ELSE
C
C         INEQUALITY BOUND
C
          IF(XVEC(I).LT.XLWR(I)-CONTOL) THEN
            TERM = XLWR(I)-XVEC(I)
          ELSEIF(XVEC(I).GT.XUPR(I)+CONTOL) THEN
            TERM = XUPR(I)-XVEC(I)
          ELSE
            TERM = ZERO
          ENDIF
C
        ENDIF
C
        WORK(MCON+I) = -WORK(MCON+I) - TWO*TERM
C
      ENDDO
C
C             ITERATION PRINT
C
      SRPHI = SQRT(PHI*CONTOL)
      SNORM = ZERO
      MACTIV = 0
      DO I=1,NDIM
        IF (ISTATV(I).GT.0)  MACTIV = MACTIV + 1
        SNORM = SNORM + SVEC(I)**2
      ENDDO
      SNORM = SQRT(SNORM)
      DO I=1,MCON
        IF (ISTATC(I).GE.1 .AND. ISTATC(I).LE.3)  MACTIV = MACTIV + 1
      ENDDO
      NDOF = NDIM - MACTIV
C
      IF(IPCFEZ.GE.10) THEN
        CALL HHCPRS(CHITER,' ','-')
        WRITE(IPU,1001) CHITER
          WRITE(LFTVAL,'(SP,1PG14.6)') SRPHI
          WRITE(RITVAL,'(SP,1PG14.6)') CNDNUM
          CALL DBLSTR('Constraint Error',LFTVAL,
     $               CNDLBL,RITVAL,IPU)
          WRITE(LFTVAL,'(I8)') MACTIV
          WRITE(RITVAL,'(I8)') NDOF
          CALL DBLSTR('Active Constraints',LFTVAL,
     $               'Degrees of Freedom',RITVAL,IPU)
          WRITE(LFTVAL,'(I8)') NQPITR
          WRITE(RITVAL,'(I8)') NUMKTF
          CALL DBLSTR('QP Iterations',LFTVAL,
     $               'Matrix Factorizations',RITVAL,IPU)
      ELSEIF(IPCFEZ.GT.0) THEN
        CALL WRTURS(NPRFLG,CHITER,NQPITR,NDOF,ALFA,SNORM,SRPHI,NUMKTF,
     $    XDUM,SRPHI,XDUM,XDUM,CNDNUM,IPU) 
      ENDIF
C
C       ----------------------------------------------------------------
C
C         SAVE ITERATION LOG QUANTITIES
C
      ITP    = IT    
      KEYMP  = 0 
      PPGNOR = ZERO
      PFNOM  = SRPHI  
      PERREQ = ZERO
      PERRIN = ZERO
      PFMRT  = SRPHI  
      PPENMA = ZERO
      PDIAGN = ZERO
      PCNDNU = CNDNUM
      PEIGMI = ZERO
      PEIGMA = ZERO
      MACTIP = MACTIV
      NDOFP  = NDOF  
      NQPITP = NQPITR
      NUMKTP = NUMKTF
C
C       ----------------------------------------------------------------
C
C
      IF(IT.EQ.1) STEPL = SNORM
C
C ----------------------------------------------------------------------
C ----------------- CONVERGENCE TESTS ----------------------------------
C ----------------------------------------------------------------------
C
C           AT LEAST ONE CONSTRAINT IS VIOLATED.
C           TERMINATE IF STEP IS TOO SMALL.
C
      IF (.NOT.CONTST.AND.IT.GT.1) THEN
C
        DO I = 1,MCON
            CHNG = ABS(CVEC(I) - COLD(I)) - MAX(SMALL,CONTOL)
            IF(CHNG.GT.ZERO) GO TO 160
        ENDDO
C
C           ALL CONSTRAINT CHANGES ARE SMALL.
C           CHECK FOR SMALL STEP LENGTHS.
C
        XNRM = ZERO
        DO I=1,NDIM
          XNRM = XNRM + XVEC(I)**2
        ENDDO
        XNRM = SQRT(XNRM)
        OLDSTP = STEPL/MAX(ZEROMN,XNRM)
        STPTST = ABS(ALFA)*SNORM/MAX(ZEROMN,XNRM)
        IF(STPTST.LT.SMALL.AND.OLDSTP.LT.SMALL) THEN
C
C           TWO SUCCESSIVE STEPS ARE SMALL AND AT LEAST ONE
C           CONSTRAINT IS VIOLATED -- TERMINATE
C
          ISTERM = 2
          GO TO 200
C
        ENDIF
C
      ENDIF
C
 160  CONTINUE
C
C           CHECK FOR MAXIMUM ITERATION COUNT.
C
      IF (IT.GT.NITMAX) THEN
        ISTERM = 3
        GO TO 200
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------- ONE DIMENSIONAL SEARCH -----------------------------
C ----------------------------------------------------------------------
C
C           INITIALIZATION.
C
      ITRCN = 1
      ALFA = ONE
C
C         FIRST ITERATION.  COMPUTE DF0.
C
      WDOTS = ZERO
      DO KK = 1,NDIM
        WDOTS = WDOTS + SVEC(KK)*WORK(MCON+KK)
      ENDDO
      DF0 = WDOTS/CONTOL
C
C         CHECK THAT SEARCH DIRECTION IS DOWNHILL--IF NOT USE GRADIENT
C
      IF(DF0.GT.ZERO) THEN
C
        WDOTS = ZERO
        DO KK = 1,NDIM
          SVEC(KK) = -WORK(MCON+KK)
          WDOTS = WDOTS + SVEC(KK)*WORK(MCON+KK)
        ENDDO
        DF0 = WDOTS/CONTOL
        IF(IPCFEZ.GE.10) WRITE(IPU,1010)
C
      ENDIF
C
      PHISR = PHI
C
 170  CONTINUE
C
      ABSSFZ = ABS(SFZTOL)
      CALL LNSRCH(DF0,D2F0,PHIBSR,FREDUC,PHISR,ALFA,CONTOL,ABSSFZ,SNORM,
     $    ITRCN,IT1MAX,GOODGR,IFUNER,IREVRS(1),IOFLIN,IPU,IERLIN)
C
      IF(IREVRS(1).EQ.0) GO TO 190
C
C ----------------------------------------------------------------------
C
C         FUNCTION EVALUATION.
C
C         COMPUTE NEW ESTIMATE OF VARIABLES.
C
      DO I = 1,NDIM
        XBAR(I) = XVEC(I) + ALFA*SVEC(I)
      ENDDO
      IRTRN = 2
      IREVRS(4) = 0
      IF(IPCFEZ.GE.20) IREVRS(4) = 1
C
C           RETURN FOR FUNCTION EVALUATION.
C
      RETURN
 502  CONTINUE
C
      IREVRS(4) = 0
C
C           THE FUNCTION EVALUATION IS COMPLETED.
C
      IF(IFERR.EQ.0) THEN
C
        IF(ISRFLG.EQ.1.AND.IT.EQ.1.AND.ITRCN.EQ.1.AND.LNRSTP.GE.0
     $    .AND.ALFA.EQ.ONE.AND.IERLDP.EQ.0) THEN
C
C         CHECK FOR LINEAR CONSTRAINT BEHAVIOR ON FIRST STEP 
C         OF FIRST LINE SEARCH
C
C         COMPUTE SLACK VARIABLE AND SAVE TEMPORARILY IN COLD
C
          CALL MVPSPR(1,MCON,NDIM,GMAT,IROWG,JSTRG,SVEC,COLD)
C
          DO II = 1,MCON
            COLD(II) = COLD(II) + CVEC(II)
            IF(CLWR(II)-CONTOL.LE.CBAR(II).AND.
     $         CBAR(II).LE.CUPR(II)+CONTOL.AND.
     $         ABS(CBAR(II)-COLD(II)).LT.CONTOL) LNRSTP = 1
          ENDDO
C
        ENDIF
C
C         COMPUTE CONSTRAINT ERROR AT POINT
C
        CALL FMERIT(CBAR,CLWR,CUPR,XBAR,XLWR,XUPR,MCON,NDIM,CONTOL,
     $    PHIBAR,IFERR)
C
        PHIBSR = PHIBAR
C
      ENDIF
C
      IFUNER = IFERR
      IF(IFERR.EQ.0) THEN
C
C           CHECK FOR CONSTRAINT SATISFACTION.
C
        CALL CONSAT(CLWR,CUPR,CBAR,XLWR,XUPR,XBAR,MCON,NDIM,
     $    CONTOL,CONTST)
C
        IF(CONTST) IFUNER = -1
C
      ENDIF
C
      IF(LNRSTP.EQ.1) THEN
C
C         LINEAR BEHAVIOR DETECTED, SO ACCEPT POINT REGARDLESS
C         OF CONSTRAINT ERROR AND CONTINUE SEARCH
C
        FERROR = SQRT(CONTOL*PHIBAR)
        ERRCHG = FERROR - SQRT(CONTOL*PHI)
        IF(IOFLAG.GE.10) WRITE(IPU,1012) FERROR,ERRCHG
C
        IF(IPCFEZ.GE.10) WRITE(IPU,1011)
        LNRSTP = 0
        IREVRS(1) = 0
        GO TO 190
C
      ENDIF
C
C ----------------------------------------------------------------------
C
C         FUNCTION EVALUATED. REENTER ONE DIMENSIONAL SEARCH
C
      GO TO 170
C
 190  CONTINUE
C
      IF(IFERR.EQ.-100) THEN
C
C             MAX. NO. OF. FUNCTION EVALS.
C
        ISTERM = 8
        GO TO 200
C
      ELSEIF(IERLIN.EQ.-1) THEN
C
C             UPHILL SEARCH DIRECTION
C
        ISTERM = 9
        GO TO 200
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------- UPDATE INFORMATION ---------------------------------
C ----------------------------------------------------------------------
C
C             UPDATE POINT, XVEC.
C
      DO KK = 1,NDIM
        XVEC(KK) = XBAR(KK)
      ENDDO
C
C         UPDATE CONSTRAINTS--STORE NEW VALUE IN CVEC, OLD IN COLD
C
      DO KK = 1,MCON
        COLD(KK) = CVEC(KK)
        CVEC(KK) = CBAR(KK)
      ENDDO
C
C             UPDATE CONSTRAINT ERROR, PHI, STEP LENGTH, STEPL, AND
C             ITERATION COUNTER, IT.
C
      PHI = PHIBAR
C
      IF(ISRFLG.EQ.2.AND.ALFA.EQ.ONE) THEN
C
        IF(IPCFEZ.GE.10) WRITE(IPU,1004)
        ISTART = 1
        ISRFLG = 1
C
      ENDIF
C
C
      STEPL = ABS(ALFA)*SNORM
      IT = IT + 1
C
C ----------------------------------------------------------------------
C ----------------- END OF ITERATION -----------------------------------
C ----------------------------------------------------------------------
C
      GO TO 110
C
C ----------------------------------------------------------------------
C ----------------- ALGORITHM TERMINATION PROCESSING -------------------
C ----------------------------------------------------------------------
 200  CONTINUE
C
C         GET CONSTRAINT VALUES BACK IN SYNC AT FINAL POINT
C
      DO KK = 1,MCON
        CBAR(KK) = CVEC(KK)
      ENDDO
C
      RETURN
C
 1001 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'-------------------------------
     $------ Iteration ',A6,'----------------------------------',
     $   T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'COLD START QP ALGORITHM',
     $    T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'RELAXATION STRATEGY FAILED',
     $    T106,'*')
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'SWITCH TO LEAST DISTANCE PROGRA
     $M',T106,'*')
 1005 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'SWITCH TO RELAXATION STRATEGY',
     $    T106,'*')
 1006 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'QP ALGORITHM FAILED---IER =',
     $    I5,T106,'*')
 1007 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'+++++ SINGULAR JACOBIAN MATRIX'
     $  ,T106,'*')
 1008 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'+++++ ILL-CONDITIONED JACOBIAN 
     $MATRIX',T106,'*')
 1009 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'+++++ LOCALLY INFEASIBLE LINEAR
     $ CONSTRAINTS',T106,'*')
 1010 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'+++++ USE GRADIENT DIRECTION',
     $    T106,'*')
 1011 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'LINEAR CONSTRAINT SATISFACTION 
     $STEP',T106,'*')
 1012 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'...............Error =',
     $    1PG16.8,5X,'Change =',1PG16.8,T106,'*')
C
      END
