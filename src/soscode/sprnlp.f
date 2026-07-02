


      SUBROUTINE SPRNLP( IRVCOM ,IREVRS ,XBAR   ,XLWR   ,XUPR   
     $          ,ISTATV ,VECNU  ,NDIM   ,FBAR   ,RESVEC ,MAXRES 
     $          ,NRES   ,DELF   ,RMAT   ,IROWR  ,JCOLR  ,NONZR  
     $          ,HMAT   ,IROWH  ,JSTRH  ,NONZH  ,CBAR   ,CLWR     
     $          ,CUPR   ,ISTATC ,MAXCON ,MCON   ,VECLAM ,GMAT      
     $          ,IROWG  ,JCOLG  ,NONZG  ,IFERR  ,NFEVAL ,HOLD   
     $          ,NHOLD  ,IHOLD  ,NIHOLD ,NEEDED ,IERNLP )
C
C ======================================================================
C     SPRNLP===>sprnlp   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C *** PURPOSE...NLPSPR ALGORITHM INTERFACE
C
C        SPRNLP CONTROLS EXECUTION OF THE "NON-LINEAR PROGRAMMING,
C        SPARSE" (NLPSPR) ALGORITHM.  IT IS THE MAIN
C        USER INTERFACE FOR INPUT OF ALGORITHM QUANTITIES TO NLPSPR. 
C        A REVERSE COMMUNICATION STRUCTURE IS USED TO OBTAIN THE
C        NECESSARY FUNCTION EVALUATIONS.  
C        
C *** CALLING ARGUMENTS
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C ======================================================================
C *** ARGUMENTS TO CONTROL REVERSE COMMUNICATION 
C ======================================================================
C
C     IRVCOM     I/O     CONTINUATION CONTROL FLAG.  THIS ARGUMENT
C                        MUST BE INITIALIZED.  THEREAFTER IT IS SET
C                        BY THE ALGORITHM AND MUST NOT BE CHANGED
C                        BY THE CALLING PROGRAM.
C
C                        = -1      INITIALIZATION PASS
C                        = 0       TERMINATION (NORMAL OR ABNORMAL)
C                        = +1      PERFORM AN "EVALUATION" AND THEN 
C                        REENTER SPRNLP.  OPERATIONS TO BE PERFORMED
C                        ARE DEFINED BY THE OUTPUT VARIABLES
C                        IREVRS(1), IREVRS(2), IREVRS(3), AND IREVRS(4) 
C                        BELOW.  IREVRS(5) CONTAINS INFORMATION WHICH
C                        MAY BE VALUABLE TO THE INTERFACE BUT DOES
C                        NOT NECESSARILY REQUIRE ANY ACTION ON THE
C                        PART OF THE USER
C     IREVRS(1)   O      FUNCTION EVALUATION REQUEST (EVALUATE
C                        OBJECTIVE FUNCTION AND ALL CONSTRAINTS
C                        AT CURRENT POINT---XBAR)
C                        = 0  FUNCTION EVALUATION NOT REQUESTED;
C                             <<< WARNING: DO NOT CHANGE THE CURRENT 
C                             FUNCTION AND CONSTRAINT VALUES >>>
C                        = 1  FUNCTION EVALUATION REQUESTED
C     IREVRS(2)   O      GRADIENT EVALUATION REQUEST
C                        = 0  GRADIENT EVALUATION NOT REQUESTED
C                             <<< WARNING: DO NOT CHANGE THE CURRENT 
C                             GRADIENT AND JACOBIAN >>>
C                        = 1  APPROXIMATE (E.G. FORWARD DIFFERENCE) 
C                             GRADIENT REQUESTED
C                        = 2  ACCURATE (E.G.CENTRAL DIFFERENCE) 
C                             GRADIENT REQUESTED
C     IREVRS(3)   O      HESSIAN EVALUATION REQUEST
C                        = 0  HESSIAN EVALUATION NOT REQUESTED
C                             <<< WARNING: DO NOT CHANGE THE CURRENT 
C                             HESSIAN >>>
C                        = 1  HESSIAN DIAGONAL REQUESTED (AND GET 
C                             READY FOR FULL HESSIAN)
C                        = 2  FULL HESSIAN EVALUATION REQUESTED
C     IREVRS(4)   O      SYSTEM PRINT (OUTPUT) REQUEST
C                        = 0  SYSTEM PRINT NOT REQUESTED
C                        = 1  SYSTEM PRINT REQUESTED
C     IREVRS(5)   O      ALGORITHM INFORMATION FLAG
C                        = 0  CALL WAS MADE FROM THE INTERFACE PORTION 
C                             OF THE ALGORITHM
C                        = 1  CALL WAS MADE FROM THE OPTIMIZATION PORTION
C                             OF THE ALGORITHM
C                        = 2  CALL WAS MADE FROM THE OPTIMIZATION PORTION 
C                             OF THE ALGORITHM REQUESTING HESSIAN RESET
C                        = 3  CALL WAS MADE FROM THE OPTIMIZATION PORTION 
C                             OF THE ALGORITHM REQUESTING OPTIMAL 
C                             PERTURBATION SIZE ADJUSTMENT
C                        = 4  CALL WAS MADE FROM THE FEASIBLITY PORTION 
C                             OF THE ALGORITHM
C
C =======================================================================
C *** INDEPENDENT VARIABLE DATA 
C =======================================================================
C
C     XBAR       I/O     CURRENT VARIABLE VALUES (NDIM)
C                        INITIAL GUESS MUST BE INPUT WHEN IRVCOM=-1.
C                        FINAL POINT IS OUTPUT WHEN IRVCOM = 0.
C                        INTERMEDIATE ITERATION POINTS ARE
C                        OUTPUT WHEN IRVCOM = 1.
C     XLWR       I       VARIABLE LOWER BOUNDS (NDIM)
C     XUPR       I       VARIABLE UPPER BOUNDS (NDIM)
C     ISTATV     I/O     INTEGER VARIABLE STATUS (NDIM)
C                        = 0  --- FREE VARIABLE 
C                        = 1  --- FIXED ON LOWER BOUND
C                        = 2  --- FIXED ON UPPER BOUND
C                        = 3  --- FIXED PERMANENTLY 
C     VECNU      I/O     BOUND MULTIPLIERS (NDIM)
C     NDIM       I       NUMBER OF VARIABLES
C
C =======================================================================
C *** OBJECTIVE FUNCTION DATA 
C =======================================================================
C
C     FBAR       I/O     OBJECTIVE FUNCTION EVALUATED AT XBAR
C                        FBAR IS INPUT WHEN NRES=0, OUTPUT OTHERWISE
C     RESVEC     I       RESIDUAL VECTOR FOR LEAST SQUARES OBJECTIVE 
C                        EVALUATED AT XBAR (NRES)
C     MAXRES     I       RESIDUAL ARRAY DIMENSION
C     NRES       I       NUMBER OF RESIDUALS
C
C           ----- GRADIENT VECTOR 
C
C     DELF       I/O     GRADIENT OF FBAR AT XBAR (NDIM)
C                        DELF IS INPUT WHEN NRES=0, OUTPUT OTHERWISE
C
C           ----- RESIDUAL JACOBIAN MATRIX (NRES .GT. 0)
C
C     RMAT       I       RESIDUAL DERIVATIVES AT XBAR (NONZR)
C     IROWR      I       ROW INDICES OF JACOBIAN NONZEROS (NONZR)
C     JCOLR      I       COLUMN INDICES OF JACOBIAN NONZEROS (NONZR)
C     NONZR      I       NUMBER OF JACOBIAN NONZEROS
C
C           ----- HESSIAN MATRIX
C
C     HMAT       I       NONZERO ELEMENTS OF HESSIAN MATRIX OF THE
C                        LAGRANGIAN FUNCTION (NONZH)
C                        SINCE THE HESSIAN IS SYMMETRIC, ONLY THE
C                        LOWER TRIANGULAR PART IS INPUT
C     IROWH      I       ROW INDICES OF HESSIAN NONZEROS (NONZH)
C     JSTRH      I       COLUMN START INDICES OF NONZEROS (NDIM+1)
C     NONZH      I       NUMBER OF NONZERO HESSIAN ELEMENTS
C
C     ==================================================================
C
C        NOTE >> THE HESSIAN OF THE LAGRANGIAN (L) IS A FUNCTION OF VECLAM
C
C        L(XBAR,VECLAM) = FBAR - VECLAM(1)*CBAR(1) - ... 
C                              - VECLAM(MCON)*CBAR(MCON)
C
C
C     EXAMPLE:           ( NDIM = 3, NONZH = 4 = [JSTRH(NDIM+1)-1] )
C
C               | 1.1       (3.1)|      HMAT(*) = (1.1, 3.1, 2.2, 3.3)
C     HESSIAN = |       2.2      |      IROWH(*) = (1, 3, 2, 3)
C               | 3.1        3.3 |      JSTRH(*) = (1, 3, 4, 5)
C
C     ==================================================================
C
C
C =======================================================================
C *** CONSTRAINT DATA 
C =======================================================================
C
C     CBAR       I       CONSTRAINTS EVALUATED AT XBAR (MCON)
C     CLWR       I       CONSTRAINT LOWER BOUNDS (MCON)
C     CUPR       I       CONSTRAINT UPPER BOUNDS (MCON)
C     ISTATC     I/O     INTEGER CONSTRAINT STATUS (MCON)
C                        = 0  --- FREE (INACTIVE) INEQUALITY
C                        = 1  --- FIXED ON LOWER BOUND
C                        = 2  --- FIXED ON UPPER BOUND
C                        = 3  --- EQUALITY
C                        = 4  --- IGNORED CONSTRAINT
C     MAXCON     I       CONSTRAINT ARRAY DIMENSION
C     MCON       I       NUMBER OF CONSTRAINTS
C     VECLAM     I/O     LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MCON)
C
C           ----- JACOBIAN MATRIX
C
C     GMAT       I       CONSTRAINT DERIVATIVES AT XBAR (NONZG)
C     IROWG      I       ROW INDICES OF JACOBIAN NONZEROS (NONZG)
C     JCOLG      I       COLUMN INDICES OF JACOBIAN NONZEROS (NONZG)
C     NONZG      I       NUMBER OF JACOBIAN NONZEROS
C
C     ==================================================================
C
C     EXAMPLE:           ( NDIM = 5, MCON = 3, NONZG = 4)
C
C                |        3.1    5.1|     GMAT(*) = (2.1, 3.1, 3.3, 5.1)
C     JACOBIAN = | 2.1              |     IROWG(*) = (2, 1, 3, 1)
C                |        3.3       |     JCOLG(*) = (1, 3, 3, 5)
C
C     ==================================================================
C
C ======================================================================
C *** ALGORITHM CONTROL DATA 
C ======================================================================
C
C     IFERR      I       FUNCTION EVALUATION ERROR FLAG
C                        = 1        WHEN FBAR CANNOT BE EVALUATED
C                        = 0        OTHERWISE
C     NFEVAL     I       NUMBER OF FUNCTION EVALUATIONS
C
C ======================================================================
C *** WORKING STORAGE 
C ======================================================================
C
C     HOLD       I       HOLD ARRAY (WORKING STORAGE NOT DESTROYED)
C     NHOLD      I       DIMENSION OF HOLD ARRAY; (WHEN IERNLP = -127,
C                        -131 NEEDED CONTAINS THE REQUIRED STORAGE)
C     IHOLD      I       INTEGER HOLD ARRAY (WORKING STORAGE NOT 
C                        DESTROYED)
C     NIHOLD     I       DIMENSION OF IHOLD ARRAY; (WHEN IERNLP = -128,
C                        -132 NEEDED CONTAINS THE STORAGE NEEDED)
C     NEEDED     O       STORAGE REQUIRED WHEN IERNLP = -127, -128,
C                        -131,-132.
C
C ======================================================================
C
C *** OPTIONAL INPUTS SET BY CALLS TO INSNLP.  INPUT IS NOT NEEDED
C                        UNLESS DEFAULT VALUES ARE INAPPRORIATE.
C                        ALL OPTIONAL INPUTS ARE CONTAINED IN ONE
C                        OF THREE COMMONS: NPSPRR; NPSPRI; NPSPRC
C
C ======================================================================
C ======================================================================
C *** CONVERGENCE TOLERANCES 
C ======================================================================
C
C     CONTOL     I       CONSTRAINT TOLERANCE
C     OBJTOL     I       OBJECTIVE FUNCTION TOLERANCE
C     PGDTOL     I       PROJECTED GRADIENT TOLERANCE
C     SLPTOL     I       LINE SEARCH SLOPE TOLERANCE
C
C     ==================================================================
C
C     CONVERGENCE REQUIRES:
C
C        (A)   MAX ABSOLUTE ERROR IN ACTIVE CONSTRAINTS AND BOUNDS 
C              .LT. CONTOL              
C        (B)   MAX ABSOLUTE ERROR IN KT CONDITIONS .LT. PGDTOL*|DELF|
C        (C)   |FBAR - FMIN| .LT. OBJTOL
C        (D)   |STEPLENGTH| .LT. 
C                        SQRT( OBJTOL/(ONE + |FBAR|) )*(ONE + XNORM)
C        (E)   CORRECT SIGN FOR ALL LAGRANGE MULTIPLIERS
C
C     ==================================================================
C
C ======================================================================
C *** ALGORITHM CONTROL DATA 
C ======================================================================
C
C     MAXNFE     I       MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C
C           *****     ALGORITHM CONTROL PARAMETERS     *****
C
C     ALGOPT     I       CHARACTER VARIABLE SPECIFYING ALGORITHM 
C                        OPTION
C                        =  FM     FIND FEASIBLE POINT THEN MINIMIZE
C                        =  FME    FIND FEASIBLE POINT THEN MINIMIZE 
C                                  WITH EQUALITIES BINDING
C                        =  M      MINIMIZE FROM THE INITIAL POINT
C                        =  F      FIND FEASIBLE POINT ONLY
C     KTOPTN     I       CHARACTER VARIABLE SPECIFYING KT FACTORIZATION
C                        =  SMALL  ACTIVE CONSTRAINTS ONLY
C                        =  LARGE  ALL EQUALITY AND INEQUALITIES
C     QPOPTN     I       CHARACTER VARIABLE SPECIFYING QP ALGORITM
C                        =  SPARSE SCHUR-COMPLEMENT QP
C                        =  DENSE  DENSE NULLSPACE QR
C
C     ==================================================================
C
C     IN GENERAL ALL INPUTS AND OUTPUTS DESCRIBED ASSUME A NORMAL (COLD) 
C     START.  
C
C     ==================================================================
C
C
C     IT1MAX     I       MAXIMUM NUMBER OF STEPS IN LINE SEARCH
C     IOFLAG     I       OUTPUT CONTROL FLAG
C     IPUNLP     I       OUTPUT DEVICE NUMBER
C     NITMAX     I       MAXIMUM NUMBER OF ITERATIONS, I.E. THE NUMBER
C                        OF TIMES THE QP IS CALLED FOR MINIMIZATION 
C                        AND FEASIBILITY SEARCH
C
C       OPTIONAL, DIAGNOSTIC, AND OUTPUT CONTROL FLAGS
C                                                
C                                                                 DEFAULT
C     IOFLIN     I       LINE SEARCH DIAGNOSTIC OUTPUT FLAG            -1
C     IOFMFR     I       MULTIFRONTAL OUTPUT FLAG                       0  
C     IOFSHR     I       SCHUR-QP OUTPUT FLAG                           0  
C     IOFSRC     I       SRCHFZ OUTPUT FLAG                             0
C     IPUDRF     I       I/O UNIT FOR SCHUR-QP DUMP FROM SRCHDR         0
C     ITDRQP     I       DUMP SCHUR-QP ON SRCHDR ITERATION NO.         -1
C     IPUFZF     I       I/O UNIT FOR SCHUR-QP DUMP FROM SRCHFZ         0
C     ITFZQP     I       DUMP SCHUR-QP ON SRCHFZ ITERATION NO.         -1
C     IPUSTF     I       I/O UNIT FOR SCHUR-QP DUMP FROM QPSTRT         0
C     IPUMF1     I       MULTIFRONTAL SQFILE I/O UNIT IN SRCHDR        11
C     IPUMF2     I       MULTIFRONTAL WAFIL1 I/O UNIT IN SRCHDR        12
C     IPUMF3     I       MULTIFRONTAL WAFIL2 I/O UNIT IN SRCHDR        13
C     IPUMF4     I       MULTIFRONTAL SQFILE I/O UNIT IN L.D.P.        14
C     IPUMF5     I       MULTIFRONTAL WAFIL1 I/O UNIT IN L.D.P.        15
C     IPUMF6     I       MULTIFRONTAL WAFIL2 I/O UNIT IN L.D.P.        16  
C     SFZTOL     I       LINE SEARCH SLOPE TOLERANCE IN SRCHFZ        .01
C
C ======================================================================
C *** OUTPUT ARGUMENTS SET BY ALGORITHM BEFORE TERMINATION (IRVCOM = 0)
C ======================================================================
C
C     IERNLP     O       SUCCESS/ERROR CODE
C
C          --------------
C
C          0     NORMAL TERMINATION
C
C          --------------
C                FATAL (USUALLY INPUT) ERRORS 
C          --------------
C
C          -101  (MCON.LT.0)
C          -102  (NDIM.LT.1)
C          -103  (MCON.GT.MAXCON) 
C          -104  (|IHESHN|.GT.3) 
C          -105  (NITMAX.LT.MAX(NITMIN,1))
C          -106  (SLPTOL.GE.1) .OR. (SLPTOL.LE.ONEEM5)
C          -107  (IT1MAX.LT.1)
C          -108  (IOFLAG.LT.0) .OR.(IOFLAG.GT.30)
C          -109  INVALID INPUT FOR ALGOPT
C          -110  (OBJTOL.LE.ONEEP1*HDMCON(5))
C          -111  (PGDTOL.LT.SQRT(HDMCON(5))).OR.(PGDTOL.GT.ONEEM2)
C          -112  CONTOL.LT.SQRT(HDMCON(5))
C          -113  (ISTATC.LT.0).OR.(ISTATC.GT.4)
C          -114  (NONZG.LE.0) .OR. (NONZG .GT. NDIM*MCON)
C          -115  (JCOLG(*).LE.0) .OR. (JCOLG(*).GT.NDIM)
C          -116  (IROWG(*).LE.0) .OR. (IROWG(*).GT.MCON)
C          -117  (CUPR.LT.CLWR)
C          -118  (CUPR.EQ.CLWR).AND.(ISTATC.NE.3) .OR.
C                (CUPR.NE.CLWR).AND.(ISTATC.EQ.3)
C          -119  (CUPR.NE.CLWR).AND.((|CUPR-CLWR|.LT.CONTOL)
C          -120  (ISTATV.LT.0).OR.(ISTATV.GT.3)
C          -121  (JSTRH(1).NE.1).OR.(JSTRH(I).GE.JSTRH(I+1))
C          -122  (NZHDIM.LE.0) .OR. (NZHDIM.GT.NDNSH)
C                .OR. (NZHDIM.NE.NONZH)      
C                WHERE NZHDIM = JSTRH(NDIM+1)-1, 
C                AND NDNSH = NDIM*(NDIM+1)/2
C          -123  (IROWH(*) IS INCORRECT)
C          -124  (XUPR.LT.XLWR)
C          -125  (XUPR.EQ.XLWR).AND.(ISTATV.NE.3) .OR.
C                (XUPR.NE.XLWR).AND.(ISTATV.EQ.3)
C          -126  (XUPR.NE.XLWR).AND.((|XUPR-XLWR|.LT.CONTOL)
C          -127  REAL HOLD ARRAY TOO SMALL
C          -128  INTEGER HOLD ARRAY TOO SMALL
C          -129  FUNCTION ERROR AT INITIAL POINT OR DURING GRADIENT EVAL.
C          -130  (MACTIV.GT.NDIM)
C          -131  INSUFFICIENT REAL STORAGE DETECTED 
C          -132  INSUFFICIENT INTEGER STORAGE DETECTED 
C          -133  SINGULAR JACOBIAN ON SUCCESSIVE ITERATIONS
C          -134  INSNLP INPUT ERROR
C          -135  (IFERR.LT.0) .OR. (IFERR.GT.1)
C          -136  (|SFZTOL|.GE.1) .OR. (|SFZTOL|.LE.ONEEM5)
C          -137  CONFLICT BETWEEN USER AND MULTIFRONTAL FILE NUMBER
C          -138  NRES = 0, AND NEWTON .NE. 0, 1
C          -139  (NRES.LT.0)
C          -140  (NRES.GT.MAXRES) 
C          -141  (NONZR.LE.0) .OR. (NONZR .GT. NDIM*NRES)
C          -142  (JCOLR(*).LE.0) .OR. (JCOLR(*).GT.NDIM)
C          -143  (IROWR(*).LE.0) .OR. (IROWR(*).GT.NRES)
C          -144  A ROW OF THE RESIDUAL JACOBIAN IS IDENTICALLY ZERO
C          -145  IROWH(JSTRH(I)).NE.I FOR SOME I
C          -146  ALGOPT.EQ.LLSQ .AND. NRES.LE.0
C          -147  UNEXPECTED ERROR
C          -148  (|IRVCOM|.GT.1) 
C          -149  INCORRECT VALUE FOR QPOPTN
C          -150  QPOPTN = DENSE AND ALGOPT = FME OR NRES .GT. 0
C          -151  CONSTRAINT VIOLATES ITS BOUNDS AND CANNOT BE CHANGED
C          -152  INCORRECT VALUE FOR KTOPTN
C          -153  I/O ERROR (INSUFFICIENT DISK SPACE)
C          -154  TOLKTC < 1
C          -155  TOLPVT < 0 OR TOLPVT > .5
C          -156  TOLFIL < 0 
C          -157  0 < BIGCON < CONTOL
C          -158  ITDRQP > 0 and IPUDRF.LE.0 or ITFZQP > 0 and IPUFZF.LE.0
C          -159  NRES > 0, AND NEWTON .NE. 0, 1, 2
C          -800  USER EXTERNAL KILL
C
C          --------------
C                WARNING ERRORS---EVALUATION BY USER SUGGESTED
C          --------------
C
C          +101  WEAK SOLUTION FOUND (MULTIPLIERS NEAR ZERO)
C          +102  (MEQUAL.EQ.NDIM).AND.ALGOPT.NE.F___
C          +103  MAX. NO. OF CONSECUTIVE FUNCTION ERRORS
C          +104  MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C          +105  SMALL STEP TERM FROM EQCMIN; SUBOPTIMAL FEASIBLE
C                POINT FOUND
C          +106  MAX. NO. OF ITER. IN EQCMIN
C          +107  MAX. NO. OF ITER. IN SRCHFZ
C          +108  FEASIBLE POINT NOT FOUND
C          +109  MAX. NO. OF INTERVAL HALVES IN LINE SEARCH
C          +110  EITHER MAX HESSIAN DIAGONAL OR VIOLATED SLOPE CONDITION;
C                SUBOPTIMAL FEASIBLE POINT FOUND  
C          +111  PROJ. GRAD. CALCULATION FAILED
C          +112  QPSTRT FAILED TO COMPUTE MULTIPLIERS
C          +113  SUBOPTIMAL FEASIBLE POINT FOUND
C          +114  SCHUR-QP FAILED WITH UNEXPECTED ERROR
C          +115  CONTOL.GT.OBJTOL
C          +116  UPHILL DIRECTION DETECTED IN LINE SEARCH
C          +117  REDUCED OBJECTIVE FUNCTION IS LINEAR
C          +118  CONSTRAINTS IGNORED
C          +119  TERMINATE AFTER DIAGNOSTIC LINE SEARCH
C          +120  RECURSIVE HESSIAN ESTIMATE AND NEWTON.NE.0
C          +122  TERMINATE AFTER USER ABORT (IFERR = -1000)
C          +123  CONTOL .LE. BIGCON .LE. 1.D-4
C
C
C     ==================================================================
C
C         NOTE: ALGORITHM PERFORMANCE STATISTICS ARE CURRENTLY ACCUMULATED
C               USING THE CLOCK ROUTINES CLKXXX, AND NUMBERS 1,3,4,5,AND 6.
C               THE VALUES ACCUMULATED ARE DEFINED IN THE COMMENTS BELOW.
C
C     ==================================================================
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C *** DESCRIPTION...
C
C        SPRNLP IS THE ALGORITHM INTERFACE FOR NLPSPR.  SPRNLP
C        PERFORMS THOSE TASKS REQUIRED FOR A STAND ALONE
C        OPTIMIZATION PACKAGE.
C        THESE FUNCTIONS ARE BROKEN INTO THE
C        FOLLOWING SEQUENCES.
C
C        STORAGE POINTER SEQUENCE 
C
C                  INPUT ERROR CHECKING IS PERFORMED
C                  STORAGE POINTERS USED TO DIVIDE THE WORK
C                  ARRAY ARE SET.  
C
C        FIRST ENTRANCE SEQUENCE 
C
C                  THE FUNCTION AND GRADIENT FLAGS ARE SET.
C                  FIRST ENTRY AND SYSTEM PRINT FLAGS ARE SET
C                  WE BRANCH TO THE RETURN
C
C        ERROR CHECKING SEQUENCE
C
C                  DIAGNOSTIC
C                  OUTPUT IS PERFORMED IF NECESSARY.  MAXIMUM NUMBER
C                  OF FUNCTION ERRORS AND EVALUATIONS IS CHECKED.
C
C        INITIAL SUMMARY PRINT (SUMRY = .TRUE.)
C
C                  THE INITIAL SUMMARY OF NLPSPR DATA IS OUTPUT.
C
C        ALGORITHM CALLING SEQUENCE
C
C                  CALL THE OPTIMIZATION ALGORITHM, NLPSPR.
C
C        FUNCTION EVALUATION SEQUENCE 
C
C                  BRANCH TO THE RETURN FOR A FUNCTION EVALUATION.
C
C        TERMINATION SEQUENCE 
C
C                  RETURN
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
c
C ======================================================================
      double precision, allocatable, dimension(:) :: cvec 
      double precision, allocatable, dimension(:) :: cold 
      double precision, allocatable, dimension(:) :: pgrd 
      double precision, allocatable, dimension(:) :: svec 
      double precision, allocatable, dimension(:) :: xold 
      double precision, allocatable, dimension(:) :: yvec 
      double precision, allocatable, dimension(:) :: ybar 
c
      integer, allocatable, dimension(:) :: irwg 
      integer, allocatable, dimension(:) :: jstr 
      integer, allocatable, dimension(:) :: iprmg 
      integer, allocatable, dimension(:) :: iprmc 
      integer, allocatable, dimension(:) :: isto 
      integer, allocatable, dimension(:) :: irwr 
      integer, allocatable, dimension(:) :: jsrr 
      integer, allocatable, dimension(:) :: iprmr 
      integer, allocatable, dimension(:) :: iscr 
C ======================================================================
c
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,ONEEP1=1.0D1,
     $      ONEEM1=1.0D-1,ONEEM2=1.0D-2,ONEEM4=1.0D-4,ONEEM5=1.0D-5,
     $      POINT5=5.0D-1)
C
      DIMENSION     HOLD(NHOLD),CBAR(MAXCON),GMAT(NONZG),DELF(NDIM),
     &      VECLAM(MAXCON),VECNU(NDIM),XBAR(NDIM),IROWG(NONZG),  
     &      JCOLG(NONZG),RESVEC(MAXRES) ,RMAT(NONZR), IROWR(NONZR), 
     &      JCOLR(NONZR),HMAT(NONZH), IROWH(NONZH) ,JSTRH(NDIM+1),
     &      ISTATV(NDIM) ,XLWR(NDIM) ,XUPR(NDIM),ISTATC(MAXCON) ,
     &      CLWR(MAXCON) ,CUPR(MAXCON)
      DIMENSION IREVRS(5)
      PARAMETER (LENBIG=5)
      DIMENSION BIGELM(LENBIG),IROWBG(LENBIG),JCOLBG(LENBIG)
      DIMENSION IHOLD(NIHOLD)
      CHARACTER(LEN=20) ITITLE
      CHARACTER(LEN=60) CODEV, LIBRV
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
      COMMON /DUALGS/ INDUAL
C
      LOGICAL  OPUN,PRMUTE,SUMRY,SPRSCK
C
** error messages:
      CHARACTER(LEN=8)  SUBNAM
      PARAMETER  (NERMSG=83)
      DIMENSION  IERMSG(NERMSG)
      CHARACTER(LEN=100)  FMTMSG(NERMSG),ERRMSG(1)
C
      DATA  (IERMSG(I),I=1,NERMSG) /
     & -101, -102, -103, -104, -105, -106, -107, -108, -109, -110,
     & -111, -112, -113, -114, -115, -116, -117, -118, -119, -120,
     & -121, -122, -123, -124, -125, -126, -127, -128, -129, -130,
     & -131, -132, -133, -134, -135, -136, -137, -138, -139, -140,
     & -141, -142, -143, -144, -145, -146, -147, -148, -149, -150,
     & -151, -152, -153, -154, -155, -156, -157, +101, +102, +103,
     & +104, +105, +106, +107, +108, +109, +110, +111, +112, +113,  
     & +114, +115, +116, +117, +118, +119, +120, +121, +122, +123,
     & -158, -800, -159 /
C
      DATA  (FMTMSG(I),I=1,10)/
     &'(T1,"(MCON.LT.0)")',
     &'(T1,"(NDIM.LT.1)")',
     &'(T1,"(MCON.GT.MAXCON)")', 
     &'(T1,"(|IHESHN|.GT.3)")', 
     &'(T1,"(NITMAX.LT.MAX(NITMIN,1))")',
     &'(T1,"(SLPTOL.GE.1) .OR. (SLPTOL.LE.ONEEM5)")',
     &'(T1,"(IT1MAX.LT.1)")',
     &'(T1,"(IOFLAG.LT.0) .OR.(IOFLAG.GT.30)")',
     &'(T1,"INVALID INPUT FOR ALGOPT")',
     &'(T1,"(OBJTOL.LE.ONEEP1*HDMCON(5))")'/
C
      DATA  (FMTMSG(I),I=11,20)/
     &'(T1,"(PGDTOL.LT.SQRT(HDMCON(5))).OR.(PGDTOL.GT.ONEEM2)")',
     &'(T1,"CONTOL.LT.SQRT(HDMCON(5))")',
     &'(T1,"(ISTATC.LT.0).OR.(ISTATC.GT.4)")',
     &'(T1,"(NONZG.LE.0) .OR. (NONZG .GT. NDIM*MCON)")',
     &'(T1,"(JCOLG(*).LE.0) .OR. (JCOLG(*).GT.NDIM)")',
     &'(T1,"(IROWG(*).LE.0) .OR. (IROWG(*).GT.MCON)")',
     &'(T1,"(CUPR.LT.CLWR) Constraint:",I6)',
     &'(T1,"((CUPR.EQ.CLWR).AND.(ISTATC.NE.3)).OR.((CUPR.NE.CLWR).AND.(I
     &STATC.EQ.3)) Constraint:",I6)',
     &'(T1,"(CUPR.NE.CLWR).AND.(|CUPR-CLWR|.LT.CONTOL) Constraint:",I6)'
     &,
     &'(T1,"(ISTATV.LT.0).OR.(ISTATV.GT.3)")'/
C
      DATA  (FMTMSG(I),I=21,30)/
     &'(T1,"(JSTRH(1).NE.1).OR.(JSTRH(I).GE.JSTRH(I+1)) Column:",I6)',
     &'(T1,"(NZHDIM.LE.0) .OR. (NZHDIM.GT.NDNSH) .OR. (NZHDIM.NE.NONZH)"
     &)',
     &'(T1,"(IROWH(*) IS INCORRECT)")',
     &'(T1,"(XUPR.LT.XLWR) Variable:",I6)',
     &'(T1,"((XUPR.EQ.XLWR).AND.(ISTATV.NE.3)).OR.((XUPR.NE.XLWR).AND.(I
     &STATV.EQ.3)) Variable:",I6)',
     &'(T1,"(XUPR.NE.XLWR).AND.(|XUPR-XLWR|.LT.CONTOL) Variable:",I6)',
     &'(T1,"REAL HOLD ARRAY TOO SMALL")',
     &'(T1,"INTEGER HOLD ARRAY TOO SMALL")',
     &'(T1,"FUNCTION ERROR AT INITIAL POINT OR DURING GRADIENT EVAL.")',
     &'(T1,"(MACTIV.GT.NDIM)")'/
C
      DATA  (FMTMSG(I),I=31,40)/
     &'(T1,"INSUFFICIENT REAL STORAGE")',
     &'(T1,"INSUFFICIENT INTEGER STORAGE")',
     &'(T1,"SINGULAR JACOBIAN ON SUCCESSIVE ITERATIONS")',
     &'(T1,"INSNLP INPUT ERROR")',
     &'(T1,"(IFERR.LT.0) .OR. (IFERR.GT.1)")',
     &'(T1,"(|SFZTOL|.GE.1) .OR. (|SFZTOL|.LE.ONEEM5)")',
     &'(T1,"CONFLICT BETWEEN USER AND MULTIFRONTAL FILE NUMBER")',
     &'(T1,"(NRES .EQ. 0) .AND. (NEWTON.LT.0) .OR. (NEWTON.GT.1)")',
     &'(T1,"(NRES.LT.0)")',
     &'(T1,"(NRES.GT.MAXRES)")'/
C
      DATA  (FMTMSG(I),I=41,50)/
     &'(T1,"(NONZR.LE.0) .OR. (NONZR .GT. NDIM*NRES)")',
     &'(T1,"(JCOLR(*).LE.0) .OR. (JCOLR(*).GT.NDIM)")',
     &'(T1,"(IROWR(*).LE.0) .OR. (IROWR(*).GT.NRES)")',
     &'(T1,"A ROW OF THE RESIDUAL JACOBIAN IS IDENTICALLY ZERO")',
     &'(T1,"IROWH(JSTRH(I)).NE.I FOR SOME I")',
     &'(T1,"ALGOPT.EQ.LLSQ .AND. NRES.LE.0")',
     &'(T1,"UNEXPECTED ERROR")',
     &'(T1,"|IRVCOM| .GT. 1")',
     &'(T1,"INVALID INPUT FOR QPOPTN")',
     &'(T1,"QPOPTN.EQ.DENSE AND (ALGOPT.EQ.FME .OR. NRES.GT.0 .OR. SPARS
     &E JACOBIAN/HESSIAN)")'/
C
      DATA  (FMTMSG(I),I=51,60)/
     &'(T1,"INCONSISTENT Constraint:",I6)',
     &'(T1,"INVALID INPUT FOR KTOPTN")',
     &'(T1,"I/O ERROR (INSUFFICIENT DISK SPACE)")',
     &'(T1,"TOLKTC < 1")',
     &'(T1,"TOLPVT < 0 OR TOLPVT > .5")',
     &'(T1,"TOLFIL < 0")',
     &'(T1,"0 < BIGCON < CONTOL")',
     &'(T1,"WEAK SOLUTION FOUND (MULTIPLIERS NEAR ZERO)")',
     &'(T1,"(MEQUAL.EQ.NDIM).AND.ALGOPT.NE.F___")',
     &'(T1,"MAX. NO. OF CONSECUTIVE FUNCTION ERRORS")'/
C
      DATA  (FMTMSG(I),I=61,70)/
     &'(T1,"MAXIMUM NUMBER OF FUNCTION EVALUATIONS")',
     &'(T1,"SMALL STEP TERMINATION; SUBOPTIMAL FEASIBLE POINT FOUND")',
     &'(T1,"MAX. NO. OF ITER. IN OPTIMIZATION PHASE")',
     &'(T1,"MAX. NO. OF ITER. IN FEASIBILITY PHASE")',
     &'(T1,"FEASIBLE POINT NOT FOUND")',
     &'(T1,"MAX. NO. OF INTERVAL HALVES IN LINE SEARCH")',
     &'(T1,"EITHER MAX HESSIAN DIAGONAL OR VIOLATED SLOPE CONDITION; SUB
     &OPTIMAL FEASIBLE POINT FOUND")',
     &'(T1,"PROJ. GRAD. CALCULATION FAILED")',
     &'(T1,"LAGRANGE MULTIPLIERS NOT COMPUTED; DEGENERATE CONSTRAINTS")'
     &,
     &'(T1,"SUBOPTIMAL FEASIBLE POINT FOUND")'/
C
      DATA  (FMTMSG(I),I=71,NERMSG)/
     &'(T1,"SCHUR-QP FAILED WITH UNEXPECTED ERROR")',
     &'(T1,"CONTOL.GT.OBJTOL")',
     &'(T1,"UPHILL DIRECTION DETECTED IN LINE SEARCH")',
     &'(T1,"REDUCED OBJECTIVE FUNCTION IS LINEAR")',
     &'(T1,"CONSTRAINTS IGNORED")',
     &'(T1,"TERMINATE AFTER DIAGNOSTIC LINE SEARCH")',
     &'(T1,"RECURSIVE HESSIAN ESTIMATE AND NEWTON.NE.0")',
     &'(T1,"TERMINATE AFTER POSTOPTIMALITY ANALYSIS")',
     &'(T1,"TERMINATE AFTER USER ABORT")',
     &'(T1,"CONTOL .LE. BIGCON .LE. 1.D-4")',
     &'(T1,"INVALID INPUTS FOR (IPUDRF,ITDRQP) or (IPUFZF,ITFZQP)")',
     &'(T1,"USER EXTERNAL KILL")',
     &'(T1,"(NRES .GT. 0) .AND. (NEWTON.LT.0) .OR. (NEWTON.GT.2)")'/
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      IF(IRVCOM.LT.0) THEN
C
C             (>>>first pass begins)
C
C
C         INITIALIZE THE TOTAL TIME FOR SOLVING THE NLP ON CLOCK 1
C
        CALL CLKSET(1)
        CALL CLKBEG(1)
C
C         INITIALIZE THE TOTAL TIME OUTSIDE THE NLP ON CLOCK 3
C
        CALL CLKSET(3)
C
C         INITIALIZE THE TOTAL TIME IN SCHUR-QP ON CLOCK 4
C
        CALL CLKSET(4)
C
C         INITIALIZE THE MAX TIME IN SCHUR-QP ON CLOCK 5
C
        CALL CLKSET(-5)
C
C         INITIALIZE THE TOTAL TIME FOR K-T FACTORIZATIONS ON CLOCK 6
C
        CALL CLKSET(6)
C
C         INITIALIZE THE TOTAL TIME FOR MATRIX-VECTOR PRODUCTS IN QPGETD
C         ON CLOCK 11
C         INITIALIZE THE TOTAL TIME FOR MATRIX-VECTOR PRODUCTS IN LPCORE
C         ON CLOCK 12
C         INITIALIZE CLOCKS FOR DENSE QP (QPOPT) ON CLOCKS 13-18
C         INITIALIZE CLOCKS FOR DENSE QP (QPCORE) ON CLOCKS 19-26
C         INITIALIZE CLOCKS FOR DENSE QP (QPCORE) ON CLOCKS 27-31
C
        DO II = 11,31
          CALL CLKSET(II)
        ENDDO
C
C         CLOCK NUMBER ASSIGNMENT 
C
C     CLOCK           PROGRAM
C
C     1               SNLP
C     2                    SOCX
C     3               SNLP
C     4               SNLP
C     5               SNLP
C     6               SNLP
C     7                    SOCX
C     8                    SOCX,DENSTS
C     9                    SOCX,DENSTS
C     10              total run 
C     11              SNLP
C     12              SNLP
C     13              SNLP -- dual use sparse or dense QP
C     :               :
C     31              SNLP
C     32                   SOCX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- STORAGE POINTER SEQUENCE ---------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C             INPUT ERROR CHECKING IS PERFORMED THEN
C             STORAGE POINTERS WHICH DIVIDE THE WORK ARRAY INTO
C             ITS SUBARRAYS ARE SET.  THE STORAGE ALLOCATION FLAG
C             IS TURNED OFF.
C
        IF(NRES.EQ.0) THEN
          SUBNAM='SNLPMN  '
        ELSE
          SUBNAM='SNLPLS  '
        ENDIF
C
C             SET SUCCESS/ERROR FLAG.
C
        IERNLP = 0
        IMESSG = 0
        MODE = 0
        CALL NPESET(SUBNAM,IER,ERRMSG,-1)
C
C             SET DEFAULTS IF THEY HAVE NOT BEEN SET
C
        IF(JHSN74(1).EQ.-1) THEN
          CALL INSNLP('SPARSE DEFAULT')
        ENDIF
        CALL HHERRS(1,IPUNLP)
C
        IF(INNPER.NE.0) THEN
          IERNLP = -134
          GO TO 270
        ENDIF
C
        IF(IRVCOM.LT.-1) THEN
          IERNLP = -148
          GO TO 270
        ENDIF
C
C             DEFINE INTERNAL PRINT CONTROL IOFLAG FROM INPUT AND IPUNLP.
C
        IF(IPUNLP.LE.0) IOFLAG = 0
C
C             CHECK NLPSPR INPUT QUANTITIES.
C
        IF((MCON.LT.0)) THEN
          IERNLP = -101
          GO TO 270
        ENDIF
        IF(NDIM.LT.1) THEN
          IERNLP = -102
          GO TO 270
        ENDIF
        IF(MCON.GT.MAXCON) THEN
          IERNLP = -103
          GO TO 270
        ENDIF
        IF(ABS(IHESHN).GT.3) THEN
          IERNLP = -104
          GO TO 270
        ENDIF
        IF(NITMAX.LT.MAX(NITMIN,1)) THEN
          IERNLP = -105
          GO TO 270
        ENDIF
        IF(IT1MAX.LT.1) THEN
          IERNLP = -107
          GO TO 270
        ENDIF
        IF((IOFLAG.LT.0) .OR.(IOFLAG.GT.30)) THEN
          IERNLP = -108
          GO TO 270
        ENDIF
        IF(IOFLIN.LT.0) IOFLIN = IOFLAG
        IF(NRES.LT.0) THEN
          IERNLP = -139
          GO TO 270
        ENDIF
        IF(NRES.GT.MAXRES) THEN
          IERNLP = -140
          GO TO 270
        ENDIF
C
C             TURN OFF BCSLIB ERROR MESSAGES.
C
        IF(IOFLAG.EQ.0) CALL HHERPT(0)
C
C         CONSTRUCT INTEGER FLAG MINOPT FROM INPUT CHARACTER VARIABLE ALGOPT
C
        IF(ALGOPT(1:4).EQ.'FM  ') THEN
          MINOPT = 1
        ELSEIF(ALGOPT(1:4).EQ.'FME ') THEN
          MINOPT = 2
        ELSEIF(ALGOPT(1:4).EQ.'M   ') THEN
          MINOPT = 3
        ELSEIF(ALGOPT(1:4).EQ.'F   ') THEN
          MINOPT = 4
        ELSEIF(ALGOPT(1:4).EQ.'LLSQ') THEN
          MINOPT = 5
        ELSE
          IERNLP = -109
          GO TO 270
        ENDIF
C  
        IF(ALGOPT.EQ.'LLSQ'.AND.NRES.LE.0) THEN
          IERNLP = -146
          GO TO 270
        ENDIF
C
        IF(QPOPTN.NE.'SPARSE'.AND.QPOPTN.NE.'DENSE ') THEN
          IERNLP = -149
          GO TO 270
        ENDIF
C
        IF(QPOPTN.EQ.'DENSE '.AND.ALGOPT(1:4).EQ.'FME ') THEN
          IERNLP = -150
          GO TO 270
        ENDIF
        IF(QPOPTN.EQ.'DENSE '.AND.NRES.GT.0) THEN
          IERNLP = -150
          GO TO 270
        ENDIF
C
        IF(KTOPTN.NE.'SMALL '.AND.KTOPTN.NE.'LARGE ') THEN
          IERNLP = -152
          GO TO 270
        ENDIF
C
        IF(TOLKTC.LT.ONE) THEN
          IERNLP = -154
          GO TO 270
        ENDIF
C
        IF(TOLPVT.LT.ZERO.OR.TOLPVT.GT.POINT5) THEN
          IERNLP = -155
          GO TO 270
        ENDIF
C
        IF(TOLFIL.LT.ZERO) THEN
          IERNLP = -156
          GO TO 270
        ENDIF
C
        IF( (ITDRQP.GT.0.AND.IPUDRF.LE.0) .OR. 
     $    (ITFZQP.GT.0.AND.IPUFZF.LE.0) ) THEN
          IERNLP = -158
          GO TO 270
        ENDIF
C
C             SPARSE MATRIX CHECK FLAG
C
        SPRSCK = QPOPTN.EQ.'SPARSE'
C
C             COMPUTE MACHINE PRECISION ZERO.
C
        ZEROMN = HDMCON(5)
C
C             COMPUTE SQUARE ROOT OF MACHINE PRECISION
C
        ZEROOT = SQRT(ZEROMN)
C
C             COMPUTE LARGEST MACHINE CONSTANT.
C
        BIGNUM = HDMCON(2)
C
C             COMPUTE SQUARE ROOT OF BIGNUM.
C
        BGROOT = SQRT(BIGNUM)
C
C             BIG BOUND VALUE
C
        BIGBND = ONEEM2/ZEROMN
C
C             BIG CONDITION NUMBER
C
        BIGCND = (ONE/ZEROMN)**.8
C
        IF(SLPTOL.GE.ONE.OR.SLPTOL.LE.ONEEM5) THEN
          IERNLP = -106
          GO TO 270
        ENDIF
        IF(OBJTOL.LE.ONEEP1*ZEROMN) THEN
          IERNLP = -110
          GO TO 270
        ENDIF
        IF((PGDTOL.LT.ZEROOT).OR.(PGDTOL.GT.ONEEM2)) THEN
          IERNLP = -111
          GO TO 270
        ENDIF
        IF(ABS(SFZTOL).GE.ONE.OR.ABS(SFZTOL).LE.ONEEM5) THEN
          IERNLP = -136
          GO TO 270
        ENDIF
C
C        CHECK CONSTRAINT AND VARIABLE INPUTS FOR CONSISTENCY
C
        IF(NONZG.LE.0) THEN
          IERNLP = -114
          GO TO 270
        ENDIF
C
        IF(MCON.GT.0) THEN
C
          IF(CONTOL.LT.ZEROOT) THEN
            IERNLP = -112
            GO TO 270
          ENDIF
C
          IF(BIGCON.GT.ZERO.AND.BIGCON.LT.CONTOL) THEN
            IERNLP = -157
            GO TO 270
          ELSEIF(BIGCON.GE.CONTOL.AND.BIGCON.LE.ONEEM4) THEN
            IERNLP = +123
            CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
          ENDIF
C
          IF(CONTOL.GT.OBJTOL) THEN
            IERNLP = +115
            CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
          ENDIF
C
          ICONBD = 0
          DO I=1,MCON
            IF (ISTATC(I).GT.4 .OR. ISTATC(I).LT.0)  ICONBD = ICONBD + 1
          ENDDO
          IF(ICONBD.GT.0) THEN
            IERNLP = -113
            GO TO 270
          ENDIF
C
          XNDIM = NDIM
          XMCON = MCON
          XNONZG = NONZG
          IF(XNONZG .GT. XNDIM*XMCON) THEN
            IERNLP = -114
            GO TO 270
          ENDIF
C
          XNDIM = NDIM
          XNDNSH = XNDIM*(XNDIM+ONE)/TWO
          XNONZH = NONZH
          IF(QPOPTN.EQ.'DENSE '.AND.(XNONZG.LT.XNDIM*XMCON
     $    .OR.XNONZH.LT.XNDNSH)) THEN
            IERNLP = -150
            GO TO 270
          ENDIF
C
          IF(SPRSCK) THEN
            JCOLMN = 1
            JCOLMX = NDIM
            DO I=1,NONZG
              JCOLMN = MIN(JCOLMN,JCOLG(I))
              JCOLMX = MAX(JCOLMX,JCOLG(I))
            ENDDO
            IF(JCOLMN.LE.0 .OR. JCOLMX.GT.NDIM) THEN
              IERNLP = -115
              GO TO 270
            ENDIF
C
            IROWMN = 1
            IROWMX = MCON
            DO I=1,NONZG
              IROWMN = MIN(IROWMN,IROWG(I))
              IROWMX = MAX(IROWMX,IROWG(I))
            ENDDO
            IF(IROWMN.LE.0 .OR. IROWMX.GT.MCON) THEN
              IERNLP = -116
              GO TO 270
            ENDIF
          ENDIF
C
          DO I = 1,MCON
C
            IF(CUPR(I).LT.CLWR(I)) THEN
              IMESSG = I
              IERNLP = -117
              GO TO 270
            ENDIF
C
            IF(CUPR(I).EQ.CLWR(I)) THEN
              IF(ISTATC(I).LT.3.OR.ISTATC(I).GT.4) THEN
                IMESSG = I
                IERNLP = -118
                GO TO 270
              ENDIF
            ELSEIF(ISTATC(I).EQ.3) THEN
              IMESSG = I
              IERNLP = -118
              GO TO 270
            ELSE
              IF(ABS(CUPR(I)-CLWR(I)).LT.CONTOL) THEN
                IMESSG = I
                IERNLP = -119
                GO TO 270
              ENDIF
            ENDIF
C
          ENDDO
C  
        ENDIF
C
        IVARBD = 0
        DO I=1,NDIM
          IF (ISTATV(I).GT.3 .OR. ISTATV(I).LT.0)  IVARBD = IVARBD + 1
        ENDDO
        IF(IVARBD.GT.0) THEN
          IERNLP = -120
          GO TO 270
        ENDIF
C
        IF(NONZH.LE.0) THEN
          IERNLP = -122
          GO TO 270
        ENDIF
C
        IF(MINOPT.LT.4) THEN
C
          IF(SPRSCK) THEN
C
            XNDIM = NDIM
            NZHDIM = JSTRH(NDIM+1)-1
            XNZHDM = NZHDIM
            XNDNSH = XNDIM*(XNDIM+ONE)/TWO
            IF(JSTRH(1).NE.1) THEN
              IMESSG = 1
              IERNLP = -121
              GO TO 270
            ENDIF
C
            IF(NZHDIM.LT.NDIM.OR.XNZHDM.GT.XNDNSH.OR.NZHDIM.NE.NONZH) 
     $        THEN
              IERNLP = -122
              GO TO 270
            ENDIF
C
            IROWMN = 1
            IROWMX = NDIM
            DO I=1,NONZH
              IROWMN = MIN(IROWMN,IROWH(I))
              IROWMX = MAX(IROWMX,IROWH(I))
            ENDDO
            IF(IROWMN.LE.0 .OR. IROWMX.GT.NDIM) THEN
              IERNLP = -123
              GO TO 270
            ENDIF
          ENDIF
C
        ELSEIF(MINOPT.EQ.5) THEN
C
          IF(NONZH.NE.NDIM) THEN
            IERNLP = -122
            GO TO 270
          ELSE
            JSTRH(NDIM+1) = NDIM + 1
          ENDIF
C
        ENDIF
C
        DO I = 1,NDIM
C
          IF(SPRSCK) THEN
            IF(MINOPT.LT.4) THEN
              IF(JSTRH(I).GE.JSTRH(I+1)) THEN
                IMESSG = I
                IERNLP = -121
                GO TO 270
              ENDIF
              IF(IROWH(JSTRH(I)).NE.I) THEN
                IERNLP = -145
                GO TO 270
              ENDIF
C
            ELSEIF(MINOPT.EQ.5) THEN
C
C         SET UP THE HESSIAN QUANTITIES FOR LINEAR LEAST SQUARES CASE
C
              JSTRH(I) = I
              IROWH(I) = I
              HMAT(I) = ZERO
C
            ENDIF
          ENDIF
C
          IF(XUPR(I).LT.XLWR(I)) THEN
            IMESSG = I
            IERNLP = -124
            GO TO 270
          ENDIF
C
          IF(XUPR(I).EQ.XLWR(I)) THEN
            IF(ISTATV(I).NE.3) THEN
              IMESSG = I
              IERNLP = -125
              GO TO 270
            ELSE
              XBAR(I) = XLWR(I)
            ENDIF
          ELSEIF(ISTATV(I).EQ.3) THEN
            IMESSG = I
            IERNLP = -125
            GO TO 270
          ELSE
            IF(ABS(XUPR(I)-XLWR(I)).LT.CONTOL) THEN
              IMESSG = I
              IERNLP = -126
              GO TO 270
            ENDIF
          ENDIF
C
        ENDDO
C
        IF(SPRSCK) THEN
C
C         CHECK FOR MONOTONICITY OF IROWH WITHIN A COLUMN
C
          DO JCOL = 1,NDIM
            IROWCK = JCOL
            DO I = JSTRH(JCOL),JSTRH(JCOL+1)-1
              IF(IROWH(I).LT.IROWCK.OR.IROWH(I).GT.NDIM) THEN
                IERNLP = -123
                GO TO 270
              ENDIF
              IROWCK = IROWH(I)
            ENDDO
          ENDDO
        ENDIF
C
C        CHECK RESIDUAL INPUTS FOR CONSISTENCY
C
        IF(NONZR.LT.0.AND.ALGOPT.NE.'F') THEN
          IERNLP = -141
          GO TO 270
        ENDIF
C
        IF(NRES.GT.0) THEN
C
          XNONZR = NONZR
          XNDIM = NDIM
          XNRES = NRES
          IF(XNONZR .GT. XNDIM*XNRES) THEN
            IERNLP = -141
            GO TO 270
          ENDIF
C
          JCOLMN = 1
          JCOLMX = NDIM
          DO I=1,NONZR
            JCOLMN = MIN(JCOLMN,JCOLR(I))
            JCOLMX = MAX(JCOLMX,JCOLR(I))
          ENDDO
          IF(JCOLMN.LE.0 .OR. JCOLMX.GT.NDIM) THEN
            IERNLP = -142
            GO TO 270
          ENDIF
C
          IROWMN = 1
          IROWMX = NRES
          DO I=1,NONZR
            IROWMN = MIN(IROWMN,IROWR(I))
            IROWMX = MAX(IROWMX,IROWR(I))
          ENDDO
          IF(IROWMN.LE.0 .OR. IROWMX.GT.NRES) THEN
            IERNLP = -143
            GO TO 270
          ENDIF
C  
        ENDIF
C
C         SPECIAL SETTINGS FOR NORMAL MATRIX TREATMENT OF LEAST SQUARES 
C         PROBLEMS.   CURRENTLY IF NORMAL > 0, THE NORMAL MATRIX WILL
C         BE FORMED, INSTEAD OF THE SPARSE TABLEAU FORMAT.  THE NEWTON
C         FLAG IS ALSO USED TO FORCE CALCULATION OF THE RESIDUAL HESSIAN
C         FOR EACH ITERATION.   THE HESSIAN H, IS STORED IN THE LOCATION
C         OF V TO SAVE STORAGE.   IF V IS SAVED (USING ADDITIONAL STORAGE),
C         THE NEWTON FLAG DOES NOT NEED TO BE SET TO 1.
C
C         DEFINE NRESNP = # OF RESIDUALS FOR NLP PROBLEM
C                       = NRES FOR SPARSE TABLEAU FORMAT (NORMAL=0)
C                       = 0    FOR NORMAL MATRIX FORMAT (NORMAL.NE.0)
C
        NRESNP = NRES
        IF(NRES.GT.0.AND.NORMAL.NE.0) THEN
          NRESNP = 0
          NEWTON = 1
        ENDIF
C
C             CHECK STORAGE ALLOCATION FOR ALGORITHM DEPENDENT ARRAYS
C             I.E. THE WORK ARRAY HOLD.
C
        NEQL = 0
        DO I=1,NDIM
          IF (ISTATV(I).EQ.3)  NEQL = NEQL + 1
        ENDDO
        MINEQL = 0
        DO I=1,MCON
          IF (ISTATC(I).EQ.3) THEN
            NEQL = NEQL + 1
          ELSEIF (ISTATC(I).GE.0 .AND. ISTATC(I).LE.2) THEN
            MINEQL = MINEQL + 1
          ENDIF
        ENDDO
C
C         CHECK THAT ALGORITHM OPTION IS CONSISTENT FOR SEARCHS
C
        IF(NDIM.EQ.NEQL.AND.MINOPT.LT.4.AND.SPRSCK) THEN
          IERNLP = +102
          CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
          MINOPT = 4
        ENDIF
C
        MACTVB = MIN(NDIM,MCON)
C
C         COMPUTE LOWER BOUNDS ON THE WORK ARRAY STORAGE FOR THE
C         VARIOUS ALGORITHM ROUTINES
C
C         NOTE THE CURRENT STORAGE SCHEME SPLITS THE WORK ARRAY
C         INTO TWO PARTS THAT ARE USED SIMULTANEOUSLY BY PRJGRD
C         AND CONELM--HENCE THE FACTOR OF 2 IN THE DEFINITION OF LNHOLD
C         
C         COMPUTE ESTIMATE FOR THE LENGTH OF THE REAL AND INTEGER WORK
C         ARRAY.  THE CALCULATION IS ORGANIZED BY COMPUTING THE 
C         ESTIMATES FOR EACH CALLING SEQUENCE THAT REQUIRES 
C         A CALL TO THE MULTIFRONTAL CODE.
C
C         (1)  LSDPQP ---> SRCHFZ ---> NLPSPR ---> SPRNLP
C
        NZHLDP = NDIM
        NZGLDP = NONZG
        IF(SPRSCK) THEN
          CALL SQPSTR(MCON,MINEQL,NDIM,NZHLDP,NZGLDP,IWKLDP,IRWLDP)
        ELSE
          IWKLDP = 0
          IRWLDP = 11*NDIM + 9*MCON + 2*NDIM**2
        ENDIF
        IWKLDP = IWKLDP + 2*NDIM + MCON + 1
        IRWLDP = IRWLDP + 6*NDIM + 3*MCON 
C
C         (2)  LSQRQP ---> SRCHDR ---> EQCMIN ---> NLPSPR ---> SPRNLP
C
        IF(NRES.GT.0.AND.NORMAL.EQ.0) THEN
          NDIMA = NDIM + NRES
          MCONA = MCON + NRES
          NZHLSQ = NONZH + NRES
          NZGLSQ = NONZG+NONZR+NRES
          CALL SQPSTR(MCONA,MINEQL,NDIMA,NZHLSQ,NZGLSQ,IWKLSQ,IRWLSQ)
          IWKLSQ = IWKLSQ + 2*NDIMA + MCONA + 1 + 2*NRES + NONZH + NONZG
     $         + NONZR + MAX(NONZH + NONZG + NRES,NDIMA +1)
          IWKLSQ = IWKLSQ + MCON + NRES
          IRWLSQ = IRWLSQ + 3*NDIMA + MCONA + 2*NRES + NONZH + NONZG 
     $         + NONZR
          IRWLSQ = IRWLSQ + 4*NDIM + 2*MCON + 5*NRES
        ELSE
          IWKLSQ = 0
          IRWLSQ = 0
        ENDIF
C
C         (3)  QPSTRT ---> NLPSPR ---> SPRNLP
C
        NZHQPS = NDIM
        NZGQPS = NONZG
        IF(SPRSCK) THEN
          CALL SQPSTR(MCON,MINEQL,NDIM,NZHQPS,NZGQPS,IWKQPS,IRWQPS)
        ELSE
          IWKQPS = 0
          IRWQPS = 11*NDIM + 9*MCON + 2*NDIM**2
        ENDIF
        IWKQPS = IWKQPS + 2*NDIM + MCON + 1
        IRWQPS = IRWQPS + 4*NDIM + 2*MCON 
C
C         (4)  RELXQP ---> SRCHFZ ---> NLPSPR ---> SPRNLP
C
        NDIMV = NDIM + MCON
        NZHRLX = NDIMV
        NZGRLX = NONZG + MCON
        IF(SPRSCK) THEN
          CALL SQPSTR(MCON,MINEQL,NDIMV,NZHRLX,NZGRLX,IWKRLX,IRWRLX)
          IWKRLX = IWKRLX + 3*NDIMV + 3*MCON + 1 + NONZG 
     $       + MAX(NZGRLX,NDIMV+1)
          IRWRLX = IRWRLX + 7*NDIMV + 4*MCON + NONZG
        ELSE
          IWKRLX = 0
          IRWRLX = 0
        ENDIF
C
C         (5)  SRCHDR ---> EQCMIN ---> NLPSPR ---> SPRNLP
C
        IF(NRES.EQ.0) THEN
          NZHSDR = NONZH
          NZGSDR = NONZG
          IF(SPRSCK) THEN
            CALL SQPSTR(MCON,MINEQL,NDIM,NZHSDR,NZGSDR,IWKSDR,IRWSDR)
          ELSE
            IWKSDR = 0
            IRWSDR = 11*NDIM + 9*MCON + 2*NDIM**2
          ENDIF
          IWKSDR = IWKSDR + MCON 
          IRWSDR = IRWSDR + 4*NDIM + 2*MCON 
        ELSE
          IWKSDR = 0
          IRWSDR = 0
        ENDIF
C
C         CONSTRUCT THE FINAL ESTIMATE (LARGEST OF ALL BRANCHES)
C
        LNIHLD = MAX(IWKLDP,IWKLSQ,IWKQPS,IWKRLX,IWKSDR)
        LNHOLD = MAX(IRWLDP,IRWLSQ,IRWQPS,IRWRLX,IRWSDR)
        IF(MINOPT.EQ.2) THEN
          LNIHLD = 2*LNIHLD
          LNHOLD = 2*LNHOLD
        ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C        ALLOCATION FOR THE REAL ARRAY
c
C        CVEC     VALUES OF CONSTRAINTS
C        COLD     OLD VALUES OF CONSTRAINTS
C        PGRD     PROJECTED GRADIENT OF F
C        SVEC     SEARCH DIRECTION VECTOR
C        XOLD     OLD ESTIMATE OF X
C        YVEC     TRANSFORMED (INTERNAL) VARIABLES
C        YBAR     TRANSFORMED (INTERNAL) VARIABLES
C
        NRESSZ = NRES
        IF(NORMAL.NE.0) NRESSZ = 0
        IF(BIGCON.GT.ZERO.AND.MINOPT.NE.2) THEN
          LNISCR = MAX(100,2*NITMAX)
        ELSE
          LNISCR = 200
        ENDIF
c
C ======================================================================
c
        if(allocated(cvec)) deallocate(cvec)
        if(allocated(cold)) deallocate(cold)
        if(allocated(pgrd)) deallocate(pgrd)
        if(allocated(svec)) deallocate(svec)
        if(allocated(xold)) deallocate(xold)
        if(allocated(yvec)) deallocate(yvec)
        if(allocated(ybar)) deallocate(ybar)
        if(allocated(irwg)) deallocate(irwg)
        if(allocated(jstr)) deallocate(jstr)
        if(allocated(iprmg)) deallocate(iprmg)
        if(allocated(iprmc)) deallocate(iprmc)
        if(allocated(isto)) deallocate(isto)
        if(allocated(irwr)) deallocate(irwr)
        if(allocated(jsrr)) deallocate(jsrr)
        if(allocated(iprmr)) deallocate(iprmr)
        if(allocated(iscr)) deallocate(iscr)

        allocate(cvec(1:maxcon))
        allocate(cold(1:maxcon))
        allocate(pgrd(1:ndim))
        allocate(svec(1:ndim+maxcon+1))
        allocate(xold(1:ndim))
        allocate(yvec(1:ndim))
        allocate(ybar(1:ndim))
c
c ======================================================================
c
C         CHECK LNHOLD AGAINST INPUT
C
        NEEDED = 0
        IF(LNHOLD.GT.NHOLD) THEN
          NEEDED = LNHOLD
C
C           IF THIS IS A LEAST SQUARES PROBLEM INCLUDE GRADIENT
C           OFFSET FROM SNLPLS
C
          IF(NRES.GT.0.AND.NORMAL.EQ.0) NEEDED = NEEDED + NDIM
          IERNLP = -127
          GO TO 270
        ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C             ALLOCATION FOR THE INTEGER ARRAY
C
C        IRWG    NONZERO ROWS IN JACOBIAN.
C        JSTR    NONZERO COLUMNS IN JACOBIAN.
C        IPRMG   PERMUTATION TO EXTERNAL ORDER.
C        IPRMC   PERMUTATION OF CONSTRAINT ROWS
C        ISTO    OLD VARIABLE/CONSTRAINT STATUS
C        IRWR    NONZERO ROWS IN RESIDUAL JACOBIAN
C        JSRR    NONZERO COLUMNS IN RESIDUAL JACOBIAN 
C        IPRMR   RESIDUAL JACOBIAN PERMUTATION
C
        LNJCOL = MAX(NONZG,NDIM+1)
        LNRCOL = MAX(NONZR,NDIM+1)
C
        IF(SPRSCK) THEN
c
C ======================================================================
c
          allocate(irwg(1:nonzg))
          allocate(jstr(1:lnjcol))
          allocate(iprmg(1:nonzg))
          allocate(iprmc(1:maxcon))
          allocate(isto(1:ndim+maxcon))
          allocate(irwr(1:nonzr))
          allocate(jsrr(1:lnrcol))
          allocate(iprmr(1:lnrcol))
          allocate(iscr(1:lniscr))
c
C ======================================================================
c
        ELSE
c
C ======================================================================
c
          allocate(irwg(1:1))
          allocate(jstr(1:1))
          allocate(iprmg(1:1))
          allocate(iprmc(1:maxcon))
          allocate(isto(1:ndim+maxcon))
          allocate(irwr(1:1))
          allocate(jsrr(1:1))
          allocate(iprmr(1:1))
          allocate(iscr(1:lniscr))
c
C ======================================================================
c
        ENDIF
C      
C         COMPUTE LENGTH FOR INTEGER HOLD ARRAY EXCLUDING WORK ARRAY
C
        IF (LNIHLD.GT.NIHOLD) THEN
          NEEDED = LNIHLD
          IERNLP = -128
          GO TO 270
        ENDIF
C
        NZERO = 0
        IF(MCON.GT.0) THEN
C
C         CHECK THE INPUT BOUNDS FOR LARGE VALUES; IF THEY ARE
C         HUGE RESET ISTATC AND PRINT A WARNING
C
          NINFBD = 0
          DO I = 1,MCON
            IF(CUPR(I).GE.BIGBND.AND.CLWR(I).LE.-BIGBND
     $        .AND.ISTATC(I).NE.4) THEN
              NINFBD = NINFBD + 1
              IHOLD(NINFBD) = I
              ISTATC(I) = 4
            ENDIF
          ENDDO
C
          IF(NINFBD.GT.0) THEN
            IERNLP = +118
            CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
            IF(IOFLAG.GE.10) THEN
              WRITE(IPUNLP,1008)
              CALL INTOUT(IHOLD,NINFBD,IPUNLP)
            ENDIF
          ENDIF
C
C         CHECK FOR ZERO ROWS IN THE JACOBIAN MATRIX
C
C         COUNT THE NUMBER OF NONZEROS IN ROW I AND STORE
C         THE RESULT IN IHOLD(I)
C
          IF(SPRSCK) THEN
C
            IHOLD(1:MCON) = 0
C
            DO K = 1,NONZG
              IR = IROWG(K)
              JC = JCOLG(K)
              IF(ISTATV(JC).NE.3.AND.ISTATC(IR).NE.4) THEN
                IHOLD(IR) = IHOLD(IR) + 1
              ENDIF
            ENDDO
C
          ELSE
            IHOLD(1:MCON) = 1
          ENDIF
C
C         CHECK THAT EACH ROW HAS AT LEAST ONE ELEMENT AND SAVE THE
C         ROW NUMBER OF THE ROWS WITH NO ELEMENTS
C
          DO I = 1,MCON
            IF(IHOLD(I).LE.0.AND.ISTATC(I).NE.4) THEN
              NZERO = NZERO + 1
              IHOLD(NZERO) = I
              ISTATC(I) = 4
            ENDIF
          ENDDO
C
          IF(NZERO.GT.0) THEN
            IERNLP = +118
            CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
            IF(IOFLAG.GE.10) THEN
              WRITE(IPUNLP,1008)
              CALL INTOUT(IHOLD,NZERO,IPUNLP)
            ENDIF
          ENDIF
C
        ENDIF
C
        IF(NRES.GT.0) THEN
C
C         CHECK FOR ZERO ROWS IN THE RESIDUAL JACOBIAN MATRIX
C
C         COUNT THE NUMBER OF NONZEROS IN ROW I AND STORE
C         THE RESULT IN IHOLD(I)
C
          IHOLD(1:NRES) = 0
C  
          DO K = 1,NONZR
            IR = IROWR(K)
            IHOLD(IR) = IHOLD(IR) + 1
          ENDDO
C
          NZROW = 0
          DO K = 1,NRES
            IF(IHOLD(K).LE.0) THEN
              NZROW = NZROW + 1
              IHOLD(NZROW) = K
            ENDIF
          ENDDO
          IF(NZROW.GT.0) THEN
            IERNLP = -144
            IF(IOFLAG.GE.10) WRITE(IPUNLP,1005) 
            CALL INTOUT(IHOLD,NZROW,IPUNLP)
            GO TO 270
          ENDIF
C
        ENDIF
C
C         COPY EXTERNAL SPARSITY DEFINITION INTO LOCAL ARRAYS (I.E.
C         SAVE THE EXTERNAL DEFINITIONS)
C
        IF(SPRSCK) THEN
          IRWG(1:NONZG) = IROWG(1:NONZG)
          IF(MCON.NE.0) THEN
            JSTR(1:NONZG) = JCOLG(1:NONZG)
          ELSE
            JSTR(1:LNJCOL) = 1
          ENDIF
C
          IRWR(1:NONZR) = IROWR(1:NONZR)
          IF(NRES.NE.0) THEN
            JSRR(1:NONZR) = JCOLR(1:NONZR)
          ENDIF
        ENDIF
C
C             SET FLAG FOR INITIAL SUMMARY PRINT.
C
        SUMRY = .TRUE.
C
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- FIRST ENTRANCE SEQUENCE ----------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C         INITIALIZE ALGORITHM PERFORMANCE STATISTICS.  ALL CODE
C         IN THIS SECTION (BETWEEN +++++) CAN BE ELIMINATED IF 
C         ALGORITHM PERFORMANCE IS NOT MONITORED.  ALL INFORMATION
C         IS CONTAINED IN TWO ARRAYS; INSTAT, AND RLSTAT.  THE 
C         PARAMETERS ARE DEFINED AS FOLLOWS:
C
C         INSTAT(1) = NO. OF GRADIENT CALLS FOR CONSTRAINT ELIMINATION
C         INSTAT(2) = NO. OF GRADIENT CALLS FOR REDUCED GRADIENTS
C         INSTAT(3) = NO. OF GRADIENT CALLS FOR FEASIBLE REGION FINDER 
C         INSTAT(4) = MAX. NO. ITER. PER CONSTRAINT ELIM.
C         INSTAT(5) = TOTAL NO. ITER. FOR CONSTRAINT ELIM.
C         INSTAT(6) = TOT. NO. CONSTRAINT ELIM.
C         INSTAT(7) = MAX. NO. ITER. WASTED IN CONSTRAINT ELIM.
C         INSTAT(8) = TOTAL NO. ITER. WASTED IN CONSTRAINT ELIM.
C         INSTAT(9) = TOT. NO. WASTED CONSTRAINT ELIM.
C         INSTAT(10)= STORAGE REQUIRED IN REAL HOLD ARRAY
C         INSTAT(11)= NO. OF SCHUR-QP FAILURES
C         INSTAT(12)= NO. OF SCHUR-QP CALLS
C         INSTAT(13)= NO. OF MOST EXPENSIVE SCHUR-QP CALL
C         INSTAT(14)= MAX. NO CALLS TO LINPACK
C         INSTAT(15)= ITERATION NO. FOR MAX. NO. CALLS TO LINPACK
C         INSTAT(16)= MAX. NO UPDATES IN SCHUR-COMP
C         INSTAT(17)= ITERATION NO FOR MAX. NO. UPDATES IN SCHUR-COMP
C         INSTAT(18)= MAX. SIZE OF MATRIX IN LINPACK CALL
C         INSTAT(19)= ITERATION NO. FOR MAX. SIZE MATRIX CALL
C         INSTAT(20)= MAX. NO CALLS TO MULTIFRONTAL SOLVE
C         INSTAT(21)= ITERATION NO. FOR MAX. NO. CALLS TO MULTIFRONTAL SOLVE
C         INSTAT(22)= NO. OF FCALLS
C         INSTAT(23)= NO. OF GCALLS
C         INSTAT(24)= NO. OF HCALLS
C         INSTAT(25)= TOTAL NO. OF FUNCTION EVAL.
C         INSTAT(26)= TOTAL NO. OF QP ITERATIONS
C         INSTAT(27)= MAX.  NO. OF QP ITERATIONS ON ANY QP CALL
C         INSTAT(28)= SET TO ONE FOR OUT-OF-CORE, ZERO OTHERWISE
C         INSTAT(29)= TOTAL NO. OF QP ITERATIONS IN PHASE I (DENSE QP)
C         INSTAT(30)= STORAGE REQUIRED IN INTEGER HOLD ARRAY (if .gt. 1)
C
C         RLSTAT(1) = AV. NO. ITER/CONST. ELIM.
C         RLSTAT(2) = AV. NO. GRADIENT CALLS/CONST. ELIM.
C         RLSTAT(6) = TIME OF MOST EXPENSIVE SCHUR-QP CALL
C
        INSTAT(1:30) = 0
        IF(QPOPTN.NE.'SPARSE') THEN
          INSTAT(10) = NHOLD
          INSTAT(30) = NIHOLD
        ENDIF
        RLSTAT(1:20) = ZERO
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C             THE INITIAL VALUES OF IREVRS(1) AND IREVRS(2) ARE SET.  THE
C             SYSTEM PRINT FLAG IS SET.  THE INITIAL X VALUE IS SAVED.
C             GO TO THE RETURN AFTER THIS SEQUENCE.
C
C
C             REQUEST FUNCTION AND GRADIENT EVALUATION.
C
        IREVRS(1) = 1
        IREVRS(2) = 1
        IF(MAXNFE.LE.1) IREVRS(2) = 0
        IREVRS(3) = 0
C
C         INITIALIZE MULTIPLIER VECTORS
C
        IF(ALGOPT(1:1).NE.'M'.AND.INDUAL.GE.0) INDUAL = 0
C
        IF(INDUAL.GT.0) THEN
C
          CALL DUALIN(VECLAM,ISTATC,MAXCON,MCON,VECNU,ISTATV,NDIM)
C
        ELSE
C
          VECLAM(1:MCON) = ZERO
          VECNU(1:NDIM) = ZERO
C
        ENDIF
C
C             SET PRINT FLAG IF REQUIRED.
C
        IREVRS(4) = 0
        IF (IOFLAG.GE.10)  IREVRS(4) = 1
C
C             STORE THE INITIAL X VALUES 
C
        XOLD(1:NDIM) = XBAR(1:NDIM)
C
C             SET CONTINUATION CALL FLAG.
C
        IRVCOM = 1
C
C             INITIALIZE ALGORITHM INFORMATION FLAG
C
        IREVRS(5) = 0
C
C             RETURN FOR A FUNCTION EVALUATION.
C
        GO TO 10000
C
C         END OF INITIALIZATION PROCESSING (>>>first pass ends)
C
      ELSE
C
        IF(IRVCOM.GT.1) THEN
          IERNLP = -148
          GO TO 270
        ENDIF
C
C         CLOCK FOR EVERYTHING BUT NLP TIME
C             
        CALL CLKSUM(3)
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- ERROR CHECKING SEQUENCE ----------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C
C             BECAUSE THIS IS THE SEQUENCE WHICH IS ENTERED AFTER
C             EACH FUNCTION EVALUATION, WE HANDLE THE MAXIMUM NUMBER
C             OF FUNCTION EVALUATIONS AND THE MAXIMUM NUMBER
C             OF CONSECUTIVE FUNCTION ERRORS IN THIS SEQUENCE.
C
C             CHECK FOR HUGE VALUES IN THE COMPUTED FUNCTIONS
C
      IF(NRES.EQ.0) THEN
        HUGEF = ABS(FBAR)
      ELSE
        HUGEF = DAMAX(NRES,RESVEC,1)
      ENDIF
      IF(MCON.GT.0) THEN
        HUGEF = MAX(HUGEF,DAMAX(MCON,CBAR,1))
      ENDIF
      IF(HUGEF.GT.BIGBND) IFERR = 1
C
C             CHECK FOR FUNCTION ERROR AT INITIAL POINT.
C
      IF(SUMRY.AND.IFERR.EQ.1) THEN
        IERNLP = -129
        GO TO 270
      ENDIF
C
      IF(IREVRS(2).NE.0.AND.IOFLAG.GE.25) THEN
C
C             DIAGNOSTIC OUTPUT.  WRITE THE DERIVATIVES.
C
        WRITE(IPUNLP,1003)
        IF(NRES.EQ.0) THEN
          CALL NZROUT(DELF,NDIM,IPUNLP,'GRADIENT VECTOR')
        ELSE
          ITITLE = 'RESIDUAL JACOBIAN   '
          CALL SRG228(1,NRES,NDIM,NONZR,RMAT,IROWR,JCOLR,
     &                ITITLE,IPUNLP)
        ENDIF
        IF(MCON.GT.0) THEN
          ITITLE = 'JACOBIAN MATRIX     '
          IF(SPRSCK) THEN
            IOPT = 1
          ELSE
            IOPT = 4
          ENDIF
          CALL SRG228(IOPT,MCON,NDIM,NONZG,GMAT,IROWG,JCOLG,
     &                ITITLE,IPUNLP)
        ENDIF
C
      ENDIF
C
      IF(IREVRS(3).EQ.2.AND.IOFLAG.GE.25) THEN
C
C             DIAGNOSTIC OUTPUT.  WRITE THE HESSIAN
C
          ITITLE = 'HESSIAN MATRIX     '
          IF(SPRSCK) THEN
            IOPT = 3
          ELSE
            IOPT = 5
          ENDIF
          CALL SRG228(IOPT,NDIM,NDIM,NONZH,HMAT,IROWH,JSTRH,
     &                ITITLE,IPUNLP)
C
      ENDIF
C
C             TEST FOR USER ABORT 
C
      IF(IFERR.EQ.-1000) THEN
C
        IRVCOM = 0
        IERNLP = +122
        GO TO 270
C
      ENDIF
C
C             CHECK FOR INCORRECT VALUE OF FUNCTION ERROR FLAG
C
      IF(IFERR.LT.0.OR.IFERR.GT.1) THEN
        IRVCOM = 0
        IERNLP = -135
        GO TO 270
      ENDIF
C
C             CHECK FOR MAXIMUM CONSECUTIVE FUNCTION ERRORS.
C
      IF(IREVRS(1).EQ.1.OR.IREVRS(2).GE.1) THEN
C
        IF (IFERR.EQ.1) THEN
          NOFERR = NOFERR + 1
          IF (NOFERR.GE.5) THEN
              IRVCOM = 0
              IERNLP = +103
              GO TO 270
          ENDIF
        ELSE
          NOFERR = 0
        ENDIF
C
      ENDIF
C
C             TEST FOR MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
C
      IF (NFEVAL.GE.MAXNFE) THEN
C
C          SET FUNCTION ERROR FLAG TO THE SPECIAL VALUE (-100)
C
        IFERR = -100
C
      ENDIF
C
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- INITIAL SUMMARY PRINT (SUMRY = .TRUE.) -------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C
C             THIS IS THE FIRST RETURN TO THE ALGORITHM WITH AN
C             EVALUATED FUNCTION.  
C
C             WRITE OUT THE INITIAL SUMMARY PRINT.
C
      IF(.NOT. SUMRY) GO TO 240
C
C             CHECK FOR MAX FUNCTION EVALUATIONS AT FIRST POINT
C
      IF(IFERR.EQ.-100) THEN
        IERNLP = +104
        GO TO 270
      ENDIF
C
      SUMRY = .FALSE.
C
C             SET IREVRS(1), IREVRS(2), AND IREVRS(4).
C
      IREVRS(1) = 0
      IREVRS(2) = 0
      IF (IOFLAG.GE.10) IREVRS(4) = 0
C
C ======================================================================
C ================= INITIAL JACOBIAN ANALYSIS ==========================
C ======================================================================
C
C         WHEN A ROW OF THE JACOBIAN HAS NO NONZEROS THE CONSTRAINT
C         IS IGNORED.  CHECK THAT THE VALUE IS CONSISTENT WITH BOUNDS
C
      zeroloop: DO II = 1,NZERO
        ICZERO = IHOLD(II)
        if(istatc(iczero).eq.4) cycle zeroloop
        IF(CUPR(ICZERO).LT.BIGBND.OR.CLWR(ICZERO).GT.-BIGBND) THEN
          IF(CBAR(ICZERO).GT.CUPR(ICZERO)
     $       .OR.CBAR(ICZERO).LT.CLWR(ICZERO)) THEN
            IERNLP = -151
            IMESSG = ICZERO
            GO TO 270
          ENDIF
        ENDIF
      ENDDO zeroloop
C
C         COMPUTE THE LARGEST ELEMENTS IN THE JACOBIAN
C
      BIGELM(1) = ZERO
      IF(SPRSCK.AND.MCON.GT.0) CALL BIGGLM(GMAT,IRWG,
     $    JSTR,NONZG,ISTATC,MCON,BIGELM,IROWBG,JCOLBG,
     $    LENBIG,IHOLD,IERNLP)
C
      IF(IERNLP.LT.0) GO TO 270
C
C         CONSTRUCT THE CONSTRAINT REORDERING PERMUTATIONS
C
      CALL SHUFUL( MAXCON ,MCON    ,MEQUAL , MIGNOR, MTOTAL  
     $     ,NDIM   ,NONZG   ,NONZGT ,CBAR   ,CLWR       
     $     ,CUPR   ,CVEC,VECLAM ,GMAT   ,IRWG
     $     ,JSTR   ,ISTATC ,IPRMC   
     $     ,IPRMG ,IHOLD ,IERNLP ,PRMUTE)
C
      IF(IERNLP.LT.0) GO TO 270
C
      IF(.NOT.PRMUTE) THEN
C
C         PERMUTATION ARRAYS ARE NOT NEEDED--DEALLOCATE STORAGE
C
c
C ======================================================================
c
        if(allocated(iprmg)) deallocate(iprmg)
        if(allocated(iprmc)) deallocate(iprmc)
c
C ======================================================================
        IF(SPRSCK) THEN
C ======================================================================
c
          allocate(iprmg(1:1))
          allocate(iprmc(1:1))
c
C ======================================================================
c
          IRWR(1:NONZR) = IROWR(1:NONZR)
          IF(NRES.NE.0) THEN
            JSRR(1:NONZR) = JCOLR(1:NONZR)
          ENDIF
c
        ELSE
C ======================================================================
c
          allocate(iprmg(1:1))
          allocate(iprmc(1:1))
c
C ======================================================================
        ENDIF
C
      ENDIF
C
C ======================================================================
C ================= END INITIAL JACOBIAN ANALYSIS ======================
C ======================================================================
C
C ======================================================================
C ================= INITIAL RESIDUAL JACOBIAN ANALYSIS =================
C ======================================================================
C
      IF(NRES.GT.0) THEN
C
C         CONSTRUCT PERMUTATION TO TAKE RESIDUAL JACOBIAN FROM
C         EXTERNAL TO INTERNAL ORDER.  FIRST INITIALIZE MASTER 
C         JACOBIAN PERMUTATION RECORD
C
C
C         CONSTRUCT COLUMN STORAGE SCHEME FOR RMAT (NOTE:
C         THE CALCULATION OF THE DENSE MATRIX LOCATION BELOW
C         MUST BE ABLE TO COMPUTE AN INTEGER OF SIZE NDIM*NRES.
C         THIS IS DONE IN DOUBLE PRECISION TO AVOID INTEGER WORD
C         SIZE RESTRICTIONS.)
C
        DO II = 1,NONZR
          IR = IRWR(II)
          JC = JSRR(II)
          HOLD(II) = DBLE(IR) + DBLE(JC-1)*DBLE(NRES)
        ENDDO
C
C         CONSTRUCT PERMUTATION TO COLUMN ORDER
C
        CALL HDSRTN(HOLD,NONZR,0,0,IPRMR,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          GO TO 270
        ENDIF
C
C         REORDER RMAT,IROWR, AND JCOLR INTO COLUMN ORDER
C
        CALL HJPRMX(IRWR,NONZR,IPRMR,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          GO TO 270
        ENDIF
        CALL HJPRMX(JSRR,NONZR,IPRMR,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          GO TO 270
        ENDIF
        CALL HDPRMX(RMAT,NONZR,IPRMR,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          GO TO 270
        ENDIF
C
C         CONVERT THE NONZR ELEMENTS OF THE JACOBIAN TO COLUMN FORMAT 
C
        CALL STORCH(1,NONZR,NDIM,JSRR,IHOLD,LNRCOL,IERCH)
        IF(IERCH.NE.0) THEN
          IERNLP = -147
          GO TO 270
        ENDIF
C
C         COMPUTE LEAST SQUARES OBJECTIVE VALUE AT INITIAL POINT
C
        FBAR = DOT_PRODUCT(RESVEC(1:NRES),RESVEC(1:NRES))/TWO
C
C         COMPUTE GRADIENT VECTOR DELF = (RESVEC**T)(RMAT) 
C
        CALL MVPSPR(11,NDIM,NRES,RMAT,IRWR,JSRR,
     $                  RESVEC,DELF)
C
      ENDIF
C
C ======================================================================
C ================= END INITIAL RESIDUAL JACOBIAN ANALYSIS =============
C ======================================================================
C
C ======================================================================
C ================= INITIAL DIAGNOSTIC ANALYSIS ========================
C ======================================================================
C
C         CALL DIAGNOSTIC OUTPUT PROCEDURE
C
      LENWRK = NDIM + MAX(NDIM,MCON,NRES)
      LENIWK = 2*MAX(NONZR,NONZG,NDIM)
      CALL DGNOST( XBAR   ,XLWR   ,XUPR  ,NDIM  
     $          ,NRES   ,DELF   ,RMAT   ,IRWR ,JSRR
     $          ,NONZR  ,IROWH  ,JSTRH  ,NONZH       
     $          ,MCON   ,MTOTAL ,IPRMC ,GMAT   ,IRWG
     $          ,JSTR ,NONZG  ,NONZGT ,hold  ,LENWRK 
     $          ,IHOLD  ,LENIWK ,BIGELM,IROWBG,JCOLBG,LENBIG
     $          ,MINOPT, IERNLP ,PRMUTE)
C
      IF(IERNLP.LT.0) GO TO 270
C
      MACTIV = 0
      DO I=1,NDIM
        IF (ISTATV(I).NE.0)  MACTIV = MACTIV + 1
      ENDDO
      DO I=1,MTOTAL
        IF (ISTATC(I).NE.0)  MACTIV = MACTIV + 1
      ENDDO
      IF(MACTIV.GT.NDIM.AND.SPRSCK) THEN
        IERNLP = -130
        GO TO 270
      ENDIF
C
      OPUN = .FALSE.
      IF(SPRSCK) THEN
C
C         CHECK FOR CONFLICT BETWEEN USER AND MULTIFRONTAL FILES
C
        IF(IPUMF1.GE.0) INQUIRE(IPUMF1,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
        IF(IPUMF2.GE.0) INQUIRE(IPUMF2,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
        IF(IPUMF3.GE.0) INQUIRE(IPUMF3,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
        IF(IPUMF4.GE.0) INQUIRE(IPUMF4,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
        IF(IPUMF5.GE.0) INQUIRE(IPUMF5,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
        IF(IPUMF6.GE.0) INQUIRE(IPUMF6,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
        IF(IPUMF7.GE.0) INQUIRE(IPUMF7,OPENED=OPUN)
        IF(OPUN) THEN
          IERNLP = -137
          GO TO 270
        ENDIF
C
      ENDIF
C
      IF(IPUDRF.GT.0) INQUIRE(IPUDRF,OPENED=OPUN)
      IF(OPUN) THEN
        IERNLP = -137
        GO TO 270
      ENDIF
      IF(IPUFZF.GT.0) INQUIRE(IPUFZF,OPENED=OPUN)
      IF(OPUN) THEN
        IERNLP = -137
        GO TO 270
      ENDIF
      IF(IPUSTF.GT.0) INQUIRE(IPUSTF,OPENED=OPUN)
      IF(OPUN) THEN
        IERNLP = -137
        GO TO 270
      ENDIF
C
C         CHECK FOR VALID NEWTON OPTION FLAG
C
      IF(NRES.EQ.0) THEN
        IF(NEWTON.LT.0.OR.NEWTON.GT.1) THEN
          IERNLP = -138
          GO TO 270
        ENDIF
      ELSEIF(NRES.GT.0) THEN
        IF(NEWTON.LT.0.OR.NEWTON.GT.2) THEN
          IERNLP = -159
          GO TO 270
        ENDIF
      ENDIF
C
      IF(IHESHN.NE.0.AND.NEWTON.EQ.0) THEN
        IERNLP = +120
        CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
      ENDIF
C
C ======================================================================
C ================= END OF INITIAL DIAGNOSTIC ANALYSIS =================
C ======================================================================
C
C             WRITE THE INITIAL SUMMARY PRINT
C
      IF(IOFLAG.GE.10) THEN
C
        WRITE(IPUNLP,1002) ALGNAM
        IF(IOFLAG.GE.20) THEN
          CALL SOSVER( CODEV, LIBRV )
          WRITE(IPUNLP,1010) CODEV
        ENDIF
        WRITE(IPUNLP,1003)
        IF(NRES.GT.0) THEN
           WRITE(IPUNLP,1009) NRES
        ELSE
           WRITE(IPUNLP,1006)
        ENDIF
C
        HOLD(1:NDIM) = XUPR(1:NDIM) - XLWR(1:NDIM)
        ISTV0 = 0
        ISTV1 = 0
        ISTV2 = 0
        ISTV3 = 0
        NBNDVR = 0
        DO I=1,NDIM
          IF (ISTATV(I).EQ.0) THEN
            ISTV0 = ISTV0 + 1
          ELSEIF (ISTATV(I).EQ.1) THEN
            ISTV1 = ISTV1 + 1
          ELSEIF (ISTATV(I).EQ.2) THEN
            ISTV2 = ISTV2 + 1
          ELSEIF (ISTATV(I).EQ.3) THEN
            ISTV3 = ISTV3 + 1
          ENDIF
          IF (HOLD(I).LT.TWO*BIGBND)  NBNDVR = NBNDVR + 1 
        ENDDO
        ISTC0 = 0
        ISTC1 = 0
        ISTC2 = 0
        ISTC3 = 0
        DO I=1,MCON
          IF (ISTATC(I).EQ.0) THEN
            ISTC0 = ISTC0 + 1
          ELSEIF (ISTATC(I).EQ.1) THEN
            ISTC1 = ISTC1 + 1
          ELSEIF (ISTATC(I).EQ.2) THEN
            ISTC2 = ISTC2 + 1
          ELSEIF (ISTATC(I).EQ.3) THEN
            ISTC3 = ISTC3 + 1
          ENDIF
        ENDDO
        ISTCN3 = MCON - ISTC3
        WRITE(IPUNLP,1007) MTOTAL,NDIM,
     $      ISTC3,ISTV3,
     $      ISTCN3-MIGNOR,NBNDVR,
     $      ISTC0,ISTV0,
     $      ISTC1,ISTV1,
     $      ISTC2,ISTV2,
     $      MACTIV,NDIM-MACTIV
C
C           ADJUST CHARACTER STRING ALGOPT TO RIGHT
C
        CALL HHADJF(ALGOPT,' ',' ','R',NSHIFT,IER)
C
        IF(IOFLAG.EQ.10) THEN
          CALL INSNLP('SUMMARY')
        ELSEIF(IOFLAG.LT.20) THEN
          CALL INSNLP('OPTION')
        ELSE 
          CALL INSNLP('FULL OPTION')
        ENDIF
C
C           ADJUST CHARACTER STRING ALGOPT TO LEFT
C
        CALL HHADJF(ALGOPT,' ',' ','L',NSHIFT,IER)
C
        WRITE(IPUNLP,1003)
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- ALGORITHM CALLING SEQUENCE -------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
  240 CONTINUE
C
      IF(IMAX(4,IREVRS,1).GE.1) THEN
C
        IF(MCON.GT.0.AND.PRMUTE) THEN
C
C         RESTORE ALL EXTERNAL QUANTITIES TO INTERNAL ORDER
C         EQUALITY CONSTRAINTS FIRST
C
          CALL HJPRMX(ISTATC,MCON,IPRMC,IERP)
          CALL HDPRMX(CBAR,MCON,IPRMC,IERP)
          CALL HDPRMX(CLWR,MCON,IPRMC,IERP)
          CALL HDPRMX(CUPR,MCON,IPRMC,IERP)
          CALL HDPRMX(VECLAM,MCON,IPRMC,IERP)
          IF(IREVRS(2).GE.1.OR.(IREVRS(3).GE.1.AND.IHESHN.NE.0)
     $                     .OR.JACPRM.EQ.1) THEN
            IF(SPRSCK) THEN
              CALL HDPRMX(GMAT,NONZG,IPRMG,IERP)
            ELSE
              IRWNDX = 1 - MAXCON
              DO JCLNDX = 1,NDIM
                IRWNDX = IRWNDX + MAXCON
                CALL HDPRMX(GMAT(IRWNDX),MCON,IPRMC,IERP)
              ENDDO
            ENDIF
          ENDIF
C
        ENDIF
C
        IF(NRES.GT.0) THEN
C
          IF(IREVRS(1).GE.1) THEN
C
C           COMPUTE LEAST SQUARES OBJECTIVE VALUE
C
            FBAR = DOT_PRODUCT(RESVEC(1:NRES),RESVEC(1:NRES))/TWO
C
          ENDIF
c
          IF(IREVRS(2).GE.1) THEN
C
C           RESTORE ALL EXTERNAL QUANTITIES TO INTERNAL ORDER
C
            CALL HDPRMX(RMAT,NONZR,IPRMR,IERP)
C
C           COMPUTE GRADIENT VECTOR DELF = (RESVEC**T)(RMAT)
C
            CALL MVPSPR(11,NDIM,NRES,RMAT,IRWR,JSRR,
     $                  RESVEC,DELF)
C
          ENDIF
C
          IF(IREVRS(3).GE.2.AND.NORMAL.NE.0) THEN
C
            CALL NORMAT(RMAT,IRWR,JSRR,NONZR,NRES,
     $        HMAT,IROWH,JSTRH,NONZH,NDIM,HOLD,LNHOLD)
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
C             TURN OFF BCSLIB ERROR MESSAGES.
C
      IF(IOFLAG.LT.30) CALL HHERPT(0)
C
C             CALL THE OPTIMIZATION ALGORITHM NLPSPR
C
      CALL NLPSPR(   GMAT  ,CLWR ,CUPR ,CBAR  ,DELF ,ISTATC
     $   ,ISTATV ,XBAR ,XLWR ,XUPR ,FBAR  ,MTOTAL ,MAXCON     
     $   ,NDIM ,NONZH, NONZGT ,IRVCOM ,IREVRS ,IFERR ,IERNLP 
     $   ,ISCR ,LNISCR, CVEC   
     $   ,COLD ,PGRD ,MAXRES, RESVEC ,NRESNP ,RMAT
     $   ,IRWR ,JSRR ,NONZR,  HMAT ,IROWH ,JSTRH
     $   ,SVEC ,XOLD ,VECLAM,VECNU  ,YVEC    
     $   ,YBAR ,MINOPT ,IRWG  ,JSTR 
     $   ,IPRMC ,MCON  ,ISTO  ,HOLD  ,NHOLD  
     $   ,IHOLD ,NIHOLD, NEEDED ,PRMUTE)
C
C             TURN ON BCSLIB ERROR MESSAGES.
C
      IF(IOFLAG.LT.30.AND.IOFLAG.GT.0) CALL HHERPT(1)
C
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- FUNCTION EVALUATION SEQUENCE -----------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C
      IF(IMAX(4,IREVRS,1).GE.1) THEN
C
        IF(MCON.GT.0.AND.PRMUTE) THEN
C
C         RESTORE ALL INTERNAL QUANTITIES TO EXTERNAL ORDER
C
          CALL HJPRMY(ISTATC,MCON,IPRMC,IERP)
          CALL HDPRMY(CBAR,MCON,IPRMC,IERP)
          CALL HDPRMY(CLWR,MCON,IPRMC,IERP)
          CALL HDPRMY(CUPR,MCON,IPRMC,IERP)
          CALL HDPRMY(VECLAM,MCON,IPRMC,IERP)
          IF(IREVRS(2).GE.1.OR.(IREVRS(3).GE.1.AND.IHESHN.NE.0)
     $                     .OR.JACPRM.EQ.1) THEN
            IF(SPRSCK) THEN
              CALL HDPRMY(GMAT,NONZG,IPRMG,IERP)
            ELSE 
              IRWNDX = 1 - MAXCON
              DO JCLNDX = 1,NDIM
                IRWNDX = IRWNDX + MAXCON
                CALL HDPRMY(GMAT(IRWNDX),MCON,IPRMC,IERP)
              ENDDO
            ENDIF
          ENDIF
C
        ENDIF
C
        IF(NRES.GT.0) THEN
          IF(IREVRS(2).GE.1) THEN
C
            CALL HDPRMY(RMAT,NONZR,IPRMR,IERP)
C
          ENDIF
        ENDIF
C
      ENDIF
C
C             BRANCH TO THE RETURN UPON COMPLETION OF THIS SECTION
C
      IF(IRVCOM.NE.0) GO TO 10000
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- TERMINATION SEQUENCE (IRVCOM = 0) ------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
  270 CONTINUE
C
      IF((IERNLP.GE.0.OR.IERNLP.EQ.-131.OR.IERNLP.EQ.-132
     $    .OR.IERNLP.EQ.-133) .AND. .NOT.SUMRY ) THEN
C
        IF(MCON.GT.0.AND.PRMUTE) THEN
C
C         RESTORE ALL INTERNAL QUANTITIES TO EXTERNAL ORDER
C
          CALL HJPRMY(ISTATC,MCON,IPRMC,IERP)
          CALL HDPRMY(CBAR,MCON,IPRMC,IERP)
          CALL HDPRMY(CLWR,MCON,IPRMC,IERP)
          CALL HDPRMY(CUPR,MCON,IPRMC,IERP)
          CALL HDPRMY(VECLAM,MCON,IPRMC,IERP)
          IF(SPRSCK) THEN
            CALL HDPRMY(GMAT,NONZG,IPRMG,IERP)
          ELSE  
            IRWNDX = 1 - MAXCON
            DO JCLNDX = 1,NDIM
              IRWNDX = IRWNDX + MAXCON
              CALL HDPRMY(GMAT(IRWNDX),MCON,IPRMC,IERP)
            ENDDO
          ENDIF
C
        ENDIF
C
        IF(NRES.GT.0) THEN
C
          CALL HDPRMY(RMAT,NONZR,IPRMR,IERP)
C
        ENDIF
C
      ENDIF
C
C         MARK THE FINAL VALUES OF ISTATC AND ISTATV IF THEY
C         APPEAR TO BE INCORRECT
C
      IF(IERNLP.NE.-101.AND.IERNLP.NE.-103.AND.IERNLP.NE.-113 .AND.
     $  .NOT.SUMRY) THEN
C
        DO I = 1,MCON
          IF(ISTATC(I).EQ.1) THEN
            IF(CBAR(I).LT.(CLWR(I)-CONTOL) ) ISTATC(I) = 10
          ELSEIF(ISTATC(I).EQ.2) THEN
            IF(CBAR(I).GT.(CUPR(I)+CONTOL) ) ISTATC(I) = 20
          ELSEIF(ISTATC(I).EQ.3) THEN
            IF(CBAR(I).LT.(CLWR(I)-CONTOL) ) ISTATC(I) = 30
            IF(CBAR(I).GT.(CUPR(I)+CONTOL) ) ISTATC(I) = 30
          ENDIF
        ENDDO
C
      ENDIF
C
      IF(IERNLP.NE.-102.AND.IERNLP.NE.-120) THEN
C
        DO I = 1,NDIM
          IF(ISTATV(I).EQ.1) THEN
            IF(XBAR(I).LT.(XLWR(I)-CONTOL) ) ISTATV(I) = 10
          ELSEIF(ISTATV(I).EQ.2) THEN
            IF(XBAR(I).GT.(XUPR(I)+CONTOL) ) ISTATV(I) = 20
          ELSEIF(ISTATV(I).EQ.3) THEN
            IF(XBAR(I).LT.(XLWR(I)-CONTOL) ) ISTATV(I) = 30
            IF(XBAR(I).GT.(XUPR(I)+CONTOL) ) ISTATV(I) = 30
          ENDIF
        ENDDO
C
      ENDIF
C
      IF (IERNLP.NE.0) THEN
        IRVCOM = 0
        IF (IERNLP.GT.0) THEN
C
C         ---WARNING ERRORS
C
          MODE = 0
C
        ELSE
C
C         ---FATAL ERRORS
C
C         INPUT ARGUMENT ERROR
C
          MODE = 1
C
C         STORAGE ERROR
C
          IF(IERNLP.EQ.-127.OR.IERNLP.EQ.-128) MODE = 2
C
          IF(IERNLP.EQ.-131) THEN
C
C         INSUFFICIENT REAL STORAGE FOR SRCHFZ OR EQCMIN
C
            MODE = 2
            NEEDED = LNHOLD + NEEDED
C
C         IF THIS IS A LEAST SQUARES PROBLEM INCLUDE GRADIENT
C         OFFSET FROM SNLPLS
C
            IF(NRES.GT.0) NEEDED = NEEDED + NDIM
C
          ENDIF
C
          IF(IERNLP.EQ.-132) THEN
C
C         INSUFFICIENT INTEGER STORAGE FOR SRCHFZ OR EQCMIN
C
            MODE = 2
            NEEDED = LNIHLD + NEEDED
C
          ENDIF
C
C         FATAL PROCESS ERROR
C
          IF(IERNLP.EQ.-133.OR.IERNLP.EQ.-134) THEN
            MODE = 3
          ENDIF
C
        ENDIF
C
        IF(IERNLP.EQ.+104) THEN
C
C          RESTORE FUNCTION ERROR FLAG 
C
          IFERR = 0
C
        ENDIF
C
        IF (IOFLAG.GE.10) THEN
          IF(IERNLP.EQ.+101) THEN
             WRITE(IPUNLP,1004)
          ELSE
            WRITE(IPUNLP,1001)
          ENDIF
        ENDIF
C
        CALL HHERR(MODE,SUBNAM,IERNLP,NEEDED)
C
        DO I=1,NERMSG
          IF (IERMSG(I).EQ.IERNLP) THEN
            MESNUM = I
            EXIT
           ENDIF
        ENDDO
        IERNPE = IERMSG(MESNUM)
        IF(IMESSG.EQ.0) THEN
          WRITE(ERRMSG,FMTMSG(MESNUM))
        ELSE
          WRITE(ERRMSG,FMTMSG(MESNUM)) IMESSG
        ENDIF
        CALL NPESET(SUBNAM,IERNPE,ERRMSG,1)
C
      ELSE
C
        IF (IOFLAG.GE.10)  WRITE(IPUNLP,1004)
C
      ENDIF
C
      IF(SPRSCK) THEN
C
C         CLOSE ALL MULTIFRONTAL SCRATCH FILES
C
        IF(IPUMF1.GT.0) CLOSE(IPUMF1)
        IF(IPUMF2.GT.0) CLOSE(IPUMF2)
        IF(IPUMF3.GT.0) CLOSE(IPUMF3)
        IF(IPUMF4.GT.0) CLOSE(IPUMF4)
        IF(IPUMF5.GT.0) CLOSE(IPUMF5)
        IF(IPUMF6.GT.0) CLOSE(IPUMF6)
        IF(IPUMF7.GT.0) CLOSE(IPUMF7)
C
      ENDIF
C
C        COMPUTE THE TOTAL NLP SOLUTION TIME
C
      CALL CLKSUM(1)
C
C        DISPLAY ALGORITHM PERFORMANCE STATISTICS
C
      INSTAT(25) = NFEVAL
      INSTAT(10) = INSTAT(10) + 1
      INSTAT(30) = INSTAT(30) + 1
      IF(IERNLP.GE.0) CALL STATSM
c
C ======================================================================
c
      if(allocated(cvec)) deallocate(cvec)
      if(allocated(cold)) deallocate(cold)
      if(allocated(pgrd)) deallocate(pgrd)
      if(allocated(svec)) deallocate(svec)
      if(allocated(xold)) deallocate(xold)
      if(allocated(yvec)) deallocate(yvec)
      if(allocated(ybar)) deallocate(ybar)
c
      if(allocated(irwg)) deallocate(irwg)
      if(allocated(jstr)) deallocate(jstr)
      if(allocated(iprmg)) deallocate(iprmg)
      if(allocated(iprmc)) deallocate(iprmc)
      if(allocated(isto)) deallocate(isto)
      if(allocated(irwr)) deallocate(irwr)
      if(allocated(jsrr)) deallocate(jsrr)
      if(allocated(iprmr)) deallocate(iprmr)
      if(allocated(iscr)) deallocate(iscr)
c
C ======================================================================
C
C     ------------------------------------------------------------------
C
10000 CONTINUE
C
C         START CLOCK FOR EVERYTHING BUT NLP TIME
C
      IF(IRVCOM.NE.0) CALL CLKBEG(3)
C
      IF(IOFLAG.LT.10) IREVRS(4) = 0
C
      IF(IREVRS(1).EQ.1) INSTAT(22) = INSTAT(22) + 1
      IF(IREVRS(2).GE.1)
     $    INSTAT(23) = INSTAT(23) + 1
      IF(IREVRS(3).EQ.2) INSTAT(24) = INSTAT(24) + 1
C
      RETURN
C
C     ------------------------------------------------------------------
C             FORMAT STATEMENTS                                         
C     ------------------------------------------------------------------
C
C
 1001 FORMAT(T3,'******************************************* ABNORMAL TE
     *RMINATION ***************************************')
C
 1002 FORMAT(2X,104('*')/T3,'*',T106,'*'/T3,'*',T35,'.....OPTIMIZATION O
     $PERATOR ',A6,'.....',T106,'*')
C
 1003 FORMAT(T3,'*',T106,'*'/2X,104('*'))
C
 1004 FORMAT(T3,'*********************************************** CONVERG
     $ENCE ********************************************')
C
 1005 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING ROWS OF THE
     $RESIDUAL JACOBIAN MAY BE ZERO',T106,'*'
     $   /T3,'*',T106,'*')
C
 1006 FORMAT(T3,'*',T106,'*')
C
 1007 FORMAT(T3,'*',T106,'*'
     $  /T3,'*',T11,85('-'),T106,'*',
     $  /T3,'*',T11,'Constraints',T41,'=',I6,T54,'|'
     $  ,T59,'Variables',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,'  Equalities',T41,'=',I6,T54,'|'
     $  ,T59,'  Equalities',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,'  Inequalities',T41,'=',I6,T54,'|'
     $  ,T59,'  Bounds',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,'    Inactive',T41,'=',I6,T54,'|'
     $  ,T59,'    Free',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,'    Fixed on Lower Bound',T41,'=',I6,T54,'|'
     $  ,T59,'    Fixed on Lower Bound',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,'    Fixed on Upper Bound',T41,'=',I6,T54,'|'
     $  ,T59,'    Fixed on Upper Bound',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,85('-'),T106,'*',
     $  /T3,'*',T11,'Number of Active Constraints',T41,'=',I6,T54,'|'
     $  ,T59,'Number of Degrees of Freedom',T89,'=',I6,T106,'*'
     $  /T3,'*',T11,85('-'),T106,'*')
C
 1008 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING CONSTRAINTS 
     $WILL BE IGNORED',T106,'*'
     $   /T3,'*',T106,'*')
C
 1009 FORMAT(T3,'*',T106,'*'
     $  /T3,'*',T11,85('-'),T106,'*',/T3,'*',T106,'*',
     $  /T3,'*',T11,'Least Squares Objective Function',T59,'Number of Re
     $siduals',T89,'=',I6,T106,'*')
C
 1010 FORMAT(T3,'*',T24,A60,T106,'*')
C
      END
