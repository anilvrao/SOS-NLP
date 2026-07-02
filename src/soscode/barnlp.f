      SUBROUTINE BARNLP( IRVCOM ,IREVRS ,XBAR   ,XLWR   ,XUPR   
     $          ,ISTATV ,VECNU  ,NDIM   ,FBAR   ,RESVEC ,MAXRES 
     $          ,NRES   ,DELF   ,RMAT   ,IROWR  ,JCOLR  ,NONZR  
     $          ,HMAT   ,IROWH  ,JSTRH  ,NONZH  ,CBAR   ,CLWR     
     $          ,CUPR   ,ISTATC ,MAXCON ,MCON   ,VECLAM ,GMAT      
     $          ,IROWG  ,JCOLG  ,NONZG  ,IFERR  ,NFEVAL ,HOLD   
     $          ,NHOLD  ,IHOLD  ,NIHOLD ,NEEDED ,IERNLP )
C
C ======================================================================
C     BARNLP===>barnlp   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C *** PURPOSE...BARYER ALGORITHM INTERFACE
C
C        BARNLP CONTROLS EXECUTION OF THE INTERIOR POINT (LOG-BARRIER)
C        SPARSE NON-LINEAR PROGRAMMING ALGORITHM.  IT IS THE MAIN
C        USER INTERFACE FOR INPUT OF ALGORITHM QUANTITIES.
C
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
C                        REENTER BARNLP.  OPERATIONS TO BE PERFORMED
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
C                             OF THE ALGORITHM REQUESTING OPTIMAL PERTURBATION
C                             SIZE ADJUSTMENT
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
C     IOFSHR     I       BARRIER NLP OUTPUT FLAG                        0  
C     IOFSRC     I       SRCHFZ OUTPUT FLAG                             0
C     IPUDRF     I       I/O UNIT FOR BARRIER NLP DUMP FROM SRCHDR      0
C     ITDRQP     I       DUMP BARRIER NLP ON SRCHDR ITERATION NO.      -1
C     IPUFZF     I       I/O UNIT FOR BARRIER NLP DUMP FROM SRCHFZ      0
C     ITFZQP     I       DUMP BARRIER NLP ON SRCHFZ ITERATION NO.      -1
C     IPUSTF     I       I/O UNIT FOR BARRIER NLP DUMP FROM QPSTRT      0
C     IPUMF1     I       MULTIFRONTAL SQFILE I/O UNIT IN SRCHDR        11
C     IPUMF2     I       MULTIFRONTAL WAFIL1 I/O UNIT IN SRCHDR        12
C     IPUMF3     I       MULTIFRONTAL WAFIL2 I/O UNIT IN SRCHDR        13
C     IPUMF4     I       MULTIFRONTAL SQFILE I/O UNIT IN L.D.P.        14
C     IPUMF5     I       MULTIFRONTAL WAFIL1 I/O UNIT IN L.D.P.        15
C     IPUMF6     I       MULTIFRONTAL WAFIL2 I/O UNIT IN L.D.P.        16  
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
C          -103  (INCORRECT VALUE FOR QPOPTN) 
C          -104  (|IHESHN|.GT.3) 
C          -105  (NITMAX.LT.MAX(NITMIN,1))
C         <-106> <currently not used>
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
C          -123  (IROWH(*).LE.0) .OR. (IROWH(*).GT.NDIM)
C          -124  (XUPR.LT.XLWR)
C          -125  (XUPR.EQ.XLWR).AND.(ISTATV.NE.3) .OR.
C                (XUPR.NE.XLWR).AND.(ISTATV.EQ.3)
C          -126  (XUPR.NE.XLWR).AND.((|XUPR-XLWR|.LT.CONTOL)
C          -127  REAL HOLD ARRAY TOO SMALL
C          -128  INTEGER HOLD ARRAY TOO SMALL
C          -129  FUNCTION ERROR AT INITIAL POINT OR DURING GRADIENT EVAL.
C          -130  BIGCON .LT. CONTOL
C          -131  INSUFFICIENT REAL STORAGE DETECTED 
C          -132  INSUFFICIENT INTEGER REAL STORAGE DETECTED 
C          -133  SINGULAR JACOBIAN ON SUCCESSIVE ITERATIONS
C          -134  INSNLP INPUT ERROR
C          -135  (IFERR.LT.0) .OR. (IFERR.GT.1)
C          -136  FEATOL.LT.CONTOL
C          -137  CONFLICT BETWEEN USER AND MULTIFRONTAL FILE NUMBER
C          -138  NEWTON .NE. 0, 1, OR 2
C          -139  (NRES.LT.0)
C         <-140> <currently not used>
C          -141  (NONZR.LE.0) .OR. (NONZR .GT. NDIM*NRES)
C          -142  (JCOLR(*).LE.0) .OR. (JCOLR(*).GT.NDIM)
C          -143  (IROWR(*).LE.0) .OR. (IROWR(*).GT.NRES)
C          -144  A ROW OF THE RESIDUAL JACOBIAN IS IDENTICALLY ZERO
C          -145  IROWH(JSTRH(I)).NE.I FOR SOME I
C          -146  PMULWR < ZEROMN
C          -147  UNEXPECTED ERROR
C          -148  (|IRVCOM|.GT.1) 
C          -149  PTHTOL < ZEROOT
C          -150  RHOLWR < ZEROMN 
C          -151  CONSTRAINT VIOLATES ITS BOUNDS AND CANNOT BE CHANGED
C          -152  IMAXMU < 1
C          -153  I/O ERROR (INSUFFICIENT DISK SPACE)
C          -154  TOLKTC < 1
C          -155  TOLPVT < 0 OR TOLPVT > .5
C          -156  IRELAX < 0 OR IRELAX > 2
C          -157  |MUCALC| < 1 OR |MUCALC| > 4
C          -158  MXQPIT < 1 
C         <-159> <currently not used>
C          -160  TIMSOS .LE. 0
C
C          --------------
C                WARNING ERRORS---EVALUATION BY USER SUGGESTED
C          --------------
C
C          +101  WEAK SOLUTION FOUND (RELAXATION REQUIRED AND/OR MULTIPLIERS NEAR ZERO)
C          +102  (MEQUAL.EQ.NFREE).AND.ALGOPT.NE.F___
C          +103  MAX. NO. OF CONSECUTIVE FUNCTION ERRORS
C          +104  MAXIMUM NUMBER OF FUNCTION EVALUATIONS
C          +105  SMALL STEP TERMINATION;  SUBOPTIMAL FEASIBLE
C                POINT FOUND
C          +106  MAXIMUM NUMBER OF ITERATIONS
C         <+107> <currently not used>
C          +108  FEASIBLE POINT NOT FOUND
C          +109  MAX. NO. OF INTERVAL HALVES IN LINE SEARCH
C          +110  EITHER MAX HESSIAN DIAGONAL OR VIOLATED SLOPE CONDITION;
C                SUBOPTIMAL FEASIBLE POINT FOUND  
C          +111  PROJ. GRAD. CALCULATION FAILED
C          +112  QPSTRT FAILED TO COMPUTE MULTIPLIERS
C          +113  SUBOPTIMAL FEASIBLE POINT FOUND
C          +114  BARRIER NLP FAILED WITH UNEXPECTED ERROR
C          +115  CONTOL.GT.OBJTOL
C          +116  UPHILL DIRECTION DETECTED IN LINE SEARCH
C          +117  REDUCED OBJECTIVE FUNCTION IS LINEAR
C          +118  CONSTRAINTS IGNORED
C          +119  TERMINATE AFTER DIAGNOSTIC LINE SEARCH
C          +120  RECURSIVE HESSIAN ESTIMATE AND NEWTON.NE.0
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
C        BARNLP IS THE ALGORITHM INTERFACE FOR SPRBAR.  
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
C                  THE INITIAL SUMMARY OF BARYER DATA IS OUTPUT.
C
C        ALGORITHM CALLING SEQUENCE
C
C                  CALL THE OPTIMIZATION ALGORITHM, BARYER.
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
C ======================================================================
c
      double precision, allocatable, dimension(:) :: ALWI 
      double precision, allocatable, dimension(:) :: AUPI 
      double precision, allocatable, dimension(:) :: BMAT 
      double precision, allocatable, dimension(:) :: BVEC 
      double precision, allocatable, dimension(:) :: CMAT 
      double precision, allocatable, dimension(:) :: CVEC 
      double precision, allocatable, dimension(:) :: ETAV 
      double precision, allocatable, dimension(:) :: GVEC 
      double precision, allocatable, dimension(:) :: RSKR 
      double precision, allocatable, dimension(:) :: VLAM 
      double precision, allocatable, dimension(:) :: WMAT 
      double precision, allocatable, dimension(:) :: XLWI 
      double precision, allocatable, dimension(:) :: XUPI 
      double precision, allocatable, dimension(:) :: YVEC 
c
      integer, allocatable, dimension(:) :: IBND
      integer, allocatable, dimension(:) :: IPRMC
      integer, allocatable, dimension(:) :: IPRMG
      integer, allocatable, dimension(:) :: IPRMH 
      integer, allocatable, dimension(:) :: IPRMR 
      integer, allocatable, dimension(:) :: IPRMX
      integer, allocatable, dimension(:) :: IRWB 
      integer, allocatable, dimension(:) :: IRWC
      integer, allocatable, dimension(:) :: IRWH
      integer, allocatable, dimension(:) :: IRWR
      integer, allocatable, dimension(:) :: IRWW
      integer, allocatable, dimension(:) :: ISKR 
      integer, allocatable, dimension(:) :: JCLB 
      integer, allocatable, dimension(:) :: JCLC
      integer, allocatable, dimension(:) :: JSRR
      integer, allocatable, dimension(:) :: JSTH
      integer, allocatable, dimension(:) :: JSTW 
c
C ======================================================================
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
      DIMENSION  IHOLD(NIHOLD)
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
C
      LOGICAL  OPUN,SUMRY,RELAX,FZMODE
C
** error messages:
      CHARACTER(LEN=8)  SUBNAM
      PARAMETER  (NERMSG=78)
      DIMENSION  IERMSG(NERMSG)
      CHARACTER(LEN=100)  FMTMSG(NERMSG),ERRMSG(1)
C
      DATA  (IERMSG(I),I=1,NERMSG) /
     & -101, -102, -103, -104, -105, -106, -107, -108, -109, -110,
     & -111, -112, -113, -114, -115, -116, -117, -118, -119, -120,
     & -121, -122, -123, -124, -125, -126, -127, -128, -129, -130,
     & -131, -132, -133, -134, -135, -136, -137, -138, -139, -140,
     & -141, -142, -143, -144, -145, -146, -147, -148, -149, -150,
     & -151, -152, -153, -154, -155, -156, -157, -158, +101, +102, 
     & +103, +104, +105, +106, +107, +108, +109, +110, +111, +112,  
     & +113, +114, +115, +116, +117, +118, +119, +120/
C
      DATA  (FMTMSG(I),I=1,10)/
     &'(T1,"(MCON.LT.0)")',
     &'(T1,"(NDIM.LT.1)")',
     &'(T1,"(INCORRECT VALUE FOR QPOPTN)")', 
     &'(T1,"(|IHESHN|.GT.3)")', 
     &'(T1,"(NITMAX.LT.1)")',
C            <CURRENTLY NOT USED>
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
     &'(T1,"(IROWH(*).LE.0) .OR. (IROWH(*).GT.NDIM)")',
     &'(T1,"(XUPR.LT.XLWR) Variable:",I6)',
     &'(T1,"((XUPR.EQ.XLWR).AND.(ISTATV.NE.3)).OR.((XUPR.NE.XLWR).AND.(I
     &STATV.EQ.3)) Variable:",I6)',
     &'(T1,"(XUPR.NE.XLWR).AND.(|XUPR-XLWR|.LT.CONTOL) Variable:",I6)',
     &'(T1,"REAL HOLD ARRAY TOO SMALL")',
     &'(T1,"INTEGER HOLD ARRAY TOO SMALL")',
     &'(T1,"FUNCTION ERROR AT INITIAL POINT OR DURING GRADIENT EVAL.")',
     &'(T1,"BIGCON.LT.CONTOL")'/
C
      DATA  (FMTMSG(I),I=31,40)/
     &'(T1,"INSUFFICIENT REAL STORAGE")',
     &'(T1,"INSUFFICIENT INTEGER STORAGE")',
     &'(T1,"SINGULAR JACOBIAN ON SUCCESSIVE ITERATIONS")',
     &'(T1,"INSNLP INPUT ERROR")',
     &'(T1,"(IFERR.LT.0) .OR. (IFERR.GT.1)")',
     &'(T1,"FEATOL.LT.CONTOL")',
     &'(T1,"CONFLICT BETWEEN USER AND MULTIFRONTAL FILE NUMBER")',
     &'(T1,"(NEWTON.LT.0) .OR. (NEWTON.GT.2)")',
     &'(T1,"(NRES.LT.0)")',
C            <CURRENTLY NOT USED>
     &'(T1,"(NRES.GT.MAXRES)")'/
C
      DATA  (FMTMSG(I),I=41,50)/
     &'(T1,"(NONZR.LE.0) .OR. (NONZR .GT. NDIM*NRES)")',
     &'(T1,"(JCOLR(*).LE.0) .OR. (JCOLR(*).GT.NDIM)")',
     &'(T1,"(IROWR(*).LE.0) .OR. (IROWR(*).GT.NRES)")',
     &'(T1,"A ROW OF THE RESIDUAL JACOBIAN IS IDENTICALLY ZERO")',
     &'(T1,"IROWH(JSTRH(I)).NE.I FOR SOME I")',
     &'(T1,"PMULWR < ZEROMN")',
     &'(T1,"UNEXPECTED ERROR")',
     &'(T1,"|IRVCOM| .GT. 1")',
     &'(T1,"PTHTOL < ZEROOT")',
     &'(T1,"RHOLWR < ZEROMN")'/
C
      DATA  (FMTMSG(I),I=51,60)/
     &'(T1,"INCONSISTENT Constraint:",I6)',
     &'(T1,"IMAXMU < 1")',
     &'(T1,"I/O ERROR (INSUFFICIENT DISK SPACE)")',
     &'(T1,"TOLKTC < 1")',
     &'(T1,"TOLPVT < 0 OR TOLPVT > .5")',
     &'(T1,"(IRELAX < 0) .OR. (IRELAX > 2)")',
     &'(T1,"(|MUCALC| < 1) .OR. (|MUCALC| > 4)")',
     &'(T1,"MXQPIT < 1")',
     &'(T1,"WEAK SOLUTION FOUND (RELAXATION REQUIRED AND/OR MULTIPLIERS 
     &NEAR ZERO)")',
     &'(T1,"(MEQUAL.EQ.NFREE).AND.ALGOPT.NE.F___")'/
C
      DATA  (FMTMSG(I),I=61,70)/
     &'(T1,"MAX. NO. OF CONSECUTIVE FUNCTION ERRORS")',
     &'(T1,"MAXIMUM NUMBER OF FUNCTION EVALUATIONS")',
     &'(T1,"SMALL STEP TERMINATION; SUBOPTIMAL FEASIBLE POINT FOUND")',
     &'(T1,"MAXIMUM NUMBER OF ITERATIONS")',
     &'(T1,"MAX. NO. OF ITER. IN FEASIBILITY PHASE")',
     &'(T1,"FEASIBLE POINT NOT FOUND")',
     &'(T1,"MAX. NO. OF INTERVAL HALVES IN LINE SEARCH")',
     &'(T1,"EITHER MAX HESSIAN DIAGONAL OR VIOLATED SLOPE CONDITION; SUB
     &OPTIMAL FEASIBLE POINT FOUND")',
     &'(T1,"PROJ. GRAD. CALCULATION FAILED")',
     &'(T1,"LAGRANGE MULTIPLIERS NOT COMPUTED; DEGENERATE CONSTRAINTS")'
     &/
C
      DATA  (FMTMSG(I),I=71,NERMSG)/
     &'(T1,"SUBOPTIMAL FEASIBLE POINT FOUND")',
     &'(T1,"BARRIER NLP FAILED WITH UNEXPECTED ERROR")',
     &'(T1,"CONTOL.GT.OBJTOL")',
     &'(T1,"UPHILL DIRECTION DETECTED IN LINE SEARCH")',
     &'(T1,"REDUCED OBJECTIVE FUNCTION IS LINEAR")',
     &'(T1,"CONSTRAINTS IGNORED")',
     &'(T1,"TERMINATE AFTER DIAGNOSTIC LINE SEARCH")',
     &'(T1,"RECURSIVE HESSIAN ESTIMATE AND NEWTON.NE.0")'/
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
C         INITIALIZE THE TOTAL TIME IN BRSRCH ON CLOCK 4
C
      CALL CLKSET(4)
C
C         INITIALIZE THE MAX TIME IN BRSRCH ON CLOCK 5
C
      CALL CLKSET(-5)
C
C         INITIALIZE THE TOTAL TIME FOR K-T FACTORIZATIONS ON CLOCK 6
C
      CALL CLKSET(6)
C
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
C     7               
C     8                    SOCX,DENSTS
C     9                    SOCX,DENSTS
C     10              
C     11              SNLP
C     :               :
C     31              SNLP
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
        SUBNAM='SBRNLP  '
      ELSE
        SUBNAM='SBRLSQ  '
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
        CALL INSNLP('BARRIER DEFAULT')
      ENDIF
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
      CALL HHERRS(1,IPUNLP)
C
C             CHECK BARYER INPUT QUANTITIES.
C
      IF(MCON.LT.0) THEN
        IERNLP = -101
        GO TO 270
      ENDIF
      IF(NDIM.LT.1) THEN
        IERNLP = -102
        GO TO 270
      ENDIF
      IF(QPOPTN.NE.'SPARSE') THEN
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
      IF(IMAXMU.LT.1) THEN
        IERNLP = -152
        GO TO 270
      ENDIF
      IF(IRELAX.LT.0.OR.IRELAX.GT.2) THEN
        IERNLP = -156
        GO TO 270
      ENDIF
      IF(ABS(MUCALC).LT.1.OR.ABS(MUCALC).GT.4) THEN
        IERNLP = -157
        GO TO 270
      ENDIF
      IF(MXQPIT.LT.1) THEN
        IERNLP = -158
        GO TO 270
      ENDIF
C
C             TURN OFF BCSLIB ERROR MESSAGES.
C
      IF(IOFLAG.EQ.0) CALL HHERPT(0)
C
C         CHECK ALGORITHM OPTION 
C
      select case(ALGOPT(1:4))
      case('FM  ','M   ','F   ')
      case default
        IERNLP = -109
        GO TO 270
      end select
C
C         SET RELAX BASED ON INPUT (IRELAX)
C
      select case(irelax)
      case(0,1)
        RELAX = .FALSE.
      case(2)
        RELAX = .TRUE.
      case default
        IERNLP = -156
        GO TO 270
      end select
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
      IF(OBJTOL.LE.ONEEP1*ZEROMN) THEN
        IERNLP = -110
        GO TO 270
      ENDIF
      IF((PGDTOL.LT.ZEROOT).OR.(PGDTOL.GT.ONEEM2)) THEN
        IERNLP = -111
        GO TO 270
      ENDIF
      IF(PMULWR.LT.ZEROMN) THEN
        IERNLP = -146
        GO TO 270
      ENDIF
      IF(PTHTOL.LT.ZEROOT) THEN
        IERNLP = -149
        GO TO 270
      ENDIF
      IF(RHOLWR.LT.ZEROMN) THEN
        IERNLP = -150
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
        IF(BIGCON.LT.CONTOL) THEN
          IERNLP = -130
          GO TO 270
        ENDIF
C
        IF(FEATOL.LT.CONTOL) THEN
          IERNLP = -136
          GO TO 270
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
        JCOLMX = NDIM
        JCOLMN = 1
        DO I=1,NONZG
          JCOLMX = MAX(JCOLMX,JCOLG(I))
          JCOLMN = MIN(JCOLMN,JCOLG(I))
        ENDDO
        IF(JCOLMN.LE.0 .OR. JCOLMX.GT.NDIM) THEN
          IERNLP = -115
          GO TO 270
        ENDIF
C
        IROWMX = MCON
        IROWMN = 1
        DO I=1,NONZG
          IROWMX = MAX(IROWMX,IROWG(I))
          IROWMN = MIN(IROWMN,IROWG(I))
        ENDDO
        IF(IROWMN.LE.0 .OR. IROWMX.GT.MCON) THEN
          IERNLP = -116
          GO TO 270
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
        enddo
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
      XNDIM = NDIM
      NZHDIM = JSTRH(NDIM+1)-1
      XNZHDM = NZHDIM
      XNDNSH = XNDIM*(XNDIM+ONE)/TWO
C
      IF(JSTRH(1).NE.1) THEN
        IMESSG = 1
        IERNLP = -121
        GO TO 270
      ENDIF
C
      IF(NZHDIM.LT.NDIM.OR.XNZHDM.GT.XNDNSH.OR.NZHDIM.NE.NONZH) THEN
        IERNLP = -122
        GO TO 270
      ENDIF
C
      IROWMX = NDIM
      IROWMN = 1
      DO I=1,NONZH
        IROWMX = MAX(IROWMX,IROWH(I))
        IROWMN = MIN(IROWMN,IROWH(I))
      ENDDO
      IF(IROWMN.LE.0 .OR. IROWMX.GT.NDIM) THEN
        IERNLP = -123
        GO TO 270
      ENDIF
C
      DO I = 1,NDIM
C
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
        IF(ISTATV(I).NE.3) THEN
C
C         MAKE THE INITIAL GUESS FEASIBLE WITH RESPECT TO THE BOUNDS
C
          FEASTP = MAX(FEATOL,FEATOL*ABS(XBAR(I)))
          FEASVL = MIN(FEASTP,POINT5*ABS(XUPR(I)-XLWR(I)))
          XBAR(I) = MAX(XLWR(I)+FEASVL,MIN(XBAR(I),XUPR(I)-FEASVL))
C
        ENDIF
C
      enddo
C
C        CHECK RESIDUAL INPUTS FOR CONSISTENCY
C
      IF(NONZR.LE.0) THEN
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
        JCOLMX = NDIM
        JCOLMN = 1
        DO I=1,NONZR
          JCOLMX = MAX(JCOLMX,JCOLR(I))
          JCOLMN = MIN(JCOLMN,JCOLR(I))
        ENDDO
        IF(JCOLMN.LE.0 .OR. JCOLMX.GT.NDIM) THEN
          IERNLP = -142
          GO TO 270
        ENDIF
C
        IROWMX = NRES
        IROWMN = 1
        DO I=1,NONZR
          IROWMX = MAX(IROWMX,IROWR(I))
          IROWMN = MIN(IROWMN,IROWR(I))
        ENDDO
        IF(IROWMN.LE.0 .OR. IROWMX.GT.NRES) THEN
          IERNLP = -143
          GO TO 270
        ENDIF
C  
      ENDIF
C
      IF(NIHOLD.LT.NONZG+MCON) THEN
        NEEDED = NONZG + MCON + 1
        IERNLP = -128
        GO TO 270
      ENDIF
C
C         TEMPORARILY SET LCIWRK = NONZG+1 FOR DIAGNOSTICS
C
      LCIWRK = NONZG + 1
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
            IHOLD(LCIWRK+NINFBD-1) = I
            ISTATC(I) = 4
          ENDIF
        enddo
C
        IF(NINFBD.GT.0) THEN
          IERNLP = +118
          CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
          IF(IOFLAG.GE.10) THEN
            WRITE(IPUNLP,1008)
            CALL INTOUT(IHOLD(LCIWRK),NINFBD,IPUNLP)
          ENDIF
        ENDIF
C
C         CHECK FOR ZERO ROWS IN THE JACOBIAN MATRIX
C
C         COUNT THE NUMBER OF NONZEROS IN ROW I AND STORE
C         THE RESULT IN IHOLD(I)
C
        ihold(1:mcon) = 0
C
        DO K = 1,NONZG
          IR = IROWG(K)
          JC = JCOLG(K)
          IF(ISTATV(JC).NE.3.AND.ISTATC(IR).NE.4) THEN
            IHOLD(IR) = IHOLD(IR) + 1
          ENDIF
        enddo
C
C         CHECK THAT EACH ROW HAS AT LEAST ONE ELEMENT AND SAVE THE
C         ROW NUMBER OF THE ROWS WITH NO ELEMENTS
C
        DO I = 1,MCON
          IF(IHOLD(I).LE.0.AND.ISTATC(I).NE.4) THEN
            NZERO = NZERO + 1
            IHOLD(LCIWRK+NZERO-1) = I
            ISTATC(I) = 4
          ENDIF
        enddo
C
        IF(NZERO.GT.0) THEN
          IERNLP = +118
          CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
          IF(IOFLAG.GE.10) THEN
            WRITE(IPUNLP,1008)
            CALL INTOUT(IHOLD(LCIWRK),NZERO,IPUNLP)
          ENDIF
        ENDIF
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C             COMPUTE DIMENSION PARAMETERS
C
C        -----NUMBER OF FIXED VARIABLES
      NFIX = 0
C        -----NUMBER OF FREE VARIABLES
      NFREE = 0
C        -----NUMBER OF REAL LOWER BOUNDS
      NXLWR = 0
C        -----NUMBER OF REAL UPPER BOUNDS
      NXUPR = 0
C
      DO I = 1,NDIM
        IF(ISTATV(I).EQ.3) THEN
          NFIX = NFIX + 1
        ELSE
          NFREE = NFREE + 1
          IF(XLWR(I).GT.-BIGBND) NXLWR = NXLWR + 1
          IF(XUPR(I).LT.BIGBND) NXUPR = NXUPR + 1
        ENDIF
      enddo
C
C        -----NUMBER OF EQUALITY CONSTRAINTS
      MEQUAL = 0 
C        -----NUMBER OF IGNORED CONSTRAINTS
      MIGNOR = 0
C        -----NUMBER OF INEQUALITY CONSTRAINTS
      MINEQL = 0 
C        -----NUMBER OF REAL LOWER CONSTRAINT BOUNDS
      NALWR = 0
C        -----NUMBER OF REAL UPPER CONSTRAINT BOUNDS
      NAUPR = 0
      DO I = 1,MCON
        IF(ISTATC(I).EQ.4) THEN
          MIGNOR = MIGNOR + 1
        ELSEIF(ISTATC(I).EQ.3) THEN
          MEQUAL = MEQUAL + 1
        ELSE
          MINEQL = MINEQL + 1
          IF(CLWR(I).GT.-BIGBND) NALWR = NALWR + 1
          IF(CUPR(I).LT.BIGBND) NAUPR = NAUPR + 1
        ENDIF
      enddo
C
C        -----NUMBER OF JACOBIAN NONZEROS REMOVED BY FIXED VARIABLES
      NZRMVC = 0
      IF(MCON.GT.0) THEN
        DO K = 1,NONZG
          JC = JCOLG(K)
          IF(ISTATV(JC).EQ.3) NZRMVC = NZRMVC + 1
        enddo
      ENDIF
C
C        -----NUMBER OF HESSIAN NONZEROS REMOVED BY FIXED VARIABLES
      NZRMVH = 0
      DO J = 1,NDIM
        DO I = JSTRH(J),JSTRH(J+1)-1
          IR = IROWH(I)
          IF(ISTATV(IR).EQ.3.OR.ISTATV(J).EQ.3) NZRMVH = NZRMVH + 1
        enddo
      enddo
C        -----NUMBER OF EQUALITY CONSTRAINTS PLUS FIXED VARIABLES
      NEQL = MEQUAL + NFIX
C        -----NUMBER OF EQUALITY CONSTRAINTS 
      MSUBE = MCON - MIGNOR
C        -----NUMBER OF SLACK VARIABLES (RELAXATION) 
      NSLKR = MINEQL + 2*MEQUAL + NALWR + NAUPR
C        -----NUMBER OF VARIABLES (RELAXATION) 
      NVARR = NFREE + NSLKR
C        -----NUMBER OF NONZEROS IN CMAT MATRIX (RELAXATION) 
      NONZCR = NONZG - NZRMVC + MINEQL + 2*MEQUAL
      NONZCR = MAX(NONZCR,NONZG,1)
C        -----NUMBER OF BOUND INEQUALITY CONSTRAINTS (RELAXATION) 
      MSUBBR = NXLWR + NXUPR + NALWR + NAUPR + 2*MEQUAL + NALWR + NAUPR
C        -----NUMBER OF NONZEROS IN BMAT (RELAXATION) 
      NONZBR = NXLWR + NXUPR + 3*(NALWR + NAUPR) + 2*MEQUAL 
      NONZBR = MAX(NONZBR,1)
C        -----NUMBER OF NONZEROS IN WMAT MATRIX (RELAXATION) 
      NONZWR = NONZH + MAX(0,NSLKR-NZRMVH)
C        -----NUMBER OF SLACK VARIABLES 
      NSLK = MINEQL
C        -----NUMBER OF VARIABLES 
      NVARN = NFREE + NSLK
C        -----NUMBER OF NONZEROS IN CMAT MATRIX 
      NONZCN = NONZG - NZRMVC + NSLK
      NONZCN = MAX(NONZCN,1)
C        -----NUMBER OF BOUND INEQUALITY CONSTRAINTS 
      MSUBBN = NXLWR + NXUPR + NALWR + NAUPR
C        -----NUMBER OF NONZEROS IN BMAT 
      NONZBN = NXLWR + NXUPR + NALWR + NAUPR
      NONZBN = MAX(NONZBN,1)
C        -----NUMBER OF NONZEROS IN WMAT MATRIX 
      NONZWN = NONZH - NZRMVH + NSLK
      NONZWN = MAX(NONZWN,1)
C
C        -----MAX BOUND DIMENSION
      MAXBND = MAX(MSUBBR,1)
C        -----MAX VARIABLE DIMENSION
      NVARMX = MAX(NVARR,NDIM)
C        -----LENGTH OF REAL SCRATCH ARRAY
      LNRSKR = MAX(NONZH,NONZG,NONZBR,NONZR,MAX(MCON,NRES,NDIM),
     $         2*(NDIM + MAXCON))
C        -----LENGTH OF INTEGER SCRATCH ARRAY
      LNISKR = MAX(2*(NDIM+MCON),NONZBR)
C        -----MAX DIMENSION FOR HESSIAN MATRIX
      LNHCOL = MAX(NONZH,NDIM+1)
C        -----MAX DIMENSION FOR JACOBIAN MATRIX
      LNJCOL = MAX(NONZG,NDIM+1)
C        -----MAX DIMENSION FOR RESIDUAL JACOBIAN MATRIX
      LNRCOL = MAX(NONZR,NDIM+1)
C        -----MAX DIMENSION FOR RESIDUAL JACOBIAN MATRIX COLUMN ARRAY
      LNJSRR = MAX(NONZR,NVARR+1)
C        -----MAX FILTER SIZE
      MXFLTR = IMAXMU*IT1MAX
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C             STORAGE ALLOCATION FOR THE REAL ARRAY
C             ALLOCATION BASED ON RELAXATION MODE
C
C ======================================================================
c
      if(allocated(alwi)) deallocate(alwi)
      if(allocated(aupi)) deallocate(aupi)
      if(allocated(bmat)) deallocate(bmat)
      if(allocated(bvec)) deallocate(bvec)
      if(allocated(cmat)) deallocate(cmat)
      if(allocated(cvec)) deallocate(cvec)
      if(allocated(etav)) deallocate(etav)
      if(allocated(gvec)) deallocate(gvec)
      if(allocated(rskr)) deallocate(rskr)
      if(allocated(vlam)) deallocate(vlam)
      if(allocated(wmat)) deallocate(wmat)
      if(allocated(xlwi)) deallocate(xlwi)
      if(allocated(xupi)) deallocate(xupi)
      if(allocated(yvec)) deallocate(yvec)
      if(allocated(ibnd)) deallocate(ibnd)
      if(allocated(iprmc)) deallocate(iprmc)
      if(allocated(iprmg)) deallocate(iprmg)
      if(allocated(iprmh)) deallocate(iprmh)
      if(allocated(iprmr)) deallocate(iprmr)
      if(allocated(iprmx)) deallocate(iprmx)
      if(allocated(irwb)) deallocate(irwb)
      if(allocated(irwc)) deallocate(irwc)
      if(allocated(irwh)) deallocate(irwh)
      if(allocated(irwr)) deallocate(irwr)
      if(allocated(irww)) deallocate(irww)
      if(allocated(iskr)) deallocate(iskr)
      if(allocated(jclb)) deallocate(jclb)
      if(allocated(jclc)) deallocate(jclc)
      if(allocated(jsrr)) deallocate(jsrr)
      if(allocated(jsth)) deallocate(jsth)
      if(allocated(jstw)) deallocate(jstw)

      allocate(alwi(1:maxcon))
      allocate(aupi(1:maxcon))
      allocate(bmat(1:nonzbr))
      allocate(bvec(1:maxbnd))
      allocate(cmat(1:nonzcr))
      allocate(cvec(1:maxcon))
      allocate(etav(1:maxcon))
      allocate(gvec(1:nvarmx))
      allocate(rskr(1:lnrskr))
      allocate(vlam(1:maxbnd))
      allocate(wmat(1:nonzwr))
      allocate(xlwi(1:ndim))
      allocate(xupi(1:ndim))
      allocate(yvec(1:nvarmx))
c
c ======================================================================
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C             STORAGE ALLOCATION FOR THE INTEGER ARRAY
c
      allocate(IBND(1:MAXBND))
      allocate(IPRMC(1:maxcon))
      allocate(IPRMG(1:NONZG))
      allocate(IPRMH(1:LNHCOL))
      allocate(IPRMR(1:LNRCOL))
      allocate(IPRMX(1:NDIM))
      allocate(IRWB(1:NONZBR))
      allocate(IRWC(1:NONZCR))
      allocate(IRWH(1:LNHCOL))
      allocate(IRWR(1:NONZR))
      allocate(IRWW(1:NONZWR))
      allocate(ISKR(1:LNISKR))
      allocate(JCLB(1:NVARR+1))
      allocate(JCLC(1:NVARR+1))
      allocate(JSRR(1:LNJSRR))
      allocate(JSTH(1:NDIM+1))
      allocate(JSTW(1:NVARR+1))
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      IF(NRES.GT.0) THEN
C
C         CHECK FOR ZERO ROWS IN THE RESIDUAL JACOBIAN MATRIX
C
C         COUNT THE NUMBER OF NONZEROS IN ROW I AND STORE
C         THE RESULT IN IHOLD(I)
C
        ihold(1:nres) = 0
C  
        DO K = 1,NONZR
          IR = IROWR(K)
          IHOLD(IR) = IHOLD(IR) + 1
        enddo
C
        NZROW = 0
        DO K = 1,NRES
          IF(IHOLD(K).LE.0) THEN
            NZROW = NZROW + 1
            IHOLD(NZROW) = K
          ENDIF
        enddo
        IF(NZROW.GT.0) THEN
          IERNLP = -144
          IF(IOFLAG.GE.10) WRITE(IPUNLP,1005) 
          CALL INTOUT(IHOLD,NZROW,IPUNLP)
          GO TO 270
        ENDIF
C
      ENDIF
C
      irwr(1:nonzr) = irowr(1:nonzr)
      IF(NRES.NE.0) jsrr(1:nonzr) = jcolr(1:nonzr)
C
      FZMODE = ALGOPT(1:1).EQ.'F'.OR.(NFREE.EQ.MEQUAL.AND.IRELAX.EQ.0)
      IF(ALGOPT(1:1).EQ.'F'.AND.NFREE.EQ.MEQUAL.AND.IRELAX.EQ.2) 
     $   FZMODE = .FALSE.
C
C         CHECK THAT ALGORITHM OPTION IS CONSISTENT FOR SEARCHS
C
      IF(NFREE.EQ.MEQUAL.AND.ALGOPT(1:1).NE.'F') THEN
        FZMODE = .TRUE.
        ALGOPT = 'F   '
        IERNLP = +102
        CALL HHERR(0,SUBNAM,IERNLP,NEEDED)
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
C         INSTAT(10)= NO. OF ELEMENTS NEEDED IN REAL HOLD ARRAY
C         INSTAT(11)= NO. OF KKT SOLUTION FAILURES
C         INSTAT(12)= NO. OF KKT SYSTEM CALLS
C         INSTAT(13)= NLP ITERATION NUMBER FOR MOST EXPENSIVE KKT SOLUTION
C         INSTAT(22)= NO. OF FCALLS
C         INSTAT(23)= NO. OF GCALLS
C         INSTAT(24)= NO. OF HCALLS
C         INSTAT(25)= TOTAL NO. OF FUNCTION EVAL.
C         INSTAT(28)= SET TO ONE FOR OUT-OF-CORE, ZERO OTHERWISE
C         INSTAT(30)= NO. OF ELEMENTS NEEDED IN INTEGER HOLD ARRAY
C
C
      instat(1:30) = 0
      rlstat(1:20) = zero
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
      IF(MUCALC.NE.4) THEN
C
C         INITIALIZE MULTIPLIER VECTOR
C
        veclam(1:mcon) = zero
C
C         INITIALIZE BOUND MULTIPLIER VECTOR
C
        vecnu(1:ndim) = zero
C
      ENDIF
C
C             SET PRINT FLAG IF REQUIRED.
C
      IREVRS(4) = 0
      IF (IOFLAG.GE.10)  IREVRS(4) = 1
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
        HUGEF = ZERO
        DO I=1,NRES
          HUGEF = MAX(HUGEF,ABS(RESVEC(I)))
        ENDDO
      ENDIF
      DO I=1,MCON
        HUGEF = MAX(HUGEF,ABS(CBAR(I)))
      ENDDO
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
          IOPT = 1
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
        IOPT = 3
        CALL SRG228(IOPT,NDIM,NDIM,NONZH,HMAT,IROWH,JSTRH,
     &                ITITLE,IPUNLP)
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
      IF(.NOT. SUMRY) GO TO 260
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
C ================= INITIAL RESIDUAL JACOBIAN ANALYSIS =================
C ======================================================================
C
      IF(NRES.GT.0) THEN
C
C         COMPUTE LEAST SQUARES OBJECTIVE VALUE AT INITIAL POINT
C
        FBAR = DOT_PRODUCT(RESVEC(1:NRES),RESVEC(1:NRES))/TWO
C
C         COMPUTE GRADIENT VECTOR DELF = (RMAT**T)(RESVEC) 
C
        CALL MVPSRC(11,NDIM,NRES,NONZR,RMAT,IROWR,JCOLR,RESVEC,DELF)
C
        CALL SHUFLR( ISTATV ,NDIM   ,NVARR ,NRES  ,IPRMX ,RMAT
     $          ,IRWR ,JCOLR  ,NONZR  ,JSRR   
     $          ,IPRMR ,ISKR ,RSKR  ,IERNLP )

        IF(IERNLP.LT.0) GO TO 270
C
      ENDIF
C
C ======================================================================
C ================= END INITIAL RESIDUAL JACOBIAN ANALYSIS =============
C ======================================================================
C
C ======================================================================
C ================= INITIAL JACOBIAN ANALYSIS ==========================
C ======================================================================
C
C
C         COMPUTE THE LARGEST ELEMENTS IN THE JACOBIAN
C
C
      IF(LENBIG.GT.NIHOLD) THEN
        NEEDED = LENBIG + 2
        IERNLP = -128
        GO TO 270
      ENDIF
c
      BIGELM(1) = ZERO
      IF(MCON.GT.0) CALL BIGGLM(GMAT,IROWG,
     $    JCOLG,NONZG,ISTATC,MCON,BIGELM,IROWBG,JCOLBG,
     $    LENBIG,IHOLD,IERNLP)
C
      IF(IERNLP.LT.0) GO TO 270
C
C         CONSTRUCT THE EXTERNAL TO INTERNAL REORDERING PERMUTATIONS
C
      CALL SHUFLE( IPRMC  ,IPRMG     
     $     ,IPRMX  ,IPRMH ,MAXCON ,MCON  ,MEQUAL ,MIGNOR    
     $     ,CBAR   ,CLWR   ,CUPR   
     $     ,VECLAM ,ISTATC ,GMAT   ,IROWG  ,JCOLG  ,NONZG   
     $     ,DELF   ,XBAR   ,XLWR ,XUPR ,VECNU  ,ISTATV ,NDIM   ,NFREE  
     $     ,HMAT   ,IROWH  ,JSTRH  ,NONZH  ,FBAR  ,CVEC ,MSUBE  
     $     ,ALWI  ,AUPI    ,ETAV  
     $     ,BMAT   ,IRWB  ,JCLB  ,NONZBR   
     $     ,CMAT   ,IRWC  ,JCLC  ,NONZCR   
     $     ,GVEC   ,YVEC   ,NVARR ,NVARMX  ,XLWI 
     $     ,XUPI  ,BVEC   ,MSUBBR ,MAXBND  ,NXUPR ,NXLWR   
     $     ,NAUPR ,NALWR  ,IBND ,VLAM   
     $     ,WMAT  ,IRWW ,JSTW,NONZWR,FNLPBR
     $     ,ISKR ,LNISKR ,RSKR  ,LNRSKR ,IERNLP 
     $     ,IFERR ,RELAX  ,FZMODE)
C
      IF(IERNLP.LT.0) GO TO 270
C
C
C ======================================================================
C ================= END INITIAL JACOBIAN ANALYSIS ======================
C ======================================================================
C
C
C ======================================================================
C ================= INITIAL DIAGNOSTIC ANALYSIS ========================
C ======================================================================
C
C         CALL DIAGNOSTIC OUTPUT PROCEDURE
C
      NEDIWK = 2*MAX(NONZR,NONZG,NDIM)
C
      IF(NEDIWK.GT.NIHOLD) THEN
        NEEDED = NEDIWK + 2
        IERNLP = -128
        GO TO 270
      ENDIF
C
      CALL DGNSTO( XBAR   ,XLWR   ,XUPR  ,NFREE
     $          ,NRES   ,DELF   ,RMAT   ,IRWR ,JSRR
     $          ,NONZR   ,IRWW  ,JSTW  
     $          ,NONZWN  ,MCON   ,MSUBE ,IPRMC ,CMAT     
     $          ,IRWC  ,JCLC  ,NONZCN ,RSKR   
     $          ,LNRSKR ,IHOLD  ,NIHOLD ,BIGELM,IROWBG,JCOLBG
     $          ,LENBIG , IERNLP )
C
      IF(IERNLP.LT.0) GO TO 270
C
      OPUN = .FALSE.
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
C
C         CHECK FOR VALID NEWTON OPTION FLAG
C
      IF(NEWTON.LT.0.OR.NEWTON.GT.2) THEN
        IERNLP = -138
        GO TO 270
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
        MSUBX = NXLWR + NXUPR
        MSUBS = NALWR + NAUPR
        IF(RELAX) THEN
          NVARP = NVARR
          NSLKP = NSLKR
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
C
C           ADJUST CHARACTER STRING ALGOPT TO RIGHT
C
        CALL HHADJF(ALGOPT,' ',' ','R',NSHIFT,IER)
C
        IF(IOFLAG.EQ.10) THEN
          CALL INSNLP('BARRIER SUMMARY')
        ELSEIF(IOFLAG.LT.20) THEN
          CALL INSNLP('BARRIER OPTION')
        ELSE 
          CALL INSNLP('BARRIER FULL OPTION')
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
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------- ALGORITHM CALLING SEQUENCE -------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
  260 CONTINUE
C
      IF(IMAX(4,IREVRS,1).GE.1) THEN
C
        IF(FZMODE) THEN
          MODE = 1
        ELSE
          MODE = -1
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
          IF(IREVRS(2).GE.1) THEN
C
C           COMPUTE GRADIENT VECTOR DELF = (RMAT**T)(RESVEC)
C
            CALL MVPSRC(11,NDIM,NRES,NONZR,RMAT,IROWR,JCOLR,RESVEC,DELF)
C
C           RESTORE ALL EXTERNAL QUANTITIES TO INTERNAL ORDER
C
            CALL HDPRMX(RMAT,NONZR,IPRMR,IERP)
C
          ENDIF
C
        ENDIF
C
        CALL EXTINT(  MODE    ,RELAX   ,IREVRS ,IPRMC 
     $        ,IPRMG  ,IPRMX  ,IPRMH  ,MAXCON     
     $        ,MCON   ,MEQUAL ,MIGNOR ,CBAR   ,VECLAM ,GMAT       
     $        ,IROWG  ,JCOLG  ,NONZG  ,DELF   ,XBAR      
     $        ,NDIM   ,NFREE  ,HMAT   ,IROWH  ,JSTRH  ,NONZH  ,FBAR        
     $        ,CVEC   ,MSUBE  ,ALWI   ,AUPI   
     $        ,ETAV   ,CMAT   ,IRWC  
     $        ,JCLC  ,NONZCR  ,GVEC   ,YVEC     
     $        ,NVARR   ,XLWI   ,XUPI   ,BVEC     
     $        ,MSUBBR  ,NXUPR  ,NXLWR  ,NAUPR  ,NALWR  ,IBND 
     $        ,VLAM   ,WMAT   ,IRWW    
     $        ,JSTW  ,NONZWR  ,FNLPBR,RSKR ,IFERR 
     $        ,IERNLP)
C
        IF(IERNLP.NE.0)  GO TO 270
C
      ENDIF
C
C             TURN OFF BCSLIB ERROR MESSAGES.
C
      IF(IOFLAG.LT.30) CALL HHERPT(0)
C
C             CALL THE OPTIMIZATION ALGORITHM BARYER
C
      CALL  BARYER( FNLPBR   ,GVEC ,YVEC ,NVARN ,NVARR
     $   ,NFIX  ,NFREE ,NXLWR, NXUPR ,IFERR    ,MAXRES   ,NRES  ,RMAT   
     $   ,IRWR,JSRR,NONZR ,CVEC ,MSUBE      
     $   ,MIGNOR,MEQUAL,MINEQL ,NALWR, NAUPR ,MAXCON ,ETAV     
     $   ,CMAT ,IRWC  ,JCLC  ,NONZCN ,NONZCR        
     $   ,BVEC ,MSUBBN ,MSUBBR ,MAXBND ,VLAM       
     $   ,BMAT ,IRWB  ,JCLB  ,NONZBN, NONZBR        
     $   ,WMAT ,IRWW  ,JSTW  ,NONZWN, NONZWR    
     $   ,IHOLD,NIHOLD ,HOLD   ,NHOLD  ,NEEDED ,IREVRS 
     $   ,IRVCOM ,RELAX,FZMODE ,IERNLP)
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
        IF(FZMODE) THEN
          MODE = 2
        ELSE
          MODE = -2
        ENDIF
C
        CALL EXTINT(  MODE    ,RELAX  ,IREVRS ,IPRMC  
     $        ,IPRMG  ,IPRMX  ,IPRMH  ,MAXCON     
     $        ,MCON   ,MEQUAL ,MIGNOR ,CBAR   ,VECLAM ,GMAT       
     $        ,IROWG  ,JCOLG  ,NONZG  ,DELF   ,XBAR      
     $        ,NDIM   ,NFREE  ,HMAT   ,IROWH  ,JSTRH  ,NONZH  ,FBAR        
     $        ,CVEC   ,MSUBE  ,ALWI   ,AUPI    
     $        ,ETAV   ,CMAT   ,IRWC  
     $        ,JCLC  ,NONZCR  ,GVEC   ,YVEC     
     $        ,NVARR   ,XLWI   ,XUPI   ,BVEC     
     $        ,MSUBBR  ,NXUPR  ,NXLWR  ,NAUPR  ,NALWR  ,IBND 
     $        ,VLAM   ,WMAT   ,IRWW    
     $        ,JSTW  ,NONZWR  ,FNLPBR,RSKR ,IFERR 
     $        ,IEREXI)
C
        IF(IEREXI.NE.0) IERNLP = IEREXI
C
        IF(IERNLP.LT.0)  GO TO 270
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
      IF((IERNLP.GE.0.OR.IERNLP.EQ.-133) .AND. .NOT.SUMRY ) THEN
C
        IF(FZMODE) THEN
          MODE = 2
        ELSE
          MODE = -2
        ENDIF
C
        CALL EXTINT(  MODE    ,RELAX  ,IREVRS ,IPRMC    
     $        ,IPRMG  ,IPRMX  ,IPRMH  ,MAXCON    
     $        ,MCON   ,MEQUAL ,MIGNOR ,CBAR   ,VECLAM ,GMAT        
     $        ,IROWG  ,JCOLG  ,NONZG  ,DELF   ,XBAR     
     $        ,NDIM   ,NFREE  ,HMAT   ,IROWH  ,JSTRH  ,NONZH  ,FBAR        
     $        ,CVEC   ,MSUBE  ,ALWI   ,AUPI    
     $        ,ETAV   ,CMAT   ,IRWC  
     $        ,JCLC  ,NONZCR  ,GVEC   ,YVEC     
     $        ,NVARR   ,XLWI   ,XUPI   ,BVEC     
     $        ,MSUBBR  ,NXUPR  ,NXLWR  ,NAUPR  ,NALWR  ,IBND 
     $        ,VLAM   ,WMAT   ,IRWW    
     $        ,JSTW  ,NONZWR  ,FNLPBR,RSKR ,IFERR 
     $        ,IEREXT)
C
        IF(IEREXT.NE.0)  THEN
          IF(IERNLP.EQ.0) IERNLP = IEREXT
          GO TO 10000
        ENDIF
C
C         CONSTRUCT ACTIVE SET INFORMATION
C
        CALL  CRSOVR( MAXCON ,MCON   ,CBAR   ,CLWR   ,CUPR   ,VECLAM    
     $         ,ISTATC ,GMAT   ,IROWG  ,JCOLG  ,NONZG  ,DELF   ,XBAR     
     $         ,XLWR   ,XUPR   ,VECNU  ,ISTATV ,NDIM )  
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
      IF(IERNLP.NE.-101.AND.IERNLP.NE.-103.AND.IERNLP.NE.-113) THEN
C
        DO I = 1,MCON
          IF(ISTATC(I).EQ.3) THEN
            IF(CBAR(I).LT.(CLWR(I)-CONTOL) .OR. 
     $         CBAR(I).GT.(CUPR(I)+CONTOL) ) ISTATC(I) = 30
          ELSEIF(ISTATC(I).EQ.1) THEN
            IF(CBAR(I).LT.(CLWR(I)-CONTOL) ) ISTATC(I) = 10
          ELSEIF(ISTATC(I).EQ.2) THEN
            IF(CBAR(I).GT.(CUPR(I)+CONTOL) ) ISTATC(I) = 20
          ENDIF
        enddo
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
        enddo
C
      ENDIF
C
      IF (IERNLP.NE.0) THEN
        IRVCOM = 0
        IF (IERNLP.GT.0) THEN
C
C           ---WARNING ERRORS
C
            MODE = 0
C
        ELSE
C
C           ---FATAL ERRORS
C
C           INPUT ARGUMENT ERROR
C
            MODE = 1
C
C           STORAGE ERROR
C
            IF(IERNLP.EQ.-127.OR.IERNLP.EQ.-128) MODE = 2
C
            IF(IERNLP.EQ.-127.OR.IERNLP.EQ.-131) THEN
C
C           INSUFFICIENT REAL STORAGE
C
              MODE = 2
C
C           IF THIS IS A LEAST SQUARES PROBLEM INCLUDE GRADIENT
C           OFFSET FROM SNLPLS
C
              IF(NRES.GT.0) NEEDED = NEEDED + NDIM
C
            ENDIF
C
            IF(IERNLP.EQ.-132) THEN
C
C           INSUFFICIENT INTEGER STORAGE FOR SRCHFZ OR EQCMIN
C
              MODE = 2
C
            ENDIF
C
C           FATAL PROCESS ERROR
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
C         CLOSE ALL MULTIFRONTAL SCRATCH FILES
C
      IF(IPUMF1.GT.0) CLOSE(IPUMF1)
      IF(IPUMF2.GT.0) CLOSE(IPUMF2)
      IF(IPUMF3.GT.0) CLOSE(IPUMF3)
      IF(IPUMF4.GT.0) CLOSE(IPUMF4)
      IF(IPUMF5.GT.0) CLOSE(IPUMF5)
      IF(IPUMF6.GT.0) CLOSE(IPUMF6)
C
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
      IF(IERNLP.GE.0) CALL STATBR
c
C ======================================================================
c
      if(allocated(ALWI)) deallocate(ALWI)
      if(allocated(AUPI)) deallocate(AUPI)
      if(allocated(BMAT)) deallocate(BMAT)
      if(allocated(BVEC)) deallocate(BVEC)
      if(allocated(CMAT)) deallocate(CMAT)
      if(allocated(CVEC)) deallocate(CVEC)
      if(allocated(ETAV)) deallocate(ETAV)
      if(allocated(GVEC)) deallocate(GVEC)
      if(allocated(RSKR)) deallocate(RSKR)
      if(allocated(VLAM)) deallocate(VLAM)
      if(allocated(WMAT)) deallocate(WMAT)
      if(allocated(XLWI)) deallocate(XLWI)
      if(allocated(XUPI)) deallocate(XUPI)
      if(allocated(YVEC)) deallocate(YVEC)
c
      if(allocated(IBND)) deallocate(IBND)
      if(allocated(IPRMC)) deallocate(IPRMC)
      if(allocated(IPRMG)) deallocate(IPRMG)
      if(allocated(IPRMH)) deallocate(IPRMH)
      if(allocated(IPRMR)) deallocate(IPRMR)
      if(allocated(IPRMX)) deallocate(IPRMX)
      if(allocated(IRWB)) deallocate(IRWB)
      if(allocated(IRWC)) deallocate(IRWC)
      if(allocated(IRWH)) deallocate(IRWH)
      if(allocated(IRWR)) deallocate(IRWR)
      if(allocated(IRWW)) deallocate(IRWW)
      if(allocated(ISKR)) deallocate(ISKR)
      if(allocated(JCLB)) deallocate(JCLB)
      if(allocated(JCLC)) deallocate(JCLC)
      if(allocated(JSRR)) deallocate(JSRR)
      if(allocated(JSTH)) deallocate(JSTH)
      if(allocated(JSTW)) deallocate(JSTW)
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
 1001 FORMAT(T3,'*********************************************** ABNORMA
     *L TERMINATION ***********************************')
C
 1002 FORMAT(2X,104('*')/T3,'*',T106,'*'/T3,'*',T34,'.....OPTIMIZATION O
     *PERATOR ',A6,'.....',T106,'*')
C
 1003 FORMAT(T3,'*',T106,'*',/2X,104('*'))
C
 1004 FORMAT(T3,'*********************************************** CONVERG
     *ENCE ********************************************')
C
 1005 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING ROWS OF THE
     $RESIDUAL JACOBIAN MAY BE ZERO',T106,'*',
     $   /T3,'*',T106,'*')
C
 1006 FORMAT(T3,'*',T106,'*')
C
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
C
 1008 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING CONSTRAINTS
     $ WILL BE IGNORED',T106,'*',
     $   /T3,'*',T106,'*')
C
 1009 FORMAT(T3,'*',T106,'*',
     $  /T3,'*',T11,85('-'),T106,'*',/T3,'*',T106,'*',
     $  /T3,'*',T11,'Least Squares Objective Function',T59,'Number of Re
     $siduals',T89,'=',I6,T106,'*')
C
 1010 FORMAT(T3,'*',T24,A60,T106,'*')
C
      END
