
      SUBROUTINE SHURQP(NDIM,NCON,HMAT,IROWH,JSTRH,NZHDIM,
     $    CVCT,GZERO,QUAD,AMAT,IROW,JCOLST,NZADIM,BUPR,BLWR,CONMLT,
     $    ISTATC,XVEC,XUPR,XLWR,VARMLT,ISTATV,WORK,NWORK,IWORK,NIWORK,
     $    NEEDED,BIGBND,IPU,IPC,ISTART,CNDNUM,IER)
C
C ======================================================================
C     SHURQP===>shurqp   J.T. BETTS
C ======================================================================
C
C     ==================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C         PURPOSE:
C
C             MINIMIZE THE QUADRATIC OBJECTIVE 
C
C               QUAD = .5*(XVEC**T)*HMAT*XVEC + (CVCT**T)*XVEC
C
C             SUBJECT TO THE NCON CONSTRAINTS
C
C               BLWR .LE. AMAT*XVEC .LE. BUPR
C
C             AND THE NDIM BOUNDS
C
C               XLWR .LE. XVEC .LE. XUPR
C
C             EQUALITY CONSTRAINTS ARE IMPOSED BY SETTING BLWR = BUPR.
C             BOTH THE HESSIAN HMAT AND THE JACOBIAN AMAT ARE SPARSE
C             MATRICES.  THE ALGORITHM IMPLEMENTS THE SCHUR-COMPLEMENT 
C             QUADRATIC PROGRAMMING METHOD, PROPOSED BY GILL, MURRAY,
C             SAUNDERS, AND WRIGHT.
C
C             REF.  "A Schur-Complement Method for Sparse Quadratic 
C                   Programming,"  Philip E. Gill, Walter Murray,
C                   Michael A. Saunders and Margaret H. Wright,  Tech
C                   Report SOL 87-12, October, 1987, Systems Optimization
C                   Laboratory, Dept. of Operations Research, Stanford
C                   University, Stanford, Calilfornia 94305-4022.
C
C
C         ARGUMENTS:
C
C             ----    Problem Data
C
C             NDIM    I    NUMBER OF VARIABLES
C             NCON    I    NUMBER OF CONSTRAINTS
C
C             ----    Objective Function Data
C
C             HMAT    I    HESSIAN MATRIX (NZHDIM)
C             IROWH   I    INTEGER ROW INDEX FOR HESSIAN (NZHDIM) 
C             JSTRH   I    COLUMN START ARRAY FOR HESSIAN (NDIM+1)
C             NZHDIM  I    DIMENSION OF NONZERO HESSIAN ELEMENTS
C                          NZHDIM.GE.NONZH = JSTRH(NDIM+1)-1 
C             CVCT    I    OBJECTIVE FUNCTION LINEAR TERM (NDIM)
C             GZERO   O    GRADIENT OF OBJECTIVE AT INITIAL POINT (NDIM)
C             QUAD    O    OPTIMAL OBJECTIVE FUNCTION VALUE
C
C             ----    Constraint Data (Dummy arguments when NCON = 0)
C
C             AMAT    I    JACOBIAN MATRIX (NZADIM)
C             IROW    I    INTEGER ROW INDEX FOR JACOBIAN (NZADIM)
C             JCOLST  I    COLUMN START ARRAY FOR JACOBIAN (MAX(NZERO,NDIM+1))
C             NZADIM  I    DIMENSION OF NONZERO JACOBIAN ELEMENTS
C                          NZADIM.GE.NZERO = JCOLST(NDIM+1)-1
C                          WHEN ISTART .LT. 2
C             BUPR    I    CONSTRAINT UPPER BOUND VECTOR (NCON)
C             BLWR    I    CONSTRAINT LOWER BOUND VECTOR (NCON)
C             CONMLT  O    LAGRANGE MULTIPLIERS FOR CONSTRAINTS (NCON)
C             ISTATC  I/O  INTEGER CONSTRAINT STATUS (NCON)
C                          = 0  --- FREE (INACTIVE) INEQUALITY
C                          = 1  --- FIXED ON LOWER BOUND
C                          = 2  --- FIXED ON UPPER BOUND
C                          = 3  --- EQUALITY
C                          = 4  --- IGNORED CONSTRAINT
C                          = -1 --- VIOLATES LOWER BOUND
C                          = -2 --- VIOLATES UPPER BOUND
C
C             ----    Independent Variable Data
C
C             XVEC    I    INITIAL GUESS FOR INDEPENDENT VARIABLES (NDIM)
C                     O    FINAL VALUE FOR INDEPENDENT VARIABLES (NDIM)
C             XUPR    I    UPPER BOUND FOR INDEPENDENT VARIABLES (NDIM)
C             XLWR    I    LOWER BOUND FOR INDEPENDENT VARIABLES (NDIM)
C             VARMLT  O    LAGRANGE MULTIPLIERS FOR VARIABLES (NDIM)
C             ISTATV  I/O  INTEGER VARIABLE STATUS (NDIM)
C                          = 0  --- FREE VARIABLE 
C                          = 1  --- FIXED ON LOWER BOUND
C                          = 2  --- FIXED ON UPPER BOUND
C                          = 3  --- FIXED PERMANENTLY 
C                          = 4  --- IGNORED BOUND
C                          = -1 --- VIOLATES LOWER BOUND
C                          = -2 --- VIOLATES UPPER BOUND
C
C             ----    Algorithm Related Data
C
C             WORK    I    WORK ARRAY (NWORK)
C             NWORK   I    LENGTH OF WORK
C             IWORK   I    INTEGER WORK ARRAY (NIWORK)
C             NIWORK  I    LENGTH OF IWORK
C             NEEDED  O    STORAGE REQUIRED WHEN EITHER NWORK OR NIWORK IS TOO SMALL
C             BIGBND  I    BIG BOUND VALUE
C             IPU     I    OUTPUT UNIT NO.
C             IPC     I    OUTPUT CONTROL FLAG
C                          = 0 --- NO PRINT
C                          IPC .GT. 100 --- DEBUG PRINT
C             ISTART  I    INTEGER START OPTION FLAG
C                          = -1 --- HOT START
C                          =  0 --- WARM START
C                          =  1 --- COLD START; JACOBIAN IN COLUMN FORMAT
C                          =  2 --- COLD START; JACOBIAN IN TRIPLE FORMAT
C             CNDNUM  O    CONDITION NUMBER OF K-T MATRIX
C             IER     O    INTEGER ERROR RETURN FLAG
C                          = 0    --- SUCCESS
C     -------------------------------------------------------------------
C            .GT. 0 --- WARNING ERRORS
C                          = 1     --- WEAK SOLUTION
C                          = 2     --- SOLUTION WITH RELAXED CONSTRAINT TOLERANCE
C                          = 3     --- RHO .GT. RHOMAX: INFEASIBLE CONSTRAINTS
C                          = 1017  --- MAXIMUM NUMBER OF KT REFACTORIZATIONS
C     -------------------------------------------------------------------
C            .LE. 0 --- FATAL (USUALLY INPUT) ERRORS
C                          = -1001 --- NDIM .LE. 0
C                          = -1002 --- NCON .LT. 0
C                          = -1003 --- WARM/HOT START WITHOUT COLD START FIRST
C                          = -1004 --- JSTRH(1) .LE.0 OR JSTRH(I) .GT. JSTRH(I+1)
C                          = -1005 --- NONZH .LE. 0 OR NONZH .GT. N(N+1)/2 OR
C                                      NONZH .GT. NZHDIM
C                          = -1006 --- IROWH(I) .LE. IDIAG OR IROWH(I) .GT. NDIM
C                                      WHERE IDIAG = ROW INDEX OF DIAGONAL.
C                                      HESSIAN MUST BE LOWER TRIAGULAR.
C                          = -1007 --- JCOLST(1) .LE.0 OR JCOLST(I) .GT. JCOLST(I+1)
C                          = -1008 --- NZERO .LE. 0 OR NZERO .GT. NDIM*NCON
C                                      OR NZERO .GT. NZADIM
C                          = -1009 --- IROW(I) .LE. 0 OR IROW(I) .GT. NCON
C                          = -1010 --- BUPR(i) .LT. BLWR(i) OR
C                                      BUPR(i) .GE. BIGBND AND BLWR(i) .LE. -BIGBND
C                          = -1011 --- XUPR(i) .LT. XLWR(i) 
C                          = -1012 --- MEQUAL + NFIXVR .GT. NDIM
C                          = -1013 --- IPC .LT. 0 AND IPU .GT. 0
C                          = -1014 --- INSUFFICIENT REAL WORK STORAGE
C                          = -1016 --- ZERO ROW IN CONSTRAINT MATRIX
C                                      WHEN IPC .GT. 100 ZERO ROW INDICES PRINTED
C                          = -1017 --- ERROR IN CNVRTR
C                          = -1018 --- ERROR IN DCNVRT
C                          = -1019 --- INSUFFICIENT INTEGER WORK STORAGE
C                          = -1020 --- INCONSISTENT LINEAR CONSTRAINTS (NEQNS > NEQNBD)
C                          = -1021 --- INCORRECT VALUE FOR KTOPTN; DEFAULTS NEED TO BE SET
C                          = -1103 --- SINGULAR KT MATRIX (MULTIFRONTAL IER = -503)
C                          = -1104 --- INERTIA OF KT MATRIX IS INCORRECT
C                          = -1105 --- OTHER MULTIFRONTAL ERRORS
C                          = -1106 --- CONDITION NUMBER OF KT MATRIX IS TOO LARGE
C                                      (OR IRETCD = 12 OR 14)
C                          = -1107 --- INCORRECT INERTIA AFTER ACTIVE SET CHANGE
C                                      AFTER WARM START (IRETCD=11).
C                          = -1108 --- EXCESSIVE FILL DURING NUMERIC FACTORIZATION
C                                      PRESUMABLY BECAUSE OF ILL-CONDITIONING
C                          = -1109 --- INCORRECT INERTIA FOR SCHUR-COMPLEMENT
C                                      WITH RELAXED FEASIBILITY TOLERANCE
C                          = -1111 --- I/O ERROR (INSUFFICIENT DISK SPACE)
C                          = -1999 --- EXTERNAL USER KILL
C
C
C     ===================================================================
C
      DIMENSION HMAT(NZHDIM),IROWH(NZHDIM),JSTRH(NDIM+1),CVCT(NDIM),
     $    GZERO(NDIM),AMAT(NZADIM),IROW(NZADIM),JCOLST(*),
     $    BUPR(*),BLWR(*),XVEC(NDIM),XUPR(NDIM),
     $    XLWR(NDIM),WORK(NWORK),IWORK(NIWORK),CONMLT(*),VARMLT(NDIM),
     $    ISTATC(*),ISTATV(NDIM)
C
C
C     KTOPTN DETERMINES WHETHER TO USE OLD OR NEW SPARSE QP. THE NEW
C     SPARSE QP HAS ONLY ACTIVE CONSTRAINTS IN THE KKT MATRIX.
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO = 2.D0)
C
C     ===================================================================
C
      IER = 0
C
C         START SCHUR-QP TIMING CLOCKS
C
      CALL CLKBEG(4)
      CALL CLKBEG(5)
C
      IF(ISTART.GT.0) THEN
C
C         INPUT CHECKING
C
        IF(NDIM.LE.0) THEN
          IER = -1001
          GO TO 10000
        ENDIF
        IF(NCON.LT.0) THEN
          IER = -1002
          GO TO 10000
        ENDIF
C
        IF(JSTRH(1).LE.0) THEN
          IER = -1004
          GO TO 10000
        ENDIF
        DO I = 1,NDIM
          IF(JSTRH(I).GT.JSTRH(I+1)) THEN
            IER = -1004
            GO TO 10000
          ENDIF
        enddo
C
        XNDIM = NDIM
        NONZH = JSTRH(NDIM+1)-1
        XNONZH = NONZH
        XNDNSH = XNDIM*(XNDIM+ONE)/TWO
        IF(NONZH.LE.0.OR.XNONZH.GT.XNDNSH.OR.NONZH.GT.NZHDIM) THEN
          IER = -1005
          GO TO 10000
        ENDIF
C
        DO I = 1,NDIM
          JROW = JSTRH(I)
          IROWMN = IROWH(JROW)
          IROWMX = IROWH(JROW)
          DO J=JROW+1,JSTRH(I+1)-1
            IROWMN = MIN(IROWMN,IROWH(J))
            IROWMX = MAX(IROWMX,IROWH(J))
          ENDDO
          IF(IROWMN.LT.I .OR. IROWMX.GT.NDIM) THEN
            IER = -1006
            GO TO 10000
          ENDIF
        enddo
C
        IF(NCON.GT.0) THEN
          IF(ISTART.EQ.1) THEN
            IF(JCOLST(1).LE.0) THEN
              IER = -1007
              GO TO 10000
            ENDIF
            DO I = 1,NDIM
              IF(JCOLST(I).GT.JCOLST(I+1)) THEN
                IER = -1007
                GO TO 10000
              ENDIF
            enddo
C
C           NOTE THE NUMBER OF NONZEROS IN THE JACOBIAN, NZERO IS
C           .GE. NONZA (THE NUMBER OF NONZEROS IN THE NONIGNORED
C           PART OF THE JACOBIAN)
C
            NZERO = JCOLST(NDIM+1)-1
            XNDIM = NDIM
            XNCON = NCON
            XNZERO = NZERO
            XNDNSA = XNDIM*XNCON
            IF(NZERO.LE.0.OR.XNZERO.GT.XNDNSA.OR.NZERO.GT.NZADIM) THEN
              IER = -1008
              GO TO 10000
            ENDIF
          ELSE
            NZERO = NZADIM
            XNDIM = NDIM
            XNCON = NCON
            XNZERO = NZERO
            XNDNSA = XNDIM*XNCON
            IF(NZERO.LE.0.OR.XNZERO.GT.XNDNSA) THEN
              IER = -1008
              GO TO 10000
            ENDIF
            JCOLMN = 1
            JCOLMX = NDIM
            DO I=1,NZERO
              JCOLMN = MIN(JCOLMN,JCOLST(I))
              JCOLMX = MAX(JCOLMX,JCOLST(I))
            ENDDO
            IF(JCOLMN.LE.0 .OR. JCOLMX.GT.NDIM) THEN
              IER = -1007
              GO TO 10000
            ENDIF
          ENDIF
C
          IROWMN = 1
          IROWMX = NCON
          DO I=1,NZERO
            IROWMN = MIN(IROWMN,IROW(I))
            IROWMX = MAX(IROWMX,IROW(I))
          ENDDO
          IF(IROWMN.LE.0 .OR. IROWMX.GT.NCON) THEN
            IER = -1009
            GO TO 10000
          ENDIF
        ENDIF
C
        MEQUAL = 0
        MIGNOR = 0
C
C         LOOK FOR EQUALITY CONSTRAINTS AND CHECK BOUNDS
C 
        DO I = 1,NCON
C
          IF(ISTATC(I).NE.4) THEN
            IF(BUPR(I).LT.BLWR(I) .OR. 
     $      (BUPR(I).GE.BIGBND .AND. BLWR(I).LE.-BIGBND) ) THEN
C
C         CONSTRAINT BOUNDS ARE INCORRECT
C
              IER = -1010
              GO TO 10000
C
            ELSEIF(BUPR(I).EQ.BLWR(I)) THEN
C
C         CONSTRAINT IS AN EQUALITY
C
              MEQUAL = MEQUAL + 1
C
            ENDIF
C
          ELSE
C
C          CONSTRAINT IS IGNORED
C
            MIGNOR = MIGNOR + 1
C
          ENDIF
C
        enddo
C
        NFIXVR = 0
C
C         CHECK VARIABLE BOUNDS
C
        DO I = 1,NDIM
C
          IF(XUPR(I).LT.XLWR(I)) THEN
C
C         VARIABLE BOUNDS ARE INCORRECT
C
            IER = -1011
            GO TO 10000
C
          ELSEIF(XUPR(I).EQ.XLWR(I)) THEN
C
C         VARIABLE BOUND IS FIXED 
C
            NFIXVR = NFIXVR + 1
C
          ENDIF
C
        enddo
C
        IF(MEQUAL+NFIXVR.GT.NDIM) THEN
          IER = -1012
          GO TO 10000
        ENDIF
C
C         NUMBER OF INEQUALITIES
C
        MINEQL = NCON - MEQUAL - MIGNOR
C
C         NUMBER OF NONIGNORED CONSTRAINTS = MCON
C         NOTE:  NCON = TOTAL NO. OF CONSTRAINTS .GE. MCON
C
        MCON = MEQUAL + MINEQL
C
C         LENGTH OF STANDARD VARIABLE ARRAYS
C
        NSTD = MINEQL + 1 + NDIM
C
        IF(IPC.LT.0.AND.IPU.GT.0) THEN
          IER = -1013
          GO TO 10000
        ENDIF
C
C         SET UP WORK ARRAY POINTERS 
C
C         --- INTEGER WORK ARRAY 
C             (BEGIN IN IWORK(1))
C
        LCWRKI = 1
        LNWRKI = MAX(2*NCON,NZERO,NDIM)
        LCIPRM = LCWRKI + LNWRKI
        LCIPRC = LCIPRM + NZERO
        LCIVST = LCIPRC + MCON
        LCIWRK = LCIVST + (MCON + 1 + NDIM) 
        LNIWRK = NIWORK - LCIWRK + 1
C
        IF(LNIWRK.LE.0) THEN
          NEEDED = LCIWRK + 1
          IER = -1019
          GO TO 10000
        ENDIF
C
C         --- REAL WORK ARRAY (BEGIN IN WORK(1))
C
        LCXSTD = 1
        LCBSTD = LCXSTD + NSTD
        LCBLST = LCBSTD + MCON
        LCBUST = LCBLST + NSTD
        LCSFEZ = LCBUST + NSTD
        LCRWRK = LCSFEZ + (MCON + 1) 
        LNRWRK = NWORK - LCRWRK + 1
C
        IF(LNRWRK.LT.MAX(1,MCON)) THEN
          NEEDED = LCRWRK + MAX(1,MCON)
          IER = -1014
          GO TO 10000
        ENDIF
C
      ELSEIF(LCWRKI.NE.1) THEN
C
C         INVALID SEQUENCE:  WARM OR HOT START BEFORE A COLD START
C                            ISTART NOT PROPERLY INITIALIZED
C
        IER = -1003
        GO TO 10000
C
      ENDIF
C
C         CONVERT PROBLEM TO INTERNAL FORMAT -- EQUALITIES FIRST, 
C         THEN INEQUALITIES, THEN IGNORED CONSTRAINTS
C
      IF(NCON.GT.0) THEN
C
        LNJCOL = MAX(NZERO,NDIM+1)
        CALL CNVRTR(MAX(0,ISTART),NDIM,NCON,NZERO,AMAT,IROW,JCOLST,
     $    LNJCOL,BUPR,BLWR,ISTATC,IWORK(LCWRKI),LNWRKI,NONZA,
     $    IWORK(LCIPRM),IWORK(LCIPRC),IERCNV) 
C
        IF(IERCNV.NE.0) THEN
          IER = -1017
          GO TO 10000
        ENDIF
C
      ENDIF
C
C         CONSTRAINT TOLERANCE
C
      TYTTOL = HDMCON(5)**.7
      BIGTOL = SQRT(TYTTOL)
      TOLCON = TYTTOL
C
 160  CONTINUE
      LNNNZA = MAX(1,NONZA)
C
C         INITIALIZE VARIABLES AND PUT QP INTO STANDARD FORM
C
      CALL INSHUR(AMAT,IROW,JCOLST,LNNNZA,XVEC,XLWR,XUPR,
     $    BLWR,BUPR,WORK(LCRWRK),BIGBND,TOLCON,MCON,MAX(1,MCON),
     $    MEQUAL,NDIM,NSTD,WORK(LCXSTD),WORK(LCBSTD),WORK(LCBLST),
     $    WORK(LCBUST),WORK(LCSFEZ),IWORK(LCIVST),MAX(0,ISTART))
C
C
C$PDF
C     NOW THAT INSHUR HAS TESTED SFEZ AGAINST TOLCON, SCALE SFEZ
C     TO THE NEAREST DECADE SO THAT IT IS ON THE ORDER OF UNITY. 
C     THE FEASIBILITY VARIABLE AND ITS BOUNDS MUST ALSO BE SCALED
C     TO THE INVERSE OF THE SCALING USED FOR SFEZ.
C     THE ABOVE SCALING IS USED TO KEEP THE KT MATRIX WELL SCALED.
C
      RESNRM = DAMAX(MCON,WORK(LCSFEZ),1)
C$PDF IF ( RESNRM .GT. ZERO ) THEN
C
C$PDF   IEXP = INT( LOG10(RESNRM) )
C$PDF   SFZSCL = ONE / (10.0D0**IEXP)
C$PDF   XISCL = ONE/SFZSCL
C
C$PDF   WORK(LCXSTD+MINEQL) = XISCL*WORK(LCXSTD+MINEQL)
C$PDF   WORK(LCBLST+MINEQL) = XISCL*WORK(LCBLST+MINEQL)
C$PDF   WORK(LCBUST+MINEQL) = XISCL*WORK(LCBUST+MINEQL)
C
C$PDF ENDIF
C      
C$PDF
C         CHECK THE JACOBIAN MATRIX FOR A ZERO ROW
C
      IF(MCON.GT.0) CALL ROWCHK(IROW,JCOLST,IWORK(LCIVST),
     $    IWORK(LCWRKI),NDIM,MCON,LNNNZA,IPC,IER)
C
      IF(IER.NE.0) GO TO 10000
C
C     ===================================================================
C
C         CALL THE ALGORITHM QPSALG/QPSALN TO SOLVE THE QUADRATIC PROGRAM
C         IN STANDARD FORM
C
      IF ( KTOPTN .EQ. 'SMALL ' ) THEN
        CALL QPSALN(NDIM,MCON,MEQUAL,HMAT,IROWH,JSTRH,NONZH,
     $      CVCT,GZERO,AMAT,IROW,JCOLST,LNNNZA,WORK(LCSFEZ),
     $      WORK(LCXSTD),WORK(LCBUST),WORK(LCBLST),IWORK(LCIVST),
     $      WORK(LCRWRK),LNRWRK,IWORK(LCIWRK),LNIWRK,NEEDED,
     $      IPU,IPC,ISTART,CONMLT,VARMLT,ISTATC,ISTATV,CNDNUM,IER)
C
      ELSEIF ( KTOPTN .EQ. 'LARGE ' ) THEN
        CALL QPSALG(NDIM,MCON,MEQUAL,HMAT,IROWH,JSTRH,NONZH,
     $      CVCT,GZERO,AMAT,IROW,JCOLST,LNNNZA,WORK(LCSFEZ),
     $      WORK(LCXSTD),WORK(LCBUST),WORK(LCBLST),IWORK(LCIVST),
     $      WORK(LCRWRK),LNRWRK,IWORK(LCIWRK),LNIWRK,NEEDED,
     $      IPU,IPC,ISTART,CONMLT,VARMLT,ISTATC,ISTATV,CNDNUM,IER)
      ELSE
C
C         KTOPTN IS INVALID:  SET DEFAULTS
C
        IER = -1021
        GO TO 10000
      ENDIF
C
C         IF ALGORITHM WAS SUCCESSFUL LOAD RESULTS
C
      IF(IER.EQ.0.OR.IER.EQ.1) THEN
C
C         COMPUTE OPTIMAL OBJECTIVE FUNCTION (use xvec as work array)
C
        CALL SQFORM(HMAT,IROWH,JSTRH,CVCT,WORK(LCXSTD+MINEQL+1),
     $      XVEC,NDIM,QUAD)
C
C         LOAD REAL VARIABLES
C
        XVEC(1:NDIM) = WORK(LCXSTD+MINEQL+1:LCXSTD+MINEQL+NDIM)
C
C         RESET ERROR RETURN IF RELAXED TOLERANCES WERE REQUIRED
C
        IF(TOLCON.EQ.BIGTOL) IER = 2
C
      ELSE
C
        IF(IER.EQ.-1102) THEN
          NEEDED = NEEDED + LCRWRK - 1
          IER = -1015
          GO TO 10000
        ENDIF
C
        IF(IER.EQ.-1110) THEN
          NEEDED = NEEDED + LCIWRK - 1
          IER = -1019
          GO TO 10000
        ENDIF
C
        IF(IER.EQ.-1999) THEN
          GO TO 10000
        ENDIF
C
        IF(ISTART.EQ.0.AND.(IER.EQ.-1103.OR.IER.EQ.-1104
     $    .OR.IER.EQ.-1106.OR.IER.EQ.-1108.OR.IER.EQ.1017)) THEN
C
C         INITIAL KT MATRIX WAS SINGULAR FOR A WARM START
C         RESET TO A COLD START AND TRY AGAIN
C
          ISTART = 1
          IER = 0
          GO TO 160
C
        ENDIF
C
        IF(IER.EQ.-1109) THEN
          IF(TOLCON.EQ.TYTTOL.AND.RESNRM.LT.BIGTOL) THEN
            TOLCON = BIGTOL
            GO TO 160
          ENDIF
        ENDIF
C
      ENDIF
C
      IF(NCON.GT.0) THEN
C
        CALL DCNVRT(NDIM,NCON,MEQUAL,MINEQL,MIGNOR,NZERO,AMAT,
     $    IROW,JCOLST,LNJCOL,BUPR,BLWR,CONMLT,IWORK(LCWRKI),
     $    LNWRKI,IWORK(LCIPRM),IWORK(LCIPRC),ISTATC,IERDCN)
C
        IF(IERDCN.NE.0) THEN
          IER = -1018
          GO TO 10000
        ENDIF
C
      ENDIF
C
10000 CONTINUE
C
C         COMBINE REAL AND INTEGER STORAGE ERRORS INTO A SINGLE ERROR
C
      IF(IER.EQ.-1014.OR.IER.EQ.-1015.OR.IER.EQ.-1101.OR.IER.EQ.-1102)
     $     THEN
        IER = -1014

      ELSEIF(IER.EQ.-1019.OR.IER.EQ.-1110) THEN
        IER = -1019
      ENDIF
C
C         CALL ERROR HANDLER IF NECESSARY
C
      IF(IER.NE.0) THEN
        IF(IER.GT.0) THEN
          MODE = 0
        ELSEIF(IER.EQ.-1014.OR.IER.EQ.-1015) THEN
          MODE = 2
        ELSEIF(IER.EQ.-1019) THEN
          MODE = 2
        ELSEIF(IER.GT.-91) THEN
          MODE = 1
        ELSE
          MODE = 3
        ENDIF
C
        CALL HHERR(MODE,'SHURQP  ',IER,NEEDED)
C
      ENDIF
C
C         STOP SCHUR-QP TIMING CLOCKS
C
      CALL CLKSUM(4)
      CALL CLKMAX(5)
C
C     ===================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
