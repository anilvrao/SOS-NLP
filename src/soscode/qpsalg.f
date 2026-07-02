
      SUBROUTINE QPSALG(NDIM,MCON,MEQUAL,HMAT,IROWH,JSTRH,
     $    NONZH,CVCT,GZERO,AMAT,IROW,JCOLST,NONZA,SFEZ,XSTD,
     $    BUPSTD,BLWSTD,IVSTAT,WORK,NWORK,IWORK,NIWORK,NEEDED,
     $    IPU,IPC,ISTART,CONMLT,VARMLT,ISTATC,ISTATV,CNDNUM,IER)
C
C ======================================================================
C     QPSALG===>qpsalg   J.T. BETTS
C ======================================================================
C
C     ===================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:
C
C             SOLVE THE QUADRATIC PROGRAM STATED IN STANDARD FORM.
C             MINIMIZE THE QUADRATIC OBJECTIVE 
C
C               QUAD = .5*(XSTD**T)*HMATST*XSTD + (CVCTST**T)*XSTD
C
C             SUBJECT TO THE MCON CONSTRAINTS
C
C               AMATST*XSTD = BSTD 
C
C             AND THE NSTD BOUNDS
C
C               BLWSTD .LE. XSTD .LE. BUPSTD
C
C             THE STANDARD FORM PROBLEM IS CONSTRUCTED FROM THE ORIGINAL
C             PROBLEM AS FOLLOWS:
C
C             ---     Independent variables
C
C             NSTD = No. of Standard var. = MINEQL + 1 + NDIM
C
C                     |slack vars.    |   | xs |
C             XSTD =  |artificial var.| = | xi |
C                     |real vars.     |   | xr |
C
C                       |  ls  |                   |  us  |
C             BLWSTD =  |  0   |         BUPSTD =  |  1   |
C                       | XLWR |                   | XUPR |
C
C             where the slack bounds ls and us are constructed from
C             the real variables and inequality constraints.
C
C             ---     Constraints
C
C             AMATST =  | 0  : se : Ae | = | I' : s : Ar |
C                       | Is : si : Ai |
C
C             BSTD =    | be | = b
C                       | bi |
C
C             where "i" and "e" denote inequality and equalities.  The
C             matrix Ae is (MEQUAL X NDIM), Ai is (MINEQL X NDIM), Is
C             is an identity matrix of size (MINEQL X MINEQL).  The 
C             feasiblity vector s, and right hand side vector b, are 
C             both partitioned into equality and inequality parts.  The
C             constraints are then; 
C
C             I'(xs) + s(xi) + Ar(xr) = b.
C
C             ---     Objective Function 
C
C                      | 0  0  0  |            |  0  |
C             HMATST = | 0  1  0  |   CVCTST = | rho |
C                      | 0  0  Hr |            | cr  |
C
C             where Hr = HMAT, cr = CVCT, chi = CHIFEZ and rho = RHOFEZ.
C             The objective function is;
C
C             f = .5(xr**T)Hr(xr) + (cr**T)xr + (rho + .5*chi*xi)xi
C
C     ===================================================================
C
C
C         INPUT:
C
C             ----    Problem Data
C
C             NDIM    NUMBER OF REAL VARIABLES
C             MCON    NUMBER OF CONSTRAINTS
C             MEQUAL  NUMBER OF EQUALITY CONSTRAINTS
C                     NSTD = MINEQL + 1 + NDIM
C                          = (MCON-MEQUAL) + 1 + NDIM
C
C             ----    Objective Function Data
C
C             HMAT    HESSIAN MATRIX (NONZH)
C             IROWH   INTEGER ROW INDEX FOR HESSIAN (NONZH) 
C             JSTRH   COLUMN START ARRAY FOR HESSIAN (NDIM+1)
C             NONZH   NO. OF NONZERO HESSIAN ELEMENTS = JSTRH(NDIM+1)-1 
C             CVCT    OBJECTIVE FUNCTION LINEAR TERM (NDIM)
C             GZERO   GRADIENT OF OBJECTIVE AT INITIAL POINT (NDIM)
C
C             ----    Constraint Data (Dummy arguments when MCON = 0)
C
C             AMAT    JACOBIAN MATRIX (NONZA)
C             IROW    INTEGER ROW INDEX FOR JACOBIAN (NONZA)
C             JCOLST  COLUMN START ARRAY FOR JACOBIAN (NDIM+1)
C             NONZA   NO. OF NONZERO JACOBIAN ELEMENTS = JCOLST(NDIM+1)-1
C             SFEZ    CONSTRAINT FEASIBILITY RESIDUAL VECTOR (MCON)
C
C             ----    Independent Variable Data
C
C             XSTD    INITIAL GUESS FOR INDEPENDENT VARIABLES (NSTD)
C             BUPSTD  UPPER BOUND FOR INDEPENDENT VARIABLES (NSTD)
C             BLWSTD  LOWER BOUND FOR INDEPENDENT VARIABLES (NSTD)
C             IVSTAT  INTEGER VARIABLE STATUS VECTOR (MCON+1+NDIM)
C
C             ----    Algorithm Related Data
C
C             WORK    WORK ARRAY (NWORK)
C             NWORK   LENGTH OF WORK
C             IWORK   INTEGER WORK ARRAY (NIWORK)
C             NIWORK  LENGTH OF IWORK
C             IPU     OUTPUT UNIT NO.
C             IPC     OUTPUT CONTROL FLAG
C             ISTART  INTEGER START OPTION
C                     = -1 --- HOT START
C                     =  0 --- WARM START
C                     =  1 --- COLD START
C
C         OUTPUT:
C
C             XSTD    FINAL VALUE FOR INDEPENDENT VARIABLES (NSTD)
C             CONMLT  LAGRANGE MULTIPLIERS FOR REAL CONSTRAINTS (MCON)
C             VARMLT  LAGRANGE MULTIPLIERS FOR REAL VARIABLES (NDIM)
C             ISTATC  INTEGER CONSTRAINT STATUS (MCON)
C             ISTATV  INTEGER VARIABLE STATUS (NDIM)
C             CNDNUM  CONDITION NUMBER OF K-T MATRIX
C             NEEDED  STORAGE REQUIRED WHEN NWORK OR NIWORK IS TOO SMALL
C             IER     INTEGER ERROR RETURN FLAG
C                     = 0    --- SUCCESS
C                     .NE. 0 --- ERROR
C
C
C     ===================================================================
C
      DIMENSION HMAT(NONZH),IROWH(NONZH),JSTRH(NDIM+1),CVCT(NDIM),
     $    GZERO(NDIM),AMAT(NONZA),IROW(NONZA),JCOLST(NDIM+1),
     $    SFEZ(*),XSTD(MCON-MEQUAL+1+NDIM),
     $    BUPSTD(MCON-MEQUAL+1+NDIM),BLWSTD(MCON-MEQUAL+1+NDIM),
     $    IVSTAT(MCON+1+NDIM),WORK(NWORK),IWORK(NIWORK),CONMLT(*),
     $    VARMLT(NDIM),ISTATC(*),ISTATV(NDIM)
      COMMON /PERCOM/ ISQPER(20)
C
      COMMON /ITEREF/ MAXREF,MAXRFN,IREFIN
C
      LOGICAL GOTRHO, GOTCHI
C
      PARAMETER (ZERO=0.D0, ONE=1.D0)
C
C     ===================================================================
C
C         ALGORITHM CONTROL DATA
C
C         --- MAX NO. OF SCHUR-COMPLEMENT UPDATES ALLOWED
C
      MAXNM = MAX(NDIM,MCON)
      KMAX = MIN(100,2*MAXNM)
C
C         --- MAX NO. OF KT FACTORIZATIONS ALLOWED, AND COUNTER
C
      MAXFAC = MAX(5, NINT( DBLE(10*MAXNM) / DBLE(KMAX) ) )
      KTFACT = 1
C
C         --- MAX NO. OF PUNT MODE KT FACTORIZATIONS ALLOWED, AND COUNTER
C
      MXPUNT = 5
      KTPUNT = 1
C
C         --- TOLERANCE ON RECIPROCAL CONDITION NUMBER OF SCHUR-COMPLEMENT
C
      CNDTOL = 1.0D3*HDMCON(6)
C
C         --- ERROR RETURN FLAG
C
      IER = 0
C
C     ===================================================================
C
      MINEQL = MCON - MEQUAL
      NSTD = MINEQL + 1 + NDIM
      NEQNBD = NSTD + MCON
      NONZA = JCOLST(NDIM+1)-1
      NONZA = MAX(1,NONZA)
      NONZH = JSTRH(NDIM+1)-1 
      XKMAX = KMAX
      XNEQNB = NEQNBD
      XNONZH = NONZH
      XNDIM = NDIM 
      XNSQR = XNDIM**2
      XNONZA = NONZA
      XMCON = MAX(1,MCON)
      XNDMC = XNDIM*XMCON
      NONZU = XKMAX*XNEQNB*(XNONZH/XNSQR + 
     $    XNONZA/XNDMC)*1.1D0 + MIN(MCON,NDIM)
      NONZU = MAX(1,NONZU)
C
C         DEFINE WORK ARRAY LOCATIONS
C
C         --- INTEGER ARRAYS
C
      LCIFRE = 1
      LCIRWU = LCIFRE + NSTD
      LCJSTU = LCIRWU + NONZU
      LCIFRK = LCJSTU + (KMAX+1)
      LCITSH = LCIFRK + NSTD
      LCINSH = LCITSH + KMAX
      LCIRUV = LCINSH + KMAX
      LCINDX = LCIRUV + NEQNBD 
      LCIDRP = LCINDX + NSTD
      LCKPVT = LCIDRP + NSTD
      LCITMP = LCKPVT + KMAX
      LCIWKT = LCITMP + NEQNBD
      LNIWKT = NIWORK - LCIWKT + 1 
C
C         CHECK STORAGE SO FAR
C
      IF(LNIWKT.LE.0) THEN
        IER = -1110 
        NEEDED = LCIWKT 
        GO TO 10000
      ENDIF
C
C         --- REAL ARRAYS
C
      LCWORK = 1
      LCGBAR = LCWORK + NEQNBD
      LCPVEC = LCGBAR + NDIM
      LCPIVC = LCPVEC + NSTD
      LCALAM = LCPIVC + MCON
      LCKINF = LCALAM + NSTD
      LCUMAT = LCKINF + NEQNBD
      LCUVEC = LCUMAT + NONZU
      LCSHRC = LCUVEC + NEQNBD
      LCSHSV = LCSHRC + KMAX**2
      LCWRKS = LCSHSV + KMAX**2
      LCDRPL = LCWRKS + KMAX
      LCWSHR = LCDRPL + NSTD
      LCZSHR = LCWSHR + KMAX
      LCYSHR = LCZSHR + KMAX
      LCLMSH = LCYSHR + NEQNBD
      LCQSTP = LCLMSH + NSTD
      LCATPI = LCQSTP + NSTD
      LCUSHR = LCATPI + NDIM
      LCVSHR = LCUSHR + NEQNBD
      LCTSHR = LCVSHR + KMAX
      LCXBAR = LCTSHR + KMAX
      LCRHSK = LCXBAR + NSTD
      LCSOLN = LCRHSK + NEQNBD
      LCRWKT = LCSOLN + NEQNBD
      LNRWKT = NWORK - LCRWKT + 1
C
C         CHECK STORAGE SO FAR
C
      IF(LNRWKT.LE.0) THEN
        IER = -1102 
        NEEDED = LCRWKT 
        GO TO 10000
      ENDIF
C
C         INITIALIZE FEASIBLITY LINEAR AND QUADRATIC TERM PENALTY
C         WEIGHTS AND GOTRHO AND GOTCHI
C
      GOTRHO = .FALSE.
      RHOFEZ = ONE
      GOTCHI = .FALSE.
      CHIFEZ = ONE
C
 110  CONTINUE
C
C         COMPUTE NUMBER OF FREE SLACK VARIABLES
C
      NSLK = 0
      DO I=MEQUAL+1,MCON
        IF (IVSTAT(I).EQ.0)  NSLK = NSLK + 1
      ENDDO
C
C         COMPUTE NUMBER OF FEASIBILITY SLACK VARIABLES
C
      IF(IVSTAT(MCON+1).EQ.0) THEN
        NFEZ = 1
      ELSE
        NFEZ = 0
      ENDIF
C
C         COMPUTE NUMBER OF FREE REAL VARIABLES
C
      NVAR = 0
      DO I=1,NDIM
        IF (IVSTAT(MCON+1+I).EQ.0)  NVAR = NVAR + 1
      ENDDO
C
C         DEFINE NUMBER OF EQUATIONS IN AUGMENTED SYSTEM
C
      NEQNS = NSLK + NFEZ + NVAR + MCON 
C
C         DEFINE QUANTITIES FOR THE FIRST STEP
C
      CALL FRSTEP(AMAT,IROW,JCOLST,NONZA,MCON,MEQUAL,MINEQL,NSLK,
     $    NFEZ,NVAR,NEQNS,NDIM,XSTD,HMAT,IROWH,JSTRH,NONZH,
     $    CVCT,GZERO,WORK(LCWORK),NEQNBD,
     $    SFEZ,RHOFEZ,
     $    CHIFEZ,IVSTAT,WORK(LCPVEC),WORK(LCPIVC),WORK(LCRHSK),
     $    WORK(LCSOLN),WORK(LCRWKT),
     $    LNRWKT,NEEDED,IWORK(LCIFRE),IWORK(LCITMP),CNDNUM,
     $    IPU,IPC,ISTART,IER)
C
      IF(IER.NE.0) THEN
C
        IF(IER.EQ.-1102) THEN
C
C         MORE STORAGE NEEDED FOR KT SOLVE
C
          NEEDED = NEEDED + LCRWKT - 1
C
        ENDIF
        GO TO 10000
C
      ENDIF
C
C     ===================================================================
C
C         PERFORM ADDITIONAL STEPS USING SCHUR-COMPLEMENT ALGORITHM
C
      NFREE0 = NEQNS - MCON
C
C         STORE FIRST KT SOLUTION 
C
      WORK(LCKINF:LCKINF+NFREE0-1) = -WORK(LCPVEC:LCPVEC+NFREE0-1)
      WORK(LCKINF+NFREE0:LCKINF+NFREE0+MCON-1) =
     &              WORK(LCPIVC:LCPIVC+MCON-1)
C
      NRWSIZ = MAX(1,MCON)
C
      ISQPER(1) = ISQPER(1) + 1
C
      CALL SHURDV(MCON,NDIM,MINEQL,MEQUAL,NSTD,KMAX,NEQNS,IPC,
     $    IPU,NFREE0,IWORK(LCIFRE),NONZA,IROW,JCOLST,AMAT,NONZH,
     $    IROWH,JSTRH,HMAT,LNRWKT,WORK(LCRWKT),
     $    IVSTAT,CVCT,BLWSTD,
     $    BUPSTD,XSTD,GZERO,WORK(LCRHSK),WORK(LCKINF),SFEZ,NONZU,
     $    NRWSIZ,IWORK(LCIRWU),IWORK(LCJSTU),WORK(LCUMAT),
     $    IWORK(LCIFRK),IWORK(LCITSH),IWORK(LCINSH),IWORK(LCIRUV),
     $    WORK(LCUVEC),WORK(LCSHRC),WORK(LCSHSV),IWORK(LCKPVT),
     $    WORK(LCWRKS),WORK(LCWORK),IWORK(LCINDX),IWORK(LCIDRP), 
     $    WORK(LCDRPL),WORK(LCWSHR),WORK(LCZSHR),
     $    WORK(LCYSHR),WORK(LCLMSH),WORK(LCPVEC),WORK(LCQSTP),
     $    WORK(LCATPI),WORK(LCUSHR),WORK(LCVSHR),WORK(LCTSHR),
     $    RHOFEZ,GOTRHO,CHIFEZ,GOTCHI,WORK(LCGBAR),
     $    WORK(LCPIVC),WORK(LCALAM),
     $    WORK(LCXBAR),NEEDUM,CNDTOL,IRETCD,IERSHR)
C
C         CHECK ERROR FLAG
C
      IF(IERSHR.NE.0) THEN
C
C         ERROR RETURN FROM SHURDV
C
        IER = IERSHR
C
        GO TO 10000
C
      ELSE
C
        IF ( IRETCD .EQ. 3 ) THEN
          IER = IRETCD 
          GO TO 10000
        ENDIF
C
C         NORMAL RETURN FROM SHURDV -- CHECK RETURN CODE
C
        IF(IRETCD.LE.1) THEN
C
C         SOLUTION OBTAINED -- SAVE DATA AND TERMINATE
C
          XSTD(1:NSTD) = WORK(LCXBAR:LCXBAR+NSTD-1)
C
C         SAVE MULTIPLIERS AND STATUS FLAGS
C
          DO I = 1,MCON
C
            II = I - MEQUAL
C
C         CHECK CONSTRAINT FLAG
C
            IF(IVSTAT(I).EQ.-1) THEN
C
C         EQUALITY CONSTRAINT
C
              CONMLT(I) = WORK(LCPIVC+I-1)
              ISTATC(I) = 3
C
            ELSEIF(IVSTAT(I).EQ.0) THEN
C
C         INACTIVE INEQUALITY CONSTRAINT
C
              CONMLT(I) = ZERO
              ISTATC(I) = 0
C
            ELSEIF(IVSTAT(I).EQ.1) THEN
C
C         SLACK FIXED AT LOWER BOUND 
C
              IF(BLWSTD(II).LE.ZERO) THEN
C
C         INEQUALITY FIXED AT UPPER BOUND
C
                CONMLT(I) = WORK(LCPIVC+I-1)
                ISTATC(I) = 2
C
              ELSE
C
C         INEQUALITY FIXED AT LOWER BOUND
C
                CONMLT(I) = WORK(LCPIVC+I-1)
                ISTATC(I) = 1
C
              ENDIF
C
            ELSEIF(IVSTAT(I).EQ.2) THEN
C
C         SLACK FIXED AT UPPER BOUND 
C
              IF(BLWSTD(II).LE.ZERO) THEN
C
C         INEQUALITY FIXED AT LOWER BOUND
C
                CONMLT(I) = WORK(LCPIVC+I-1)
                ISTATC(I) = 1
C
              ELSE
C
C         INEQUALITY FIXED AT UPPER BOUND
C
                CONMLT(I) = WORK(LCPIVC+I-1)
                ISTATC(I) = 2
C
              ENDIF
C
            ELSEIF(IVSTAT(I).EQ.3) THEN
C
C         IGNORED CONSTRAINT
C
              CONMLT(I) = HDMCON(1)
              ISTATC(I) = 4
C
            ELSE 
C
C         INCORRECT VALUE FOR IVSTAT -- TERMINATE WITH ERROR 
C
              WRITE(*,*) 'INCORRECT IVSTAT'
              RETURN
C
            ENDIF
C
          enddo
C
          DO I = 1,NDIM 
C
            IV = MCON + 1 + I
C
C         CHECK VARIABLE FLAG
C
            IF(IVSTAT(IV).EQ.0) THEN
C 
C         FREE VARIABLE 
C
              VARMLT(I) = ZERO
              ISTATV(I) = 0
C
            ELSEIF(IVSTAT(IV).EQ.1) THEN
C
C         REAL VARIABLE FIXED AT LOWER BOUND 
C
              VARMLT(I) = WORK(LCALAM+MINEQL+I)
              ISTATV(I) = 1
C
            ELSEIF(IVSTAT(IV).EQ.2) THEN
C
C         REAL VARIABLE FIXED AT UPPER BOUND 
C
              VARMLT(I) = WORK(LCALAM+MINEQL+I)
              ISTATV(I) = 2
C
            ELSEIF(IVSTAT(IV).EQ.3) THEN
C
C         REAL VARIABLE FIXED PERMANENTLY 
C
              VARMLT(I) = WORK(LCALAM+MINEQL+I)
              ISTATV(I) = 3
C
            ELSE 
C
C         INCORRECT VALUE FOR IVSTAT -- TERMINATE WITH ERROR 
C
              WRITE(*,*) 'INCORRECT IVSTAT'
              RETURN
C
            ENDIF
C
          enddo
C
          IER = IRETCD
C
        ELSEIF(IRETCD.EQ.11) THEN
C
C         INERTIA WENT BAD AFTER A WARM START
C
          IER = -1107
C
        ELSEIF(IRETCD.EQ.12) THEN
C
C         INITIAL KT MATRIX IS TOO ILL-CONDITIONED FOR SCHUR COMPLEMENT
C         UPDATES
C
          IER = -1106
C
        ELSEIF(IRETCD.EQ.13) THEN
C
C         ILL-CONDITIONING WITH FEASIBILITY VARIABLE NOT ON A BOUND
C
          IER = -1109
C
        ELSEIF(IRETCD.EQ.14) THEN
C
C         SOLVE FAILED BECAUSE OF ILL-CONDITIONING
C
          IER = -1106
C
        ELSE
C
C         SOLUTION NOT FOUND -- NEW KT FACTORIZATION REQUIRED
C         SAVE CURRENT POINT IN XSTD
C
          XSTD(1:NSTD) = WORK(LCXBAR:LCXBAR+NSTD-1)
C
          IF(KTFACT.LT.MAXFAC.AND.IRETCD.EQ.4) THEN
C
C         INCREMENT THE NUMBER OF KT FACTORIZATIONS
C
            KTFACT = KTFACT + 1
C
            GO TO 110
C
          ELSEIF(KTPUNT.LT.MXPUNT.AND.IRETCD.NE.4) THEN
C
C         ACTIVATE ITERATIVE REFINEMENT 
C
            IF(IRETCD.EQ.7) IREFIN = 2
            IF(IRETCD.EQ.5) THEN
              IREFIN = 2
              CNDTOL = CNDTOL/1.0D3
            ENDIF
C
C         INCREMENT THE NUMBER OF KT FACTORIZATIONS
C
            KTPUNT = KTPUNT + 1
C
            GO TO 110
C
          ELSE
C
C         MAX. NUMBER OF REFACTORIZATIONS EXCEEDED
C
            IER = 1017
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
C     ===================================================================
C
10000 CONTINUE
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
