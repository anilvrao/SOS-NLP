
      SUBROUTINE SHURUP( NEQNS, IPC, IPU, K, KMAX, U, V, SIGMA, SOLVTL,
     $                   DPFSUP, KPVT, SHRCSV, SHURC, MXNZUM,
     $                   IROWU, JSTRU, UMAT, WRK, NWKKT, WORKKT,
     $                   T, GAMMA, CHI, GOTCHI, IER )
C
C ======================================================================
C     SHURUP===>shurup   J.T. BETTS
C ======================================================================
C
C     THIS ROUTINE UPDATES THE SCHUR-COMPLEMENT MATRIX
C                         T    -1
C            SHURC = V - U (K0)  U, WHERE
C     U IS AN NEQNS BY K MATRIX AND V IS A K BY K MATRIX.
C     (K IS THE TOTAL NUMBER OF SCHUR-COMPLEMENT UPDATES
C     PRIOR TO THE CURRENT ONE. NEQNS IS THE SIZE OF THE
C     KT SYSTEM AT THE MOST RECENT FACTORIZATION.)
C             _                        _
C     THE NEW U = | U  u | AND THE NEW V = | V        v  | .
C                 |      |                 |(v)TR  SIGMA |
C                   _____
C     THUS, THE NEW SHURC = | SHURC       t   |  ,
C                           | (t)TR     GAMMA | 
C
C                    T    -1                       T    -1
C     WHERE t = v - U (K0)  u AND GAMMA = SIGMA - u (K0)  u .
C
C
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
C
      LOGICAL DPFSUP, GOTCHI
C
      INTEGER NEQNS, NRHS, NWKKT, NEEDS, JER
      INTEGER I, J, IPC, IPU, IER, IOPT, K
      INTEGER KMAX, MXNZUM, KPVT(KMAX)
      INTEGER IROWU(MXNZUM), JSTRU(KMAX+1) 
C
      DOUBLE PRECISION CHI, GAMMA, SIGMA, UTKIVU, TTCIVT
      DOUBLE PRECISION WRK(NEQNS), U(NEQNS), WORKKT(NWKKT)
      DOUBLE PRECISION UMAT(MXNZUM), T(KMAX), V(KMAX)
      DOUBLE PRECISION SHRCSV(KMAX,KMAX), SHURC(KMAX,KMAX)
      DOUBLE PRECISION DELTA, SOLVTL
      COMMON /PERCOM/ ISQPER(20)
      INTEGER ISQPER
C
C
      IF ( IPC .GT. 0 ) THEN
        WRITE(IPU,*) ' AFTER UPDATE: K+1 = ', K+1
      ENDIF
C
C
      IER = 0
C     
C                 -1 
C     WRK <-- (K0)  * u
C
      WRK(1:NEQNS) = U(1:NEQNS)
C
      ISQPER(10) = ISQPER(10) + 1
C
C         SOLVE LINEAR SYSTEM
C
      CALL REFITR(U,WRK,NEQNS,WORKKT,NWKKT,NEEDS,JER)
C
      IF ( JER .NE. 0 ) THEN
C       ALLOW EXCEEDING REFITR'S SOLUTION ERROR GOAL, BUT DON'T
C       ACCEPT AN ERROR GREATER THAN SOLVTL 
csolvtl        IF ( JER .NE. -701 .OR. DELTA .GT. SOLVTL ) THEN
        IF ( JER .NE. -701  ) THEN
          IER = 1
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' SHURUP: ERROR, JER IN HDSLSL=', JER
          ENDIF
          GO TO 999
        ENDIF
      ENDIF
C
C       
      IF ( K .GT. 0 ) THEN
C
C              T          T      -1
C       t <-- U * WRK == U * (K0)  * u
C
        IOPT = 11
        CALL MVPSPR( IOPT, K, NEQNS, UMAT, IROWU, JSTRU, WRK, T )
C
C                           T     -1
C       t <-- v - t == v - U *(K0)  * u
C
        DO I = 1, K
          T(I) = V(I) - T(I)
        enddo
C
      ENDIF
C
C                 T      -1
C     UTKIVU <-- u * (K0)  * u
C
      UTKIVU =  DOT_PRODUCT(U(1:NEQNS),WRK(1:NEQNS))
C
C
      IF ( DPFSUP ) THEN
C
C       THIS UPDATE IS A RESULT OF DROPPING THE FEASIBILITY VAR. SPECIAL
C       CARE MUST BE TAKEN SO THAT THE FEASIBILITY VAR WILL GET SMALLER
C       THAN ITS BOUND (1.0) WHEN RHO IS MADE LARGE MADE ENOUGH.
C       THAT IS, IT MUST BE MADE POSSIBLE FOR A NEGATIVE STEP IN THE
C       FEASIBILITY VAR TO REACH A QP STATIONARY POINT AS RHO IS MADE
C       LARGE ENOUGH.
C
        IF ( K .GT. 0 ) THEN
C                                                T
C         SOLVE C*WRK = t FOR WRK, I.E. SOLVE RDR *WRK = t FOR WRK.
C
          WRK(1:K) = T(1:K)
C
          CALL DSISL( SHURC, KMAX, K, KPVT, WRK )
C
C                     T    T -1      T
C         TTCIVT <-- t (LDL )  t == t *WRK
C
          TTCIVT = DOT_PRODUCT(T(1:K),WRK(1:K))
        ELSE
          TTCIVT = ZERO
        ENDIF
C                                               T 
C       ANALYSIS OF UPDATE AND SOLVE FOR THE LDL  FACTORIZATION
C       OF THE SCHUR COMPLEMENT INDICATES THAT A SUFFICIENCY CONDITION
C       FOR THE FEASIBILITY VARIABLE TO DECREASE FROM ITS UPPER
C       BOUND, AS RHO IS MADE LARGE ENOUGH, IS TO REQUIRE THAT
C                T    -1     T    T -1 
C       SIGMA - u (K0)  u - t (LDL )  t  > 0 .
C       THUS, WE NEED SIGMA > UTKIVU + TTCIVT.
C       NOTE: THE ABOVE ALSO KEEPS THE MATRIX INERTIA CORRECT.
C
        CHI = MAX( ONE, ONE + UTKIVU + TTCIVT )
        GOTCHI = .TRUE.
C
C       REDEFINE SIGMA TO REFLECT THE NEW VALUE OF CHI
C       NOTE: THE FEASIBILTY TERM IN THE OBJECTIVE FUNCTION
C             IS RHO*XI + 0.5*CHI*(XI**2)
C
        SIGMA = CHI
C
      ENDIF
C
C                         T      -1
C     GAMMA <--  SIGMA - u * (K0)  * u == SIGMA - UTKOIVU
C
      GAMMA = SIGMA - UTKIVU
C
C
C     ...CONSTRUCT SHURC AND UPDATE SHRCSV.
C
C     THE LEADING K BY K SUBMATRIX OF SHURC <-- SHRCSV
C
      DO J = 1, K
        SHURC(1:K,J) = SHRCSV(1:K,J)
      enddo
C
C
C     SHURC    <-- T IN ROW AND COL K+1 AND SHURC(K+1,K+1) = GAMMA
C     SHRCSV <-- T IN ROW AND COL K+1 AND SHRCSV(K+1,K+1) = GAMMA
C
      IF ( K .GT. 0 ) THEN
        SHURC(1:K,K+1) = T(1:K)
        SHURC(K+1,1:K) = T(1:K)
C
        SHRCSV(1:K,K+1) = T(1:K)
        SHRCSV(K+1,1:K) = T(1:K)
      ENDIF
C
      SHURC(K+1,K+1) = GAMMA
      SHRCSV(K+1,K+1) = GAMMA
C
      IF ( IPC .GT. 100 ) THEN
        CALL MATPRN( SHURC, KMAX, K+1, K+1, 'SHURC', 5, 6 )
      ENDIF
C
 999  CONTINUE
C
      RETURN
      END
