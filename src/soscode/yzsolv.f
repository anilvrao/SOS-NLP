      SUBROUTINE YZSOLV( NEQNS, K, IPC, IPU, KMAX, SOLVTL, SHURC,
     $                   KPVT, MXNZUM, IROWU, JSTRU, UMAT,
     $                   AKIVF0, F0, W, NWKKT, WORKKT,  
     $                   WRK, WRKBIG, DRPFES, INSHUR, RHO,  
     $                   RHOMAX, MINEQL, Z, Y, IER                  )
C
C ======================================================================
C     YZSOLV===>yzsolv   J.T. BETTS
C ======================================================================
C
C                                     T          
C     ROUTINE TO SOLVE SHURC*z = w - U AKIVF0 FOR z  (WHERE THE
C                                                       T
C     SCHUR-COMPLEMENT MATRIX SHURC IS FACTORED INTO RDR )
C     AND THEN SOLVE K0*y = f0 - U*z FOR y. THE VECTORS y AND z
C     ARE THEN UNSCRAMBLED BY SUBROUTINE DCODYZ TO OBTAIN
C     THE STEP q FROM X0, THE LAGRANGE MULTIPLIERS PI FOR
C     THE GENERAL LINEAR CONSTRAINTS AND THE SUBSET OF LAGRANGE
C     MULTIPLIERS IN (ALMBDA)FX CORRESPONDING TO VARIABLES SET
C     TO THEIR BOUNDS BY SCHUR-COMPLEMENT UPDATES.
C
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D0 )
C
      LOGICAL DRPFES
C
      INTEGER IER, IOPT, K, KMAX, I, NEQNS
      INTEGER IPC, IPU, FEASDX, JER
      INTEGER NWKKT, MXNZUM, MINEQL, NEEDS
      INTEGER KPVT(KMAX), INSHUR(KMAX)
      INTEGER IROWU(MXNZUM), JSTRU(KMAX+1)
C
      DOUBLE PRECISION WRK(KMAX), AKIVF0(NEQNS), UMAT(MXNZUM)
      DOUBLE PRECISION W(KMAX), SHURC(KMAX,KMAX)
      DOUBLE PRECISION Z(KMAX), F0(NEQNS), WRKBIG(NEQNS)
      DOUBLE PRECISION Y(NEQNS), WORKKT(NWKKT)
      DOUBLE PRECISION RHO, RHOMAX
      DOUBLE PRECISION DELTA, SOLVTL
C
      INTEGER ISQPER
      COMMON /PERCOM/ ISQPER(20)
C
C
      IER = 0
C
C                  T         T 
C     --- SOLVE RDR z = w - U AKIVF0
C
C              T
C     WRK <-- U AKIVF0
C
      IOPT = 11
      CALL MVPSPR( IOPT, K, NEQNS, UMAT, IROWU, JSTRU,
     $             AKIVF0, WRK )
C
C                             T
C     WRK <-- w - WRK == w - U AKIVF0
C
      DO I = 1, K
        WRK(I) = W(I) - WRK(I)
      enddo
C
C
C                       T
C     z <-- WRK == w - U AKIVF0
C     
      Z(1:K) = WRK(1:K)
C
C                  T          T
C     NOW SOLVE RDR *z = w - U AKIVF0 FOR z.
C
      IF ( DRPFES ) THEN
C       THE FEASIBILITY VAR WAS DROPPED AS PART OF THE NEW
C       SCHUR-COMPLEMENT UPDATE, RHO MUST BE MADE POSITIVE ENOUGH SUCH
C       THAT THE FEASIBILITY VAR GOES LOWER FROM IT UPPER BOUND (1.0).
C       THAT IS, THE QP STATIONARY POINT MUST BE LESS THAN 1.0 FOR XI.
C
        FEASDX = 0
        DO I = 1, K
          IF ( INSHUR(I) .EQ. MINEQL + 1 ) THEN
            FEASDX = I
          ENDIF
        enddo
C
        IF ( FEASDX .EQ. 0 ) THEN
          IER = 2
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' YZSOLV: ERROR, EXPECTED THE FEAS VAR',
     $                   ' AS AN UPDATE'
          ENDIF
          GO TO 999
        ENDIF
C
C       START OF RHO INCREASE LOOP.
 35     CONTINUE
C
          CALL DSISL( SHURC, KMAX, K, KPVT, Z )
C
          IF ( Z(FEASDX) .LE. ZERO ) THEN
C           RHO WILL BECOME 100 TIMES ITS VALUE, ACCOUNT FOR THIS.
C
            IF ( IPC .GT. 0 ) THEN
              WRITE(IPU,*) ' RHO INCREASED FOR NEGATIVE FEAS. STEP  -- B
     $AD SCALING OR CONDITIONING'
            ENDIF
            WRK(FEASDX) = WRK(FEASDX) + 99.0D0*RHO
            W(FEASDX) = W(FEASDX) +  99.0D0*RHO
            RHO = 1.0D2*RHO
C
            IF ( RHO .LE. RHOMAX ) THEN
C             SOLVE AGAIN WITH THE LARGER RHO.
C
              Z(1:K) = WRK(1:K)
              GO TO 35
            ELSE
             IER = 3
             IF ( IPC .GT. 0 ) THEN
               WRITE(IPU,*) ' YZSOLV: CANNOT GET NEG STEP ON FEAS VAR'
               WRITE(IPU,*) ' WHEN LEAVING ITS UB, FOR ANY RHO < RHOMAX'
             ENDIF
             GO TO 999
            ENDIF
          ENDIF
C       END OF RHO INCREASE LOOP. 
C
      ELSE
        CALL DSISL( SHURC, KMAX, K, KPVT, Z )
      ENDIF
C
C
C     ---SOLVE K0*Y = f0 - U*z FOR y
C
C     Y <-- U*z
C
      IOPT = 1
      CALL MVPSPR( IOPT, NEQNS, K, UMAT, IROWU, JSTRU, Z, Y )
C
C     y <-- f0 - y == f0 - U*z
C
      DO I = 1, NEQNS
        Y(I) = F0(I) - Y(I)
      enddo
C
C               -1           -1 
C     y <-- (K0)  * y == (K0)  (f0 - U*z)
C
      ISQPER(10) = ISQPER(10) + 1
      WRKBIG(1:NEQNS) = Y(1:NEQNS)
C
C         SOLVE LINEAR SYSTEM
C
      CALL REFITR(WRKBIG,Y,NEQNS,WORKKT,NWKKT,NEEDS,JER)
C
C
      IF ( JER .NE. 0 ) THEN
C       ALLOW EXCEEDING REFITR'S SOLUTION ERROR GOAL, BUT DON'T
C       ACCEPT AN ERROR GREATER THAN SOLVTL 
csolvtl        IF ( JER .NE. -701 .OR. DELTA .GT. SOLVTL ) THEN
        IF ( JER .NE. -701  ) THEN
          IER = 1
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' YZSOLV: ERROR, JER IN HDSLSL=', JER
          ENDIF
          GO TO 999
        ENDIF
      ENDIF
C
C
 999  CONTINUE
C
      RETURN
      END
