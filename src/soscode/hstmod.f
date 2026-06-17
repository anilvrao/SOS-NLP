
      SUBROUTINE HSTMOD( K, NX, NEQNS, NDIM, KMAX, MINEQL, IPC, IPU,
     .                   GOTRHO, RHO, CHI, X0, IFREE0, BL, BU, INSHUR, 
     .                   ITPSHR, F0, GRAD0, W, AKIVF0, SOLVTL, NWKKT, 
     .                   WORKKT, IER )
C
C ======================================================================
C     HSTMOD===>hstmod   J.T. BETTS
C ======================================================================
C
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D0 )
C
      LOGICAL GOTRHO
C
      INTEGER I, KK, NX, JVAR, IPC, IPU
      INTEGER NRHS, JER, IER, NDIM, NVFRE0
      INTEGER MINEQL, K, KMAX, NEQNS, NWKKT, NEEDS
C
      INTEGER IFREE0(NX), INSHUR(KMAX), ITPSHR(KMAX)
C
      DOUBLE PRECISION F0(NEQNS), GRAD0(NDIM), RHO, CHI
      DOUBLE PRECISION AKIVF0(NEQNS), X0(NX), BL(NX), BU(NX)
      DOUBLE PRECISION WORKKT(NWKKT)
      DOUBLE PRECISION W(KMAX)
      DOUBLE PRECISION DELTA, SOLVTL
      COMMON /PERCOM/ ISQPER(20)
      INTEGER ISQPER
C
C
      IER = 0
C
C     REDEFINE F0 BASED ON THE NEW GRAD0 AT THE NEW X0.
C
      DO I = 1, NX
C
C       NVFRE0 = { <0, IF VAR I WAS FIXED AT PRIOR KT FACTORIZATION.
C                 { J>0, IF VAR I WAS J-TH FREE VARIABLE
C                        AT PRIOR KT FACTORIZATION.
C
        NVFRE0 = IFREE0(I)
C
        IF ( I .GT. MINEQL + 1 .AND. NVFRE0 .GT. 0 ) THEN
C        
C         AN ORIGINAL PROBLEM VARIABLE.
C
          F0(NVFRE0) = GRAD0( I - (MINEQL+1) )
C
        ELSEIF ( I .EQ. MINEQL + 1 .AND. NVFRE0 .GT. 0 ) THEN
C
C         FEASIBILITY VARIABLE.
C
          IF ( .NOT. GOTRHO ) THEN
            IER = 1
            IF ( IPC .GT. 0 ) THEN
              WRITE(IPU,*) ' HSTMOD:ERROR, GOTRHO=F BUT FEAS VAR IN KT'
            ENDIF
          ENDIF
C
          F0(NVFRE0) = RHO + CHI*X0(MINEQL+1)
C
        ENDIF
C
      enddo
C
C                         -1
C     REDFINE AKIVF0 == K0  F0.
C
      ISQPER(10) = ISQPER(10) + 1
C
C         SOLVE LINEAR SYSTEM
C
      CALL REFITR(F0,AKIVF0,NEQNS,WORKKT,NWKKT,NEEDS,JER)
C
      IF ( JER .NE. 0 ) THEN
C       ALLOW EXCEEDING REFITR'S SOLUTION ERROR GOAL, BUT DON'T
C       ACCEPT AN ERROR GREATER THAN SOLVTL 
csolvtl        IF ( JER .NE. -701 .OR. DELTA .GT. SOLVTL ) THEN
csolvtl    check solvtl,  ???
        IF ( JER .NE. -701  ) THEN
          IER = 2
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' HSTMOD: ERROR, JER IN HDSLSL=', JER
          ENDIF
        ENDIF
      ENDIF
C
C
C     REDEFINE W BASED ON THE NEW X0 AND THE NEW
C     GRAD0 AT THE NEW X0.
C
      DO KK = 1, K
        JVAR = ABS(INSHUR(KK))
C
        IF ( ITPSHR(KK) .GT. 0 ) THEN
C
C         A VARIABLE FREED AT X0.
C
          IF ( JVAR .GT. MINEQL + 1 ) THEN
C        
C           AN ORIGINAL PROBLEM VARIABLE FREED AT UPDATE KK.
C
            W(KK) = GRAD0( JVAR - (MINEQL+1) )
C 
          ELSEIF ( JVAR .EQ. MINEQL + 1 ) THEN
C
C           FEASIBILITY VARIABLE FREED AT UPDATE KK.
C
            W(KK) = RHO + CHI*X0(MINEQL+1)
C
          ENDIF
C
        ELSEIF( ITPSHR(KK) .EQ. -1 ) THEN
C
C         A VARIABLE SET TO ITS LOWER BOUND AT UPDATE KK.
C
          W(KK) = X0(JVAR) - BL(JVAR)
C
        ELSEIF( ITPSHR(KK) .EQ. -2 ) THEN
C
C         A VARIABLE SET TO ITS LOWER BOUND AT UPDATE KK.
C
          W(KK) = X0(JVAR) - BU(JVAR)
C
        ENDIF

      enddo
C
C
      RETURN
      END
