      SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      DOUBLE PRECISION X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
C
C
C     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X
C     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY DSVDC.
C
C         LDX       INTEGER.
C                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C         N         INTEGER.
C                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C         P         INTEGER.
C                   P IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C         LDU       INTEGER.
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
C                   (SEE BELOW).
C
C         WORK      DOUBLE PRECISION(N).
C                   WORK IS A SCRATCH ARRAY.
C
C         JOB       INTEGER.
C                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
C                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
C                   WITH THE FOLLOWING MEANING
C
C                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
C                                  VECTORS.
C                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
C                                  IN U.
C                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
C                                  VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         DOUBLE PRECISION(P).
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF
C                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2
C                                   THEN K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WITH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
C                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
C                   VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL DROT
C     BLAS DSWAP,DROTG
C     FORTRAN MOD
C
C     INTERNAL VARIABLES
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     *        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      DOUBLE PRECISION T
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,SCALE,SHIFT,SL,SM,SN,
     *                 SMM1,T1,TEST,ZTEST,FACT
      LOGICAL WANTU,WANTV
C
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN(N-1,P)
      NRT = MAX(0,MIN(P-2,N))
      LU = MAX(NCT,NRT)
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .LE. NCT) THEN
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            CNORM = 0.0D0
            DO I=L,N
              CNORM = CNORM + X(I,L)**2
            ENDDO
            S(L) = SQRT(CNORM)
            IF (S(L) .NE. 0.0D0) THEN
               IF (X(L,L) .NE. 0.0D0) S(L) = SIGN(S(L),X(L,L))
               FACT = 1.0D0/S(L)
               X(L:N,L) = FACT*X(L:N,L)
               X(L,L) = 1.0D0 + X(L,L)
            ENDIF
            S(L) = -S(L)
         ENDIF
         DO J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .NE. 0.0D0) THEN
C
C              APPLY THE TRANSFORMATION.
C
               T = -DOT_PRODUCT(X(L:N,L),X(L:N,J))/X(L,L)
               DO I=L,N
                 X(I,J) = X(I,J) + T*X(I,L)
               ENDDO
            ENDIF
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
         ENDDO
         IF (WANTU .AND. L .LE. NCT) THEN
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO I = L, N
               U(I,L) = X(I,L)
            ENDDO
         ENDIF
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            CNORM = 0.0D0
            DO I=LP1,P
              CNORM = CNORM + E(I)**2
            ENDDO
            E(L) = SQRT(CNORM)
            IF (E(L) .NE. 0.0D0) THEN
               IF (E(LP1) .NE. 0.0D0) E(L) = SIGN(E(L),E(LP1))
               FACT = 1.0D0/E(L)
               E(LP1:P) = FACT*E(LP1:P)
               E(LP1) = 1.0D0 + E(LP1)
            ENDIF
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               WORK(LP1:N) = 0.0D0
               DO J = LP1, P
                  DO I=LP1,N
                    WORK(I) = WORK(I) + E(J)*X(I,J)
                  ENDDO
               ENDDO
               DO J = LP1, P
                  T = -E(J)/E(LP1)
                  DO I=LP1,N
                    X(I,J) = X(I,J) + T*WORK(I)
                  ENDDO
               ENDDO
  120       CONTINUE
            IF (WANTV) THEN
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO I = LP1, P
                  V(I,L) = E(I)
               ENDDO
            ENDIF
  150    CONTINUE
  160 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
C
C     IF REQUIRED, GENERATE U.
C
      IF (WANTU) THEN
         DO J = NCTP1, NCU
            DO I = 1, N
               U(I,J) = 0.0D0
            ENDDO
            U(J,J) = 1.0D0
         ENDDO
         DO LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) THEN
               U(1:N,L) = 0.0D0
               U(L,L) = 1.0D0
               CYCLE
            ENDIF
            LP1 = L + 1
            DO J = LP1, NCU
               T = -DOT_PRODUCT(U(L:N,L),U(L:N,J))/U(L,L)
               DO I=L,N
                 U(I,J) = U(I,J) + T*U(I,L)
               ENDDO
            ENDDO
            U(L:N,L) = -U(L:N,L)
            U(L,L) = 1.0D0 + U(L,L)
            LM1 = L - 1
            IF (LM1 .GE. 1)  U(1:LM1,L) = 0.0D0
         ENDDO
      ENDIF
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (WANTV) THEN
         DO LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO J = LP1, P
                  T = -DOT_PRODUCT(V(LP1:P,L),V(LP1:P,J))/V(LP1,L)
                  DO I=LP1,P
                    V(I,J) = V(I,J) + T*V(I,L)
                  ENDDO
               ENDDO
  320       CONTINUE
            V(1:P,L) = 0.0D0
            V(L,L) = 1.0D0
         ENDDO
      ENDIF
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      MAXIT = 30 * M
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .GE. MAXIT) THEN
            INFO = M
C     ......EXIT
            GO TO 620
         ENDIF
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) EXIT
            TEST = ABS(S(L)) + ABS(S(L+1))
            ZTEST = TEST + ABS(E(L))
            IF (ZTEST .EQ. TEST) THEN
               E(L) = 0.0D0
C        ......EXIT
               EXIT
            ENDIF
         ENDDO
         IF (L .EQ. M - 1) THEN
            KASE = 4
            GO TO 480
         ENDIF
            LP1 = L + 1
            MP1 = M + 1
            DO LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) EXIT
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + ABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + ABS(E(LS-1))
               ZTEST = TEST + ABS(S(LS))
               IF (ZTEST .EQ. TEST) THEN
                  S(LS) = 0.0D0
C           ......EXIT
                  EXIT
               ENDIF
            ENDDO
            IF (LS .EQ. L) THEN
               KASE = 3
            ELSEIF (LS .EQ. M) THEN
               KASE = 1
            ELSE
               KASE = 2
               L = LS
            ENDIF
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         IF (KASE.EQ.2) THEN
           GO TO 520
         ELSEIF (KASE.EQ.3) THEN
           GO TO 540
         ELSEIF (KASE.EQ.4) THEN
           GO TO 570
         ENDIF
C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .NE. L) THEN
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
               ENDIF
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
            ENDDO
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO K = L, M
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
            ENDDO
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = MAX(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),
     *                    ABS(S(L)),ABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .NE. 0.0D0 .OR. C .NE. 0.0D0) THEN
               SHIFT = SQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
            ENDIF
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     *            CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
            ENDDO
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .LT. 0.0D0) THEN
               S(L) = -S(L)
               IF (WANTV) V(1:P,L) = -V(1:P,L)
            ENDIF
C
C           ORDER THE SINGULAR VALUE.
C
  590       CONTINUE
            IF (L .EQ. MM) GO TO 600
C           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     *            CALL DSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     *            CALL DSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END
