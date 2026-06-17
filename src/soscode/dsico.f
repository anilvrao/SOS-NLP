      SUBROUTINE DSICO(A,LDA,N,KPVT,RCOND,Z)
      INTEGER LDA,N,KPVT(*)
      DOUBLE PRECISION A(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
C     DSICO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION
C     WITH SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE
C     MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DSIFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DSICO BY DSISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSICO BY DSISL.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSICO BY DSIDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSICO BY DSIDI.
C     TO COMPUTE  INERTIA(A), FOLLOW DSICO BY DSIDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DSIFA
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,EK,T
      DOUBLE PRECISION ANORM,S,YNORM
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS
C
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      DO J = 1, N
         S = 0.0D0
         DO I=1,J
           S = S + ABS(A(I,J))
         ENDDO
         Z(J) = S
         JM1 = J - 1
         DO I = 1, JM1
            Z(I) = Z(I) + ABS(A(I,J))
         ENDDO
      ENDDO
      ANORM = 0.0D0
      DO J = 1, N
         ANORM = MAX(ANORM,Z(J))
      ENDDO
C
C     FACTOR
C
      CALL DSIFA(A,LDA,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = 1.0D0
      Z(1:N) = 0.0D0
      K = N
   60 CONTINUE
      IF (K .NE. 0) THEN
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = ABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .NE. KPS) THEN
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
         ENDIF
         IF (Z(K) .NE. 0.0D0) EK = SIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         DO I=1,K-KS
           Z(I) = Z(I) + Z(K)*A(I,K)
         ENDDO
         IF (KS .NE. 1) THEN
            IF (Z(K-1) .NE. 0.0D0) EK = SIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            DO I=1,K-KS
              Z(I) = Z(I) + Z(K-1)*A(I,K-1)
            ENDDO
         ENDIF
         IF (KS .NE. 2) THEN
            IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
               S = ABS(A(K,K))/ABS(Z(K))
               Z(1:N) = S*Z(1:N)
               EK = S*EK
            ENDIF
            IF (A(K,K) .NE. 0.0D0) THEN
              Z(K) = Z(K)/A(K,K)
            ELSE
              Z(K) = 1.0D0
            ENDIF
         ELSE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
         ENDIF
         K = K - KS
         GO TO 60
      ENDIF
      S = 0.0D0
      DO J=1,N
        S = S + ABS(Z(J))
      ENDDO
      IF (S.GT.0.0D0)  Z(1:N) = Z(1:N)/S
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
  130 CONTINUE
      IF (K .LE. N) THEN
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .NE. 1) THEN
            Z(K) = Z(K) + DOT_PRODUCT(A(1:K-1,K),Z(1:K-1))
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + DOT_PRODUCT(A(1:K-1,K+1),Z(1:K-1))
            KP = ABS(KPVT(K))
            IF (KP .NE. K) THEN
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
            ENDIF
         ENDIF
         K = K + KS
         GO TO 130
      ENDIF
      S = 0.0D0
      DO J=1,N
        S = S + ABS(Z(J))
      ENDDO
      IF (S.GT.0.0D0)  Z(1:N) = Z(1:N)/S
C
      YNORM = 1.0D0
C
C     SOLVE U*D*V = Y
C
      K = N
  170 CONTINUE
      IF (K .NE. 0) THEN
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .NE. KS) THEN
            KP = ABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .NE. KPS) THEN
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
            ENDIF
            DO I=1,K-KS
              Z(I) = Z(I) + Z(K)*A(I,K)
            ENDDO
            IF (KS .EQ. 2) THEN
              DO I=1,K-KS
                Z(I) = Z(I) + Z(K-1)*A(I,K-1)
              ENDDO
            ENDIF
         ENDIF
         IF (KS .NE. 2) THEN
            IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
               S = ABS(A(K,K))/ABS(Z(K))
               Z(1:N) = S*Z(1:N)
               YNORM = S*YNORM
            ENDIF
            IF (A(K,K) .NE. 0.0D0) THEN
              Z(K) = Z(K)/A(K,K)
            ELSE
              Z(K) = 1.0D0
            ENDIF
         ELSE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
         ENDIF
         K = K - KS
         GO TO 170
      ENDIF
      S = 0.0D0
      DO J=1,N
        S = S + ABS(Z(J))
      ENDDO
      IF (S.GT.0.0D0)  Z(1:N) = Z(1:N)/S
      YNORM = YNORM/S
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
  240 CONTINUE
      IF (K .LE. N) THEN
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .NE. 1) THEN
            Z(K) = Z(K) + DOT_PRODUCT(A(1:K-1,K),Z(1:K-1))
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + DOT_PRODUCT(A(1:K-1,K+1),Z(1:K-1))
            KP = ABS(KPVT(K))
            IF (KP .NE. K) THEN
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
            ENDIF
         ENDIF
         K = K + KS
          GO TO 240
      ENDIF
C     MAKE ZNORM = 1.0
      S = 0.0D0
      DO J=1,N
        S = S + ABS(Z(J))
      ENDDO
      IF (S.GT.0.0D0)  Z(1:N) = Z(1:N)/S
      YNORM = YNORM/S
C
      IF (ANORM .NE. 0.0D0) THEN
        RCOND = YNORM/ANORM
      ELSE
        RCOND = 0.0D0
      ENDIF
      RETURN
      END
