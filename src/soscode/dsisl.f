      SUBROUTINE DSISL(A,LDA,N,KPVT,B)
      INTEGER LDA,N,KPVT(*)
      DOUBLE PRECISION A(LDA,*),B(*)
C
C     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY DSIFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA,N)
C                THE OUTPUT FROM DSIFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DSIFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  DSICO  HAS SET RCOND .EQ. 0.0
C        OR  DSIFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DSIFA(A,LDA,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DSISL(A,LDA,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     INTERNAL VARIABLES.
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,TEMP
      INTEGER K,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
   10 CONTINUE
      IF (K .NE. 0) THEN
         IF (KPVT(K) .GE. 0) THEN
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .NE. 1) THEN
               KP = KPVT(K)
               IF (KP .NE. K) THEN
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
               ENDIF
C
C              APPLY THE TRANSFORMATION.
C
               DO I=1,K-1
                 B(I) = B(I) + B(K)*A(I,K)
               ENDDO
            ENDIF
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/A(K,K)
            K = K - 1
            GO TO 10
         ENDIF
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .NE. 2) THEN
               KP = ABS(KPVT(K))
               IF (KP .NE. K - 1) THEN
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
               ENDIF
C
C              APPLY THE TRANSFORMATION.
C
               DO I=1,K-2
                 B(I) = B(I) + B(K)*A(I,K) + B(K-1)*A(I,K-1)
               ENDDO
            ENDIF
C
C           APPLY D INVERSE.
C
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/A(K-1,K)
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            GO TO 10
      ENDIF
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
   90 CONTINUE
      IF (K .LE. N) THEN
         IF (KPVT(K) .GE. 0) THEN
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .NE. 1) THEN
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DOT_PRODUCT(A(1:K-1,K),B(1:K-1))
               KP = KPVT(K)
               IF (KP .NE. K) THEN
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
               ENDIF
            ENDIF
            K = K + 1
         ELSE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .NE. 1) THEN
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DOT_PRODUCT(A(1:K-1,K),B(1:K-1))
               B(K+1) = B(K+1) + DOT_PRODUCT(A(1:K-1,K+1),B(1:K-1))
               KP = ABS(KPVT(K))
               IF (KP .NE. K) THEN
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
               ENDIF
            ENDIF
            K = K + 2
         ENDIF
         GO TO 90
      ENDIF
C
      RETURN
      END
