      SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)
      INTEGER LDA,N,KPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
C     DSIFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION
C     WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW DSIFA BY DSISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSIFA BY DSISL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSIFA BY DSIDI.
C     TO COMPUTE  INERTIA(A) , FOLLOW DSIFA BY DSIDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSIFA BY DSIDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA,N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
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
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DSWAP
C     BLAS2 DSYR
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP
      LOGICAL SWAP
C
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + SQRT(17.0D0))/8.0D0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .LE. 1) THEN
            KPVT(1) = 1
            IF (A(1,1) .EQ. 0.0D0) INFO = 1
C     ......EXIT
            GO TO 200
         ENDIF
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         ABSAKK = ABS(A(K,K))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = 1
         COLMAX = ABS(A(1,K))
         DO J=2,K-1
           IF (ABS(A(J,K)).GT.COLMAX) THEN
             IMAX = J
             COLMAX = ABS(A(J,K))
           ENDIF
         ENDDO
         IF (ABSAKK .GE. ALPHA*COLMAX) THEN
            KSTEP = 1
            SWAP = .FALSE.
            GO TO 90
         ENDIF
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            DO J = IMAXP1, K
               ROWMAX = MAX(ROWMAX,ABS(A(IMAX,J)))
            ENDDO
            DO J=1,IMAX-1
              ROWMAX = MAX(ROWMAX,ABS(A(J,IMAX)))
            ENDDO
            IF (ABS(A(IMAX,IMAX)) .GE. ALPHA*ROWMAX) THEN
               KSTEP = 1
               SWAP = .TRUE.
            ELSEIF (ABSAKK .GE. ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
               KSTEP = 1
               SWAP = .FALSE.
            ELSE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
            ENDIF
   90    CONTINUE
         IF (MAX(ABSAKK,COLMAX) .EQ. 0.0D0) THEN
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
            GO TO 190
         ENDIF
         IF (KSTEP .NE. 2) THEN
C
C           1 X 1 PIVOT BLOCK.
C
            IF (SWAP) THEN
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
               ENDDO
            ENDIF
C
C           PERFORM THE ELIMINATION.
C
            MULK = -1.0D0 / A (K, K)
            CALL DSYR ( 'UPPER', KM1, MULK, A (1, K), 1, A, LDA )
            A(1:KM1,K) = MULK*A(1:KM1,K)
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         ELSE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (SWAP) THEN
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
               ENDDO
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
            ENDIF
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .NE. 0) THEN
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0D0 - AK*AKM1
               DO JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  DO I=1,J
                    A(I,J) = A(I,J) + T*A(I,K)
                  ENDDO
                  T = MULKM1
                  DO I=1,J
                    A(I,J) = A(I,J) + T*A(I,K-1)
                  ENDDO
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
               ENDDO
            ENDIF
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
         ENDIF
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
