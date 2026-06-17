
      SUBROUTINE MVPSRC(IOPT,M,N,NZ,A,IROW,JCOL,X,Y)
C
C     PURPOSE:
C        Compute sparse matrix - vector product.
C
C     INPUT:
C        IOPT      Integer option code.
C                   = 1    Y = A*X
C                   .NE.1  Y = (A**T)*X
C        M         Number of elements in Y.
C        N         Number of elements in X.
C        NZ        Number of nonzero elements in A.
C        A         M x N matrix stored as a vector of length NZ. 
C                  For IOPT .NE. 1, A**T has dimension M x N.
C        IROW      Integer array of length NZ, containing row index
C                  of nonzero element in A.
C        JCOL      Integer array of length NZ, containing column
C                  index of nonzero element in A.
C        X         N-vector.
C
C     OUTPUT:
C        Y         M-vector.
C
C     *******************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      DIMENSION A(NZ),X(N),Y(M),IROW(NZ),JCOL(NZ)
C
C     *******************************************************
C
      IF(IOPT.EQ.1) THEN
C
C             Multiply by the matrix A.
C
        Y(1:M) = ZERO
        DO K = 1,NZ
          I = IROW(K)
          J = JCOL(K)
          Y(I) = Y(I) + A(K)*X(J)
        ENDDO
C
      ELSE
C
C             Multiply by the matrix (A**T).
C
        Y(1:M) = ZERO
        DO K = 1,NZ
          I = IROW(K)
          J = JCOL(K)
          Y(J) = Y(J) + A(K)*X(I)
        ENDDO
C
      ENDIF
C
      RETURN
      END
