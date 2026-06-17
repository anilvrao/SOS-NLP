
      SUBROUTINE MVPSPR(IOPT,M,N,A,IROW,JCOLST,X,Y)
C
C ======================================================================
C     MVPSPR===>mvpspr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  COMPUTE SPARSE/DENSE MATRIX - VECTOR PRODUCT
C
C         INPUT:
C
C            IOPT   INTEGER OPTION CODE
C                   = 1  Y = A*X
C                   = 11 Y = (A**T)*X
C            M      NUMBER OF ELEMENTS IN Y
C            N      NUMBER OF ELEMENTS IN X
C            A      M X N MATRIX STORED AS A VECTOR OF LENGTH NZ 
C                   FOR IOPT = 11, A**T HAS DIMENSION M X N.
C            IROW   INTEGER ARRAY OF LENGTH NZ, CONTAINING ROW INDEX
C                   OF NONZERO ELEMENT IN A 
C            JCOLST INTEGER ARRAY OF LENGTH (N+1), CONTAINING COLUMN
C                   START (NZ = JCOLST(N+1)-1) 
C            X      N-VECTOR
C
C         OUTPUT:
C
C            Y      M-VECTOR
C
C     *******************************************************
C
      PARAMETER (ZERO=0.0D0)
      DIMENSION A(*),X(*),Y(*),IROW(*),JCOLST(*)
C
      INCLUDE '../commons/NLPSPR.CMN'
C
C     *******************************************************
C
      IF(M.EQ.0) RETURN
      Y(1:M) = ZERO
      IF(N.EQ.0) RETURN
C
      IF(QPOPTN.EQ.'SPARSE') THEN
C
        IF(IOPT.EQ.1) THEN
C
C         MULTIPLY BY THE MATRIX A
C
C           LOOP OVER COLUMNS
C
          DO J = 1,N
            DO K=JCOLST(J),JCOLST(J+1)-1
              IR = IROW(K)
              Y(IR) = Y(IR) + A(K)*X(J)
            ENDDO
          ENDDO
C
        ELSE
C
C         MULTIPLY BY THE MATRIX (A**T)
C
          DO J = 1,M
            YJ = ZERO
            DO K=JCOLST(J),JCOLST(J+1)-1
              IR = IROW(K)
              YJ = YJ + A(K)*X(IR)
            ENDDO
            Y(J) = YJ
          ENDDO
C
        ENDIF
C
      ELSE
C
        MROWA = JCOLST(1)
        CALL MVPROD(IOPT,M,N,A,MROWA,X,Y)
C
      ENDIF
C
      RETURN
      END
