

      SUBROUTINE MVPROD(IOPT,M,N,A,MROWA,X,Y)
C
C ======================================================================
C     MVPROD===>mvprod   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  COMPUTE DENSE MATRIX - VECTOR PRODUCT
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
C            MROWA  LEADING DIMENSION OF A
C            X      N-VECTOR
C
C         OUTPUT:
C
C            Y      M-VECTOR
C
C     *******************************************************
C
      PARAMETER (ZERO=0.0D0)
      DIMENSION A(MROWA,*),X(*),Y(*)
C
C     *******************************************************
C
      IF(M.EQ.0) RETURN
      Y(1:M) = ZERO
      IF(N.EQ.0) RETURN
C
      IF(IOPT.EQ.1) THEN
C
C         MULTIPLY BY THE MATRIX A
C
        DO IROW = 1,M
          Y(IROW) = ZERO
          DO JCOL = 1,N
            Y(IROW) = Y(IROW) + A(IROW,JCOL)*X(JCOL)
          ENDDO
        ENDDO
C
      ELSE
C
C         MULTIPLY BY THE MATRIX (A**T)
C
        DO IROW = 1,M
          Y(IROW) = ZERO
          DO JCOL = 1,N
            Y(IROW) = Y(IROW) + A(JCOL,IROW)*X(JCOL)
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
