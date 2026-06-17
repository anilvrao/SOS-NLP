C
C
C
      SUBROUTINE GETOBJ( NDIM, X, MXNZHM, IROWH,
     .                   JSTRH, HMAT, CVCT, WRK, OBJ   )
C
C ======================================================================
C     GETOBJ===>getobj   J.T. BETTS
C ======================================================================
C
C     ROUTINE TO COMPUTE THE OBJECTIVE FUNCTION OF THE "ORIGINAL"
C     QP PROBLEM
C
      INTEGER NDIM, MXNZHM
      INTEGER IROWH(MXNZHM), JSTRH(NDIM+1)
C
      DOUBLE PRECISION HMAT(MXNZHM), X(NDIM)
      DOUBLE PRECISION WRK(NDIM), CVCT(NDIM)
      DOUBLE PRECISION OBJ
C
C
C     WRK <-- HESSIAN*X
C     (NOTE: X IN GETOBJ REPRESENTS THE ORIGINAL PROBLEM
C            VARIABLES; I.E. NO SLACKS AND ARTIFICIAL VARS.)
C
      CALL SYMMVP( HMAT, IROWH, JSTRH, NDIM, X, WRK )
C
C                    T        T          T                T
C     OBJ <-- 1/2 WRK X + CVCT X == 1/2 X HESSIAN*X + CVCT X
C
      OBJ = 0.5D0*DOT_PRODUCT(WRK(1:NDIM),X(1:NDIM))
     .          + DOT_PRODUCT(CVCT(1:NDIM),X(1:NDIM))
C
C
      RETURN
      END
