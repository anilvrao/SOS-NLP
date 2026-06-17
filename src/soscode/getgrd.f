C
C
C
      SUBROUTINE GETGRD( NDIM, NX, X, MXNZHM, IROWH,
     .                    JSTRH, HMAT, CVCT, GRAD      )
C
C ======================================================================
C     GETGRD===>getgrd   J.T. BETTS
C ======================================================================
C
C     ROUTINE TO COMPUTE THE GRADIENT OF THE "ORIGINAL"
C     QP PROBLEM
C
      INTEGER I, NDIM, NX, MXNZHM
      INTEGER IROWH(MXNZHM), JSTRH(NDIM+1)
C
      DOUBLE PRECISION HMAT(MXNZHM), X(NX)
      DOUBLE PRECISION GRAD(NDIM), CVCT(NDIM)
      INTEGER ISQPER
      COMMON /PERCOM/ ISQPER(20)
      ISQPER(3) = ISQPER(3) + 1
C
C
C     GRAD <-- HESSIAN*X
C     (NOTE: X IN GETGRD REPRESENTS THE ORIGINAL PROBLEM
C            VARIABLES; I.E. NO SLACKS AND ARTIFICIAL VARS.)
C
      CALL SYMMVP( HMAT, IROWH, JSTRH, NDIM, X, GRAD )
C
C     GRAD <-- GRAD + CVCT == HESSIAN*X + CVCT
C
      DO I = 1, NDIM
        GRAD(I) = GRAD(I) + CVCT(I)
      enddo
C
C
      RETURN
      END
