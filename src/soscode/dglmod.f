
      SUBROUTINE DGLMOD(NDIM,HMAT,JSTRH,DIAGMD)
C
C ======================================================================
C     DGLMOD===>diagmd   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  MODIFY THE DIAGONAL OF THE HESSIAN MATRIX
C                   BY A MULTIPLE OF THE IDENTITY
C
C
      DIMENSION HMAT(*),JSTRH(NDIM+1)
C
      INCLUDE '../commons/NLPSPR.CMN'
C
      IF  ( QPOPTN.NE.'SPARSE' )  THEN
        NDIM2 = 2*NDIM
        DO I = 1,NDIM
          JSTRHD = I - ((I-1)*(I-NDIM2))/2
          HMAT(JSTRHD) = HMAT(JSTRHD) + DIAGMD
        ENDDO
      ELSE
        DO I = 1,NDIM
          HMAT(JSTRH(I)) = HMAT(JSTRH(I)) + DIAGMD
        ENDDO
      ENDIF
C
      RETURN
      END
