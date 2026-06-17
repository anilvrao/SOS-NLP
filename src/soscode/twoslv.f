
      SUBROUTINE TWOSLV(AMAT,BVEC,XVEC,IER)
C
C ======================================================================
C     TWOSLV===>twoslv   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       SUBROUTINE TO SOLVE A 2 X 2 LINEAR SYSTEM USING GAUSSIAN
C       ELIMINATION WITH COMPLETE PIVOTING
C
C       INPUT
C
C         AMAT  COEFFICIENT MATRIX (2 X 2)
C         BVEC  RIGHT HAND SIDE VECTOR (2)
C
C       OUTPUT
C
C         XVEC  SOLUTION VECTOR (2)
C         IER   ERROR RETURN FLAG
C               = 0   NORMAL RETURN
C               = -1  COEFFICIENT MATRIX IS ZERO
C               = -2  COEFFICIENT MATRIX IS SINGULAR
C
      DIMENSION AMAT(2,2),BVEC(2),XVEC(2)
C
      PARAMETER (ZERO = 0.0D0)
C
      IER = 0
C
C         COMPUTE THE PIVOT ELEMENT
C
      APIV = ZERO
      DO I = 1,2
        DO J = 1,2
          ABSAMT = ABS(AMAT(I,J))
          IF(ABSAMT.GT.APIV) THEN
            APIV = ABSAMT
            MAXROW = I
            MAXCOL = J
          ENDIF
        ENDDO
      ENDDO
C
C         CHECK THAT A NONZERO PIVOT EXISTS
C
      IF(APIV.EQ.ZERO) THEN
        IER = -1
        RETURN
      ENDIF
C
C         SELECT THE APPROPRIATE CASE
C
      IF(MAXROW.EQ.1) THEN
C
        IF(MAXCOL.EQ.1) THEN
C
C         CASE 1:  AMAT(1,1) IS PIVOT
C
          TERM1 = AMAT(2,1) / AMAT(1,1)
          TERM2 = AMAT(2,2) - AMAT(1,2)*TERM1
C
          IF(TERM2.EQ.ZERO) THEN
            IER = -2
            RETURN
          ENDIF
C
          RHSTRM = BVEC(2) - BVEC(1)*TERM1
          XVEC(2) = RHSTRM/TERM2
          XVEC(1) = (BVEC(1) - AMAT(1,2)*XVEC(2))/AMAT(1,1)
C
        ELSE
C
C         CASE 2:  AMAT(1,2) IS PIVOT
C
          TERM1 = AMAT(2,2) / AMAT(1,2)
          TERM2 = AMAT(2,1) - AMAT(1,1)*TERM1
C
          IF(TERM2.EQ.ZERO) THEN
            IER = -2
            RETURN
          ENDIF
C
          RHSTRM = BVEC(2) - BVEC(1)*TERM1
          XVEC(1) = RHSTRM/TERM2
          XVEC(2) = (BVEC(1) - AMAT(1,1)*XVEC(1))/AMAT(1,2)
C
        ENDIF
C
      ELSE
C
        IF(MAXCOL.EQ.1) THEN
C
C         CASE 3:  AMAT(2,1) IS PIVOT
C
          TERM1 = AMAT(1,1) / AMAT(2,1)
          TERM2 = AMAT(1,2) - AMAT(2,2)*TERM1
C
          IF(TERM2.EQ.ZERO) THEN
            IER = -2
            RETURN
          ENDIF
C
          RHSTRM = BVEC(1) - BVEC(2)*TERM1
          XVEC(2) = RHSTRM/TERM2
          XVEC(1) = (BVEC(2) - AMAT(2,2)*XVEC(2))/AMAT(2,1)
C
        ELSE
C
C         CASE 4:  AMAT(2,2) IS PIVOT
C
          TERM1 = AMAT(1,2) / AMAT(2,2)
          TERM2 = AMAT(1,1) - AMAT(2,1)*TERM1
C
          IF(TERM2.EQ.ZERO) THEN
            IER = -2
            RETURN
          ENDIF
C
          RHSTRM = BVEC(1) - BVEC(2)*TERM1
          XVEC(1) = RHSTRM/TERM2
          XVEC(2) = (BVEC(2) - AMAT(2,1)*XVEC(1))/AMAT(2,2)
C
        ENDIF
C
      ENDIF
C
      RETURN
      END
