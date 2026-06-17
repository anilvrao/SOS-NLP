

      SUBROUTINE SYMMVP(HMAT,IROWH,JSTRH,NCOL,XVEC,YVEC)
C
C ======================================================================
C     SYMMVP===>symmvp   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:
C
C             COMPUTE THE SYMMETRIC MATRIX VECTOR PRODUCT 
C
C                 YVEC = HMAT*XVEC
C
C             WHERE HMAT IS A SYMMETRIC NCOL X NCOL MATRIX.
C
      DIMENSION HMAT(*),IROWH(*),JSTRH(*),XVEC(*),YVEC(*)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
C         ACCUMULATE THE PRODUCT OF HMAT*XVEC IN YVEC
C
      YVEC(1:NCOL) = ZERO
C
      DO K = 1,NCOL
C
        IF(K.EQ.IROWH(JSTRH(K))) THEN
C
C         ADD THE DIAGONAL ELEMENT CONTRIBUTION
C
          YVEC(K) = YVEC(K) + HMAT(JSTRH(K))*XVEC(K)
C  
          JSTRHK = JSTRH(K) + 1
C
        ELSE
C
          JSTRHK = JSTRH(K)
C
        ENDIF        
C
C         COMPUTE CONTRIBUTION FOR THE REMAINDER OF THE ELEMENTS IN
C         COLUMN K OF HMAT
C
        DO I =  JSTRHK,JSTRH(K+1)-1
C
          IROWK = IROWH(I)
          JCOLK = K
          HMATK = HMAT(I)
C
          YVEC(IROWK) = YVEC(IROWK) + HMATK*XVEC(JCOLK)
          YVEC(JCOLK) = YVEC(JCOLK) + HMATK*XVEC(IROWK)
C
        ENDDO
      ENDDO
C
      RETURN
C
      END
