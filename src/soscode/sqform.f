




      SUBROUTINE SQFORM(HMAT,IROWH,JSTRH,CVEC,XVEC,WORK,NVAR,QUAD)
C
C ======================================================================
C     SQFORM===>sqform   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE: EVALUATE THE QUADRATIC FORM
C
C         QUAD = .5*(XVEC**T)HMAT*XVEC + (CVEC**T)XVEC
C
C         INPUT:
C
C           HMAT   HESSIAN MATRIX NONZEROS (LOWER TRIANGLE) 
C           IROWH  INTEGER ROW INDEX ARRAY
C           JSTRH  INTEGER COLUMN START ARRAY (NVAR+1)
C           CVEC   OBJECTIVE FUNCTION LINEAR TERM (NVAR)
C           XVEC   VARIABLES (NVAR)
C           WORK   WORK ARRAY (NVAR)
C           NVAR   NUMBER OF VARIABLES
C
C         OUTPUT:
C 
C           QUAD   VALUE OF THE QUADRATIC FORM
C
      DIMENSION HMAT(*),IROWH(*),JSTRH(NVAR+1),CVEC(NVAR),
     $    XVEC(NVAR),WORK(NVAR)
C
      PARAMETER (ZERO=0.0D0,POINT5=5.0D-1)
C
      CALL SYMMVP(HMAT,IROWH,JSTRH,NVAR,XVEC,WORK)
      QUAD = ZERO
      DO I = 1,NVAR
        WORK(I) = POINT5*WORK(I) + CVEC(I)
        QUAD = QUAD + WORK(I)*XVEC(I)
      ENDDO
C
      RETURN
      END
