

      SUBROUTINE NLPSTA(IPF,IPU,IER)
C
C ======================================================================
C     NLPSTA===>nlpsta   J.T. BETTS
C ======================================================================
C
C     PURPOSE:
C        Set/print error status for Sparse NLP routines.
C
C     INPUT:
C        IPF       Print level flag.
C        IPU       FORTRAN output unit number.
C
C     OUTPUT:
C        IER       Error flag.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION  (A-H,O-Z)
C
C  Commons:
C
C             NLPSPR error common.
C
      PARAMETER (NERR=10)
      COMMON /NLPERR/ NUMERR,IERSNP(NERR)
      COMMON /NLPERC/ ERRSUB(NERR),ERRMSG(NERR)
      CHARACTER(LEN=8)  ERRSUB
      CHARACTER(LEN=100)  ERRMSG
C
C  Local:
C
C     ******************************************************************
C
      IER = 0
      IF (NUMERR.GT.0) IER = IERSNP(NUMERR)
C
C             Print information if required.
C
      IF (IPF.GT.0) THEN
        DO I=1,NUMERR
          WRITE(IPU,1001)  ERRSUB(I),IERSNP(I),ERRMSG(I)
        ENDDO
      ENDIF
C
      RETURN
 1001 FORMAT(/T1,'***** NLP ERROR from subroutine ',A8,', IER = ',I4,
     &       /T7,A100,:,5(/T7,A100,:))
      END
