
      SUBROUTINE NPESET(SUBNAM,IERSUB,SUBMSG,NMSG)
C
C ======================================================================
C     NPESET===>npeset   J.T. BETTS
C ======================================================================
C
C     PURPOSE:
C        Set/print error message for Sparse NLP routines.
C
C     INPUT:
C        SUBNAM    Name of subroutine where error occured.
C        IERSUB    Number of error to be returned to calling program.
C        SUBMSG    Error message to be saved/printed.
C        NMSG      Number of lines in message.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
C  Arguments:
      CHARACTER(LEN=8)  SUBNAM
      CHARACTER(LEN=100)  SUBMSG(*)
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
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
C  Local:
C
C     ******************************************************************
C
      IF (NMSG.LE.0) THEN
        NUMERR = 0
        RETURN
      ENDIF
C
      NUMERR = NUMERR + 1
      NUMERR = MIN(NUMERR,NERR)
C
C             Save off error information.
C
      IERSNP(NUMERR) = IERSUB
      ERRSUB(NUMERR) = SUBNAM
      ERRMSG(NUMERR) = SUBMSG(1)
C
C             Print information if required.
C
      IF (IOFLAG.GT.0.OR.IERSUB.EQ.-108) THEN
        WRITE(IPUNLP,1001)  SUBNAM,IERSUB,(SUBMSG(I),I=1,NMSG)
      ENDIF
C
      RETURN
 1001 FORMAT(/T1,'***** NLP ERROR from subroutine ',A8,', IER = ',I5,
     &       /T7,A100,:,5(/T7,A100,:))
      END
