
      SUBROUTINE SYMMAT(HMAT,NCOL,IPU)
C
C ======================================================================
C     SYMMAT===>symmat   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         THIS ROUTINE WRITES THE DENSE NCOL X NCOL SYMMETRIC  
C         (LOWER TRIANGULAR PACKED) MATRIX HMAT ON OUTPUT UNIT IPU. 
C
      DIMENSION  HMAT(*)
      CHARACTER(LEN=106) ROWPRT,BLANK
      DATA BLANK(1:106) / ' '/
C-----------------------------------------------------------------------
C
C
      IF (NCOL.LE.0) RETURN
      NCOL2 = 2*NCOL
C
      JUP = 0
 110  CONTINUE
C
C         DEFINE LOWER AND UPPER COLUMN INDICES
C
      JLW = JUP + 1
      JUP = MIN(JLW+4,NCOL)
C
C         COLUMN INDEX PRINT
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
C
      II = 0
      DO I=JLW,JUP
        II = II + 1
        JSTRT = 16*(II-1)+18
        JSTOP = JSTRT + 3
        WRITE(ROWPRT(JSTRT:JSTOP),'(I3)') I
      ENDDO
      WRITE(IPU,1001) ROWPRT
C
C         WRITE ARRAY
C
      DO I = 1, NCOL
        ROWPRT = BLANK
        ROWPRT(3:3) = '*'
        ROWPRT(106:106) = '*'
C
        WRITE(ROWPRT(7:10),'(I3)') I
        II = 0
        DO J=JLW,JUP
          II = II + 1
          JSTRT = 16*(II-1)+11
          JSTOP = JSTRT + 15
          IJ = I - ((J-1)*(J-NCOL2))/2
          IF(J.GT.I) CYCLE
          WRITE(ROWPRT(JSTRT:JSTOP),'(1PG16.8)') HMAT(IJ)
        ENDDO
C
        WRITE(IPU,1001) ROWPRT
C
      ENDDO
C
      IF (JUP.NE.NCOL) THEN
        WRITE (IPU,1002)
        GO TO 110
      ENDIF
C
      RETURN
 1001 FORMAT(A106)
 1002 FORMAT (T3,'*',T106,'*',/T3,'*',T8,'- - - - - - - - - - - - - - ',
     +       '- - - - - - - - - - - - - - - - - - - - - - - - - - - - ',
     +       T106,'*',/T3,'*',T106,'*')
      END
