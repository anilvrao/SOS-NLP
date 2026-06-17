

      SUBROUTINE NZROUT(RELARY,NLENG,IPU,TITLE)
C
C ======================================================================
C     NZROUT===>nzrout   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
C         PURPOSE: DISPLAY THE NONZERO VALUES OF THE REAL
C                  ARRAY RELARY
C
C         ARGUMENTS:
C
C         RELARY  REAL ARRAY (NLENG)
C
C         NLENG   LENGTH OF RELARY
C
C         IPU     OUTPUT UNIT NO.
C
C         TITLE   CHARACTER VARIABLE FOR TITLE
C
C-----------------------------------------------------------------------
C
      PARAMETER (ZERO = 0.0D0)
      DIMENSION RELARY(NLENG)
      CHARACTER(LEN=106) ROWPRT,BLANK
      CHARACTER(LEN=*)  TITLE
      DATA BLANK(1:106) / ' '/
C-----------------------------------------------------------------------
C
C
C         COUNT THE NUMBER OF NONZERO ELEMENTS IN RELARY
C
      NZERO = 0
      DO I=1,NLENG
        IF (RELARY(I).NE.ZERO)  NZERO = NZERO + 1
      ENDDO
      IF(NZERO.EQ.0) RETURN
C
C         COUNT THE NUMBER OF FULL LENGTH ROWS:
C         4 ELEMENTS PER ROW, 21 COLUMNS WIDE
C
      NROWS = NZERO/4
C
C         COUNT THE LENGTH OF THE LAST ROW
C
      LASTRW = NZERO - 4*NROWS
C
C         WRITE OUT THE TITLE
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
      LNTITL = LEN(TITLE)
      LNTITL = LNTITL + 10
      ROWPRT(11:LNTITL) = TITLE
      WRITE(IPU,9020) ROWPRT
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
      WRITE(IPU,9020) ROWPRT
C
C         WRITE OUT THE FULL LENGTH ROWS
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
C
      JSCAN = 1
      DO I = 1,NROWS
        JNZ = 0
        DO J = JSCAN,NLENG
          IF(RELARY(J).NE.ZERO) THEN
            JNZ = JNZ + 1
            JSTRT = 23*(JNZ-1)+11
            JSTOP = JSTRT + 5
            WRITE(ROWPRT(JSTRT:JSTOP),'(I6)') J
            JSTRT = JSTOP + 1
            ROWPRT(JSTRT:JSTRT) = ':'
            JSTRT = JSTRT + 1
            JSTOP = JSTRT + 15
            WRITE(ROWPRT(JSTRT:JSTOP),'(1PG16.8)') RELARY(J)
            IF(JNZ.EQ.4) EXIT
          ENDIF
        ENDDO
        JSCAN = J + 1
        WRITE(IPU,9020) ROWPRT
      ENDDO
C
      IF(LASTRW.NE.0) THEN
C
C         WRITE OUT THE LAST ROW 
C
        ROWPRT = BLANK
        ROWPRT(3:3) = '*'
        ROWPRT(106:106) = '*'
C
        JNZ = 0
        DO J = JSCAN,NLENG
          IF(RELARY(J).NE.ZERO) THEN
            JNZ = JNZ + 1
            JSTRT = 23*(JNZ-1)+11
            JSTOP = JSTRT + 5
            WRITE(ROWPRT(JSTRT:JSTOP),'(I6)') J
            JSTRT = JSTOP + 1
            ROWPRT(JSTRT:JSTRT) = ':'
            JSTRT = JSTRT + 1
            JSTOP = JSTRT + 15
            WRITE(ROWPRT(JSTRT:JSTOP),'(1PG16.8)') RELARY(J)
          ENDIF
        ENDDO
        WRITE(IPU,9020) ROWPRT
C
      ENDIF
C
      RETURN
 9020 FORMAT(A106)
      END
