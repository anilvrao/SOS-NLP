
      SUBROUTINE INTOUT(INTARY,NLENG,IPU)
C
C ======================================================================
C     INTOUT===>intout   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLEPRECISION  (A-H,O-Z)
C
C         PURPOSE DISPLAY THE INTEGER ARRAY INTARY
C
C         ARGUMENTS:
C
C         INTARY  INTEGER ARRAY (NLENG)
C
C         NLENG   LENGTH OF INTARY
C
C         IPU     OUTPUT UNIT NO.
C
C-----------------------------------------------------------------------
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
      DIMENSION INTARY(*)
      CHARACTER(LEN=106) ROWPRT,BLANK
      DATA BLANK(1:106) / ' '/
C-----------------------------------------------------------------------
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
C
C         COUNT THE NUMBER OF FULL LENGTH ROWS:
C         10 ELEMENTS PER ROW, 8 COLUMNS WIDE
C
      NROWS = NLENG/10
      IF(NROWS.GT.MAXLYN) THEN
        NROWS = MAX(0,MAXLYN)
        LASTRW = 0
        NSUPRS = NLENG - 10*NROWS
      ELSE
C
C         COUNT THE LENGTH OF THE LAST ROW
C
        LASTRW = NLENG - 10*NROWS
        NSUPRS = 0
      ENDIF
C
C         WRITE OUT THE FULL LENGTH ROWS
C
      II = 0
      DO I = 1,NROWS
        DO J = 1,10
          II = II + 1
          JSTRT = 8*(J-1)+11
          JSTOP = JSTRT + 7
          WRITE(ROWPRT(JSTRT:JSTOP),'(I8)') INTARY(II)
        ENDDO
        WRITE(IPU,9020) ROWPRT
      ENDDO
C
C         WRITE OUT THE LAST ROW 
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
C
      DO J = 1,LASTRW
        II = II + 1
        JSTRT = 8*(J-1)+11
        JSTOP = JSTRT + 7
        WRITE(ROWPRT(JSTRT:JSTOP),'(I8)') INTARY(II)
      ENDDO
      IF(LASTRW.NE.0) WRITE(IPU,9020) ROWPRT
      IF(NSUPRS.NE.0) THEN
        ROWPRT = BLANK
        ROWPRT(3:3) = '*'
        ROWPRT(106:106) = '*'
        ROWPRT(11:57) = '.....LINE LIMIT EXCEEDED, ELEMENTS SUPPRESSED:'
        WRITE(ROWPRT(58:66),'(I8)') NSUPRS
        WRITE(IPU,9020) ROWPRT
      ENDIF
C
C
      RETURN
 9020 FORMAT(A106)
      END
