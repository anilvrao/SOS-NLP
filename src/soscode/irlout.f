
      SUBROUTINE IRLOUT(INTARY,RELARY,NLENG,IPU)
C
C ======================================================================
C     IRLOUT===>irlout   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLEPRECISION  (A-H,O-Z)
C
C         PURPOSE: DISPLAY THE INTEGER AND REAL ARRAYS INTARY AND
C                  RELARY
C
C         ARGUMENTS:
C
C         INTARY  INTEGER ARRAY (NLENG)
C
C         RELARY  REAL ARRAY (NLENG)
C
C         NLENG   LENGTH OF RELARY
C
C         IPU     OUTPUT UNIT NO.
C
C-----------------------------------------------------------------------
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
      DIMENSION INTARY(*),RELARY(*)
      CHARACTER(LEN=106) ROWPRT,BLANK
      DATA BLANK(1:106) / ' '/
C-----------------------------------------------------------------------
C
      ROWPRT = BLANK
      ROWPRT(3:3) = '*'
      ROWPRT(106:106) = '*'
C
C         COUNT THE NUMBER OF FULL LENGTH ROWS:
C         5 ELEMENTS PER ROW, 16 COLUMNS WIDE
C
      NROWS = NLENG/5
      IF(NROWS.GT.MAXLYN) THEN
        NROWS = MAX(0,MAXLYN)
        LASTRW = 0
        NSUPRS = NLENG - 5*NROWS
      ELSE
C
C         COUNT THE LENGTH OF THE LAST ROW
C
        LASTRW = NLENG - 5*NROWS
        NSUPRS = 0
      ENDIF
C
C         WRITE OUT THE FULL LENGTH ROWS
C
      II = 0
      DO I = 1,NROWS
        DO J = 1,5
          II = II + 1
          INDX = INTARY(II)
          JSTRT = 17*(J-1)+12
          JSTOP = JSTRT + 5
          WRITE(ROWPRT(JSTRT:JSTOP),'(I6)') INDX
          JSTRT = JSTOP + 1
          ROWPRT(JSTRT:JSTRT) = ':'
          JSTRT = JSTRT + 1
          JSTOP = JSTRT + 9
          WRITE(ROWPRT(JSTRT:JSTOP),'(1PG10.3)') RELARY(INDX)
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
          INDX = INTARY(II)
          JSTRT = 17*(J-1)+12
          JSTOP = JSTRT + 5
          WRITE(ROWPRT(JSTRT:JSTOP),'(I6)') INDX
          JSTRT = JSTOP + 1
          ROWPRT(JSTRT:JSTRT) = ':'
          JSTRT = JSTRT + 1
          JSTOP = JSTRT + 9
          WRITE(ROWPRT(JSTRT:JSTOP),'(1PG10.3)') RELARY(INDX)
      ENDDO
      IF(LASTRW.NE.0) WRITE(IPU,9020) ROWPRT
C
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
