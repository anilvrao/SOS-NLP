
      SUBROUTINE RSTRNG(STRING,RVALUE,IPU)
C
C ======================================================================
C     RSTRNG===>rstrng   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(LEN=*) STRING
      CHARACTER(LEN=80) RITOUT,BLANK
      DATA BLANK(1:80) / ' '/
C
C         COMPUTE LENGTH OF INPUT CHARACTER STRING
C
      LNINPT = LEN(STRING)
      LNIP1 = LNINPT + 1
C
C         IF INPUT STRING IS TOO LONG EXIT
C
      IF(LNINPT.GT.80) THEN
        PRINT *,'OUTPUT STRING TOO LARGE'
        RETURN
      ENDIF
C
      RITOUT = BLANK
      RITOUT = STRING(1:LNINPT)
C
C         LOAD THE REAL VALUE
C
      WRITE(RITOUT(LNIP1:LNIP1+14),'(SP,1PG14.6)') RVALUE
C
      CALL HHADJF(RITOUT(LNIP1:LNIP1+14),' ',' ','L',NS,IERSH)
      CALL HHADJF(RITOUT(LNIP1:80),' ','.','R',NS,IERSH)
C
C         WRITE OUT THE WHOLE STRING
C
      WRITE(IPU,1001) RITOUT
C
 1001 FORMAT(T10,'|',T15,A80,T100,'|')
      RETURN
      END
