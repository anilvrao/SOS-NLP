
      SUBROUTINE BTURSE(NP,CHITER,LNFLTR,CNDNUM,ALFA,DIAGNL,CMPLMT,
     $    PENMU,ERREQL,GRDLNM,FOBJ,FMRT)
C
C
C ======================================================================
C     BTURSE===>bturse   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  WRITE THE BARNLP TERSE ITERATION SUMMARY 
C
C         INPUT:
C
C           NP      PHASE NUMBER
C                   IF NP .LT. 0, THEN PRINT TITLE, AND 
C                   SET NP .GT. 0
C                   = 1 BARRIER OPTIMIZATION PHASE
C           CHITER  CHARACTER STRING CONTAINING ITERATION NO.
C           LNFLTR  LENGTH OF FILTER
C           CNDNUM  CONDITION NUMBER OF KT SYSTEM
C           ALFA    STEP LENGTH
C           DIAGNL  LEVENBERG PARAMETER
C           CMPLMT  COMPLEMENTARITY CONDITION
C           PENMU   BARRIER PARAMETER
C           ERREQL  EQUALITY CONSTRAINT VIOLATION
C           GRDLNM  GRADIENT OF LAGRIANIAN NORM
C           FOBJ    OBJECTIVE FUNCTION 
C           FMRT    LOG-BARRIER FUNCTION 
C
      INCLUDE '../commons/NLPSPR.CMN'
C
      CHARACTER(LEN=104) RITOUT,BLANK,TITLE,PHASE
      CHARACTER(LEN=6)  CHITER
      DATA BLANK(1:104) / ' '/
C
      RITOUT = BLANK
      TITLE = BLANK
      PHASE = BLANK
C
      NPA = ABS(NP)
      IF(NPA.EQ.1) THEN
        PHASE(43:57) = 'Feasible Point'
      ELSEIF(NPA.EQ.2) THEN
        PHASE(44:56) = 'Optimization'
      ENDIF
C
      TITLE(5:8) = ' It'
      CALL HHADJF(CHITER,' ',' ','R',NSHF,IERSHF)
      IF(CHITER(6:6).NE.')') THEN
        NSHF = 1
        CALL HHADJP(CHITER,' ',' ','L',NSHF,IERSHF)
      ENDIF
      RITOUT(3:8) = CHITER
      TITLE(10:15) = 'Lnfltr'
      WRITE(RITOUT(10:14),'(I5)') LNFLTR
      TITLE(17:23) = 'KT Cond'
      WRITE(RITOUT(17:23),'(1PE7.1)') CNDNUM
      TITLE(25:32) = '  Step '
      WRITE(RITOUT(25:32),'(1PE7.1)') ALFA
      TITLE(34:41) = 'Levnbrg'
      WRITE(RITOUT(34:41),'(1PE7.1)') DIAGNL
      TITLE(43:50) = ' Cmplmt'
      WRITE(RITOUT(43:50),'(1PE7.1)') CMPLMT
      TITLE(52:59) = 'Barrier'
      WRITE(RITOUT(52:59),'(1PE7.1)') PENMU
      TITLE(61:68) = ' Violtn'
      WRITE(RITOUT(61:68),'(1PE7.1)') ERREQL
      TITLE(70:77) = ' |Grdl|'
      WRITE(RITOUT(70:77),'(1PE7.1)') GRDLNM
      TITLE(79:87) = '   Obj  '
      WRITE(RITOUT(79:87),'(1PE8.1)') FOBJ
      TITLE(89:103) = 'Log-Barrier Fn'
      WRITE(RITOUT(89:103),'(1PE14.6)') FMRT
C         
      IF(NP.LT.0) THEN
        NP = -NP
        WRITE(IPUNLP,1002) PHASE
        WRITE(IPUNLP,1001) TITLE
      ENDIF
C
      WRITE(IPUNLP,1001) RITOUT
C
 1001 FORMAT(1X,A104)
 1002 FORMAT(/1X,A104)
C 
      RETURN
      END
