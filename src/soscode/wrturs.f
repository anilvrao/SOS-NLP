
      SUBROUTINE WRTURS(NP,CHITER,NQPITR,NDOF,ALFA,PNORM,SQPMF,NKT,
     $    PGNORM,VIOL,PENLTY,DIAGNL,CNDNUM,IPU) 
C
C ======================================================================
C     WRTURS===>wrturs   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  WRITE THE NLPSPR TERSE ITERATION SUMMARY 
C
C         INPUT:
C
C           NP      PHASE NUMBER
C                   IF NP .LT. 0, THEN PRINT TITLE, AND 
C                   SET NP .GT. 0
C                   = 1 FEASIBILITY PHASE
C                   = 2 SQP OPTIMIZATION PHASE
C           CHITER  CHARACTER STRING CONTAINING ITERATION NO.
C           NQPITR  NUMBER OF QP ITERATIONS
C           NDOF    NUMBER OF DEGREES OF FREEDOM
C           ALFA    STEP LENGTH
C           PNORM   SEARCH DIRECTION MAGNITUDE
C           SQPMF   MERIT FUNCTION 
C           NKT     NUMBER OF KT FACTORIZATIONS
C           PGNORM  PROJECTED GRADIENT NORM
C           VIOL    CONSTRAINT VIOLATION
C           PENLTY  MERIT FUNCTION PENALTY WEIGHT
C           DIAGNL  LEVENBERG PARAMETER
C           CNDNUM  CONDITION NUMBER OF KT SYSTEM
C           IPU     OUTPUT UNIT NUMBER
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
        TITLE(5:8) = ' It'
        CALL HHADJF(CHITER,' ',' ','R',NSHF,IERSHF)
        IF(CHITER(6:6).NE.')') THEN
          NSHF = 1
          CALL HHADJP(CHITER,' ',' ','L',NSHF,IERSHF)
        ENDIF
        RITOUT(3:8) = CHITER
        TITLE(11:13) = 'Qit'
        WRITE(RITOUT(9:13),'(I5)') NQPITR
        TITLE(15:17) = 'Nkt'
        WRITE(RITOUT(15:17),'(I3)') NKT
        TITLE(20:24) = 'Ndof'
        WRITE(RITOUT(19:24),'(I5)') NDOF
        IF(QPOPTN.EQ.'SPARSE') THEN
          TITLE(27:34) = 'KT Cond'
        ELSE
          TITLE(27:34) = 'Cond(G)'
        ENDIF
        WRITE(RITOUT(27:34),'(1PE7.1)') CNDNUM
        TITLE(36:43) = '  Step '
        WRITE(RITOUT(36:43),'(1PE7.1)') ALFA
        TITLE(45:52) = ' Norm p'
        WRITE(RITOUT(45:52),'(1PE7.1)') PNORM
        TITLE(54:61) = ' Violtn'
        WRITE(RITOUT(54:61),'(1PE7.1)') VIOL
      ELSEIF(NPA.EQ.2) THEN
        PHASE(44:56) = 'Optimization'
        TITLE(5:8) = ' It'
        CALL HHADJF(CHITER,' ',' ','R',NSHF,IERSHF)
        IF(CHITER(6:6).NE.')'.AND.CHITER(6:6).NE.']') THEN
          NSHF = 1
          CALL HHADJP(CHITER,' ',' ','L',NSHF,IERSHF)
        ENDIF
        RITOUT(3:8) = CHITER
        TITLE(11:13) = 'Qit'
        WRITE(RITOUT(9:13),'(I5)') NQPITR
        TITLE(15:17) = 'Nkt'
        WRITE(RITOUT(15:17),'(I3)') NKT
        TITLE(20:24) = 'Ndof'
        WRITE(RITOUT(19:24),'(I5)') NDOF
        IF(QPOPTN.EQ.'SPARSE') THEN
          TITLE(27:34) = 'KT Cond'
        ELSE
          TITLE(27:34) = 'PH Cond'
        ENDIF
        WRITE(RITOUT(27:34),'(1PE7.1)') CNDNUM
        TITLE(36:43) = '  Step '
        WRITE(RITOUT(36:43),'(1PE7.1)') ALFA
        TITLE(45:52) = ' Norm p'
        WRITE(RITOUT(45:52),'(1PE7.1)') PNORM
        TITLE(54:61) = ' Violtn'
        WRITE(RITOUT(54:61),'(1PE7.1)') VIOL
        TITLE(63:70) = 'Levnbrg'
        WRITE(RITOUT(63:70),'(1PE7.1)') DIAGNL
        TITLE(72:79) = 'Penalty'
        WRITE(RITOUT(72:79),'(1PE7.1)') PENLTY
        TITLE(81:88) = 'Norm Pg'
        WRITE(RITOUT(81:88),'(1PE7.1)') PGNORM
        TITLE(90:104) = 'Merit Function'
        WRITE(RITOUT(90:104),'(1PE14.6)') SQPMF
      ENDIF
C
C         
      IF(NP.LT.0) THEN
        NP = -NP
        WRITE(IPU,1002) PHASE
        WRITE(IPU,1001) TITLE
      ENDIF
C
      WRITE(IPU,1001) RITOUT
C
 1001 FORMAT(1X,A104)
 1002 FORMAT(/1X,A104)
C 
      RETURN
      END
