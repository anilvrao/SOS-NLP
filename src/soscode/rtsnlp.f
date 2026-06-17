

      SUBROUTINE RTSNLP(STRING,IVAL,RVAL,CVAL)
C
C ======================================================================
C     RTSNLP===>RTSNLP   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  RETRIEVE VALUES IN THE NLPSPR.CMN
C                   GIVEN A CHARACTER STRING INPUT
C
C         INPUT:
C
C           STRING  CHARACTER VARIABLE CORRESPONDING TO NAME OF
C                   QUANTITY IN NLPSPR.CMN
C         OUTPUT:
C
C           IVAL    INTEGER VALUE OF "STRING" IF "STRING" IS A VALID 
C                   INTEGER VARIABLE NAME;  CLOBBER CONSTANT OTHERWISE
C
C           RVAL    REAL VALUE OF "STRING" IF "STRING" IS A VALID 
C                   REAL VARIABLE NAME;  CLOBBER CONSTANT OTHERWISE
C
C           CVAL    CHARACTER VALUE OF "STRING" IF "STRING" IS A VALID 
C                   CHARACTER VARIABLE NAME;  BLANK OTHERWISE
C
C         COMMONS:
C
C           NPSPRR  NLPSPR REAL VARIABLE INPUTS
C           NPSPRI  NLPSPR INTEGER VARIABLE INPUTS
C           NPSPRC  NLPSPR CHARACTER VARABLE INPUTS
C                   
      CHARACTER(LEN=*) STRING
C
      PARAMETER (NRSYM=11,NISYM=35,NCSYM=3)
      PARAMETER (NBRSYM=5,NBISYM=3,NBCSYM=2)
      PARAMETER (NSYMBL=NRSYM+NISYM+NCSYM)
      PARAMETER (NBSMBL=NBRSYM+NBISYM+NBCSYM)
      CHARACTER(LEN=6)  INPBUF
      CHARACTER(LEN=6)  CVAL
      CHARACTER(LEN=6)  SYMBUL(NSYMBL)
      CHARACTER(LEN=6)  BSYMBL(NBSMBL)
C
C======================================================================
C Array equivalences into include file NLPSPR.CMN
C
      INCLUDE '../commons/NLPSPR.CMN'
C
C     COMMON /NPSPRR/ ALFLWR, ...      double precision
C     COMMON /NPSPRI/ INNPER, ...      integer
C     COMMON /NPSPRC/ ALGOPT, ...      character(len=6)
C     COMMON /NPSALG/ ALGNAM           character(len=6) (no equivalence)
C
      DOUBLE PRECISION RNPSPR(NRSYM)
      INTEGER          INPSPR(NISYM+1)
      CHARACTER(LEN=6)  CNPSPR(NCSYM)
C
      EQUIVALENCE (RNPSPR(1),ALFLWR)
      EQUIVALENCE (INPSPR(1),INNPER)
      EQUIVALENCE (CNPSPR(1),ALGOPT)
C
C Array equivalences into include file BARNLP.CMN
C
      INCLUDE '../commons/BARNLP.CMN'
C
C     COMMON /BNPSPR/ BIGCON, ...      double precision
C     COMMON /BNPSPI/ IMAXMU, ...      integer
C
      DOUBLE PRECISION RBNPSP(NBRSYM)
      INTEGER          IBNPSP(NBISYM)
C
      EQUIVALENCE (RBNPSP(1),BIGCON)
      EQUIVALENCE (IBNPSP(1),IMAXMU)
C
C End of array equivalences
C======================================================================
C
C-------------------------------------------------------------
C
      DATA (SYMBUL(I),I=1,NSYMBL) /
     $    'ALFLWR', 'ALFUPR', 'CONTOL', 'EPSRLF', 'OBJTOL', 
     $    'PGDTOL', 'SLPTOL', 'SFZTOL', 'TOLFIL', 'TOLKTC', 
     $    'TOLPVT', 'IHESHN', 'IOFLAG', 'IOFLIN', 'IOFMFR',      
     $    'IOFPAT', 'IOFSHR', 'IOFSRC', 'IPUDRF', 'IPUFZF',   
     $    'IPUMF1', 'IPUMF2', 'IPUMF3', 'IPUMF4', 'IPUMF5',    
     $    'IPUMF6', 'IPUMF7', 'IPUNLP', 'IPUSTF', 'IRELAX',       
     $    'ITDRQP', 'ITFZQP', 'IT1MAX', 'JACPRM', 'LYNFNC',    
     $    'LYNOUT', 'LYNPLT', 'LYNPNT', 'LYNVAR', 'MAXLYN',    
     $    'MAXNFE', 'MNSAME', 'NEWTON', 'NITMAX', 'NITMIN', 
     $    'NORMAL', 'ALGOPT', 'KTOPTN', 'QPOPTN' /
      DATA (BSYMBL(I),I=1,NBSMBL) /
     $    'BIGCON', 'FEATOL', 'PMULWR', 'PTHTOL', 'RHOLWR',  
     $    'IMAXMU', 'MUCALC', 'MXQPIT', 'ALGOPT', 'KTOPTN' /
C
C-------------------------------------------------------------
C
C         IF INPSPR(1) = -1, THE DEFAULTS HAVE NEVER BEEN SET.  SET
C         THEM AND THEN MODIFY THEM WITH THE INPUT
C
      IF(INPSPR(1).EQ.-1) CALL SPRDFL
C
C         COMPUTE LENGTH OF INPUT CHARACTER STRING
C
      LNINPT = LEN(STRING)
C
C         SAVE DEFAULTS IN THE INTEGER, REAL, AND CHARACTER VALUES
C
      IVAL = JHMCON(1)
      RVAL = HDMCON(1)
      CVAL = '      '
C
C         IF INPUT STRING HAS INCORRECT LENGTH, EXIT
C
      IF(LNINPT.NE.6) GO TO 10000
C
      INPBUF = STRING
C
C         CAPITALIZE ALL THE LOWER CASE SYMBOLS IN THE INPUT STRING
C
      CALL CAPTUL(INPBUF,6,IERREP)
C
      IF(IERREP.NE.0) GO TO 10000
C
      ISYM = 0
      DO II = 1,NRSYM
        ISYM = ISYM + 1
        IF(INPBUF.EQ.SYMBUL(ISYM)) THEN
          RVAL = RNPSPR(II)
          GO TO 10000
        ENDIF
      ENDDO
C
      DO II = 1,NISYM
        ISYM = ISYM + 1
        IF(INPBUF.EQ.SYMBUL(ISYM)) THEN
          IVAL = INPSPR(II+1)
          GO TO 10000
        ENDIF
      ENDDO
C
      DO II = 1,NCSYM
        ISYM = ISYM + 1
        IF(INPBUF.EQ.SYMBUL(ISYM)) THEN
          CVAL = CNPSPR(II)
          GO TO 10000
        ENDIF
      ENDDO
C
      ISYM = 0
      DO II = 1,NBRSYM
        ISYM = ISYM + 1
        IF(INPBUF.EQ.BSYMBL(ISYM)) THEN
          RVAL = RBNPSP(II)
          GO TO 10000
        ENDIF
      ENDDO
C
      DO II = 1,NBISYM
        ISYM = ISYM + 1
        IF(INPBUF.EQ.BSYMBL(ISYM)) THEN
          IVAL = IBNPSP(II+1)
          GO TO 10000
        ENDIF
      ENDDO
C
C
10000 CONTINUE
C
      RETURN
      END
