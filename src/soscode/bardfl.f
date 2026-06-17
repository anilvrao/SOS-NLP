


      SUBROUTINE BARDFL
C
C ======================================================================
C     BARDFL===>bardfl   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  LOAD DEFAULT VALUES INTO THE NLPSPR.CMNMONS
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
C         DEFAULT VALUE COMMON
C
      COMMON /NPSDF1/ RDFLT(11),IDFLT(35)
      COMMON /NBSDF1/ RBDFLT(5),IBDFLT(3)
      COMMON /NPSDF2/ CDFLT(3),DFAWLT
      CHARACTER(LEN=6) CDFLT
      CHARACTER(LEN=17) DFAWLT
      LOGICAL OPUN
C
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,   
     $    ONEEM2 = 1.0D-2, ONEEM3 = 1.D-3, ONEEM5 = 1.0D-5, 
     $    ONEEM7 = 1.0D-7 ,POINT9 = 9.0D-1)
C
      DFAWLT = '(BARRIER DEFAULT)'
C
C             ALGORITHM DEFAULTS. (DEFAULT ALL SPRNLP VALUES 
C             EVEN IF THEY ARE NOT USED)
C
      ALFLWR = ZERO
      ALFUPR = ONE
      CONTOL = SQRT(HDMCON(5))
      EPSRLF = CONTOL
      OBJTOL = ONEEM7
      PGDTOL = ONEEM5
      SLPTOL = POINT9
      SFZTOL = ONEEM2
      TOLFIL = TWO
      TOLKTC = (ONE/HDMCON(5))**1.6
      TOLPVT = ONEEM3
C
      BIGCON = 1.0D2
      FEATOL = ONEEM2
      PMULWR = 1.0D-1
      PTHTOL = 1.0D1
      RHOLWR = 1.0D2
C
      RDFLT(1) = ALFLWR
      RDFLT(2) = ALFUPR
      RDFLT(3) = CONTOL
      RDFLT(4) = EPSRLF
      RDFLT(5) = OBJTOL
      RDFLT(6) = PGDTOL
      RDFLT(7) = SLPTOL
      RDFLT(8) = SFZTOL
      RDFLT(9) = TOLFIL
      RDFLT(10) = TOLKTC
      RDFLT(11) = TOLPVT
C
      RBDFLT(1) = BIGCON
      RBDFLT(2) = FEATOL
      RBDFLT(3) = PMULWR
      RBDFLT(4) = PTHTOL
      RBDFLT(5) = RHOLWR
C
      IHESHN = 0
      IOFLAG = 10
      IOFLIN = -1
      IOFMFR = 0
      IOFPAT = 0
      IOFSHR = 0
      IOFSRC = 0
      IPUDRF = 0
      IPUFZF = 0
      IPUMF1 = 11
      IPUMF2 = 12
      IPUMF3 = 13
      IPUMF4 = 14
      IPUMF5 = 15
      IPUMF6 = 16
      IPUMF7 = 17
      INQUIRE(IPUMF1,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF1)
      INQUIRE(IPUMF2,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF2)
      INQUIRE(IPUMF3,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF3)
      INQUIRE(IPUMF4,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF4)
      INQUIRE(IPUMF5,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF5)
      INQUIRE(IPUMF6,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF6)
      INQUIRE(IPUMF7,OPENED=OPUN)
      IF(OPUN) CLOSE(IPUMF7)
      IPUNLP = JHMCON(6)
      IPUSTF = 0
      IRELAX = 0
      ITDRQP = -1
      ITFZQP = -1
      IT1MAX = 20
      JACPRM = 0
      LYNFNC = 0
      LYNOUT = 0
      LYNPLT = 0
      LYNPNT = 101
      LYNVAR = 0
      MAXLYN = 5
      MAXNFE = 10000
      MNSAME = 2
      NEWTON = 0
      NITMAX = 100
      NITMIN = 0
      NORMAL = 0
C
      IMAXMU = 10
      MUCALC = 3
      MXQPIT = 1
C
      IDFLT(1) = IHESHN
      IDFLT(2) = IOFLAG
      IDFLT(3) = IOFLIN
      IDFLT(4) = IOFMFR
      IDFLT(5) = IOFPAT
      IDFLT(6) = IOFSHR
      IDFLT(7) = IOFSRC
      IDFLT(8) = IPUDRF
      IDFLT(9) = IPUFZF
      IDFLT(10) = IPUMF1
      IDFLT(11) = IPUMF2
      IDFLT(12) = IPUMF3
      IDFLT(13) = IPUMF4
      IDFLT(14) = IPUMF5
      IDFLT(15) = IPUMF6
      IDFLT(16) = IPUMF7
      IDFLT(17) = IPUNLP
      IDFLT(18) = IPUSTF
      IDFLT(19) = IRELAX
      IDFLT(20) = ITDRQP
      IDFLT(21) = ITFZQP
      IDFLT(22) = IT1MAX
      IDFLT(23) = JACPRM
      IDFLT(24) = LYNFNC
      IDFLT(25) = LYNOUT
      IDFLT(26) = LYNPLT
      IDFLT(27) = LYNPNT
      IDFLT(28) = LYNVAR
      IDFLT(29) = MAXLYN
      IDFLT(30) = MAXNFE
      IDFLT(31) = MNSAME
      IDFLT(32) = NEWTON
      IDFLT(33) = NITMAX
      IDFLT(34) = NITMIN
      IDFLT(35) = NORMAL
C
      IBDFLT(1) = IMAXMU
      IBDFLT(2) = MUCALC
      IBDFLT(3) = MXQPIT
C
      ALGOPT = 'M     '
      KTOPTN = 'SMALL '
      QPOPTN = 'SPARSE'
C
      ALGNAM = 'BARNLP'
C
      CDFLT(1) = ALGOPT
      CDFLT(2) = KTOPTN
      CDFLT(3) = QPOPTN
C
C         SET INNPER TO INDICATE DEFAULTS HAVE BEEN SET
C
      INNPER = 0
C
      RETURN
      END
