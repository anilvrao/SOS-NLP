

      SUBROUTINE STATBR
C
C ======================================================================
C     STATBR===>statbr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
C
      INSTAT(3) = INSTAT(3) + 1
C
      IF(IOFLAG.EQ.0) RETURN
C
      WRITE(IPUNLP,1003)
      WRITE(IPUNLP,1002)
      WRITE(IPUNLP,1001) ALGNAM
      WRITE(IPUNLP,1002)
C
      CALL CLKOUT(1,CPUTIM,XNTIME)
      CALL CLKOUT(3,TNONLP,XNONLP)
      CALL CLKOUT(4,TIMKKT,XNMKKT)
      CALL CLKOUT(5,TMXKKT,XINST)
      INSTAT(13) = XINST
      CALL CLKOUT(6,TMLNAL,XNKTFC)
C
      CALL RSTRNG( 'Total CPU Time',CPUTIM,IPUNLP)
C
      IF(IOFLAG.GT.10) THEN
      CALL RSTRNG( 'Total Time Outside NLP',TNONLP,IPUNLP)
      CALL RSTRNG( 'Total Time Inside NLP',CPUTIM-TNONLP,IPUNLP)
      CALL RSTRNG( 'Total Time For All KKT Solutions',TIMKKT,IPUNLP)
      CALL RSTRNG( 'Time of Longest KKT Solution',TMXKKT,IPUNLP)
      CALL RSTRNG( 'Total Time For Linear Algebra',TMLNAL,IPUNLP)
      ENDIF
C
      WRITE(IPUNLP,1002)
C
      CALL ISTRNG( 'Number of Function Calls'
     $   ,INSTAT(22),IPUNLP)
      CALL ISTRNG( 'Number of Gradient Calls'
     $   ,INSTAT(23),IPUNLP)
      CALL ISTRNG( 'Number of Hessian Calls'
     $   ,INSTAT(24),IPUNLP)
      CALL ISTRNG( 'Total Number of Function Evaluations'
     $   ,INSTAT(25),IPUNLP)
C
      IF(IOFLAG.GT.10) THEN
        WRITE(IPUNLP,1002)
      CALL ISTRNG( 'Number of KKT Solution Failures'
     $   ,INSTAT(11),IPUNLP)
      CALL ISTRNG( 'Number of KKT System Calls'
     $   ,INSTAT(12),IPUNLP)
      CALL ISTRNG( 'NLP Iteration Number for Most Expensive KKT Solution
     $'  ,INSTAT(13),IPUNLP)
      ENDIF
C
      WRITE(IPUNLP,1002)
C
      IF(INSTAT(28).EQ.0) THEN
        CALL ISTRNG( 'Storage Needed in HOLD Array'
     $   ,INSTAT(10),IPUNLP)
        IF(INSTAT(30).GT.1) CALL ISTRNG( 'Storage Needed in IHOLD Array'
     $   ,INSTAT(30),IPUNLP)
      ELSE
      CALL ISTRNG( '**** OUT OF CORE MODE; Storage Needed in HOLD Array'
     $   ,INSTAT(10),IPUNLP)
      CALL ISTRNG( '**** OUT OF CORE MODE; Storage Needed in IHOLD Array
     $'  ,INSTAT(30),IPUNLP)
      ENDIF
C
      WRITE(IPUNLP,1002)
      WRITE(IPUNLP,1003)
C
 1001 FORMAT(T10,'|',T30,A8,' ALGORITHM PERFORMANCE STATISTICS'
     $   ,T100,'|')
 1002 FORMAT(T10,'|',T100,'|')
 1003 FORMAT(T10,91('-'))
      RETURN
      END
