
      SUBROUTINE STATSM
C
C ======================================================================
C     STATSM===>statsm   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
      PARAMETER (HUNDRD = 1.0D2)
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
C
      INSTAT(3) = INSTAT(3) + 1
      RLSTAT(1) = DBLE(INSTAT(5))/DBLE(MAX(1,INSTAT(6)))
      RLSTAT(2) = DBLE(INSTAT(1))/DBLE(MAX(1,INSTAT(6)))
      TOTLGC = INSTAT(1) + INSTAT(2) + INSTAT(3)
      RLSTAT(3) = DBLE(INSTAT(1))/TOTLGC
      RLSTAT(4) = DBLE(INSTAT(2))/TOTLGC
      RLSTAT(5) = DBLE(INSTAT(3))/TOTLGC
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
      CALL CLKOUT(4,TSHRQP,XSHRQP)
      CALL CLKOUT(5,TMSRQP,XINST)
      INSTAT(13) = XINST
      CALL CLKOUT(6,TKTFAC,XNKTFC)
C
      CALL RSTRNG( 'Total CPU Time',CPUTIM,IPUNLP)
C
      IF(IOFLAG.GT.10) THEN
      CALL RSTRNG( 'Total Time Outside NLP',TNONLP,IPUNLP)
      CALL RSTRNG( 'Total Time Inside NLP',CPUTIM-TNONLP,IPUNLP)
      CALL RSTRNG( 'Total Time For All QP Subproblems',TSHRQP,IPUNLP)
      CALL RSTRNG( 'Time of Longest QP Subproblem',TMSRQP,IPUNLP)
      CALL RSTRNG( 'Total Time For Matrix Factorizations',TKTFAC,IPUNLP)
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
      CALL ISTRNG( 'Number of GCALLS For Constraint Elimination'
     $   ,INSTAT(1),IPUNLP)
      CALL ISTRNG( 'Number of GCALLS For Merit Function Minimization'
     $   ,INSTAT(2),IPUNLP)
      CALL ISTRNG( 'Number of GCALLS For Feasible Point Excluding First'
     $   ,INSTAT(3),IPUNLP)
      CALL ISTRNG( 'Maximum Number of Iterations Per Constraint Eliminat
     $ion'  ,INSTAT(4),IPUNLP)
      CALL ISTRNG( 'Total Number of Iterations For Constraint Eliminatio
     $n' ,INSTAT(5),IPUNLP)
      CALL ISTRNG( 'Total Number of Constraint Eliminations'
     $   ,INSTAT(6),IPUNLP)
      CALL ISTRNG( 'Maximum Number of Iterations Wasted in Constraint El
     $imination'  ,INSTAT(7),IPUNLP)
      CALL ISTRNG( 'Total Number of Iterations Wasted in Constraint Elim
     $ination'   ,INSTAT(8),IPUNLP)
      CALL ISTRNG( 'Total Number of Wasted Constraint Eliminations '
     $   ,INSTAT(9),IPUNLP)
      CALL RSTRNG( 'Average Number of Iterations Per Constraint Eliminat
     $ion' ,RLSTAT(1),IPUNLP)
      CALL RSTRNG( 'Average Number of GCALLS Per Constraint Elimination'
     $  ,RLSTAT(2),IPUNLP)
      CALL RSTRNG( 'Percent of GCALLS in Constraint Elimination'
     $   ,HUNDRD*RLSTAT(3),IPUNLP)
      CALL RSTRNG( 'Percent of GCALLS in Merit Function Minimization'
     $   ,HUNDRD*RLSTAT(4),IPUNLP)
      CALL RSTRNG( 'Percent of GCALLS in Feasible Point Location'
     $   ,HUNDRD*RLSTAT(5),IPUNLP)
      ENDIF
C
      IF(IOFLAG.GT.10) THEN
        WRITE(IPUNLP,1002)
        IF(QPOPTN.EQ.'SPARSE') THEN
      CALL ISTRNG( 'Number of Schur-QP Failures'
     $   ,INSTAT(11),IPUNLP)
      CALL ISTRNG( 'Number of Schur-QP Calls'
     $   ,INSTAT(12),IPUNLP)
      CALL ISTRNG( 'NLP Iteration Number for Most Expensive Schur-QP Cal
     $l' ,INSTAT(13),IPUNLP)
      CALL ISTRNG( 'Maximum Number of Calls to LINPACK'
     $   ,INSTAT(14),IPUNLP)
      CALL ISTRNG( 'NLP Iteration Number For Maximum Number Calls to LIN
     $PACK' ,INSTAT(15),IPUNLP)
      CALL ISTRNG( 'Maximum Number of Updates in Schur-Complement'
     $   ,INSTAT(16),IPUNLP)
      CALL ISTRNG( 'NLP Iteration Number For Maximum Number Updates in S
     $chur-Complement' ,INSTAT(17),IPUNLP)
      CALL ISTRNG( 'Maximum Size of Matrix in LINPACK Call'
     $   ,INSTAT(18),IPUNLP)
      CALL ISTRNG( 'NLP Iteration Number For Maximum Size Matrix Call'
     $   ,INSTAT(19),IPUNLP)
      CALL ISTRNG( 'Maximum Number of Calls to Multifrontal Solve'
     $   ,INSTAT(20),IPUNLP)
      CALL ISTRNG( 'NLP Iteration Number For Maximum Number of Multifron
     $tal Solves'   ,INSTAT(21),IPUNLP)
        ELSE
          CALL CLKOUT(11,TMVPQG,XMVPQG)
          CALL RSTRNG( 'Time for Matrix-Vector Products in QPGETD',
     $                  TMVPQG,IPUNLP)
          CALL CLKOUT(12,TMVPLP,XMVPLP)
          CALL RSTRNG( 'Time for Matrix-Vector Products in LPCORE',
     $                  TMVPLP,IPUNLP)
          CALL CLKOUT(13,TMQPOP,XMQPOP)
          CALL RSTRNG( 'Time for QPOPT Initialization',
     $                  TMQPOP,IPUNLP)
          CALL CLKOUT(14,TMQPOP,XMQPOP)
          CALL RSTRNG( 'Time for QPOPT Phase I',
     $                  TMQPOP,IPUNLP)
          CALL CLKOUT(15,TMQPOP,XMQPOP)
          CALL RSTRNG( 'Time for QPOPT Phase II',
     $                  TMQPOP,IPUNLP)
          CALL CLKOUT(16,TMQPOP,XMQPOP)
          CALL RSTRNG( 'Time for QPOPT Levenberg Modification',
     $                  TMQPOP,IPUNLP)
          CALL CLKOUT(17,TMQPOP,XMQPOP)
          CALL RSTRNG( 'Time for QPOPT Conclusion',
     $                  TMQPOP,IPUNLP)
          CALL CLKOUT(18,TMQPOP,XMQPOP)
          CALL RSTRNG( 'Time for QPOPT Diagnostics',
     $                  TMQPOP,IPUNLP)
          CALL CLKOUT(19,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE Phase II Constraint Deletion',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(20,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE Phase II Constraint Addition',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(21,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE New Search Direction',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(22,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE Distance to Nearest Constraint',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(23,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE Update of TQ Factorization',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(24,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE Check of Feasibility',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(25,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE New Gradient Calculation',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(26,TMQPCR,XMQPCR)
          CALL RSTRNG( 'Time for QPCORE Check of Indefiniteness',
     $                  TMQPCR,IPUNLP)
          CALL CLKOUT(27,TMLPCR,XMLPCR)
          CALL RSTRNG( 'Time for LPCORE Distance to Nearest Constraint',
     $                  TMLPCR,IPUNLP)
          CALL CLKOUT(28,TMLPCR,XMLPCR)
          CALL RSTRNG( 'Time for LPCORE Update of TQ Factorization',
     $                  TMLPCR,IPUNLP)
          CALL CLKOUT(29,TMLPCR,XMLPCR)
          CALL RSTRNG( 'Time for LPCORE Check of Feasibility',
     $                  TMLPCR,IPUNLP)
          CALL CLKOUT(30,TMLPCR,XMLPCR)
          CALL RSTRNG( 'Time for LPCORE New Gradient Calculation',
     $                  TMLPCR,IPUNLP)
          CALL CLKOUT(31,TMLPCR,XMLPCR)
          CALL RSTRNG( 'Time for LPCORE Constraint Deletion',
     $                  TMLPCR,IPUNLP)
          NP2ITR = INSTAT(26) - INSTAT(29)
          CALL ISTRNG( 'Number of Phase I Iterations',
     $                  INSTAT(29),IPUNLP)
          CALL ISTRNG( 'Number of Phase II Iterations',
     $                  NP2ITR,IPUNLP)
        ENDIF
      CALL ISTRNG( 'Total Number of QP Iterations'
     $   ,INSTAT(26),IPUNLP)
      CALL ISTRNG( 'Maximum Number of QP Iterations on any QP Subproblem
     $'  ,INSTAT(27),IPUNLP)
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
        instat(28) = 1
      CALL ISTRNG( '**** OUT OF CORE MODE; Storage Needed in HOLD Array'
     $   ,INSTAT(10),IPUNLP)
      CALL ISTRNG( '**** OUT OF CORE MODE; Storage Needed in IHOLD Array
     $'  ,INSTAT(30),IPUNLP)
      ENDIF
C
      WRITE(IPUNLP,1002)
      WRITE(IPUNLP,1003)
C
 1001 FORMAT(T10,'|',T30,A6,' ALGORITHM PERFORMANCE STATISTICS'
     $   ,T100,'|')
 1002 FORMAT(T10,'|',T100,'|')
 1003 FORMAT(T10,91('-'))
      RETURN
      END
