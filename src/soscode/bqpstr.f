      SUBROUTINE BQPSTR(CVEC,MSUBE,MAXCON,ETAVEC,CMAT,
     $    IROWC,JCOLC,NONZC,BVEC,MSUBB,MAXBND,BMAT,IROWB,JCOLB,
     $    NONZB,GVEC,NVAR,FKAPPA,PENMU,SCRTCH,LNSCRT)
C
C
C ======================================================================
C     BQPSTR===>bqpstr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C          PURPOSE:  THE PURPOSE OF THIS ROUTINE IS TO COMPUTE 
C                    AN ESTIMATE FOR THE BARRIER PARAMETER.
C                    CENTRAL PATH ESTIMATES lamda = mu(D_b)^(-1)e ARE USED.
C                    THE VALUE OF MU IS CHOSEN SUCH THAT
C
C                        ||g - (C^T)eta - mu(B^T)(D_b)^(-1)e||
C                        || - - - - - - - - - - - - - - - - || = K mu
C                        ||               c                 ||
C
C                    WHERE K IS AN INPUT (E.G. 1.1 X 10 = 11).
C
C         ARGUMENTS:    
C
C           CVEC    I  EQUALITY CONSTRAINTS (MAXCON)
C           MSUBE   I  NUMBER OF EQUALITY CONSTRAINTS MSUBE 
C           MAXCON  I  MAXIMUM NUMBER OF CONSTRAINTS
C           ETAVEC  I  LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C           CMAT    I  CONSTRAINT DERIVATIVES AT YVEC (NONZC)
C           IROWC   I  ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C           JCOLC   I  COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C           NONZC   I  NUMBER OF JACOBIAN NONZEROS
C           BVEC    I  BOUND INEQUALITIES (MAXBND)
C           MSUBB   I  NUMBER OF BOUNDS
C           MAXBND  I  MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C           BMAT    I  BOUND DERIVATIVES AT YVEC (NONZB)
C           IROWB   I  ROW INDICES OF BOUND JACOBIAN NONZEROS (NONZB)
C           JCOLB   I  COLUMN INDICES OF BOUND JACOBIAN NONZEROS (NVAR+1)
C           NONZB   I  NUMBER OF BOUND JACOBIAN NONZEROS
C           GVEC    I  GRADIENT AT YVEC (NVAR)
C           NVAR    I  NUMBER OF VARIABLES
C           FKAPPA  I  THE FACTOR K
C           PENMU   O  INTERIOR POINT (BARRIER) PARAMETER 
C           SCRTCH  I  REAL SCRATCH ARRAY (LNSCRT)
C           LNSCRT  I  LENGTH OF SCRATCH ARRAY (.GE. 2*NVAR+MAXBND)
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      DIMENSION CVEC(MAXCON)   ,ETAVEC(MAXCON) ,CMAT(NONZC)
     $         ,IROWC(NONZC)   ,JCOLC(NVAR+1)  ,BVEC(MAXBND)
     $         ,BMAT(NONZB)    ,IROWB(NONZB)   ,JCOLB(NVAR+1)      
     $         ,GVEC(NVAR)     ,SCRTCH(LNSCRT)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
      PENMU = ZEROOT
      IF(MSUBB.EQ.0) RETURN
C
C             STORAGE ALLOCATION FOR THE SCRATCH ARRAY
C     ---ALLOCATE THE ARRAYS (I.E. CONSTRUCT THE POINTERS)
C
      LCAVEC = 1
      LCBVEC = LCAVEC + NVAR
      LCLAMB = LCBVEC + NVAR
      LCRSCR = LCLAMB + MAXBND
      NCRSCR = LNSCRT - LCRSCR + 1
C
      IF(NCRSCR.LT.0) THEN
        NEEDED = LCRSCR - 1
        PRINT *,'LNSCRT IS TOO SMALL; NEEDED =',NEEDED
        RETURN
      ENDIF
C
C         COMPUTE (CMAT**T)(ETAVEC) AND STORE IN SCRTCH(LCAVEC)
C
      CALL MVPSPR(11,NVAR,MSUBE,CMAT,IROWC,JCOLC,ETAVEC,SCRTCH(LCAVEC))
C
C         CONSTRUCT avec = g - (C^T)eta
C
      DO I = 1,NVAR
        SCRTCH(LCAVEC+I-1) = GVEC(I) - SCRTCH(LCAVEC+I-1)
      enddo
C
C         COMPUTE vlamda = 1/bvec
C
      DO I = 1,MSUBB
        SCRTCH(LCLAMB+I-1) = ONE/BVEC(I)
      enddo
C
C         COMPUTE (bmat**T)(vlamda) AND SAVE IN SCRTCH(LCBVEC)
C
      CALL MVPSPR(11,NVAR,MSUBB,BMAT,IROWB,JCOLB,SCRTCH(LCLAMB),
     $              SCRTCH(LCBVEC))
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
      FNORM = ZERO
      DO I = 1,NVAR
        ASUBI = SCRTCH(LCAVEC+I-1)
        BSUBI = SCRTCH(LCBVEC+I-1)
        XMU1 = ASUBI/(BSUBI+FKAPPA)
        XMU2 = -ASUBI/(FKAPPA-BSUBI)
        FNORM1 = ABS(ASUBI - XMU1*BSUBI)
        FNORM2 = ABS(ASUBI - XMU2*BSUBI)
        IF(FNORM1.GT.FNORM) THEN
          FNORM = FNORM1
          PENMU = XMU1
        ENDIF
        IF(FNORM2.GT.FNORM) THEN
          FNORM = FNORM2
          PENMU = XMU2
        ENDIF
      enddo
C
      DO I = 1,MSUBE
        FNORMC = ABS(CVEC(I))
        IF(FNORMC.GT.FNORM) THEN
          FNORM = FNORMC
          PENMU = FNORMC/FKAPPA
        ENDIF
      enddo
C
      PENMU = MIN(ONE,MAX(PENMU,.1D0*ZEROOT))
C
      RETURN 
      END
