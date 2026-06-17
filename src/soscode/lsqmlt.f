
      SUBROUTINE LSQMLT(MUOPTN,NCALL,CVEC,MSUBE,MAXCON,ETAVEC,CMAT,
     $    IROWC,JCOLC,NONZC,BVEC,MSUBB,MAXBND,VLAMDA,BMAT,IROWB,JCOLB,
     $    NONZB,GVEC,NVAR,PENMU,SCRTCH,RWORK,LNRWRK,IPU,
     $    IPC,CNDNUM,GRADLN,IERLSM,NEEDED)
C
C
C ======================================================================
C     LSQMLT===>lsqmlt   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C          PURPOSE:  The purpose of this routine is to compute 
C                    least squares estimates for the Lagrange 
C                    multipliers and optionally the barrier parameter.
C                    Central path estimates lamda = mu(D_b)^(-1)e are used.
C                    The following options are implemented:
C                    (1) Compute multipliers (eta) [mu: fixed] to minimize
C
C                        ||(C^T)eta + (B^T)lamda - g||
C
C                    (2) Compute multipliers (eta,mu) to minimize
C
C                        ||(C^T)eta + mu(B^T)(D_b)^(-1)e - g||
C
C                    (3) Compute (eta,lamda,mu) to minimize
C
C                        ||g - (C^T)eta - (B^T)lamda||
C                        ||   (D_b)lamda - mu(e)    ||
C
C                    In all cases when the systems are underdetermined
C                    the minimum norm solution is computed.
C
C         ARGUMENTS:    
C
C           MUOPTN  I  INTEGER VARIABLE DEFINING OPTION
C           NCALL   IO INTEGER VARIABLE DEFINING CALL NUMBER
C                      MUST BE SET 0 BEFORE FIRST CALL, THEREAFTER
C                      SUBSEQUENT CALLS (NCALL > 1) WILL DO A SOLVE ONLY
C           CVEC    I  EQUALITY CONSTRAINTS (MAXCON)
C           MSUBE   I  NUMBER OF EQUALITY CONSTRAINTS MSUBE 
C           MAXCON  I  MAXIMUM NUMBER OF CONSTRAINTS
C           ETAVEC  O  LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C           CMAT    I  CONSTRAINT DERIVATIVES AT YVEC (NONZC)
C           IROWC   I  ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C           JCOLC   I  COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C           NONZC   I  NUMBER OF JACOBIAN NONZEROS
C           BVEC    I  BOUND INEQUALITIES (MAXBND)
C           MSUBB   I  NUMBER OF BOUNDS
C           MAXBND  I  MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C           VLAMDA  O  LAGRANGE MULTIPLIERS FOR BOUNDS (MAXBND)
C           BMAT    I  BOUND DERIVATIVES AT YVEC (NONZB)
C           IROWB   I  ROW INDICES OF BOUND JACOBIAN NONZEROS (NONZB)
C           JCOLB   I  COLUMN INDICES OF BOUND JACOBIAN NONZEROS (NVAR+1)
C           NONZB   I  NUMBER OF BOUND JACOBIAN NONZEROS
C           GVEC    I  GRADIENT AT YVEC (NVAR)
C           NVAR    I  NUMBER OF VARIABLES
C           PENMU   O  INTERIOR POINT (BARRIER) PARAMETER 
C           SCRTCH  I  REAL SCRATCH ARRAY MAX[(NVAR+2*MAXBND+MAXCON+1),2*NVAR]
C           RWORK   I  REAL WORK ARRAY (LNRWRK)
C           LNRWRK  I  LENGTH OF RWORK; LNRWRK .GE. ?
C           IPU     I  OUTPUT UNIT NO.
C           IPC     I  OUTPUT CONTROL FLAG
C           CNDNUM  O  CONDITION NUMBER OF MATRIX
C           GRADLN  O  NORM OF GRADIENT, ||(C^T)ETA + (B^T)LAMDA - G||
C           IERLSM  O  ERROR RETURN CODE
C           NEEDED  O  STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C
C
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      COMMON /ITEREF/ MAXREF,MAXRFN,IREFIN
C
      DIMENSION CVEC(MAXCON)   ,ETAVEC(MAXCON) ,CMAT(NONZC)
     $         ,IROWC(NONZC)   ,JCOLC(NVAR+1)  ,BVEC(MAXBND)
     $         ,VLAMDA(MAXBND) ,BMAT(NONZB)    ,IROWB(NONZB)
     $         ,JCOLB(NVAR+1)  ,GVEC(NVAR)     
     $         ,SCRTCH(NVAR+2*MAXBND+MAXCON+1)
     $         ,RWORK(LNRWRK)  
C
      DIMENSION INRTIA(3),TYME(6),OPCNTS(2)
C
      LOGICAL UNDER
C
      INTEGER LNZLTA
      INTEGER XDSLNI

C     -----------------------------------------------------------------

      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
C         INITIALIZE ERROR RETURN FLAG
C
      IERLSM = 0
C
C         BEGIN FACTORIZATION TIMING CLOCK
C
      CALL CLKBEG(6)
C
C         COMPUTE ROW AND COLUMN DIMENSIONS FOR A AND THE
C         TOTAL NUMBER OF NONZEROS IN THE LOWER TRIANGULAR PART OF
C         THE SYMMETRIC SYSTEM
C
      IF(MUOPTN.EQ.1) THEN
        NROWA = NVAR
        NCOLA = MSUBE 
        NZLTA = NONZC + MSUBE
      ELSEIF(MUOPTN.EQ.2) THEN
        NROWA = NVAR 
        NCOLA = MSUBE + 1
        NZLTA = NONZC + NVAR + MAX(NROWA,NCOLA)
      ELSEIF(MUOPTN.EQ.3) THEN
        NROWA = NVAR + MSUBB
        NCOLA = MSUBE + MSUBB + MIN(MSUBB,1)
        NZLTA = NONZC + NONZB + 2*MSUBB + MAX(NROWA,NCOLA)
      ENDIF
C
C         DETERMINE WHETHER PROBLEM IS OVER OR UNDERDETERMINED
C         (ASSUMING FULL RANK MATRICES)
C
      UNDER = NROWA.LE.NCOLA
C
C         COMPUTE THE NUMBER OF EQUATIONS IN THE SYMMETRIC SYSTEM
C
      NEQNS = NROWA + NCOLA
C
C     ------------------------------------------------------------
C        PREDICTED MEMORY REQUIREMENTS FOR SYMBOLIC PHASES OF 
C        BCSLIB-EXT
C
C        REQUIREMENT BELOW SHOULD ALLOW SYMBOLIC PHASES TO PROCEED
C        "IN-CORE". 
C     ------------------------------------------------------------
     
      LNZLTA = XDSLNI ( NZLTA)

      NWRKCK = 200 + 2*LNZLTA + 2*NEQNS
      NWRKCK = MAX(NWRKCK, 200+10*(NEQNS+1))
C
      IF(NWRKCK.GT.LNRWRK) THEN
        NEEDED = NWRKCK
        IERLSM = -906
        GO TO 560
      ENDIF
C
      SCRTCH(1:NEQNS) = ZERO
      IF(MUOPTN.LE.2.AND.MSUBB.GT.0) THEN
C
C         COMPUTE VLAMDA = 1/BVEC
C
        vlamda(1:msubb) = one/bvec(1:msubb)
C
C         COMPUTE (BMAT**T)(VLAMDA) AND SAVE IN SCRTCH
C
        CALL MVPSPR(11,NVAR,MSUBB,BMAT,IROWB,JCOLB,VLAMDA,
     $              SCRTCH)
C
      ENDIF
C
      IF(NCALL.GE.1) GO TO 480
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C
C         INPUT THE MATRIX STRUCTURE TO THE MULTIFRONTAL CODE
C
C
C         -- SET IN-CORE OPTION FLAG: O (IN-CORE); 1 (OUT-OF-CORE)
C
      MFRSTR = 0
C
C         -- SET MULTIFRONTAL I/O UNITS
C
      CALL XDSLIN ( NEQNS, 'SYMMETRIC', IOFMFR, IPU, IPUMF1, IPUMF2, 
     $              IPUMF3, IPUMF4, 0, IPUMF5, RWORK, LNRWRK, IER )

C
      IF(IER.EQ.-101) THEN
        IERLSM = -907
      ELSEIF(IER.EQ.-102) THEN
        IERLSM = -808
      ENDIF
      IF(IERLSM.NE.0) GO TO 600
C
C          SET PANEL SIZE TO NEQNS UNTIL A BETTER SCHEME IS FOUND
C
      CALL XDSLSP ( 'panel size', NEQNS, DUMMY, RWORK, IER )
      IF(IER.NE.0) GO TO 600
C
C          ----SMALLER TOLERANCE (FOR NO PARTICULAR REASON)
C
      PVTMLT = TOLPVT**2
      CALL XDSLSP ( 'pivot tolerance', 0, PVTMLT, RWORK, IER )
C 
      IF ( IER .NE. 0 ) GO TO 600
C
      CALL XDSLSP ( 'limit fill', 0, TOLFIL, RWORK, IER )
C 
      IF ( IER .NE. 0 ) GO TO 600
C
      CALL XDSLSP ( 'save original matrix', 0, ZERO, RWORK, IER )
C 
      IF ( IER .NE. 0 ) GO TO 600
C
      IF(UNDER) THEN
C
C ------------------------------------------------------------------
C -------------UNDERDETERMINED SYSTEM-------------------------------
C ------------------------------------------------------------------
C
C         ---- IDENTITY ON THE (1,1) BLOCK 
C
        DO I = 1,NCOLA
          IROW = I
          JCOL = I
          CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
C         ---- INPUT C^T
C
        IF(MSUBE.GT.0) THEN
          DO J = 1,NVAR
            jclcloop: DO I = JCOLC(J),JCOLC(J+1)-1
              IF(CMAT(I).EQ.ZERO) cycle jclcloop
              IRWC = IROWC(I)
              JCLC = J
              IROW = JCLC + NCOLA
              JCOL = IRWC
              CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
              IF(IERLSM.NE.0) GO TO 600
            enddo jclcloop
          enddo
        ENDIF
C
        IF(MUOPTN.EQ.2) THEN
C
C         ---- input (B^T)(D_b)^(-1)e
C
          DO I = 1,NROWA
            IROW = NCOLA + I
            JCOL = MSUBE + 1
            CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
            IF(IERLSM.NE.0) GO TO 600
          enddo
C
        ELSEIF(MUOPTN.GE.3) THEN
C
C         ---- INPUT B^T
C
          DO J = 1,NVAR
            DO I = JCOLB(J),JCOLB(J+1)-1
              IRWB = IROWB(I)
              JCLB = J
              IROW = JCLB + NCOLA
              JCOL = IRWB + MSUBE
              CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
              IF(IERLSM.NE.0) GO TO 600
            enddo
          enddo
C
C         ---- INPUT D_b and -e
C
          DO I = 1,MSUBB
C
            IROW = NVAR + NCOLA + I
            JCOL = I + MSUBE
            CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
            IF(IERLSM.NE.0) GO TO 600
C
            IF(MUOPTN.EQ.3) THEN
              IROW = NVAR + NCOLA + I
              JCOL = MSUBE + MSUBB + 1
              CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
              IF(IERLSM.NE.0) GO TO 600
            ENDIF
C
          enddo
C
        ENDIF
C
C         ---- ZERO DIAGONALS ON THE (2,2) BLOCK 
C
        DO I = 1,NROWA
          IROW = NCOLA + I
          JCOL = NCOLA + I
          CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
      ELSE
C
C ------------------------------------------------------------------
C -------------OVERDETERMINED SYSTEM--------------------------------
C ------------------------------------------------------------------
C
C         ---- IDENTITY ON THE (1,1) BLOCK 
C
        DO I = 1,NROWA
          IROW = I
          JCOL = I
          CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
C         ---- INPUT C
C
        IF(MSUBE.GT.0) THEN
          DO J = 1,NVAR
            jcloop: DO I = JCOLC(J),JCOLC(J+1)-1
              IF(CMAT(I).EQ.ZERO) cycle jcloop
              IRWC = IROWC(I)
              JCLC = J
              IROW = IRWC + NROWA
              JCOL = JCLC
              CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
              IF(IERLSM.NE.0) GO TO 600
            enddo jcloop
          enddo
        ENDIF
C
        IF(MUOPTN.EQ.2) THEN
C
C         ---- INPUT [(B^T)(D_b)^(-1)e]^T
C
          DO I = 1,NVAR
            IROW = NROWA + MSUBE + 1
            JCOL = I
            CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
            IF(IERLSM.NE.0) GO TO 600
          enddo
C
        ELSEIF(MUOPTN.GE.3) THEN
C
C         ---- INPUT B
C
          DO J = 1,NVAR
            DO I = JCOLB(J),JCOLB(J+1)-1
              IRWB = IROWB(I)
              JCLB = J
              IROW = IRWB + MSUBE + NROWA
              JCOL = JCLB
              CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
              IF(IERLSM.NE.0) GO TO 600
            enddo
          enddo
C
C         ---- INPUT D_B AND -e
C
          DO I = 1,MSUBB
C
            IROW = NROWA + MSUBE + I
            JCOL = NVAR + I
            CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
            IF(IERLSM.NE.0) GO TO 600
C
            IF(MUOPTN.EQ.3) THEN
              IROW = NROWA + MSUBE + MSUBB + 1
              JCOL = NVAR + I
              CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
              IF(IERLSM.NE.0) GO TO 600
            ENDIF
C
          enddo
C
        ENDIF
C
C         ---- ZERO DIAGONALS ON THE (2,2) BLOCK 
C
        DO I = 1,NCOLA
          IROW = NROWA + I
          JCOL = NROWA + I
          CALL XDSLI1 ( 'A', IROW, JCOL, RWORK, LNRWRK, IERLSM )
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
      ENDIF
C
C         STRUCTURAL INPUT IS COMPLETE
C
      CALL XDSLIF (  RWORK, LNRWRK, NEEDS, IERLSM  )
C
      IF(NEEDS.GT.LNRWRK) THEN
        NEEDED = NEEDS
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDED,LNRWRK
        IERLSM = -67
        GO TO 560
      ENDIF
C
      IF(IERLSM.NE.0) GO TO 600
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C         ORDER THE MATRIX
C
      MAXZER = 0
C
      CALL XDSLOR ( RWORK, LNRWRK, NEEDS, IERLSM  )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(NEEDS.GT.LNRWRK) THEN
        NEEDED = NEEDS
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDED,LNRWRK
        IERLSM = -67
        GO TO 560
      ENDIF
C
      IF(IERLSM.NE.0) GO TO 600
C
C         PERFORM SYMBOLIC FACTORIZATION
C
      CALL XDSLSF ( RWORK, LNRWRK, NEEDS, NEEDMN, IER )
C
      IF(IER.EQ.-300) THEN
        IERLSM = -822
      ELSEIF(IER.EQ.-301) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IERLSM = -924
          IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
          NEEDED = NEEDS
        ENDIF
        GO TO 560
      ELSEIF(IER.EQ.-302) THEN
        IERLSM = -925
      ELSEIF(IER.EQ.-303) THEN
        IERLSM = -626
      ENDIF
      IF(IERLSM.NE.0) GO TO 600
C
C         CHECK THAT DIMENSION OF REAL WORK ARRAY IS ADEQUATE FOR NEXT STEP
C
      IF(MFRSTR.EQ.0) THEN
        NEEDST = NEEDS
      ELSE
        NEEDST = NEEDMN
      ENDIF
C
      IF(NEEDST.GT.LNRWRK) THEN
        NEEDED = NEEDST
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDED,LNRWRK
        IERLSM = -68
        GO TO 560
      ENDIF
C
C         GET MULTIFRONTAL STATISTICS -- SAVE LENGTH NEEDED
C         FOR SYMBOLIC FACTORIZATION
C
      CALL XDSLSR(RWORK,LNSYMB,MXUSED,TYME,OPCNTS)
C
      IF(IERLSM.NE.0) GO TO 600
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C         INPUT MATRIX VALUES TO MULTIFRONTAL SOLVER
C
C
      IF(UNDER) THEN
C
C ------------------------------------------------------------------
C -------------UNDERDETERMINED SYSTEM-------------------------------
C ------------------------------------------------------------------
C
C         ---- IDENTITY ON THE (1,1) BLOCK 
C
        DO I = 1,NCOLA
          IROW = I
          JCOL = I
          VAL  = ONE
          CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
C         ---- INPUT C^T
C
        IF(MSUBE.GT.0) THEN
          DO J = 1,NVAR
            jccloop: DO I = JCOLC(J),JCOLC(J+1)-1
              IF(CMAT(I).EQ.ZERO) cycle jccloop
              IRWC = IROWC(I)
              JCLC = J
              IROW = JCLC + NCOLA
              JCOL = IRWC
              VAL = CMAT(I)
              CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
              IF(IERLSM.NE.0) GO TO 600
            enddo jccloop
          enddo
        ENDIF
C
        IF(MUOPTN.EQ.2) THEN
C
C         ---- INPUT (B^T)(D_b)^(-1)e
C
          DO I = 1,NROWA
            IROW = NCOLA + I
            JCOL = MSUBE + 1
            VAL = SCRTCH(I)
            CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
            IF(IERLSM.NE.0) GO TO 600
          enddo
C
        ELSEIF(MUOPTN.GE.3) THEN
C
C         ---- INPUT B^T
C
          DO J = 1,NVAR
            DO I = JCOLB(J),JCOLB(J+1)-1
              IRWB = IROWB(I)
              JCLB = J
              IROW = JCLB + NCOLA
              JCOL = IRWB + MSUBE
              VAL = BMAT(I)
              CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
              IF(IERLSM.NE.0) GO TO 600
            enddo
          enddo
C
C         ---- INPUT D_B AND -e
C
          DO I = 1,MSUBB
C
            IROW = NVAR + NCOLA + I
            JCOL = I + MSUBE
            VAL = BVEC(I)
            CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
            IF(IERLSM.NE.0) GO TO 600
C
            IF(MUOPTN.EQ.3) THEN
              IROW = NVAR + NCOLA + I
              JCOL = MSUBE + MSUBB + 1
              VAL = -ONE
              CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
              IF(IERLSM.NE.0) GO TO 600
            ENDIF
C
          enddo
C
        ENDIF
C
C         ---- ZERO DIAGONALS ON THE (2,2) BLOCK 
C
        DO I = 1,NROWA
          IROW = NCOLA + I
          JCOL = NCOLA + I
          VAL = ZERO
          CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
      ELSE
C
C ------------------------------------------------------------------
C -------------OVERDETERMINED SYSTEM--------------------------------
C ------------------------------------------------------------------
C
C         ---- IDENTITY ON THE (1,1) BLOCK 
C
        DO I = 1,NROWA
          IROW = I
          JCOL = I
          VAL = ONE
          CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
C         ---- INPUT C
C
        IF(MSUBE.GT.0) THEN
          DO J = 1,NVAR
            jcjloop: DO I = JCOLC(J),JCOLC(J+1)-1
              IF(CMAT(I).EQ.ZERO) cycle jcjloop
              IRWC = IROWC(I)
              JCLC = J
              IROW = IRWC + NROWA
              JCOL = JCLC
              VAL = CMAT(I)
              CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
              IF(IERLSM.NE.0) GO TO 600
            enddo jcjloop
          enddo
        ENDIF
C
        IF(MUOPTN.EQ.2) THEN
C
C         ---- INPUT [(B^T)(D_b)^(-1)e]^T
C
          DO I = 1,NVAR
            IROW = NROWA + MSUBE + 1
            JCOL = I
            VAL = SCRTCH(I)
            CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
            IF(IERLSM.NE.0) GO TO 600
          enddo
C
        ELSEIF(MUOPTN.GE.3) THEN
C
C         ---- INPUT B
C
          DO J = 1,NVAR
            DO I = JCOLB(J),JCOLB(J+1)-1
              IRWB = IROWB(I)
              JCLB = J
              IROW = IRWB + MSUBE + NROWA
              JCOL = JCLB
              VAL = BMAT(I)
              CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
              IF(IERLSM.NE.0) GO TO 600
            enddo
          enddo
C
C         ---- INPUT D_B AND -e
C
          DO I = 1,MSUBB
C
            IROW = NROWA + MSUBE + I
            JCOL = NVAR + I
            VAL = BVEC(I)
            CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
            IF(IERLSM.NE.0) GO TO 600
C
            IF(MUOPTN.EQ.3) THEN
              IROW = NROWA + MSUBE + MSUBB + 1
              JCOL = NVAR + I
              VAL = -ONE
              CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
              IF(IERLSM.NE.0) GO TO 600
            ENDIF
C
          enddo
C
        ENDIF
C
C         ---- ZERO DIAGONALS ON THE (2,2) BLOCK 
C
        DO I = 1,NCOLA
          IROW = NROWA + I
          JCOL = NROWA + I
          VAL = ZERO
          CALL XDSLV1 ('A', IROW, JCOL, VAL, RWORK, LNRWRK, IERLSM)
          IF(IERLSM.NE.0) GO TO 600
        enddo
C
      ENDIF
C
      CALL XDSLVF ( RWORK, LNRWRK, IER )
C
      IF(IER.EQ.-400) THEN
        IERLSM = -833
      ELSEIF(IER.EQ.-401) THEN
        IERLSM = -934
      ELSEIF(IER.EQ.-410) THEN
        IERLSM = -635
      ENDIF
      IF ( IERLSM .NE. 0 ) GO TO 600
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C         FACTOR THE MATRIX
C
      CALL XDSLFA ( RWORK, LNRWRK, CNDNUM, INRTIA, NEEDS, IERLSM )
      MAXRFN = 0
C
      IF(NEEDS.GT.LNRWRK) THEN
        NEEDED = NEEDS
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDED,LNRWRK
        IERLSM = -67
        GO TO 560
      ENDIF
C
C         CHECK FOR ILL-CONDITIONED OR SINGULAR KT MATRIX OR EXCESSIVE
C         FILL
C
      IF  ( CNDNUM.GT.BIGCND**2 .OR. IERLSM.EQ.-503 .OR. IERLSM.EQ.-513 
     $                          .OR. IERLSM.EQ.-509)  THEN
C
C         RESET ERROR FLAG TO THE SPECIAL VALUE (-999) AND COMPUTE ESTIMATES
C
        IF(NVAR.NE.MSUBE) THEN
C
          IERLSM = -999
C
C         CONSTRUCT INITIAL MULTIPLIERS TO SATISFY COMPLEMENTARITY 
C
          DO I = 1,MSUBE
            ETAVEC(I) = -CVEC(I)/PENMU
          enddo
C
        ELSE
C
          IERLSM = 0
          ETAVEC(1:MSUBE) = ZERO
C
        ENDIF
C
        DO I = 1,MSUBB
          VLAMDA(I) = PENMU/BVEC(I)
        enddo
C
        GO TO 560
C
      ENDIF
C
      IF(IOFMFR.GT.0) CALL XDSLPS(RWORK)
C
C         GET MULTIFRONTAL STATISTICS
C
      CALL XDSLSR(RWORK,INUSE,MXUSED,TYME,OPCNTS)
      INSTAT(10) = MAX(INSTAT(10),INUSE)
C
      IF(IERLSM.NE.0) GO TO 600
C
 480  CONTINUE
C
      NCALL = NCALL + 1
C
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C ------------------------------------------------------------------
C
C         DEFINE RIGHT HAND SIDE
C
      IF(MUOPTN.EQ.1) THEN
C
        IF(UNDER) THEN
          DO I = NVAR,1,-1
            II = NEQNS - I + 1
            SCRTCH(II) = GVEC(I) - PENMU*SCRTCH(I)
          enddo
          SCRTCH(1:NCOLA) = ZERO
        ELSE
          DO I = 1,NVAR
            SCRTCH(I) = GVEC(I) - PENMU*SCRTCH(I)
          enddo
          SCRTCH(NVAR+1:NVAR+NCOLA) = ZERO
        ENDIF
C
      ELSEIF(MUOPTN.EQ.2) THEN
        SCRTCH(1:NEQNS) = ZERO
        IF(UNDER) THEN
          SCRTCH(NCOLA+1:NCOLA+NVAR) = GVEC(1:NVAR)
        ELSE
          SCRTCH(1:NVAR) = GVEC(1:NVAR)
        ENDIF
      ELSE
        SCRTCH(1:NEQNS) = ZERO
        IF(UNDER) THEN
          SCRTCH(NCOLA+1:NCOLA+NVAR) = GVEC(1:NVAR)
          IF(MUOPTN.EQ.2) SCRTCH(NCOLA+NVAR+1:NCOLA+NVAR+MSUBB) = HATMU
        ELSE
          SCRTCH(1:NVAR) = GVEC(1:NVAR)
          IF(MUOPTN.EQ.2) SCRTCH(NVAR+1:NVAR+MSUBB) = HATMU
        ENDIF
      ENDIF
C
C         SOLVE THE SYSTEM
C
      LDRHS = NVAR+2*MAXBND+MAXCON+1
      CALL XDSLSL ( 1, SCRTCH, LDRHS, RWORK, LNRWRK,needs, IERLSM )
C
      IF(NEEDS.GT.LNRWRK) THEN
        NEEDED = NEEDS
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDED,LNRWRK
        IERLSM = -67
        GO TO 560
      ENDIF
C
      IF(IERLSM.NE.0) GO TO 600
C
C ------------------------------------------------------------------
C
C         --- STORE THE SOLUTION INTO ETAVEC, VLAMDA, AND PENMU
C
      IF(UNDER) THEN
C
        ETAVEC(1:MSUBE) = SCRTCH(1:MSUBE)
C
        IF(MSUBB.GT.0) THEN
          IF(MUOPTN.EQ.2) THEN
            PENMU = SCRTCH(MSUBE+1)
          ELSEIF(MUOPTN.EQ.3) THEN
            PENMU = SCRTCH(MSUBE+MSUBB+1)
          ENDIF
          VLAMDA(1:MSUBB) = PENMU*VLAMDA(1:MSUBB)
        ELSE
          PENMU = ZEROOT
        ENDIF
C
        IF(MUOPTN.GT.2) THEN
C
          DO I = 1,MSUBB
            VLAMDA(I) = SCRTCH(MSUBE+I)
          enddo
C
        ENDIF
C
      ELSE
C
        ETAVEC(1:MSUBE) = SCRTCH(NROWA+1:NROWA+MSUBE)
C
        IF(MSUBB.GT.0) THEN
          IF(MUOPTN.EQ.2) THEN
            PENMU = SCRTCH(NROWA+MSUBE+1)
          ELSEIF(MUOPTN.EQ.3) THEN
            PENMU = SCRTCH(NROWA+MSUBE+MSUBB+1)
          ENDIF
          VLAMDA(1:MSUBB) = PENMU*VLAMDA(1:MSUBB)
        ELSE
          PENMU = ZEROOT
        ENDIF
C
        IF(MUOPTN.GT.2) THEN
C
          DO I = 1,MSUBB
            VLAMDA(I) = SCRTCH(NROWA+MSUBE+I)
          enddo
C
        ENDIF
C
      ENDIF
C
      IF(MSUBB.GT.0) THEN
        DO I = 1,MSUBB
          VLAMDA(I) = MAX(ZEROMN,VLAMDA(I))
        enddo
      ENDIF
C
C         COMPUTE THE NORM OF ||(C^T)eta + (B^T)lamda - g||
C
C
C         COMPUTE (CMAT**T)(ETAVEC) AND STORE IN SCRTCH(NVAR+1)
C
      CALL MVPSPR(11,NVAR,MSUBE,CMAT,IROWC,JCOLC,ETAVEC,SCRTCH(NVAR+1))
C
C         COMPUTE (BMAT**T)(VLAMDA) AND SAVE IN SCRTCH
C
      CALL MVPSPR(11,NVAR,MSUBB,BMAT,IROWB,JCOLB,VLAMDA,
     $              SCRTCH)
C
      GRADLN = ZERO
      DO I = 1,NVAR
        SCRTCH(I) = GVEC(I) - SCRTCH(NVAR+I) - SCRTCH(I)
        GRADLN = MAX(GRADLN,ABS(SCRTCH(I)))
      enddo
C
 560  CONTINUE
C
C         END FACTORIZATION TIMING CLOCK
C
      CALL CLKSUM(6)
      RETURN 
C      
 1001 FORMAT(5X,'STORAGE ERROR IN LSQMLT:  NEEDED =',I9,'  LNWORK =',I9)
C
C     -----------------------------------------------------------
C     ... TERMINATE ON UNEXPECTED ERROR CONDITION FROM BCSLIB-EXT
C     -----------------------------------------------------------

 600  CONTINUE
      IERSLM = -999
      GO TO 560
C
      END
