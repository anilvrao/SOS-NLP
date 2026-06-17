
      SUBROUTINE DGNSTO( XBAR   ,XLWR   ,XUPR   ,NDIM   
     $          ,NRES   ,DELF   ,RMAT   ,IROWR  ,JCOLR  
     $          ,NONZR  ,IROWH  ,JSTRH  ,NONZH  ,MPRMC  
     $          ,MCON   ,IPRMC  ,GMAT   ,IROWG  ,JCOLG  ,NONZG        
     $          ,WORK   ,NWORK  ,IWORK  ,NIWORK ,BIGELM ,IROWBG 
     $          ,JCOLBG ,LENBIG ,IERNLP )
C
C
C ======================================================================
C     HDSNzz===>dgnsto   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C *** PURPOSE...DIAGNOSTIC OUTPUT FOR THE SPARSE NLP ALGORITHM
C
C       XBAR        CURRENT VARIABLE VALUES (NDIM)
C       XLWR        VARIABLE LOWER BOUNDS (NDIM)
C       XUPR        VARIABLE UPPER BOUNDS (NDIM)
C       NDIM        NUMBER OF FREE VARIABLES
C       NRES        NUMBER OF RESIDUALS
C       DELF        GRADIENT OF FBAR AT XBAR (NDIM)
C       RMAT        RESIDUAL DERIVATIVES AT XBAR (NONZR)
C       IROWR       ROW INDICES OF JACOBIAN NONZEROS (NONZR)
C       JCOLR       COLUMN INDICES OF JACOBIAN NONZEROS (NONZR)
C       NONZR       NUMBER OF JACOBIAN NONZEROS
C       IROWH       ROW INDICES OF HESSIAN NONZEROS (NONZH)
C       JSTRH       COLUMN START INDICES OF NONZEROS (NDIM+1)
C       NONZH       NUMBER OF NONZERO HESSIAN ELEMENTS
C       MPRMC       NUMBER OF CONSTRAINTS (including ignored)
C       MCON        NUMBER OF CONSTRAINTS EXCLUDING IGNORED
C       IPRMC       CONSTRAINT PERMUTATION ARRAY 
C       GMAT        CONSTRAINT DERIVATIVES AT XBAR (NONZG)
C       IROWG       ROW INDICES OF JACOBIAN NONZEROS (NONZG)
C       JCOLG       COLUMN INDICES OF JACOBIAN NONZEROS (NONZG)
C       NONZG       NUMBER OF NONIGNORED JACOBIAN NONZEROS
C       WORK        WORK ARRAY 
C       NWORK       DIMENSION OF WORK ARRAY 
C                   (NWORK .GT. MAX(MCON,NRES,NDIM) + NDIM)
C       IWORK       INTEGER WORK ARRAY
C       NIWORK      DIMENSION OF IWORK ARRAY
C                   (NIWORK .GT. MAX(NONZR,NONZG,2*NDIM+1) )
C       BIGELM      ARRAY OF BIG ELEMENTS IN JACOBIAN (LENBIG)
C       IROWBG      ROW INDICES OF BIG ELEMENTS (LENBIG)
C       JCOLBG      COLUMN INDICES OF BIG ELEMENTS (LENBIG)
C       LENBIG      LENGTH OF BIGELM, IROWBG, AND JCOLBG
C       IERNLP      SUCCESS/ERROR CODE SET TO -147 (UNEXPECTED ERROR)
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,ONEEP1=1.0D1,
     $      ONEEM1=1.0D-1,ONEEM2=1.0D-2,ONEEM4=1.0D-4,ONEEM5=1.0D-5)
C
      DIMENSION     WORK(NWORK),GMAT(NONZG),DELF(NDIM),
     $      XBAR(NDIM),IROWG(NONZG), 
     $      JCOLG(NONZG),RMAT(NONZR), IROWR(NONZR), 
     $      JCOLR(NONZR),IROWH(NONZH) ,JSTRH(NDIM+1),
     $      XLWR(NDIM) ,XUPR(NDIM),
     $      IWORK(NIWORK) , BIGELM(LENBIG), IROWBG(LENBIG),
     $      JCOLBG(LENBIG),IPRMC(*) 
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
C
      LOGICAL ZTEST
      CHARACTER(LEN=100) BLANK,TITLE
      DATA BLANK(1:100) / ' '/
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         DISPLAY SPARSITY PATTERN 
C
      IF((IOFLAG.GE.30.OR.IOFPAT.GE.10).AND.(QPOPTN.EQ.'SPARSE'))  THEN
        IF(IOFPAT.EQ.10) THEN
          IPUINT = -IPUNLP
        ELSE
          IF(IOFPAT.EQ.IPUMF1
     $   .OR.IOFPAT.EQ.IPUMF2.OR.IOFPAT.EQ.IPUMF3.OR.IOFPAT.EQ.IPUMF4
     $   .OR.IOFPAT.EQ.IPUMF5.OR.IOFPAT.EQ.IPUMF6)
     $    THEN 
            IPUINT = -IPUNLP 
          ELSE
            IPUINT = IOFPAT
          ENDIF
        ENDIF
C
        IF(MCON.GT.0) THEN
          TITLE = BLANK
          TITLE(34:65) = 'JACOBIAN MATRIX SPARSITY PATTERN'
          IF(MCON.GT.0) CALL SPRPAT(1,IROWG,JCOLG,
     $              NDIM,MCON,NONZG,IPUINT,TITLE)
        ENDIF
C
        IF(NDIM.GT.0) THEN
          TITLE = BLANK
          IF(NRES.GT.0) THEN
            TITLE(30:69) = 'RESIDUAL HESSIAN MATRIX SPARSITY PATTERN'
          ELSE
            TITLE(34:64) = 'HESSIAN MATRIX SPARSITY PATTERN'
          ENDIF
          CALL SPRPAT(1,IROWH,JSTRH,NDIM,NDIM,NONZH,IPUINT,TITLE)
        ENDIF
C
        IF(NRES.GT.0) THEN
          TITLE = BLANK
          TITLE(30:70) = 'RESIDUAL JACOBIAN MATRIX SPARSITY PATTERN'
          CALL SPRPAT(2,IROWR,JCOLR,NDIM,NRES,NONZR,
     $              IPUINT,TITLE)
        ENDIF
C
      ENDIF
C
C         DETERMINE A TYPICAL "RANGE" FOR THE VARIABLES
C
      NFIXVR = 0
      NZRNG = 0
      DO I = 1,NDIM
C
        IF(XUPR(I).GE.BIGBND) THEN
          IF(XLWR(I).LE.-BIGBND) THEN
C
C           CASE 1;  UNBOUNDED ABOVE, UNBOUNDED BELOW
C
            WORK(I) = ABS(XBAR(I))
C
          ELSE
C
C           CASE 2:  UNBOUNDED ABOVE, BOUNDED BELOW
C
            WORK(I) = MAX(ABS(XBAR(I)),ABS(XBAR(I)-XLWR(I))
     $                       ,ABS(XLWR(I)))
C
          ENDIF
C
        ELSE
C
          IF(XLWR(I).LE.-BIGBND) THEN
C
C           CASE 3;  BOUNDED ABOVE, UNBOUNDED BELOW
C
            WORK(I) = MAX(ABS(XBAR(I)),ABS(XBAR(I)-XUPR(I))
     $                       ,ABS(XUPR(I)))
C
          ELSE
C
C           CASE 4;  BOUNDED ABOVE AND BELOW
C
            WORK(I) = ABS(XUPR(I)-XLWR(I))
C
          ENDIF
C
        ENDIF
C
C         COUNT THE NUMBER OF FIXED VARIABLES
C
        IF(XLWR(I).EQ.XUPR(I)) NFIXVR = NFIXVR + 1
C
C         SPECIAL FOR ZERO LENGTH RANGES
C
        IF(WORK(I).LE.ZEROMN.AND.XLWR(I).NE.XUPR(I)) THEN
          NZRNG = NZRNG + 1
          IWORK(NDIM + NZRNG) = I
        ENDIF
C
      enddo
C
C         COMPUTE MEDIAN RANGE
C
      CALL HDSRTN(WORK,NDIM,0,0,IWORK,IERSRT)
      IF(IERSRT.NE.0) THEN
        IERNLP = -147
        GO TO 210
      ENDIF
C
      LCNZRO = NZRNG + NFIXVR
      IMIDL = LCNZRO + (NDIM-LCNZRO)/2
      IMIDL = MAX(1,IMIDL)
      VRMEDN = WORK(IMIDL)
C
C         CHECK FOR EXCESSIVELY LARGE AND/OR SMALL VARIABLE RANGES
C
      SMLRNG = ONEEM4
      IF(LCNZRO.EQ.NDIM) THEN
        ZTEST = .TRUE.
      ELSE
        ZTEST = WORK(1+LCNZRO).LE.SMLRNG*VRMEDN
      ENDIF
      IF((ZTEST.OR.WORK(NDIM).GE.VRMEDN/SMLRNG).AND.
     $     IOFLAG.GE.10) THEN
C
        WRITE(IPUNLP,1001)
        WRITE(IPUNLP,1020)
        WRITE(IPUNLP,1019) VRMEDN
        IF(NZRNG.NE.0) THEN
          WRITE(IPUNLP,1018) 
          CALL INTOUT(IWORK(1+NDIM),NZRNG,IPUNLP)
        ENDIF
C
C         RESTORE ORIGINAL ORDER
C
        CALL HDPRMY(WORK,NDIM,IWORK,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          GO TO 210
        ENDIF
C
C       DISPLAY LARGEST AND SMALLEST RANGES
C
        CALL MAGRYT(5,NDIM,WORK,WORK(1+NDIM),
     $    IWORK,'VARIABLE RANGE MAGNITUDES',IPUNLP)
C
      ENDIF
C
C         CHECK FOR VERY LARGE AND/OR VERY SMALL COLUMN NORMS OF AUGMENTED
C         JACOBIAN (INCLUDES GRADIENT OR RESIDUAL JACOBIAN)
C        
C         COMPUTE COLUMN NORMS OF JACOBIAN
C
      WORK(1:NDIM) = ZERO
      IF(MCON.GT.0) CALL MRNSPR(-2,MCON,NDIM,GMAT,IROWG,
     $    JCOLG,WORK)
C
C         COMPUTE COLUMN NORMS OF RESIDUAL JACOBIAN
C
      IF(NRES.GT.0) THEN
C
        WORK(NDIM+1:2*NDIM) = ZERO
        CALL MRNSPR(-2,NRES,NDIM,RMAT,IROWR,JCOLR,
     $               WORK(NDIM+1))
C
      ENDIF
C
      NSMALL = 0
      NLARGE = 0
C
      DO I = 1,NDIM
C
C         ADD CONTRIBUTION OF OBJECTIVE GRADIENT OR RESIDUAL JACOBIAN
C         TO COLUMN NORM
C
        IF(NRES.GT.0) THEN
          OBJTRM = WORK(NDIM+I)
        ELSE
          OBJTRM = ABS(DELF(I))
        ENDIF
C
        WORK(I) = MAX(WORK(I),OBJTRM)
C
C         SAVE INDEX OF LARGE AND SMALL COLUMN NORMS
C
        COLNRM = WORK(I)
        IF(COLNRM.LT.ZEROOT) THEN
          NSMALL = NSMALL + 1
          IWORK(NSMALL) = I
        ELSEIF(COLNRM.GT.BIGBND) THEN
          NLARGE = NLARGE + 1
          IWORK(NDIM+NLARGE) = I
        ENDIF
C
      enddo
C
C         WRITE INDEX OF SMALL COLUMNS
C
      IF(IOFLAG.GE.10.AND.NSMALL.NE.0) THEN
        WRITE(IPUNLP,1002) 
        CALL INTOUT(IWORK,NSMALL,IPUNLP)
      ENDIF
C
C         WRITE INDEX OF LARGE COLUMNS
C
      IF(IOFLAG.GE.10.AND.NLARGE.NE.0) THEN
        WRITE(IPUNLP,1003) 
        CALL INTOUT(IWORK(1+NDIM),NLARGE,IPUNLP)
      ENDIF
C
      NSMALC = 0
      NLARGC = 0
      IF(IOFLAG.GE.20) THEN
        SMLNRM = ONEEM1
      ELSE
        SMLNRM = ONEEM2
      ENDIF
      BIGNRM = ONE/SMLNRM
C
      DO I = 1,NDIM
C
C         SAVE INDEX OF LARGE AND SMALL COLUMN NORMS
C
        COLNRM = WORK(I)
        IF(COLNRM.LT.SMLNRM) THEN
          NSMALC = NSMALC + 1
          IWORK(NSMALC) = I
        ELSEIF(COLNRM.GT.BIGNRM) THEN
          NLARGC = NLARGC + 1
          IWORK(NDIM+NLARGC) = I
        ENDIF
C
      enddo
C
      IF(IOFLAG.GE.10.AND.(NSMALC.NE.0.OR.NLARGC.NE.0)) THEN
        WRITE(IPUNLP,1001)
        WRITE(IPUNLP,1017)
        IF(ABS(BIGELM(1)).GT.BIGNRM) THEN
          WRITE(IPUNLP,1021)
          lenloop: DO I = 1,LENBIG
            TITLE = BLANK
            IF(IROWBG(I).EQ.0) cycle lenloop
            WRITE(TITLE(1:10),'(I10)') IROWBG(I)
            CALL HHADJF(TITLE(1:10),' ',' ','L',ISHFT,IERSHF)
            WRITE(TITLE(20:30),'(I10)') JCOLBG(I)
            CALL HHADJF(TITLE(20:30),' ',' ','L',ISHFT,IERSHF)
            WRITE(TITLE(40:56),'(SP,1PG16.8)') BIGELM(I)
            CALL HHADJF(TITLE(40:56),' ',' ','L',ISHFT,IERSHF)
            WRITE(IPUNLP,1022) TITLE(1:56)
          enddo lenloop
        ENDIF
      ENDIF
C
C         WRITE INDEX OF SMALL COLUMNS
C
      IF(IOFLAG.GE.10.AND.NSMALC.NE.0) THEN
        WRITE(IPUNLP,1014) SMLNRM 
        CALL IRLOUT(IWORK,WORK,NSMALC,IPUNLP)
        WRITE(IPUNLP,1016)
      ENDIF
C
C         WRITE INDEX OF LARGE COLUMNS
C
      IF(IOFLAG.GE.10.AND.NLARGC.NE.0) THEN
        WRITE(IPUNLP,1015) BIGNRM
        CALL IRLOUT(IWORK(1+NDIM),WORK,NLARGC,IPUNLP)
        WRITE(IPUNLP,1016)
      ENDIF
C
C         CHECK FOR VERY LARGE AND/OR VERY SMALL ROW NORMS OF JACOBIAN
C          
      IF(MCON.GT.0) THEN
C
        CALL MRNSPR(-1,MCON,NDIM,GMAT,IROWG,JCOLG,WORK)
C
C         LOAD NEGATIVE VALUE INTO IGNORED ROW NORMS
C
        DO I=MCON+1,MPRMC
          WORK(I) = -ONE
        ENDDO
C
C         PERMUTE ROW NORM VECTOR TO EXTERNAL ORDER
C
         CALL HDPRMY(WORK,MPRMC,IPRMC,IERP)
         IF(IERP.NE.0) THEN
           IERNLP = -147
           GO TO 210
         ENDIF
C
      ENDIF
C
      NSMALL = 0
      NLARGE = 0
C
      DO I = 1,MPRMC
C
C         SAVE INDEX OF LARGE AND SMALL ROW NORMS
C
        ROWNRM = WORK(I)
        IF(ROWNRM.GT.-ONE.AND.ROWNRM.LT.ZEROOT) THEN
          NSMALL = NSMALL + 1
          IWORK(NSMALL) = I
        ELSEIF(ROWNRM.GT.BIGBND) THEN
          NLARGE = NLARGE + 1
          IWORK(MCON+NLARGE) = I
        ENDIF
C
      enddo
C
C         WRITE INDEX OF SMALL ROWS
C
      IF(IOFLAG.GE.10.AND.NSMALL.NE.0) THEN
        WRITE(IPUNLP,1004) 
        CALL INTOUT(IWORK,NSMALL,IPUNLP)
      ENDIF
C
C         WRITE INDEX OF LARGE ROWS
C
      IF(IOFLAG.GE.10.AND.NLARGE.NE.0) THEN
        WRITE(IPUNLP,1006) 
        CALL INTOUT(IWORK(MCON+1),NLARGE,IPUNLP)
      ENDIF
C
      NSMALR = 0
      NLARGR = 0
C
      DO I = 1,MPRMC
C
C       SAVE INDEX OF LARGE AND SMALL ROW NORMS
C
        ROWNRM = WORK(I)
        IF(ROWNRM.GT.-ONE.AND.ROWNRM.LT.SMLNRM) THEN
          NSMALR = NSMALR + 1
          IWORK(NSMALR) = I
        ELSEIF(ROWNRM.GT.BIGNRM) THEN
          NLARGR = NLARGR + 1
          IWORK(MCON+NLARGR) = I
        ENDIF
C
      enddo
C
      IF(IOFLAG.GE.10.AND.(NSMALR.NE.0.OR.NLARGR.NE.0)
     $    .AND.(NSMALC.EQ.0.AND.NLARGC.EQ.0)) THEN
          WRITE(IPUNLP,1001)
          WRITE(IPUNLP,1017)
          IF(ABS(BIGELM(1)).GT.BIGNRM) THEN
            WRITE(IPUNLP,1021)
            lnbigloop: DO I = 1,LENBIG
              TITLE = BLANK
              IF(IROWBG(I).EQ.0) cycle lnbigloop
              WRITE(TITLE(1:10),'(I10)') IROWBG(I)
              CALL HHADJF(TITLE(1:10),' ',' ','L',ISHFT,IERSHF)
              WRITE(TITLE(20:30),'(I10)') JCOLBG(I)
              CALL HHADJF(TITLE(20:30),' ',' ','L',ISHFT,IERSHF)
              WRITE(TITLE(40:56),'(SP,1PG16.8)') BIGELM(I)
              CALL HHADJF(TITLE(40:56),' ',' ','L',ISHFT,IERSHF)
              WRITE(IPUNLP,1022) TITLE(1:56)
            enddo lnbigloop
          ENDIF
      ENDIF
C
      IF(IERNLP.LT.0) GO TO 210
C
C         WRITE INDEX OF SMALL ROWS
C
      IF(IOFLAG.GE.10.AND.NSMALR.NE.0) THEN
        WRITE(IPUNLP,1008) SMLNRM 
        CALL IRLOUT(IWORK,WORK,NSMALR,IPUNLP)
        WRITE(IPUNLP,1016)
      ENDIF
C
C         WRITE INDEX OF LARGE ROWS
C
      IF(IOFLAG.GE.10.AND.NLARGR.NE.0) THEN
        WRITE(IPUNLP,1010) BIGNRM
        CALL IRLOUT(IWORK(MCON+1),WORK,NLARGR,IPUNLP)
        WRITE(IPUNLP,1016)
      ENDIF
C
      IF(NRES.GT.0) THEN
C
C         CHECK FOR VERY LARGE AND/OR VERY SMALL ROW NORMS OF RESIDUAL
C         JACOBIAN
C          
        CALL MRNSPR(-1,NRES,NDIM,RMAT,IROWR,JCOLR,WORK)
C
        NSMALL = 0
        NLARGE = 0
C
        DO I = 1,NRES
C
C         SAVE INDEX OF LARGE AND SMALL ROW NORMS
C
          ROWNRM = WORK(I)
          IF(ROWNRM.LT.ZEROOT) THEN
            NSMALL = NSMALL + 1
            IWORK(NSMALL) = I
          ELSEIF(ROWNRM.GT.BIGBND) THEN
            NLARGE = NLARGE + 1
            IWORK(NRES+NLARGE) = I
          ENDIF
C
        enddo
C
C         WRITE INDEX OF SMALL ROWS
C
        IF(IOFLAG.GE.10.AND.NSMALL.NE.0) THEN
          WRITE(IPUNLP,1005) 
          CALL INTOUT(IWORK,NSMALL,IPUNLP)
        ENDIF
C
C         WRITE INDEX OF LARGE ROWS
C
        IF(IOFLAG.GE.10.AND.NLARGE.NE.0) THEN
          WRITE(IPUNLP,1007) 
          CALL INTOUT(IWORK(NRES+1),NLARGE,IPUNLP)
        ENDIF
C
        NSMALR = 0
        NLARGR = 0
C
        DO I = 1,NRES
C
C         SAVE INDEX OF LARGE AND SMALL ROW NORMS
C
          ROWNRM = WORK(I)
          IF(ROWNRM.LT.SMLNRM) THEN
            NSMALR = NSMALR + 1
            IWORK(NSMALR) = I
          ELSEIF(ROWNRM.GT.BIGNRM) THEN
            NLARGR = NLARGR + 1
            IWORK(NRES+NLARGR) = I
          ENDIF
C
        enddo
C
        IF(IOFLAG.GE.10.AND.(NSMALR.NE.0.OR.NLARGR.NE.0)
     $    .AND.(NSMALC.EQ.0.AND.NLARGC.EQ.0)) THEN
          WRITE(IPUNLP,1001)
          WRITE(IPUNLP,1017)
          IF(ABS(BIGELM(1)).GT.BIGNRM) THEN
            WRITE(IPUNLP,1021)
            lbigloop: DO I = 1,LENBIG
              TITLE = BLANK
              IF(IROWBG(I).EQ.0) cycle lbigloop
              WRITE(TITLE(1:10),'(I10)') IROWBG(I)
              CALL HHADJF(TITLE(1:10),' ',' ','L',ISHFT,IERSHF)
              WRITE(TITLE(20:30),'(I10)') JCOLBG(I)
              CALL HHADJF(TITLE(20:30),' ',' ','L',ISHFT,IERSHF)
              WRITE(TITLE(40:56),'(SP,1PG16.8)') BIGELM(I)
              CALL HHADJF(TITLE(40:56),' ',' ','L',ISHFT,IERSHF)
              WRITE(IPUNLP,1022) TITLE(1:56)
            enddo lbigloop
          ENDIF
        ENDIF
C
        IF(IERNLP.LT.0) GO TO 210
C
C         WRITE INDEX OF SMALL ROWS
C
        IF(IOFLAG.GE.10.AND.NSMALR.NE.0) THEN
          WRITE(IPUNLP,1009) SMLNRM 
          CALL IRLOUT(IWORK,WORK,NSMALR,IPUNLP)
          WRITE(IPUNLP,1016)
        ENDIF
C
C         WRITE INDEX OF LARGE ROWS
C
        IF(IOFLAG.GE.10.AND.NLARGR.NE.0) THEN
          WRITE(IPUNLP,1011) BIGNRM
          CALL IRLOUT(IWORK(NRES+1),WORK,NLARGR,IPUNLP)
          WRITE(IPUNLP,1016)
        ENDIF
C
      ELSE
C       
C         COMPUTE ROW NORM CORRESPONDING TO OBJECTIVE GRADIENT
C
        GRDNRM = DAMAX(NDIM,DELF,1)
        IF(IOFLAG.GE.10.AND.GRDNRM.GT.BIGNRM) 
     $      WRITE(IPUNLP,1012) GRDNRM,BIGNRM
        IF(IOFLAG.GE.10.AND.GRDNRM.LT.SMLNRM) 
     $      WRITE(IPUNLP,1013) GRDNRM,SMLNRM
C
      ENDIF
C
C ======================================================================
C ================= END OF INITIAL DIAGNOSTIC ANALYSIS =================
C ======================================================================
C
 210  CONTINUE
C
      RETURN
C
C     ------------------------------------------------------------------
C             FORMAT STATEMENTS                                         
C     ------------------------------------------------------------------
C
C
C
 1001 FORMAT(T3,'*',T106,'*',/2X,104('*'))
C
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING COLUMNS OF T
     $HE AUGMENTED JACOBIAN MAY BE ZERO',T106,'*',
     $   /T3,'*',T106,'*')
C
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING COLUMNS OF T
     $HE AUGMENTED JACOBIAN MAY BE TOO LARGE',T106,'*',
     $   /T3,'*',T106,'*')
C
 1004 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING ROWS OF THE
     $JACOBIAN MAY BE ZERO',T106,'*',
     $   /T3,'*',T106,'*')
C
 1005 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING ROWS OF THE
     $RESIDUAL JACOBIAN MAY BE ZERO',T106,'*',
     $   /T3,'*',T106,'*')
C
 1006 FORMAT(T3,'*',T106,'*'/T3,1H*,T11,'.....THE FOLLOWING ROWS OF THE
     $JACOBIAN MAY BE TOO LARGE',T106,'*',
     $   /T3,'*',T106,'*')
C
 1007 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING ROWS OF THE
     $RESIDUAL JACOBIAN MAY BE TOO LARGE',T106,'*',
     $   /T3,'*',T106,'*')
C
 1008 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....ROWS OF THE JACOBIAN WITH 
     $NORM LESS THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1009 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....ROWS OF THE RESIDUAL JACOB
     $IAN WITH NORM LESS THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1010 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....ROWS OF THE JACOBIAN WITH 
     $NORM GREATER THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1011 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....ROWS OF THE RESIDUAL JACOB
     $IAN WITH NORM GREATER THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1012 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....GRADIENT NORM =',1PG10.3,
     $' GREATER THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1013 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....GRADIENT NORM =',1PG10.3,
     $' LESS THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1014 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....COLUMNS OF THE AUGMENTED J
     $ACOBIAN WITH NORM LESS THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1015 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....COLUMNS OF THE AUGMENTED J
     $ACOBIAN WITH NORM GREATER THAN',5X,1PG10.3,T106,'*',
     $   /T3,'*',T106,'*')
C
 1016 FORMAT(T3,'*',T106,'*')
C
 1017 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'SCALE INFORMATION',T106,'*',
     $  /T3,'*',T106,'*')
C
 1018 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....THE FOLLOWING VARIABLES HA
     $VE ZERO RANGE',T106,'*',/T3,'*',T106,'*')
C
 1019 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....SIGNIFICANT DISPARITY IN V
     $ARIABLE RANGES;  MEDIAN RANGE =',1PG16.8,T106,'*'/T3,'*',T106,'*')
C
 1020 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'VARIABLE SCALE INFORMATION'
     $   ,T106,'*',/T3,'*',T106,'*')
C
 1021 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'.....LARGEST JACOBIAN ELEMENTS'
     $  ,T106,'*'/T3,'*',T106,'*'/T3,'*',T16,'ROW',T35,'COLUMN',T55
     $  ,'VALUE',T106,'*',/T3,'*',T106,'*')
C
 1022 FORMAT(T3,'*',T16,A,T106,'*')
C
      END
