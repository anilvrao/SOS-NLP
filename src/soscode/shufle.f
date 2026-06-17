      SUBROUTINE SHUFLE( IPRMC  ,IPRMG  ,IPRMX  ,IPRMH  ,MAXCON  
     $          ,MCON   ,MEQUAL ,MIGNOR ,AVEC   ,ALWR   ,AUPR    
     $          ,VECLAM ,ISTATA ,AMAT   ,IROWA  ,JCOLA  ,NONZA  
     $          ,DELF   ,XVEC   ,XLWR   ,XUPR   ,VECNU  ,ISTATX   
     $          ,NDIM   ,NFREE  ,HMAT   ,IROWH  ,JCOLH  ,NONZH    
     $          ,FBAR   ,CVEC   ,MSUBE  ,ALWRI  ,AUPRI  ,ETAVEC      
     $          ,BMAT   ,IROWB  ,JCOLB  ,NONZB  ,CMAT   ,IROWC    
     $          ,JCOLC  ,NONZC  ,GVEC   ,YVEC   ,NVAR   ,NVARMX    
     $          ,XLWRI  ,XUPRI  ,BVEC   ,MSUBB  ,MAXBND ,NXUPR     
     $          ,NXLWR  ,NAUPR  ,NALWR  ,IBOUND ,VLAMDA ,WMAT     
     $          ,IROWW  ,JCOLW  ,NONZW  ,FNLP   ,IWORK  ,LNIWRK  
     $          ,RWORK  ,LNRWRK ,IERNLP ,IFERR  ,RELAX  ,FZMODE ) 
C
C
C ======================================================================
C     SHUFLE===>SHUFLE   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C *** PURPOSE
C
C        THIS ROUTINE:
C 
C        (A) CONSTRUCTS THE TRANSFORMATION FROM THE USER DATA (EXTERNAL) 
C            FORMAT TO THE BARRIER ALGORITHM (INTERNAL) FORMAT
C        (B) APPLIES THE TRANSFORMATION TO THE INITIAL POINT
C
C *** CALLING ARGUMENTS
C
C        PERMUTATION ARRAYS
C
C        IPRMC    O    INTEGER CONSTRAINT PERMUTATION ARRAY (MAXCON)
C        IPRMG    O    INTEGER JACOBIAN PERMUTATION ARRAY (NONZA)
C        IPRMX    O    INTEGER VARIABLE PERMUTATION ARRAY (NDIM)
C        IPRMH    O    INTEGER HESSIAN PERMUTATION ARRAY (NONZH)
C
C        EXTERNAL FORMAT
C
C        MAXCON   I    MAXIMUM NUMBER OF CONSTRAINTS
C        MCON     I    NUMBER OF USER CONSTRAINTS
C        MEQUAL   I    NUMBER OF USER EQUALITY CONSTRAINTS
C        MIGNOR   I    NUMBER OF USER IGNORED CONSTRAINTS
C        AVEC     I    CONSTRAINTS AT XVEC (MAXCON)
C        ALWR     I    CONSTRAINTS LOWER BOUNDS (MAXCON)
C        AUPR     I    CONSTRAINTS LOWER BOUNDS (MAXCON)
C        VECLAM   I    LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C        ISTATA   I    INTEGER CONSTRAINT STATUS (MAXCON)
C        AMAT     I    CONSTRAINT DERIVATIVES AT XVEC (NONZA)
C        IROWA    I    ROW INDICES OF JACOBIAN NONZEROS (NONZA)
C        JCOLA    I    COLUMN INDICES OF JACOBIAN NONZEROS (NONZA)
C        NONZA    I    NUMBER OF JACOBIAN NONZEROS
C        DELF     I    GRADIENT AT XVEC (NDIM)
C        XVEC     I    CURRENT POINT (NDIM)
C        XLWR     I    VARIABLE LOWER BOUNDS (NDIM)
C        XUPR     I    VARIABLE UPPER BOUNDS (NDIM)
C        VECNU    I    LAGRANGE MULTIPLIERS FOR BOUNDS (NDIM)
C        ISTATX   I    INTEGER VARIABLE STATUS (NDIM)
C        NDIM     I    NUMBER OF USER VARIABLES
C        NFREE    I    NUMBER OF FREE USER VARIABLES
C        HMAT     I    HESSIAN MATRIX AT XVEC (NONZH)
C        IROWH    I    ROW INDICES OF HESSIAN NONZEROS (NONZH)
C        JCOLH    I    COLUMN START INDICES OF HESSIAN (NDIM+1)
C        NONZH    I    NUMBER OF HESSIAN NONZEROS
C        FBAR     I    OBJECTIVE FUNCTION
C
C        INTERNAL FORMAT
C
C        CVEC     O    EQUALITY CONSTRAINTS AT YVEC (MAXCON)
C        MSUBE    O    NUMBER OF EQUALITY CONSTRAINTS MSUBE = MCON - MIGNOR
C        ALWRI    O    PERMUTED USER CONSTRAINT LOWER BOUNDS (MAXCON)
C        AUPRI    O    PERMUTED USER CONSTRAINT UPPER BOUNDS (MAXCON)
C        ETAVEC   O    LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MAXCON)
C        BMAT     O    BOUND DERIVATIVES AT YVEC (NONZB)
C        IROWB    O    ROW INDICES OF BOUND JACOBIAN NONZEROS (NONZB)
C        JCOLB    O    COLUMN INDICES OF BOUND JACOBIAN NONZEROS (NVAR+1)
C        NONZB    I    NUMBER OF BOUND JACOBIAN NONZEROS
C        CMAT     O    CONSTRAINT DERIVATIVES AT YVEC (NONZC)
C        IROWC    O    ROW INDICES OF JACOBIAN NONZEROS (NONZC)
C        JCOLC    O    COLUMN INDICES OF JACOBIAN NONZEROS (NVAR+1)
C        NONZC    I    NUMBER OF JACOBIAN NONZEROS
C        GVEC     O    GRADIENT AT YVEC (NVARMX)
C        YVEC     O    CURRENT POINT (NVARMX)
C        NVAR     I    NUMBER OF INTERNAL VARIABLES;  NVAR = NFREE + NSLK
C        NVARMX   I    MAX(NVAR,NDIM)
C        XLWRI    O    PERMUTED USER VARIABLE LOWER BOUNDS (NDIM)
C        XUPRI    O    PERMUTED USER VARIABLE UPPER BOUNDS (NDIM)
C        BVEC     O    BOUND INEQUALITIES (MAXBND)
C        MSUBB    I    NUMBER OF BOUNDS
C        MAXBND   I    MAXIMUM NUMBER OF BOUNDS MAX(MSUBB,1)
C        NXUPR    I    NUMBER OF VARIABLE UPPER BOUNDS
C        NXLWR    I    NUMBER OF VARIABLE LOWER BOUNDS
C        NAUPR    I    NUMBER OF CONSTRAINT UPPER BOUNDS
C        NALWR    I    NUMBER OF CONSTRAINT LOWER BOUNDS
C        IBOUND   O    INTEGER BOUND DEFINITION (MAXBND)
C        VLAMDA   O    LAGRANGE MULTIPLIERS FOR BOUNDS (MAXBND)
C        WMAT     O    HESSIAN MATRIX AT YVEC (NONZW)
C        IROWW    O    ROW INDICES OF HESSIAN NONZEROS (NONZW)
C        JCOLW    O    COLUMN START INDICES OF HESSIAN (NVAR+1)
C        NONZW    I    NUMBER OF HESSIAN NONZEROS
C        FNLP     O    NLP OBJECTIVE FUNCTION
C
C        WORK ARRAYS
C
C        IWORK    IO   INTEGER WORK ARRAY (LNIWRK)
C        LNIWRK   I    LENGTH OF IWORK = MAX(2*(NDIM + MCON),NONZB)
C        RWORK    IO   REAL WORK ARRAY (LNRWRK)
C        LNRWRK   I    LENGTH OF WORK = MAX(NONZH,NONZA,2*(NDIM + MCON))
C        IERNLP   O    SUCCESS/ERROR CODE
C        IFERR    O    FUNCTION ERROR FLAG
C        RELAX    I    LOGICAL FLAG: TRUE FOR RELAXATION MODE; ELSE FALSE
C        FZMODE   I    LOGICAL FLAG: TRUE FOR FEASIBILITY MODE; ELSE FALSE
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      DIMENSION IPRMC(MAXCON)  ,IPRMG(NONZA)       ,IPRMX(NDIM)
     &         ,IPRMH(NONZH)
      DIMENSION AVEC(MAXCON)   ,ALWR(MAXCON)       ,AUPR(MAXCON)  
     &         ,VECLAM(MAXCON) ,ISTATA(MAXCON)     ,AMAT(NONZA)    
     &         ,IROWA(NONZA)   ,JCOLA(NONZA)       ,DELF(NDIM)     
     &         ,XVEC(NDIM)     ,XLWR(NDIM)         ,XUPR(NDIM)
     &         ,VECNU(NDIM)    ,ISTATX(NDIM)       ,HMAT(NONZH)            
     &         ,IROWH(NONZH)   ,JCOLH(NDIM+1)      ,CVEC(MAXCON)        
     &         ,ALWRI(MAXCON)  ,AUPRI(MAXCON)      ,ETAVEC(MAXCON)      
     &         ,BMAT(NONZB)    ,IROWB(NONZB)           
     &         ,JCOLB(NVAR+1)  ,CMAT(NONZC)        ,IROWC(NONZC)            
     &         ,JCOLC(NVAR+1)  ,GVEC(NVARMX)       ,YVEC(NVARMX)               
     &         ,XUPRI(NDIM)    ,BVEC(MAXBND)       ,IBOUND(MAXBND)      
     &         ,XLWRI(NDIM)    ,VLAMDA(MAXBND)              
     &         ,WMAT(NONZW)    ,IROWW(NONZW)       ,JCOLW(NVAR+1)
     &         ,IWORK(LNIWRK)  ,RWORK(LNRWRK)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,POINT5=5.D-1)
C
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
      LOGICAL RELAX,FZMODE
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
      INCLUDE '../commons/BARNLP.CMN'
C-------------------------------------------------------------
C
      IFERR = 0
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C         COMPUTE THE DESIRED "NEW" LOCATION FOR THE VARIABLES AND 
C         TEMPORARILY SAVE IN IPRMX
C
      K1 = NFREE
      K2 = 0
      DO I = 1,NDIM
C
        IF(ISTATX(I).EQ.3) THEN
C
C         PUT FIXED VARIABLES AT END OF THE LIST
C
          K1 = K1 + 1
          IPRMX(I) = K1
C 
        ELSE
C
C         PUT FREE VARIABLES AT THE FRONT OF THE LIST
C
          K2 = K2 + 1
          IPRMX(I) = K2
C
        ENDIF
      enddo
C
C         LOOP OVER ALL HESSIAN NONZEROS AND COMPUTE NEW LOCATION NUMBER
C         REPLACE ROW INDEX WITH NEW VALUE
C
      XFREE = NFREE
      XNDIM = NDIM
      DO J = 1,NDIM
        DO I = JCOLH(J),JCOLH(J+1)-1
          IROWW(I) = IPRMX(IROWH(I))
          XROW = IROWW(I)
          JCOL = IPRMX(J)
          XCOL = JCOL
          IF(IROWW(I).GT.NFREE) XCOL = XCOL + XFREE
          RWORK(I) = XROW + (XCOL-ONE)*XNDIM
        enddo
      enddo
C
C         SORT THE COLUMN INDEX ARRAY RWORK INTO ASCENDING ORDER
C         THEREBY CONSTRUCTING THE MASTER HESSIAN PERMUTATION ARRAY
C
      CALL HDSRTN(RWORK,NONZH,0,0,IPRMH,IERP)
C
      IF(IERP.NE.0) THEN
        IERNLP = -147
        RETURN
      ENDIF
C
C         APPLY THE PERMUTATION TO THE HESSIAN ROW STORAGE ARRAY
C         CONVERTING IT TO INTERNAL ORDER.
C         NOTE:  SINCE THE HESSIAN HAS NOT YET BEEN EVALUATED 
C         HMAT AND WMAT DO NOT NEED TO BE MODIFIED.
C
      CALL HJPRMX(IROWW,NONZH,IPRMH,IERP)
      IF(IERP.NE.0) THEN
        IERNLP = -147
        RETURN
      ENDIF
C
C         CONSTRUCT THE NEW JCOLW ARRAY BY SCANNING THE ELEMENTS 
C         NOTE THIS ASSUMES EACH COLUMN HAS AT LEAST ONE NONZERO
C
      JCOL = 1
      JCOLW(JCOL) = 1
      nzhloop: DO I = 1,NONZH-1
        ROWI = IROWW(I)
        XCOLI = (RWORK(I) - ROWI)/XNDIM + ONE
        ROWP = IROWW(I+1)
        XCOLP = (RWORK(I+1) - ROWP)/XNDIM + ONE
        IF(NINT(XCOLP).GT.NINT(XCOLI)) THEN
          JCOL = JCOL + 1
          JCOLW(JCOL) = I + 1
        ENDIF
        IF(JCOL.EQ.NFREE) THEN
          NZW = I + 1
          exit nzhloop
        ELSEIF(JCOL.GT.NFREE) THEN
          NZW = I 
          exit nzhloop
         ENDIF
      enddo nzhloop
C
C         NOTE THE DEFINITION OF WMAT IS FINISHED BELOW
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      LNJCOL = MAX(NONZA,NDIM+1)
      MINEQL = 0
      MSUBE = MCON - MIGNOR
      IF(MSUBE.GT.0) THEN
C
C         COMPUTE THE DESIRED "NEW" LOCATION FOR THE CONSTRAINTS AND 
C         TEMPORARILY SAVE IN IWORK
C
        MINEQL = MSUBE - MEQUAL
C
        K1 = 0
        K2 = 0
        K3 = 0
        DO I = 1,MCON
C
          IF(ISTATA(I).EQ.3) THEN
C
C           EQUALITY CONSTRAINTS
C
            K1 = K1 + 1
            IWORK(I) = K1
C 
          ELSEIF(ISTATA(I).LT.3) THEN
C
C           INEQUALITY CONSTRAINTS
C
            K2 = K2 + 1
            IWORK(I) = MEQUAL+K2
C
          ELSE
C
C           IGNORED CONSTRAINTS
C
            K3 = K3 + 1
            IWORK(I) = MSUBE+K3
          ENDIF
C
        enddo
C
C         NOW CONSTRUCT THE CONSTRAINT PERMUTATION ARRAY IPRMC BUT LEAVE
C         IWORK IN THE ORIGINAL ORDER
C
        CALL HJSRTN(IWORK,MCON,1,0,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         LOOP OVER ALL JACOBIAN NONZEROS AND COMPUTE NEW LOCATION NUMBER
C         REPLACE ROW AND COLUMN INDICES WITH NEW VALUES
C
        XFREE = NFREE
        XMSUBE = MSUBE
        NZC = 0
        DO I = 1,NONZA
          IROWC(I) = IWORK(IROWA(I))
          XROW = IROWC(I)
          JCOL = IPRMX(JCOLA(I))
          XCOL = JCOL
          IF(IROWC(I).GT.MSUBE) THEN
            XCOL = XCOL + XFREE
          ELSE
            NZC = NZC + 1
          ENDIF
          RWORK(I) = XROW + (XCOL-ONE)*XMSUBE
        enddo
C
C         SORT THE COLUMN INDEX ARRAY RWORK INTO ASCENDING ORDER
C         THEREBY CONSTRUCTING THE MASTER JACOBIAN PERMUTATION ARRAY
C
        CALL HDSRTN(RWORK,NONZA,0,0,IPRMG,IERP)
C
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         APPLY THE PERMUTATION TO THE JACOBIAN AND ROW ARRAYS, THUS
C         CONVERTING THEM TO INTERNAL ORDER
C
        CALL HDPRMX(AMAT,NONZA,IPRMG,IERP)
        cmat(1:nzc) = amat(1:nzc)
        CALL HDPRMY(AMAT,NONZA,IPRMG,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
        CALL HJPRMX(IROWC,NONZA,IPRMG,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         CONSTRUCT THE NEW JCOLC ARRAY BY SCANNING THE ELEMENTS 
C         FIRST COUNT THE NONZEROS IN EACH COLUMN
C
        iwork(1:nfree) = 0
        DO I = 1,NONZA
          ROWI = IROWC(I)
          XCOLI = (RWORK(I) - ROWI)/XMSUBE + ONE
          JCOLI = NINT(XCOLI)
          IF(JCOLI.LE.NFREE) IWORK(JCOLI) = IWORK(JCOLI) + 1
        enddo
C
C         ADD UP NONZEROS TO FORM JCOLC
C
        jcolc(1:nvar+1) = 0
        JCOLC(1) = 1
        DO JCOL = 2,NFREE+1
          JCOLC(JCOL) = JCOLC(JCOL-1) + IWORK(JCOL-1)
        enddo
C
        NZC = JCOLC(NFREE+1) - 1
C
C         PUT -I INTO THE (2,2) BLOCK
C
        KCOL = NFREE
        DO I = 1,MINEQL
          KCOL = KCOL + 1
          NZC = NZC + 1
          CMAT(NZC) = -ONE
          IROWC(NZC) = MEQUAL + I
          JCOLC(KCOL+1) = JCOLC(KCOL) + 1
        enddo
C
C         PUT -I INTO THE (1,3) BLOCK
C
        DO I = 1,MEQUAL
          KCOL = KCOL + 1
          NZC = NZC + 1
          CMAT(NZC) = -ONE
          IROWC(NZC) = I
          JCOLC(KCOL+1) = JCOLC(KCOL) + 1
        enddo
C
C         PUT I INTO THE (1,4) BLOCK
C
        DO I = 1,MEQUAL
          KCOL = KCOL + 1
          NZC = NZC + 1
          CMAT(NZC) = ONE
          IROWC(NZC) = I
          JCOLC(KCOL+1) = JCOLC(KCOL) + 1
        enddo
C
C         FILL OUT THE REMAINING COLUMN START INDICES
C
        LENCOL = NVAR-KCOL
        DO I=2,LENCOL+1
          JCOLC(KCOL+I) = JCOLC(KCOL+1)
        ENDDO
C
C         APPLY ORDERING TO ALWR, AND AUPR
C         TO CREATE ALWRI, AND AUPRI
C
        alwri(1:mcon) = alwr(1:mcon)
        CALL HDPRMX(ALWRI,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
        aupri(1:mcon) = aupr(1:mcon)
        CALL HDPRMX(AUPRI,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         APPLY ORDERING TO AVEC TO CREATE CVEC
C
        cvec(1:mcon) = avec(1:mcon)
        CALL HDPRMX(CVEC,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         APPLY ORDERING TO VECLAM TO CREATE ETAVEC
C
        etavec(1:mcon) = veclam(1:mcon)
        CALL HDPRMX(ETAVEC,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         IF NECESSARY PERMUTE THE DIAGNOSTIC FUNCTION VALUE
C
        IF(LYNFNC.NE.0) LYNFNC = IPRMC(LYNFNC)
C
      ELSE
C
        MSUBE = 0
        DO I = 1,MCON
          IPRMC(I) = I
        enddo
        DO I = 1,NONZA
          IPRMG(I) = I
        enddo
        jcolc(1:nvar+1) = 1
C
      ENDIF
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
C
C         COMPLETE THE DEFINITION OF WMAT BY FILLING THE LOWER
C         RIGHT BLOCK WITH ZERO
C
      NSLK = NVAR - NFREE
C
      DO I = NFREE+1,NVAR
        NZW = NZW + 1
        JCOLW(I) = NZW
        IROWW(NZW) = I
        WMAT(NZW) = ZERO
      enddo
      JCOLW(NVAR+1) = NZW + 1
C
C         MODIFY THE IPRMX ARRAY---FROM HERE ON IT HAS LENGTH NFREE
C
      iwork(1:ndim) = iprmx(1:ndim)
      II = 0
      DO I = 1,NDIM
        IF(IWORK(I).LE.NFREE) THEN
          II = II + 1
          IPRMX(II) = I
        ENDIF
      enddo
c
      xlwri(1:ndim) = zero
      xupri(1:ndim) = zero
C
C         CONSTRUCT THE INTERNAL QUANTITIES CORRESPONDING TO Y
C
      DO I = 1,NFREE
C
C         CONSTRUCT THE INTERNAL VARIABLES Y
C
        YVEC(I) = XVEC(IPRMX(I))
C
C         CONSTRUCT GRADIENT WITH RESPECT TO INTERNAL VARIABLES Y
C
        GVEC(I) = DELF(IPRMX(I))
C
C         CONSTRUCT REAL VARIABLE LOWER BOUND WITH INTERNAL ORDER
C
        XLWRI(I) = XLWR(IPRMX(I))
C
C         CONSTRUCT REAL VARIABLE UPPER BOUND WITH INTERNAL ORDER
C
        XUPRI(I) = XUPR(IPRMX(I))
C
      enddo
C
      VIOLMX = ZERO
C
      IF(RELAX) THEN
C
        IF(NSLK.GT.0) THEN
          GVEC(NFREE+MINEQL+1:NFREE+NSLK) = ONE
        ENDIF
C
        LCSVEC = NFREE + 1
        LCTVEC = LCSVEC + MINEQL
        LCUVEC = LCTVEC + MEQUAL
        LCVVEC = LCUVEC + MEQUAL
        LCWVEC = LCVVEC + NALWR
C
C         INITIALIZE THE EQUALITY SLACKS
C
        DO I = 1,MEQUAL
          SLAK = CVEC(I) - ALWRI(I) 
          FEASTP = MAX(FEATOL,FEATOL*ABS(SLAK))
          IF(SLAK.GT.ZERO) THEN
            YVEC(LCTVEC+I-1) = SLAK + FEASTP
            YVEC(LCUVEC+I-1) = FEASTP
          ELSE
            YVEC(LCTVEC+I-1) = FEASTP
            YVEC(LCUVEC+I-1) = FEASTP - SLAK
          ENDIF
        enddo
C
C         INITIALIZE THE INEQUALITY SLACKS
C
        KL = 0
        KU = 0
        DO I = 1,MINEQL
          II = MEQUAL + I
          YVEC(LCSVEC+I-1) = CVEC(II)
          IF(ALWRI(II).GT.-BIGBND) THEN
            IF(AUPRI(II).LT.BIGBND) THEN
              AMALWR = CVEC(II) - ALWRI(II)
              FEASTP = MAX(FEATOL,FEATOL*ABS(AMALWR))
              KL = KL + 1
              YVEC(LCVVEC+KL-1) = MAX(FEASTP,FEASTP-AMALWR)
              AUPRMA = AUPRI(II) - CVEC(II)
              FEASTP = MAX(FEATOL,FEATOL*ABS(AUPRMA))
              KU = KU + 1
              YVEC(LCWVEC+KU-1) = MAX(FEASTP,FEASTP-AUPRMA)
            ELSE
              AMALWR = CVEC(II) - ALWRI(II)
              FEASTP = MAX(FEATOL,FEATOL*ABS(AMALWR))
              KL = KL + 1
              YVEC(LCVVEC+KL-1) = MAX(FEASTP,FEASTP-AMALWR)
            ENDIF
          ELSE
            IF(AUPRI(II).LT.BIGBND) THEN
              AUPRMA = AUPRI(II) - CVEC(II)
              FEASTP = MAX(FEATOL,FEATOL*ABS(AUPRMA))
              KU = KU + 1
              YVEC(LCWVEC+KU-1) = MAX(FEASTP,FEASTP-AUPRMA)
            ENDIF
          ENDIF
        enddo
C
C         COMPLETE CONSTRUCTION OF CVEC 
C
        DO I = 1,MEQUAL
          CVEC(I) = CVEC(I) - ALWRI(I) - YVEC(LCTVEC+I-1) 
     $                                 + YVEC(LCUVEC+I-1)
          VIOLMX = MAX(VIOLMX,ABS(CVEC(I)))
        enddo
C
        DO I = 1,MINEQL
          II = MEQUAL + I
          VIOL = MAX(ZERO,ALWRI(II)-CVEC(II),CVEC(II)-AUPRI(II))
          VIOLMX = MAX(VIOLMX,VIOL)
          CVEC(II) = CVEC(II) - YVEC(LCSVEC+I-1)
        enddo
C
      ELSE
C
        LCSVEC = NFREE + 1
        LCTVEC = LCSVEC + MINEQL
        LCUVEC = LCTVEC + MEQUAL
        LCVVEC = LCUVEC + MEQUAL
        LCWVEC = LCVVEC + NALWR
C
        IF(MCON.GT.0.AND.NVAR-LCTVEC+1.GT.0) THEN
          YVEC(LCTVEC:NVAR) = ZERO
        ENDIF
C
C         INITIALIZE THE INEQUALITY SLACKS
C
        DO I = 1,MINEQL
          II = MEQUAL + I
          FEASTP = MAX(FEATOL,FEATOL*ABS(CVEC(II)))
          IF(ALWRI(II).GT.-BIGBND) THEN
            IF(AUPRI(II).LT.BIGBND) THEN
              FEASVL = MIN(FEASTP,POINT5*ABS(AUPRI(II)-ALWRI(II)))
              YVEC(LCSVEC+I-1) = MAX(ALWRI(II)+FEASVL
     $                       ,MIN(CVEC(II),AUPRI(II)-FEASVL))
            ELSE
              YVEC(LCSVEC+I-1) = MAX(ALWRI(II)+FEASTP,CVEC(II))
            ENDIF
          ELSE
            IF(AUPRI(II).LT.BIGBND) THEN
              YVEC(LCSVEC+I-1) = MIN(CVEC(II),AUPRI(II)-FEASTP)
            ENDIF
          ENDIF
        enddo
C
C         COMPLETE CONSTRUCTION OF CVEC 
C
        DO I = 1,MEQUAL
          CVEC(I) = CVEC(I) - ALWRI(I)
          VIOLMX = MAX(VIOLMX,ABS(CVEC(I)))
        enddo
C
        DO I = 1,MINEQL
          II = MEQUAL + I
          VIOL = MAX(ZERO,ALWRI(II)-CVEC(II),CVEC(II)-AUPRI(II))
          VIOLMX = MAX(VIOLMX,VIOL)
          CVEC(II) = CVEC(II) - YVEC(LCSVEC+I-1)
        enddo
C
      ENDIF
C
C         FOR FEASIBILITY MODE CONSTRUCT OBJECTIVE AND GRADIENT 
C
      IF(FZMODE) THEN
C
C         CHECK FOR CONSTRAINT SATISFACTION 
C
        IF(VIOLMX.LT.CONTOL) IFERR = -101
C
        FNLP = POINT5*DOT_PRODUCT(CVEC(1:MSUBE),CVEC(1:MSUBE))
C
C         COMPUTE GRADIENT = (CMAT**T)(CVEC) 
C
        IF(RELAX) THEN
          NVRGVC = NVAR
        ELSE
          NVRGVC = NFREE + MINEQL
        ENDIF
C
        CALL MVPSPR(11,NVRGVC,MSUBE,CMAT,IROWC,JCOLC,CVEC,GVEC)
C
      ELSE
C
        FNLP = FBAR
        IF(NSLK.GT.0) GVEC(NFREE+1:NFREE+NSLK) = ZERO
C
      ENDIF
C
C         DEFINE THE NUMBER OF BOUNDS AND CONSTRUCT IBOUND
C
      MSBB = 0
      DO I = 1,NFREE
        IF(XLWRI(I).GT.-BIGBND) THEN
          MSBB = MSBB + 1
          IBOUND(MSBB) = I
        ENDIF
      enddo
C
      DO I = 1,NFREE
        IF(XUPRI(I).LT.BIGBND) THEN
          MSBB = MSBB + 1
          IBOUND(MSBB) = I
        ENDIF
      enddo
C
      DO I = MEQUAL+1,MEQUAL+MINEQL
        IF(ALWRI(I).GT.-BIGBND) THEN
          MSBB = MSBB + 1
          IBOUND(MSBB) = I
        ENDIF
      enddo
C
      DO I = MEQUAL+1,MEQUAL+MINEQL
        IF(AUPRI(I).LT.BIGBND) THEN
          MSBB = MSBB + 1
          IBOUND(MSBB) = I
        ENDIF
      enddo
C
      MSBB = MSBB + 2*MEQUAL + NALWR + NAUPR
C
      IF(MSBB.GT.0) THEN
C
C         MAY WANT TO REDO THIS USING VECNU AND VECLAM
C
        vlamda(1:msbb) = zero
C
        XMSUBB = MSBB
C
C         CONSTRUCT THE BOUND JACOBIAN BMAT AND ASSOCIATED IROWB, JCOLB
C
        LC1 = 1
        LC2 = LC1 + NXLWR
        LC3 = LC2 + NXUPR
        LC4 = LC3 + NALWR
        LC5 = LC4 + NAUPR
        LC6 = LC5 + 2*MEQUAL + NALWR + NAUPR
        NZB = 0
        jcolb(1:nvar+1) = 0
C
        DO JJ = LC1,LC2-1
          II = IBOUND(JJ)
          BVEC(JJ) = YVEC(II) - XLWRI(II)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
          NZB = NZB + 1
          BMAT(NZB) = ONE
          IROWB(NZB) = JJ
          JCOLB(II+1) = JCOLB(II+1) + 1
          XROW = JJ
          XCOL = II
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
        enddo
C
        DO JJ = LC2,LC3-1
          II = IBOUND(JJ)
          BVEC(JJ) = XUPRI(II) - YVEC(II)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
          NZB = NZB + 1
          BMAT(NZB) = -ONE
          IROWB(NZB) = JJ
          JCOLB(II+1) = JCOLB(II+1) + 1
          XROW = JJ
          XCOL = II
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
        enddo
C
        DO JJ = LC3,LC4-1
C
          KK = JJ - LC3 
          II = IBOUND(JJ)
          KCOL1 = LCSVEC + II - MEQUAL - 1
          KCOL2 = LCVVEC + KK 
          BVEC(JJ) = YVEC(KCOL1) + YVEC(KCOL2) - ALWRI(II)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
C
        enddo
C
        DO JJ = LC4,LC5-1
C
          KK = JJ - LC4 
          II = IBOUND(JJ)
          KCOL1 = LCSVEC + II - MEQUAL - 1
          KCOL2 = LCWVEC + KK 
          BVEC(JJ) = AUPRI(II) - YVEC(KCOL1) + YVEC(KCOL2)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
C
        enddo
C
        KCOL = LCTVEC - 1
        DO JJ = LC5,LC6-1
C
          KCOL = KCOL + 1
          BVEC(JJ) = YVEC(KCOL)
          BVEC(JJ) = MAX(ZEROMN,BVEC(JJ))
C
        enddo
C
C           CONSTRUCT THE CORRESPONDING BMAT
C
        DO JJ = LC3,LC4-1
C
          KK = JJ - LC3 
          II = IBOUND(JJ)
          KCOL1 = LCSVEC + II - MEQUAL - 1
          KCOL2 = LCVVEC + KK 
C
          NZB = NZB + 1
          BMAT(NZB) = ONE
          IROWB(NZB) = JJ
          JCOLB(KCOL1+1) = JCOLB(KCOL1+1) + 1
          XROW = JJ
          XCOL = KCOL1
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
C
          NZB = NZB + 1
          BMAT(NZB) = ONE
          IROWB(NZB) = JJ
          JCOLB(KCOL2+1) = JCOLB(KCOL2+1) + 1
          XROW = JJ
          XCOL = KCOL2
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
C
        enddo
C
        DO JJ = LC4,LC5-1
C
          KK = JJ - LC4 
          II = IBOUND(JJ)
          KCOL1 = LCSVEC + II - MEQUAL - 1
          KCOL2 = LCWVEC + KK 
C
          NZB = NZB + 1
          BMAT(NZB) = -ONE
          IROWB(NZB) = JJ
          JCOLB(KCOL1+1) = JCOLB(KCOL1+1) + 1
          XROW = JJ
          XCOL = KCOL1
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
C
          NZB = NZB + 1
          BMAT(NZB) = ONE
          IROWB(NZB) = JJ
          JCOLB(KCOL2+1) = JCOLB(KCOL2+1) + 1
          XROW = JJ
          XCOL = KCOL2
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
C
        enddo
C
        KCOL = LCTVEC - 1
        DO JJ = LC5,LC6-1
C
          KCOL = KCOL + 1
C
          NZB = NZB + 1
          BMAT(NZB) = ONE
          IROWB(NZB) = JJ
          JCOLB(KCOL+1) = JCOLB(KCOL+1) + 1
          XROW = JJ
          XCOL = KCOL
          RWORK(NZB) = XROW + (XCOL-ONE)*XMSUBB
C
        enddo
C
C
C         SORT THE COLUMN INDEX ARRAY RWORK INTO ASCENDING ORDER
C         CONSTRUCTING THE PERMUTATION ARRAY FOR BMAT
C
        CALL HDSRTN(RWORK,NZB,0,0,IWORK,IERP)
C
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         APPLY THE PERMUTATION TO THE TWO BMAT ARRAYS, THUS
C         CONVERTING THEM TO COLUMN ORDER
C
        CALL HDPRMX(BMAT,NZB,IWORK,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
        CALL HJPRMX(IROWB,NZB,IWORK,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
        JCOLB(1) = 1
        DO K = 2,NVAR+1
          JCOLB(K) = JCOLB(K-1) + JCOLB(K)
        enddo
C
      ELSE
        jcolb(1:nvar+1) = 0
      ENDIF
C
C         IF NECESSARY PERMUTE THE DIAGNOSTIC FUNCTION VALUE
C
      IF(LYNFNC.NE.0) LYNFNC = IPRMX(LYNFNC)
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      RETURN
C
      END
