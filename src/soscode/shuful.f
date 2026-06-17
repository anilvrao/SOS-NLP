

      SUBROUTINE SHUFUL( MAXCON ,MCON   ,MEQUAL ,MIGNOR ,MTOTAL  
     $          ,NDIM   ,NONZG  ,NONZGT ,CBAR   ,CLWR   ,CUPR     
     $          ,CVEC   ,VECLAM ,GMAT   ,IROWG  ,JCOLG  ,ISTATC    
     $          ,IPRMC  ,IPRMG  ,IWORK  ,IERNLP ,PRMUTE)
C
C ======================================================================
C     SHUFUL===>shuful   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C *** PURPOSE
C
C       The purpose of this routine is to reorder the constraint
C       data.  On input the constraints and associated information
C       are in the user specified (arbitrary) order.  On output
C       this routine reorders the constraint and Jacobian data
C       such that the equalities are first, followed by the inequalities,
C       and then the ignored contraints.  The Jacobian matrix is
C       reordered to correspond to the constraint ordering.  If
C       a sparse representation is used the Jacobian is permuted
C       such that the first NONZGT elements form a contiguous array
C       which contains the Jacobian of the nonignored constraints.
C       This corresponds to storing the columns of the nonignored
C       Jacobian.  The ignored constraints are stored after this
C       (i.e. beginning in location NONZGT+1).  If a dense 
C       representation of the Jacobian is used, the rows are
C       interchanged to reflect the constraint ordering, and
C       the leading row dimension is stored in the first
C       element of the array JCOLG.  The permutation information
C       is retained in two arrays; iprmc defines the constraint
C       permutation; iprmg defines the corresponding Jacobian
C       permuation.
C
C *** CALLING ARGUMENTS
C
C        MAXCON       MAXIMUM NUMBER OF CONSTRAINTS
C        MCON         NUMBER OF CONSTRAINTS
C        MEQUAL       NUMBER OF EQUALITY CONSTRAINTS
C        MIGNOR       NUMBER OF IGNORED CONSTRAINTS
C        MTOTAL       NUMBER OF NONIGNORED CONSTRAINTS
C        NDIM         NUMBER OF VARIABLES
C        NONZG        NUMBER OF JACOBIAN NONZEROS
C        NONZGT       NUMBER OF JACOBIAN NONZEROS IN NONIGNORED PORTION
C        CBAR         CONSTRAINTS EVALUATED AT XBAR (MCON)
C        CLWR         CONSTRAINT LOWER BOUNDS (MCON)
C        CUPR         CONSTRAINT UPPER BOUNDS (MCON)
C        CVEC         CONSTRAINTS AT X (MCON)
C        VECLAM       LAGRANGE MULTIPLIERS FOR CONSTRAINTS (MCON)
C        GMAT         CONSTRAINT DERIVATIVES AT XBAR (NONZG)
C        IROWG        ROW INDICES OF JACOBIAN NONZEROS (NONZG)
C        JCOLG        COLUMN INDICES OF JACOBIAN NONZEROS (NONZG)
C        ISTATC       INTEGER CONSTRAINT STATUS (MCON)
C        IPRMC        INTEGER CONSTRAINT PERMUTATION ARRAY (MCON)
C        IPRMG        INTEGER JACOBIAN PERMUTATION ARRAY (NONZG)
C        IWORK        INTEGER WORK ARRAY (NONZG)
C        IERNLP       SUCCESS/ERROR CODE
C        PRMUTE       LOGICAL FLAG--TRUE FOR CONSTRAINT/JACOBIAN
C                     PERMUTATION FROM INTERNAL TO EXTERNAL ORDER
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      DIMENSION  CBAR(MAXCON)  ,CLWR(MAXCON)    ,CUPR(MAXCON)
     &          ,CVEC(MAXCON)  ,VECLAM(MAXCON)  ,GMAT(NONZG) 
     &          ,IROWG(NONZG)  ,JCOLG(*)        ,ISTATC(MAXCON) 
     &          ,IPRMC(MAXCON) ,IPRMG(NONZG)    ,IWORK(NONZG)
      LOGICAL PRMUTE
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C
      LNJCOL = MAX(NONZG,NDIM+1)
      IF(MCON.GT.0) THEN
C
C         CONSTRUCT PERMUTATION TO PUT EQUALITY CONSTRAINTS AS THE 
C         FIRST ROWS OF THE JACOBIAN MATRIX 
C         NOTE THAT IGNORED CONSTRAINTS ARE PUT AT THE BOTTOM OF 
C         THE STACK.  LOWER LEVEL ROUTINES (NLPSPR AND BELOW) ONLY 
C         SEE THE FIRST MTOTAL CONSTRAINTS. FOR NLPSPR AND BELOW 
C         MTOTAL BECOMES MCON.
C         SIMILARLY THE CORRESPONDING NUMBER OF NONZEROS NONZGT
C         CONSTRUCTED BELOW, BECOMES NONZG FOR LOWER LEVEL ROUTINES
C
C
        MEQUAL = 0
        MIGNOR = 0
        DO I=1,MCON
          IF (ISTATC(I).EQ.3) THEN
            MEQUAL = MEQUAL + 1
          ELSEIF (ISTATC(I).EQ.4) THEN
            MIGNOR = MIGNOR + 1
          ENDIF
        ENDDO
        MTOTAL = MCON - MIGNOR
C
        K1 = 0
        K2 = 0
        K3 = 0
        DO I = 1,MCON
C
          IF(ISTATC(I).EQ.3) THEN
C
C         PUT EQUALITY CONSTRAINTS INTO ACTIVE SET
C
            K1 = K1 + 1
            IPRMC(K1) = I
C 
          ELSEIF(ISTATC(I).LT.3) THEN
C
C         PUT INEQUALITY CONSTRAINTS INTO SET
C
            K2 = K2 + 1
            IPRMC(MEQUAL+K2) = I
C
          ELSE
            K3 = K3 + 1
            IPRMC(MTOTAL+K3) = I
          ENDIF
        ENDDO
C
        IF(QPOPTN.EQ.'SPARSE') THEN
C
C         INITIALIZE MASTER JACOBIAN PERMUTATION RECORD
C
          DO I = 1,NONZG
            IPRMG(I) = I
          ENDDO
C
C         CONSTRUCT ORDERING FOR INITIAL ACTIVE SET
C         ON OUTPUT THE JACOBIAN VECTOR TRIPLE (GMAT,IROW,JSTRG)
C         IS STORED IN ROW ORDER WITH EQUALITIES IN THE FIRST
C         MEQUAL ROWS, FOLLOWED BY INEQUALITIES AND THEN IGNORED
C         CONSTRAINTS
C
          CALL PRMMAT(1,NONZG,GMAT,IROWG,JCOLG,
     $      IPRMG,IPRMC,MCON,IWORK,IERP)
          IF(IERP.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
C
C         COMPUTE THE NUMBER OF NONZEROS IN THE IGNORED PORTION
C         OF THE JACOBIAN, AND THEN THE NUMBER OF NONZEROS
C         CORRESPONDING TO THE MTOTAL CONSTRAINTS
C
          NONZGI = 0
          DO I=1,NONZG
            IF (IROWG(I).GT.MTOTAL)  NONZGI = NONZGI + 1
          ENDDO
          NONZGT = NONZG - NONZGI
C
        ELSE
C
C         CONSTRUCT THE PERMUTATION ARRAY IPRMC, BUT NOT IPRMG
C
          CALL HJSRTN(IPRMC,MCON,1,0,IWORK,IERP)
          IF(IERP.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
C
C         COMPUTE THE NUMBER OF ELEMENTS IN THE NONIGNORED JACOBIAN
C
          NONZGT = MTOTAL*NDIM
C
        ENDIF
C
C         APPLY ORDERING TO CBAR,CVEC,VECLAM AND ISTATC
C
        CALL HDPRMX(CBAR,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
        CALL HJPRMX(ISTATC,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
        CVEC(1:MCON) = CBAR(1:MCON)
        CALL HDPRMX(VECLAM,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         APPLY ORDERING TO UPPER AND LOWER BOUNDS
C
        CALL HDPRMX(CUPR,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
        CALL HDPRMX(CLWR,MCON,IPRMC,IERP)
        IF(IERP.NE.0) THEN
          IERNLP = -147
          RETURN
        ENDIF
C
C         IF NECESSARY PERMUTE THE DIAGNOSTIC FUNCTION VALUE
C
        IF(LYNFNC.NE.0) LYNFNC = IPRMC(LYNFNC)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
C
C         REORDER THE FIRST NONZGT ELEMENTS OF THE JACOBIAN BY COLUMNS
C
          IOPT = 0
          ISTBLE = 0
          IF(MTOTAL.GT.0) CALL HJSRTN(JCOLG,NONZGT,IOPT,ISTBLE,
     $     IWORK,IERP)
          IF(IERP.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
C
C         USE THE IDENTITY PERMUTATION FOR THE LAST NONZGI ELEMENTS
C
          DO I = NONZGT+1,NONZG
            IWORK(I) = I
          ENDDO
C
C         APPLY THE SORT TO THE OTHER ARRAYS
C
          CALL HJPRMX(IROWG,NONZG,IWORK,IERP)
          IF(IERP.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
          CALL HJPRMX(IPRMG,NONZG,IWORK,IERP)
          IF(IERP.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
          CALL HDPRMX(GMAT,NONZG,IWORK,IERP)
          IF(IERP.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
C
C         CONVERT THE FIRST NONZGT ELEMENTS OF THE JACOBIAN TO 
C         COLUMN FORMAT
C
          IF(MTOTAL.GT.0) THEN
            CALL STORCH(1,NONZGT,NDIM,JCOLG,IWORK,LNJCOL,IERCH)
          ELSE
            JCOLG(1:LNJCOL) = 1
          ENDIF
          IF(IERCH.NE.0) THEN
            IERNLP = -147
            RETURN
          ENDIF
        ELSE
C
          IROW = 1 - MAXCON
          DO JCOL = 1,NDIM
            IROW = IROW + MAXCON
C
C           PERMUTE COLUMN JCOL INTO THE CORRECT ORDER
C
            CALL HDPRMX(GMAT(IROW),MCON,IPRMC,IERP)
            IF(IERP.NE.0) THEN
              IERNLP = -147
              RETURN
            ENDIF
C
          ENDDO
C
C           SET JCOLG(1) TO MAXCON 
C
          JCOLG(1) = MAXCON
C
        ENDIF
C
C         MAKE SURE NONZGT IS NONZERO TO AVOID DIMENSION PROBLEMS IN
C         LOWER LEVEL ROUTINES
C
        NONZGT = MAX(1,NONZGT)
C
      ELSE
C
        NONZGT = 1
        MIGNOR = 0
        MTOTAL = 0
        JCOLG(1:LNJCOL) = 1
C
      ENDIF
C
      PRMUTE = .FALSE.
      IF(MCON.EQ.0) GO TO 170
      DO I = 1,MAXCON
        IF(IPRMC(I).NE.I) THEN
          PRMUTE = .TRUE.
          GO TO 170
        ENDIF
      ENDDO
      IF(QPOPTN.EQ.'SPARSE') THEN
        DO I = 1,NONZGT
          IF(IPRMG(I).NE.I) THEN
            PRMUTE = .TRUE.
            GO TO 170
          ENDIF
        ENDDO
      ENDIF
 170  CONTINUE
C
      RETURN
C
      END
