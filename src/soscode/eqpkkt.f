
      SUBROUTINE EQPKKT(NCON,NDIM,GMAT,IROW,JCOLST,NZGDIM,BUPR,BLWR,
     $    ISTATC,IACTIV,HMAT,IROWH,JSTRH,NZHDIM,GVEC,QUAD,PUPR,PLWR,
     $    ISTATV,PVEC,VARMLT,CONMLT,RHSKKT,NRDIM,SOLNKT,RWORK,LNRWRK,
     $    NEEDED,IFREEV,KTROW,NCALL,IPU,IPC,CNDNUM,LNSYMB,IEREKT)
C
C ======================================================================
C     EQPKKT===>eqpkkt   J.T. BETTS 
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C         PURPOSE:  SOLVE THE SET OF EQUATIONS RESULTING FROM 
C                   LINEARIZATION OF THE KUHN-TUCKER CONDITIONS.
C                   THE MULTIFRONTAL SPARSE SYMMETRIC INDEFINITE
C                   SOLVER IS APPLIED TO A SYSTEM CONSTRUCTED 
C                   FROM COLUMNS (AND ROWS) OF THE SYMMETRIC SYSTEM
C
C                       |  H        || -p | = | g  |
C                       |  Ga    0  ||  v | = | c  |
C
C                   WHERE THE SUBSCRIPT "a" REFERS TO THE ACTIVE
C                   CONSTRAINTS.
C                   THIS SYSTEM IS DERIVED FROM THE GENERAL QP:
C
C                   MINIMIZE
C
C                     q = .5(p^T)Hp + (g^T)p
C
C                   SUBJECT TO THE CONSTRAINTS
C
C                     bl < Gp < bu
C
C                   AND BOUNDS
C
C                     pl < p < pu
C
C                   THE FIXED VARIABLES ARE DEFINED BY THE INPUT ISTATV
C                   AND THE ACTIVE CONSTRAINTS ARE DEFINED BY THE INPUT ISTATC.
C
C         INPUT:    
C
C             NCON   NUMBER OF CONSTRAINTS AND NUMBER OF ROWS IN GMAT
C             NDIM   NUMBER OF VARIABLES AND NUMBER OF COLUMNS IN GMAT
C             GMAT   NONZERO ELEMENTS OF LINEAR SYSTEM (NZGDIM)
C             IROW   INTEGER ROW INDEX VECTOR (NZGDIM)
C             JCOLST INTEGER COLUMN START VECTOR (NDIM+1)
C             NZGDIM  DIMENSION OF GMAT (GE. JCOLST(NDIM+1) - 1)
C             BUPR   CONSTRAINT UPPER BOUND VECTOR (NCON)
C             BLWR   CONSTRAINT LOWER BOUND VECTOR (NCON)
C             ISTATC INTEGER CONSTRAINT STATUS (NCON)
C                    = 0  --- FREE (INACTIVE) INEQUALITY
C                    = 1  --- FIXED ON LOWER BOUND
C                    = 2  --- FIXED ON UPPER BOUND
C                    = 3  --- EQUALITY
C                    = 4  --- IGNORED CONSTRAINT
C             IACTIV INTEGER ACTIVE CONSTRAINT ROW INDEX (NCON)
C             HMAT   NONZERO ELEMENTS OF LOWER TRIANGLE OF HMAT (NZHDIM)
C             IROWH  INTEGER ROW INDEX VECTOR (NZHDIM)
C             JSTRH  INTEGER COLUMN START VECTOR (NDIM+1)
C             NZHDIM DIMENSION OF HMAT (.GE. JSTRH(NDIM+1)-1 )
C             GVEC   RIGHT HAND SIDE (GRADIENT) VECTOR (NDIM)
C             QUAD   OPTIMAL OBJECTIVE FUNCTION VALUE
C             PUPR   UPPER BOUND FOR INDEPENDENT VARIABLES (NDIM)
C             PLWR   LOWER BOUND FOR INDEPENDENT VARIABLES (NDIM)
C             ISTATV INTEGER VARIABLE STATUS (NDIM)
C                    = 0  --- FREE VARIABLE 
C                    = 1  --- FIXED ON LOWER BOUND
C                    = 2  --- FIXED ON UPPER BOUND
C                    = 3  --- FIXED PERMANENTLY 
C                    = 4  --- IGNORED BOUND
C
C             LNRWRK LENGTH OF RWORK ARRAY (SEE BELOW)
C
C                    LOWER BOUNDS FOR THE WORK ARRAYS
C                    REAL:
C                      LNRWRK > NKT + 4*NEQN + 200
C                    WHERE
C                      NKT  > MINEQL + 1 + MCON + NZHDIM + NZGDIM
C                      NEQN > MINEQL + 1 + NDIM + MCON
C
C             NEEDED STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C             KTROW  INTEGER WORK ARRAY (2*NDIM)
C             NCALL  NUMBER OF THE CALL 
C                    = 1 FIRST CALL,  DO THE FOLLOWING STEPS:
C                        (1) INPUT MATRIX STRUCTURE
C                        (2) ORDER THE MATRIX
C                        (3) SYMBOLIC FACTORIZATION
C                        (4) INPUT MATRIX VALUES
C                        (5) NUMERICAL FACTORIZATION
C                        (6) SOLVE
C                    = 2 SECOND CALL  DO STEPS 4-6
C                    = 3 THIRD CALL   DO STEP  6
C             IPU    OUTPUT UNIT NO.
C             IPC    OUTPUT CONTROL FLAG
C
C         OUTPUT:
C
C             PVEC   SEARCH DIRECTION PART OF SOLUTION TO THE 
C                    LINEAR SYSTEM (NDIM)
C             VARMLT VARIABLE MULTIPLIERS (NDIM)
C             CONMLT CONSTRAINT MULTIPLIERS (NCON)
C             RHSKKT RIGHT HAND SIDE OF THE KT SYSTEM AND SCRATCH SPACE (NRDIM)
C             NRDIM  DIMENSION OF RHSKKT; NRDIM = MAX(2*NDIM,NCON))
C             SOLNKT SOLUTION OF THE KT SYSTEM (2*NDIM)
C             RWORK  WORK ARRAY CONTAINING MULTIFRONTAL FACTORIZATION
C                    INFORMATION FOR SUBSEQUENT USE STORED AS FOLLOWS
C                    RWORK(1) -- SYMBOLIC FACTORIZATION DATA (LNSYMB)
C                    RWORK(LNSYMB+1) -- NUMERIC FACTORIZATION DATA (LNRWRK)
C             IFREEV INTEGER WORK ARRAY (NDIM)
C             CNDNUM CONDITION NUMBER OF MATRIX
C             IEREKT ERROR RETURN CODE
C
C                    0    NORMAL RETURN
C
C             ERROR RETURNS HAVE BEEN MAPPED TO CORREPOND TO SHUR-QP EQUIVALENTS
C
C                   -1     CURRENT INPUT ACTIVE SET IS NOT CORRECT
C                   -1014  REAL WORK ARRAY TOO SMALL 
C                   -1019  INTEGER WORK ARRAY TOO SMALL
C                   -1111  I/O ERROR (INSUFFICIENT DISK SPACE)
C                   -1103  SINGULAR KT
C                   -1104  WRONG INERTIA
C                   -1108  EXCESSIVE FILL
C
C             THE CONVENTION FOR ERROR RETURN NUMBERS IS AS FOLLOWS:
C             -600 > IEREKT > -700       I/O ERRORS (INSUFFICIENT DISK SPACE)
C             -700 > IEREKT > -800       ERROR REQUIRES SCHUR-QP RESPONSE
C             -800 > IEREKT > -900   --- INPUT ERRORS
C             -900 > IEREKT          --- STORAGE ERRORS REQUIRING
C                                        SCHUR-QP RESPONSE
C             THE UNITS DIGIT IS NUMBERED CHRONOLOGICALLY WITHIN EQPKKT
C
C                    -801  NCALL.LT.1 OR NCALL.GT.3
C                    -802  NCON.LT.0
C                    -803  NDIM.LE.0
C                    -1014  REAL WORK ARRAY TOO SMALL
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLIN, IER=-101)
C                    -808  NEQNS .LT. 1 (XDSLIN, IER=-102)
C                    -809  INCORRECT PATH (XDSLIC, IER=-100)
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLIC, IER=-101)
C                    -811  ILLEGAL VALUE FOR NZCOL (XDSLIC, IER=-103)
C                    -812  ILLEGAL VALUE FOR JCOL (XDSLIC, IER=-104)
C                    -813  ILLEGAL VALUE FOR JROWIN (XDSLIC, IER=-106)
C                    -1111  I/O ERROR ON WAFIL1 (XDSLIC, IER=-110)
C                    -814  INCORRECT PATH (XDSLIF, IER=-100)
C                    -1014  REAL WORK ARRAY TOO SMALL (NEEDS.GT.LNRWRK)
C                    -1111  I/O ERROR ON WAFIL1 (XDSLIF, IER=-110)
C                    -818  INCORRECT PATH (XDSLOR, IER=-200)
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLOR, IER=-201)
C                    -1111  I/O ERROR ON SQFILE (XDSLOR, IER=-202)
C                    -822  INCORRECT PATH (XDSLSF, IER=-300)
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLSF, IER=-301)
C                    -925  INTERNAL ERROR DURING SYMBOLIC FACT. (XDSLSF, IER=-302)
C                    -1111  I/O ERROR ON SQFIL2 (XDSLSF, IER=-303)
C                    -1014  REAL WORK ARRAY TOO SMALL (NEEDST.GT.LNRWRK)
C                    -829  INCORRECT PATH (XDSLVC, IER=-400)
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLVC, IER=-401)
C                    -831  JROWIN(I),JCOL PAIR NOT SPECIFIED (XDSLVC, IER=-402)
C                    -1111  I/O ERROR ON WAFILE (XDSLVC, IER=-410)
C                    -833  INCORRECT PATH (XDSLVF, IER=-400)
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLVF, IER=-401)
C                    -1111  I/O ERROR ON SQFILE, WAFILE (XDSLVF, IER=-410)
C                    -836  INCORRECT PATH (XDSLFA, IER=-500,-502)
C                    -1014  REAL WORK ARRAY TOO SMALL (XDSLFA, IER=-501,-507,-508)
C                    -1103  ZERO PIVOT; SINGULAR MATRIX (XDSLFA, IER=-503)
C                    -1111  I/O ERROR ON SQFIL1 OR SQFIL2 (XDSLFA, IER=-504)
C                    -1111  I/O ERROR ON WAFIL1 (XDSLFA, IER=-505,-511,-512)
C                    -1111  I/O ERROR ON WAFIL2 (XDSLFA, IER=-506)
C                    -1108  EXCESSIVE FILL DUE TO PIVOTING (XDSLFA, IER=-509)
C                    -1104  WRONG INERTIA OF KT MATRIX
C                    -845  INCORRECT PATH (XDSLSL, IER=-600)
C                    -1014  REAL WORK ARRAY TOO SMALL FOR SOLUTION (XDSLSL, IER=-601)
C                    -847  NRHS .LT. 1 (XDSLSL, IER=-602)
C                    -848  LDRHS .LE. NEQNS (XDSLSL, IER=-603)
C                    -1111  I/O ERROR ON WAFIL1 (XDSLSL, IER=-604,605)
C                    -750  MAX. ITERATION DURING ITERATIVE REFINEMENT (REFITR, IER=-701)
C                    -852  INCORRECT INPUT WHEN REQUESTING KT FILE OUTPUT
C                    -732  PIVOTING FAILURE -- NEED TO INCREASE PANEL SIZE 
C                          (XDSLFA, IER=-513)
C
C             LNSYMB LENGTH OF RWORK ARRAY NEEDED FOR SYMBOLIC FACTORIZATION
C             NEEDED STORAGE REQUIRED WHEN LNRWRK IS TOO SMALL
C
      INCLUDE '../commons/NLPSPR.CMN'
C
      COMMON /STATIS/ INSTAT(30),RLSTAT(20)
      COMMON /INERVL/ INREQD
      COMMON /KONSTN/ 
     *  ZEROMN  ,ZEROOT  ,BIGNUM  ,BGROOT  ,BIGBND  ,BIGCND
      DIMENSION GMAT(NZGDIM),IROW(NZGDIM),JCOLST(NDIM+1),HMAT(NZHDIM),
     $    IROWH(NZHDIM),JSTRH(NDIM+1),GVEC(NDIM),PVEC(NDIM),
     $    VARMLT(NDIM),CONMLT(*),RHSKKT(nrdim),SOLNKT(2*NDIM),
     $    RWORK(LNRWRK),IFREEV(NDIM),KTROW(2*NDIM)
      DIMENSION ISTATC(*),ISTATV(NDIM),IACTIV(*)
      DIMENSION BUPR(*),BLWR(*),PUPR(NDIM),PLWR(NDIM)
C
      DIMENSION INRTIA(3),TYME(6),OPCNTS(2)
C    
      LOGICAL OPUN,DUMPKT
      LOGICAL BADACT
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,ONEEM3=1.0D-3,TWO=2.D0)
      PARAMETER (NITREF=0)
      INTEGER XDSLNI
C
C     ===================================================================
C
C         DROP TOLERANCE...MATRIX NONZERO ELEMENT IF |A(i,j)| > drptol
C
      DRPTOL = ZERO
C
C         BEGIN K-T FACTORIZATION TIMING
C
      CALL CLKBEG(6)
      CALL CLKBEG(13)
C
C         CHECK CALL NUMBER
C
      IEREKT = 0
C
      IF(NCALL.LT.1.OR.NCALL.GT.3) THEN
        IEREKT = -801
        GO TO 280
      ENDIF
C
C         CONSTRUCT VARIABLE DIMENSIONS FROM ISTATV
C
C         NFREE:    TOTAL NUMBER OF free VARIABLES
C         CONSTRUCT free VARIABLE COLUMN NUMBER IN IFREEV
C
      PMAX = ZERO
      NFREE = 0
      DO I = 1,NDIM
        IF(ISTATV(I).EQ.0) THEN
          NFREE = NFREE + 1
          IFREEV(I) = NFREE
        ELSE
          IFREEV(I) = 0
          IF(ISTATV(I).EQ.1) THEN
            ABSLWR = ABS(PLWR(I))
            PMAX = MAX(PMAX,ABSLWR)
          ELSE
            ABSUPR = ABS(PUPR(I))
            PMAX = MAX(PMAX,ABSUPR)
          ENDIF
        ENDIF
      enddo
C
C         CHECK FIXED VARIABLE BOUNDS TO SEE THEY ARE ZERO
C
      IF(PMAX.GT.ZEROMN) THEN
        IEREKT = -1
        GO TO 280
      ENDIF
C
C         CONSTRUCT CONSTRAINT DIMENSIONS FROM ISTATC
C
C         MACTIV:   TOTAL NUMBER OF ACTIVE CONSTRAINTS
C         CONSTRUCT ACTIVE CONSTRAINT ROW NUMBER IN IACTIV
C
      MACTIV = 0
      DO I = 1,NCON
        IF(ISTATC(I).GE.1.AND.ISTATC(I).LE.3) THEN
          MACTIV = MACTIV + 1
          IACTIV(I) = MACTIV
        ELSE
          IACTIV(I) = 0
        ENDIF
      enddo
C
C     NCROW  <-- TOTAL NUMBER OF ROWS IN CONSTRAINT PORTION OF KKT
C
      NCROW = MACTIV 
C
C     NEQNS  <-- TOTAL NUMBER OF EQUATIONS IN KKT SYSTEM
C
      NEQNS = NFREE + NCROW
C
      IF(NEQNS.GT.2*NFREE) THEN
C
C         OVERDETERMINED SYSTEM
C
        IEREKT = -804
        GO TO 280
      ENDIF
C
C     (over) ESTIMATE FOR NUMBER OF NONZEROS IN LOWER TRIANGULAR PART OF
C     KT MATRIX
C      
      NONZG = JCOLST(NDIM+1)-1
      XNDIM = NDIM
      XNCON = NCON
      XNONZG = NONZG
      XNDNSA = XNDIM*XNCON
      IF(NCON.GT.0.AND.
     $  (NONZG.LE.0.OR.XNONZG.GT.XNDNSA.OR.NONZG.GT.NZGDIM)) THEN
        IEREKT = -1008
        GO TO 280
      ENDIF
C
      XNDIM = NDIM
      NONZH = JSTRH(NDIM+1)-1
      XNONZH = NONZH
      XNDNSH = XNDIM*(XNDIM+ONE)/TWO
      IF(NONZH.LE.0.OR.XNONZH.GT.XNDNSH.OR.NONZH.GT.NZHDIM) THEN
        IEREKT = -1005
        GO TO 280
      ENDIF
C
      NZLTA = NONZH + NONZG
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
      NEEDED = 0
C
      IF(NWRKCK.GT.LNRWRK) THEN
        NEEDED = NWRKCK
        IEREKT = -1014
        GO TO 280
      ENDIF
C
C         SET ITERATIVE REFINEMENT SOLVE TOLERANCE
C
      SOLVTL = ONEEM3*ZEROOT
C
C         TO DUMP THE KKT SYSTEM AND RIGHT HAND SIDE, THE FLAG IOFMFR
C         MUST BE SET NEGATIVE.   TO DUMP ON A SPECIFIC ITERATION THE
C         VALUE OF IOFMFR CAN BE SET BASED ON THE VALUES OF IREVRS.
C         SEE SOCXOC_OPTCST FOR AN EXAMPLE.
C
      DUMPKT = IOFMFR.LT.0.AND.IOFMFR.GT.-99
C
      IF (NCALL.EQ.2) THEN
        GO TO 180
      ELSEIF (NCALL.EQ.3) THEN
        GO TO 220
      ENDIF
C
 130  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 1, 2, AND 3-----------------------------------------
C ------------------------------------------------------------------
C
C         INPUT THE MATRIX STRUCTURE TO THE MULTIFRONTAL CODE
C
C
C         -- SET IN-CORE OPTION FLAG: O (IN-CORE); 1 (OUT-OF-CORE)
C
      MFRSTR = 0
C
      CALL XDSLIN ( NEQNS, 'SYMMETRIC', IOFMFR, IPU, IPUMF1, IPUMF2, 
     $              IPUMF3, IPUMF4, 0, IPUMF5, RWORK, LNRWRK, IER )
C
      IF(IER.EQ.-101) THEN
        IEREKT = -1014
      ELSEIF(IER.EQ.-102) THEN
        IEREKT = -808
      ENDIF
      IF(IEREKT.NE.0) GO TO 280
C
C          SET PANEL SIZE TO NEQNS UNTIL A BETTER SCHEME IS FOUND
C
      CALL XDSLSP ( 'panel size', NEQNS, DUMMY, RWORK, IER )
      IF(IER.NE.0) GO TO 280
C
      CALL XDSLSP ( 'pivot tolerance', 0, TOLPVT, RWORK, IER ) 
      IF ( IER .NE. 0 ) GO TO 280
C
      CALL XDSLSP ( 'limit fill', 0, TOLFIL, RWORK, IER )
      IF ( IER .NE. 0 ) GO TO 280
C
      CALL XDSLSP ( 'save original matrix', 0, ZERO, RWORK, IER )
      IF ( IER .NE. 0 ) GO TO 280
C
C         -- LOAD COLUMNS 
C
      JCOL = 0
C
      DO ICOL = 1,NDIM
C
C         LOAD REAL COLUMN ICOL INTO KT COLUMN JCOL
C
        IF(IFREEV(ICOL).GT.0) THEN
C
          JCOL = JCOL + 1
C
C         INITIALIZE THE NUMBER OF NONZERO ELEMENTS IN COLUMN JCOL
C
          NZCOL = 0
C
C         ------ HESSIAN BLOCK COLUMNS
C
C         LOAD COLUMN ICOL OF HMAT 
C
          DO I =  JSTRH(ICOL),JSTRH(ICOL+1)-1
C
            IRH = IROWH(I)
            IRWHF = IFREEV(IRH)
C
C           LOAD ROW IRH IF IT IS A FREE VARIABLE
C
            IF(IRWHF.GT.0.AND.ABS(HMAT(I)).GT.DRPTOL) THEN
C
              NZCOL = NZCOL + 1
C                
              KTROW(NZCOL) = IRWHF
C
            ENDIF
C
          enddo
C
C         ------ JACOBIAN BLOCK COLUMNS
C
C         LOAD COLUMN JCOL OF GMAT
C
          IF(NCON.GT.0) THEN
C
            DO I =  JCOLST(ICOL),JCOLST(ICOL+1)-1
C
              IR = IROW(I)
C
              IF(IACTIV(IR).NE.0.AND.ABS(GMAT(I)).GT.DRPTOL) THEN
C
                NZCOL = NZCOL + 1
C                
                KTROW(NZCOL) = IACTIV(IR) + NFREE
C
              ENDIF
C
            enddo
C
          ENDIF
C
          CALL XDSLIC ( 'A', JCOL, NZCOL, KTROW, RWORK, LNRWRK, IER )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
          select case(ier)
          case(-100)
            IEREKT = -809
          case(-101)
            IEREKT = -1014
          case(-103)
            IEREKT = -811
          case(-104)
            IEREKT = -812
          case(-106)
            IEREKT = -813
          case(-110)
            IEREKT = -1111
          end select 
          IF(IEREKT.NE.0) GO TO 280
C
        ENDIF
C
      enddo
C
C         -- LAST NCROW COLUMNS
C
      DO JCOL = NFREE+1,NFREE+NCROW
C
        NZCOL = 0
C
C         DUMMY ROW INDEX FOR NONZERO IN COLUMN JCOL 
C         (ENTER THE ROW INDEX FOR THE DIAGONAL OF THE KKT MATRIX.)
C
        KTROW(1) = JCOL
C
        CALL XDSLIC ( 'A', JCOL, NZCOL, KTROW, RWORK, LNRWRK, IER )
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
        select case(ier)
        case(-100)
          IEREKT = -809
        case(-101)
          IEREKT = -1014
        case(-103)
          IEREKT = -811
        case(-104)
          IEREKT = -812
        case(-106)
          IEREKT = -813
        case(-110)
          IEREKT = -1111
        end select 
        IF(IEREKT.NE.0) GO TO 280
C
      enddo
C
      CALL XDSLIF ( RWORK, LNRWRK, NEEDS, IER )
C
      IF(IER.EQ.-100) THEN
        IEREKT = -814
      ELSEIF(IER.EQ.-101) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
          NEEDED = NEEDS
          IEREKT = -1014
        ENDIF
      ELSEIF(IER.EQ.-110) THEN
        IEREKT = -1111
      ENDIF
      IF ( IEREKT .NE. 0 ) GO TO 280
C
C         ORDER THE MATRIX
C
      CALL XDSLOR ( RWORK, LNRWRK, NEEDS, IER )
C
      IF(IER.EQ.-200) THEN
        IEREKT = -818
      ELSEIF(IER.EQ.-201) THEN
        IEREKT = -1014
      ELSEIF(IER.EQ.-202) THEN
        IEREKT = -1111
      ENDIF
      IF(IEREKT.NE.0) GO TO 280
C
C         PERFORM SYMBOLIC FACTORIZATION
C
      CALL XDSLSF ( RWORK ,LNRWRK, NEEDS, NEEDMN, IER )
C
      IF(IER.EQ.-300) THEN
        IEREKT = -822
      ELSEIF(IER.EQ.-301) THEN
        IF(NEEDS.GT.LNRWRK) THEN
          IEREKT = -1014
          IF(IPC.GT.0) WRITE(IPU,1001) NEEDS,LNRWRK
          NEEDED = NEEDS
        ENDIF
      ELSEIF(IER.EQ.-302) THEN
        IEREKT = -925
      ELSEIF(IER.EQ.-303) THEN
        IEREKT = -1111
      ENDIF
      IF(IEREKT.NE.0) GO TO 280
C
C         CHECK DIMENSION OF REAL WORK ARRAY
C
      IF(MFRSTR.EQ.0) THEN
        NEEDST = NEEDS
      ELSE
        NEEDST = NEEDMN
      ENDIF
C
      NEEDST = INT((ONE + TOLPVT)*DBLE(NEEDST)) 
      IF(NEEDST.GT.LNRWRK) THEN
        IF(IPC.GT.0) WRITE(IPU,1001) NEEDST,LNRWRK
        NEEDED = NEEDST
        IEREKT = -1014
        GO TO 280
      ENDIF
C
C         GET MULTIFRONTAL STATISTICS -- SAVE LENGTH NEEDED
C         FOR SYMBOLIC FACTORIZATION
C
      CALL XDSLSR(RWORK,LNSYMB,MXUSED,TYME,OPCNTS)
C
 180  CONTINUE
C
C ------------------------------------------------------------------
C --------STEPS 4 AND 5---------------------------------------------
C ------------------------------------------------------------------
C
C         INPUT MATRIX VALUES TO MULTIFRONTAL SOLVER
C
C         -- LOAD COLUMNS 
C
      JCOL = 0
C
      DO ICOL = 1,NDIM
C
C         LOAD REAL COLUMN ICOL INTO KT COLUMN JCOL
C
        IF(IFREEV(ICOL).GT.0) THEN
C
          JCOL = JCOL + 1
C
C         INITIALIZE THE NUMBER OF NONZERO ELEMENTS IN COLUMN JCOL
C
          NZCOL = 0
C
C         ------ HESSIAN BLOCK COLUMNS
C
C         LOAD COLUMN ICOL OF HMAT
C
          DO I =  JSTRH(ICOL),JSTRH(ICOL+1)-1
C
            IRH = IROWH(I)
            IRWHF = IFREEV(IRH)
C
C           LOAD ROW IRH IF IT IS A FREE VARIABLE
C
            IF(IRWHF.GT.0.AND.ABS(HMAT(I)).GT.DRPTOL) THEN
C
              NZCOL = NZCOL + 1
C              
              KTROW(NZCOL) = IRWHF
              RHSKKT(NZCOL) = HMAT(I)
C
            ENDIF
C
          enddo
C
C         ------ JACOBIAN BLOCK COLUMNS
C
C         LOAD COLUMN ICOL OF GMAT
C
          IF(NCON.GT.0) THEN
C
            DO I =  JCOLST(ICOL),JCOLST(ICOL+1)-1
C
              IR = IROW(I)
C
              IF(IACTIV(IR).NE.0.AND.ABS(GMAT(I)).GT.DRPTOL) THEN
C
                NZCOL = NZCOL + 1
C              
                KTROW(NZCOL) = IACTIV(IR) + NFREE
                RHSKKT(NZCOL) = GMAT(I)
C
              ENDIF
C
            enddo
C
          ENDIF
C
          CALL XDSLVC ( 'A', JCOL, NZCOL, KTROW, RHSKKT, 
     $                RWORK, LNRWRK, IER )
C
          IF(IER.EQ.-400) THEN
            IEREKT = -829
          ELSEIF(IER.EQ.-401) THEN
            IEREKT = -1014
          ELSEIF(IER.EQ.-402) THEN
            IEREKT = -831
          ELSEIF(IER.EQ.-410) THEN
            IEREKT = -1111
          ENDIF
          IF(IEREKT.NE.0) GO TO 280
C
        ENDIF
C
      enddo
C
C         -- LAST NCROW COLUMNS ARE ZERO BLOCK; NO INPUT REQUIRED
C
      CALL XDSLVF ( RWORK ,LNRWRK, IER )
C
      IF(IER.EQ.-400) THEN
        IEREKT = -833
      ELSEIF(IER.EQ.-401) THEN
        IEREKT = -1014
      ELSEIF(IER.EQ.-410) THEN
        IEREKT = -1111
      ENDIF
      IF ( IEREKT .NE. 0 ) GO TO 280
C
C         FACTOR THE MATRIX
C
      CALL XDSLFA ( RWORK, LNRWRK, CNDNUM, INRTIA, NEEDS, IER )
C
      select case(ier)
      case(-500,-502)
        IEREKT = -836
      case(-501,-507,-508)
        IEREKT = -1014
        NEEDED = NEEDS
      case(-503)
        IEREKT = -1103
      case(-513)
        IEREKT = -732
      case(-504,-505,-506,-511,-512)
        IEREKT = -1111
      case(-509)
        IEREKT = -1108
      end select
C
C         CHECK TO SEE IF OUT-OF-CORE MODE WAS USED
C
      IF(RWORK(55).NE.ZERO) INSTAT(28) = 1
C
      IF(IEREKT.NE.0) GO TO 280
C
C         CHECK FOR ILL-CONDITIONED KT MATRIX
C
      IF(IPC.GT.0) WRITE(IPU,*) 'CNDNUM =',CNDNUM,
     $    '    INERTIA =',(INRTIA(I),I=1,3)
C
C         CHECK INERTIA OF KT MATRIX
C
      IF(INRTIA(1).NE.(NEQNS-INREQD-NCROW)
     $    .OR.INRTIA(2).NE.(NCROW+INREQD)
     $    .OR.INRTIA(3).NE.0) THEN
        IEREKT = -1104
      ENDIF
C
      IF(IOFMFR.GT.0) CALL XDSLPS(RWORK)
C
C         GET MULTIFRONTAL STATISTICS
C
      CALL XDSLSR(RWORK,INUSE,MXUSED,TYME,OPCNTS)
      INSTAT(10) = MAX(INSTAT(10),INUSE)
C
      IF(IEREKT.NE.0) GO TO 280
C
 220  CONTINUE
C
C ------------------------------------------------------------------
C --------STEP 6----------------------------------------------------
C ------------------------------------------------------------------
C
C         DEFINE RIGHT HAND SIDE
C
      IKKT = 0
      DO I = 1,NDIM
        IFR = IFREEV(I)
        IF(IFR.GT.0) THEN
          IKKT = IKKT + 1
          RHSKKT(IKKT) = GVEC(I)
        ENDIF
      enddo
C
      IF ( NCON .GT. 0 ) THEN
        DO I=1,NCON
          IF(ISTATC(I).EQ.1) THEN
            IKKT = IKKT + 1
            RHSKKT(IKKT) = -BLWR(I)           
          ELSEIF(ISTATC(I).EQ.2) THEN
            IKKT = IKKT + 1
            RHSKKT(IKKT) = -BUPR(I)            
          ELSEIF(ISTATC(I).EQ.3) THEN
            IKKT = IKKT + 1
            RHSKKT(IKKT) = -BUPR(I)           
          ENDIF
        enddo
      ENDIF
C
C         SOLVE LINEAR SYSTEM
C
      CALL REFITR(RHSKKT,SOLNKT,NEQNS,RWORK,LNRWRK,NEEDS,IER)
C
      select case(ier)
      case(-600)
        IEREKT = -845
      case(-601)
        IEREKT = -1014
        NEEDED = NEEDS
      case(-602)
        IEREKT = -847
      case(-603)
        IEREKT = -848
      case(-604,-605)
        IEREKT = -1111
      case(-701)
        IEREKT = -750
      end select
      IF(IEREKT.NE.0) GO TO 280
C
C ------------------------------------------------------------------
C
C         --- STORE THE SOLUTION INTO PVEC 
C             NOTE THE SOLUTION OF THE KKT SYSTEM IS (-PVEC)
C
      pvec(1:ndim) = zero
      DO I = 1,NDIM
        IFR = IFREEV(I)
        IF(IFR.GT.0) PVEC(I) = -SOLNKT(IFR)
      enddo
C
C     INITIALIZE ALL NCON ENTRIES OF CONMLT TO ZERO. THEN FILL IN THE
C     CONMLT VALUES CORRESPONDING TO THE ACTIVE CONSTRAINTS, USING THE
C     KKT SYSTEM SOLUTION.
C
      IF ( NCON .GT. 0 ) THEN
        DO I=1,NCON
          CONMLT(I) = ZERO
          IF(IACTIV(I).GT.0) THEN
            CONMLT(I) = SOLNKT(IACTIV(I)+NFREE)
          ENDIF
        enddo
      ENDIF
C
C     INITIALIZE ALL NDIM ENTRIES OF VARMLT TO ZERO. THEN FILL IN THE
C     VARMLT VALUES CORRESPONDING TO THE FIXED VARIABLES, USING THE
C     KKT SYSTEM SOLUTION.
C
      VARMLT(1:NDIM) = ZERO
C
      IF(NCON.GT.0) CALL MVPSPR(11,NDIM,NCON,GMAT,IROW,JCOLST,
     $              CONMLT,RHSKKT)
C
      CALL SYMMVP(HMAT,IROWH,JSTRH,NDIM,PVEC,SOLNKT)
C
      DO I = 1,NDIM
        IF(IFREEV(I).EQ.0) VARMLT(I) = GVEC(I) - RHSKKT(I) + SOLNKT(I)
      enddo
C
C         COMPUTE OPTIMAL OBJECTIVE FUNCTION (USE RHSKKT AS WORK ARRAY)
C
      CALL SQFORM(HMAT,IROWH,JSTRH,GVEC,PVEC,RHSKKT,NDIM,QUAD)
C
C         COMPUTE LINEAR CONSTRAINT ESTIMATE AND SAVE IN RHSKKT
C
      CALL MVPSPR(1,NCON,NDIM,GMAT,IROW,JCOLST,PVEC,RHSKKT)
C
C         CHECK OPTIMALITY CONDITIONS FOR LINEAR PREDICTED POINT
C
      CALL BADSET(BLWR,BUPR,RHSKKT,CONMLT,VARMLT,PLWR,PUPR,PVEC,
     $    ISTATC,ISTATV,CONTOL,NCON,NCON,NDIM,BADACT)
      IF(BADACT) THEN
        IEREKT = -1
      ENDIF
C
 280  CONTINUE
C
      IF(DUMPKT) THEN
C
C         CHECK IF |IOFMFR| IS A VALID UNIT NUMBER FOR KT OUTPUT
C
        KKTOUT = ABS(IOFMFR)
        OPUN = .FALSE.
        INQUIRE(KKTOUT,OPENED=OPUN)
        IF(OPUN) THEN
          IEREKT = -852
        ELSE
          OPEN(KKTOUT,FILE='KKTMATRIX.FIL',STATUS='UNKNOWN')
C
          CALL XDSLWM ( 'ASCII','KKT MATRIX FROM SOCX', KKTOUT, 
     $                        RWORK, LNRWRK, IER )
          IF(IER.NE.0) THEN
            IEREKT = -852
          ELSE
            WRITE(KKTOUT,'(D25.17)') (RHSKKT(II),II=1,NEQNS)
          ENDIF
        ENDIF
C
        CLOSE(KKTOUT)
C
      ENDIF
C
      IF(IEREKT.NE.0.AND.IPC.GT.0) WRITE(IPU,1002) IEREKT
C
C         END K-T FACTORIZATION TIMING
C
      CALL CLKSUM(6)
      CALL CLKSUM(13)
C      
 1001 FORMAT(5X,'EQPKKT STORAGE ERROR;  NEEDS =',I6,'  LNRWRK =',I6)
 1002 FORMAT(5X,'MULTIFRONTAL ERROR;  IEREKT =',I6)
      RETURN 
      END
