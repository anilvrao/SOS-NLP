      SUBROUTINE UVGETN( IUPTYP, JVAR, MINEQL, MEQUAL, NROW,
     $                  NDIM, IPC, IPU, NEQNS, NX, K, KMAX,
     $                  NFREE0, IFREE0, IFREKT, SFEZ,
     $                  NFRSLK, NCACTV, CON2KT, ACT2CN, INSHUR,
     $                  MXNZAM, JSTRA, IROWA, AMAT,
     $                  MXNZHM, JSTRH, IROWH, HMAT, U, NZUVEC,
     $                  UVEC, IRWUVC, V, CHI, SIGMA, IER   )
C
C ======================================================================
C     UVGETN===>uvgetn    J.T. BETTS 
C ======================================================================
C
C
C     ROUTINE WHICH GIVES VALUES TO SIGMA AND THE VECTORS
C     u AND v FOR A SCHUR-COMPLEMENT UPDATE WITH TYPE INDICATED
C     BY THE VALUE OF IUPTYP. IN ADDITION A SPARSE STORAGE VERSION
C     OF u, WITH NZUVEC NONZEROS, IS STORED IN UVEC AND IRWUVC.
C
C     NOTES: JVAR IS THE INDEX OF THE VARIABLE (AMONG ALL VARS) TO BE
C           FIXED OR FREED. 
C
C                  { -2, FIX JVAR AT UPPER BOUND  
C        IUPTYP =  { -1, FIX JVAR AT LOWER BOUND
C                  {  1, FREE JVAR
C
C        NFREE0    NUMBER OF VARIABLES (INCLUDING SLACKS AND
C                  FEASIBILITY) THAT WERE FREE WHEN THE KT SYSTEM
C                  WAS FACTORED PRIOR TO THIS CALL TO SHURDV.
C                  NOTE: NEQNS SHOULD = NFREE0 + NROW
C        
C        IFREE0    VECTOR INDICATING FREE VARIABLES AT KT
C                  FACTORIZATION. IFREE0 IS PARTITIONED
C                  AS ( MINEQL, 1, NDIM ).
C                              {  -1, VAR I WAS FIXED AT LOW BOUND
C                  IFREE0(I) = {  -2, VAR I WAS FIXED AT UPPER BOUND
C                              { J>0, VAR I WAS THE J-TH FREE VAR
C                                    AMONG THE NFREE0 TOTAL FREE VARS
C        IFREKT    SCHUR-COMPLEMENT UPDATE HISTORY FROM THE POINT
C                  OF VIEW OF THE VARIABLES.
C                              { -KK, VAR I FIXED AT UPDATE KK
C                  IFREKT(I) = {   0, VAR I NOT AFFECTED BY AN UPDATE
C                              {  KK, VAR I FREED AT UPDATE KK
C        
C        INSHUR    SCHUR-COMPLEMENT UPDATE HISTORY INDICATING WHICH
C                  VARIABLE WAS INVOLVED.
C                  INSHUR(KK) = { -J, VAR J FIXED AT UPDATE KK
C                               {  J, VAR J FREED AT UPDATE KK
C
C
C
C     NFREE0 IS THE TOTAL NUMBER OF FREE VARS AT X0
C     NFRSLK IS THE NUMBER OF FREE SLACKS AT X0
C     NFRKT IS THE NO. OF FREE VARS REPRESENTED IN THE KKT MATRIX AT X0
C           SO, NFRKT =  NFREE0 - NFRSLK
C     ALSO, IF IFREE0(JVAR) > 0 FOR JVAR > MINEQL THEN
C     JVAR IS REPRESENTED IN THE INITIAL KKT MATRIX AS THE
C     VARIABLE WITH INDEX IFREE0(JVAR) - NFRSLK.
C
C
      DOUBLE PRECISION ZERO, ONE
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
C
C
      INTEGER IUPTYP, JVAR, MINEQL, NROW, NDIM, IPC, IPU
      INTEGER IER, NEQNS, K, I, II, NFREE0, MEQUAL 
      INTEGER KMAX, NX, MXNZAM, MXNZHM
      INTEGER IIROW, IOFFST, JREAL
      INTEGER INONZ, JJCOL, NZUVEC
      INTEGER IDEX, IDXSHR
      INTEGER NFRSLK, NCACTV
      INTEGER ICON, J, IKT, KK, ISLK, JSHUR
C
      INTEGER IFREE0(NX), IFREKT(NX)
      INTEGER JSTRA(NDIM+1), IROWA(MXNZAM)
      INTEGER JSTRH(NDIM+1), IROWH(MXNZHM)
      INTEGER IRWUVC(NEQNS)
      INTEGER CON2KT(NROW+1), ACT2CN(NROW+1), INSHUR(KMAX)
C
      DOUBLE PRECISION CHI
      DOUBLE PRECISION U(NEQNS), V(KMAX), SFEZ(*)
      DOUBLE PRECISION UVEC(NEQNS), SIGMA
      DOUBLE PRECISION AMAT(MXNZAM), HMAT(MXNZHM)
C
C
      IER = 0
C
C     ERROR CHECK ON THE VALUE OF JVAR.
      IF ( .NOT. ( (  JVAR .GE. 1 ) .AND.
     $                 ( JVAR .LE. MINEQL + 1 + NDIM ) ) ) THEN
C
        IER = 1
        IF ( IPC .GT. 0 ) THEN
          WRITE(IPU,*) ' UVGETN: ERROR, JVAR VALUE NOT IN PROPER',
     $                 ' RANGE, JVAR=', JVAR
        ENDIF
        GO TO 999
      ENDIF
C
C
C     u <-- 0.0, UVEC <-- 0.0, v<-- 0.0  AND SIGMA<-- 0.0 
C
      NZUVEC = 0
      SIGMA = 0.0
      U(1:NEQNS) = ZERO
      V(1:K) = ZERO
C
C
      IF ( IUPTYP .LT. 0 ) THEN
C
C       ...JVAR IS TO BE FIXED AT A BOUND.
C
        II = IFREE0(JVAR)
C
        IF ( .NOT. ( (  II .GE. 1 ) .AND. ( II .LE. NFREE0) 
     $               .AND. (IFREKT(JVAR) .EQ. 0 )          )  ) THEN
C
C         JVAR CANNOT BE FIXED AT A BOUND BY AN UPDATE BECAUSE
C         EITHER IT WAS NOT FREE AT THE PREVIOUS KT FACTORIZATION
C         OR IT WAS FREE AND LATER UPDATED. (IN THE LATER CASE,
C         IT IS EITHER NOW FIXED OR IT SHOULD HAVE BEEN DOWNDATED.)
C
          IER = 2
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' UVGETN: FREE VAR INDICATOR ERROR WHEN',
     $                   ' FIXING VAR NO.=', JVAR
          ENDIF
C
          GO TO 999
        ENDIF
C
        IF ( II .GT. NFRSLK ) THEN
C         FIX EITHER THE FEAS VAR OR A REAL VAR. THUS, IT WAS
C         REPRESENTED IN THE KKT MATRIX AT X0 AS THE IKT-TH FREE VAR.
C
          IKT = II - NFRSLK
          
C         U <-- (E)IKT, AND V AND SIGMA REMAIN ZERO.
          U(IKT) = ONE
C
          NZUVEC = NZUVEC + 1
          IRWUVC(NZUVEC) = IKT
          UVEC(NZUVEC) = ONE
        ELSE
C         FIX A SLACK VAR. THUS, ALL FREE VARS AT X0 HAVE TO BE
C         REPRESENTED IN U & UVEC, AND NONSLACK VARS FREED BY UPDATES
C         HAVE TO BE REPRESENTED IN V.
C
C         ICON < -- CONSTRAINT NUMBER CORRESPONDING TO SLACK VAR JVAR.
C
          ICON = MEQUAL + JVAR
C
C         ...FEAS VAR
C
          IDEX = IFREE0(MINEQL+1)
          IDXSHR = IFREKT(MINEQL+1)
          IF ( IDEX .GT. 0 ) THEN
C           FEAS VAR WAS IDEX-TH FREE VAR AT X0
C
            IF ( IDXSHR .GT. 0 ) THEN
              IER = 6
              IF ( IPC .GT. 0 ) THEN
                WRITE(IPU,*) ' UVGETN: ERROR, TRYING TO',
     $                 ' FREE A FREE VARIABLE'
              ENDIF 
              GO TO 999
            ENDIF
C
            U(IDEX - NFRSLK) = SFEZ(ICON)
            NZUVEC = NZUVEC + 1
            IRWUVC(NZUVEC) = IDEX - NFRSLK
            UVEC(NZUVEC) = U(IDEX - NFRSLK)
          ELSEIF ( IDXSHR .GT. 0 ) THEN
C           FEAS VAR WAS FREED AT SCHUR-COMPLEMENT UPDATE IDXSHR
C
            V(IDXSHR) = SFEZ(ICON)
          ENDIF
C
C
C         ...REAL VARS
C
          DO J= 1, NDIM
C
C           DETERMINE WHETHER OR NOT REAL VAR J HAS A NON ZERO ENTRY
C           IN THE JACOBIAN FOR CONSTRAINT ICON.
C
            INONZ = 0
            iiloop: DO II= JSTRA(J), JSTRA(J+1) - 1
              IF ( IROWA(II) .EQ. ICON ) THEN
                INONZ = II
                exit iiloop
              ENDIF
            enddo iiloop
C
            IF ( INONZ .GT. 0 ) THEN
C
              IDEX = IFREE0(MINEQL+1+J)
              IDXSHR = IFREKT(MINEQL+1+J)
C
              IF ( IDEX .GT. 0 ) THEN
C               REAL VAR J WAS IDEX-TH FREE VAR AT X0
C
                IF ( IDXSHR .GT. 0 ) THEN
                  IER = 7
                  IF ( IPC .GT. 0 ) THEN
                    WRITE(IPU,*) ' UVGETN: ERROR, TRYING TO',
     $                     ' FREE A FREE VARIABLE'
                  ENDIF 
                  GO TO 999
                ENDIF
C
C             
                U(IDEX - NFRSLK) = AMAT(INONZ)
                NZUVEC = NZUVEC + 1
                IRWUVC(NZUVEC) = IDEX - NFRSLK
                UVEC(NZUVEC) = U(IDEX - NFRSLK)
              ELSEIF ( IDXSHR .GT. 0 ) THEN
C               REAL VAR J WAS FREED AT SCHUR-COMPLEMENT UPDATE IDXSHR
C
                V(IDXSHR) = AMAT(INONZ)
              ENDIF
            ENDIF
          enddo
        ENDIF
C
      ELSE
C       ...JVAR IS TO BE FREED.
C          NOTE: IN THIS VERSION OF THE CODE CONSTRAINTS THAT
C                WERE NOT REPRESENTED IN THE KKT MATRIX AT X0 CAN
C                SURFACE DUE TO UPDATES THAT FIX THEIR CORRESPONDING
C                SLACK VARIABLE. THUS, WHEN FREEING, VAR JVAR, THE
C                "V" VECTOR MUST INCLUDE THE JVAR COL ENTRIES OF THE
C                CONSTRAINT MATRIX FOR ALL CONSTRAINT ROWS THAT
C                SURFACE DUE TO UPDATES.
C
        IF ( .NOT. ( ( IFREE0(JVAR) .LT. 0 ) .AND. 
     $                      ( IFREKT(JVAR) .EQ. 0 ) )  ) THEN
C
          IER = 3
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' UVGETN: FREE VAR INDICATOR ERROR WHEN',
     $                   ' FREEING VAR NO.=', JVAR
          ENDIF
C
          GO TO 999
C
C         EITHER JVAR WAS NOT FIXED AT THE PREVIOUS KT FACTORIZATION
C         OR IT WAS UPDATED AND NOT DOWNDATED PROPERLY. 
C
        ENDIF
C     
        IF ( JVAR .LE. MINEQL ) THEN
C         THE (JVAR)-TH SLACK VAR IS TO BE FREED. (THIS ASSUMES THAT
C         THE SLACK VARS HAVE NO CONTRIBUTIONS TO THE HESSIAN, THE
C         ONLY NONZERO CONTRIBUTIONS TO U OR V IS THE UNIT VECTOR
C         SEGMENT OF U CORRESPONDING TO THE "SLACK IDENTITY MATRIX"
C         APPENDED TO THE end OF THE INEQUALITY CONSTRAINTS.
C
          ICON = MEQUAL + JVAR
C
          IF ( CON2KT(ICON) .LE. 0 ) THEN
            IER = 8
            IF ( IPC .GT. 0 ) THEN
              WRITE(IPU,*) ' UVGETN: NON-POSITIVE CON2KT MAP WHEN',
     $                   ' FREEING SLACK VAR NO.=', JVAR
            ENDIF
C
            GO TO 999
          ENDIF

          IKT = NFREE0 - NFRSLK + CON2KT(ICON)
          U( IKT ) = ONE
C
          NZUVEC = NZUVEC + 1
          IRWUVC(NZUVEC) = IKT
          UVEC(NZUVEC) = ONE
C
        ELSEIF ( JVAR .EQ. MINEQL + 1 ) THEN
C         THE FEASIBILITY (ARTIFICIAL) VARIABLE XI IS TO BE FREED.
C         SINCE 0.5*CHI*(XI**2) APPEARS IN THE OBJECTIVE, XI
C         CONTRIBUTES ONLY CHI TO THE APPROPRIATE DIAGONAL IN
C         THE HESSIAN PORTION OF THE KKT MATRIX.
C         THUS, THE ONLY NONZEROS IN V COME FROM SFEZ ENTRIES
C         CORRESPONDING TO SLACKS THAT WERE FREED BY SCHUR-COMPLEMENT
C         UPDATES. THE ONLY NONZERO PORTION OF U IS THE SCALED
C         RESIDUAL VECTOR SFEZ WHICH IS ENTERED IN THE "CONSTRAINT"
C         SEGMENT OF U.
C
          SIGMA = CHI
C
C
C         COMPUTE THE U VECTOR.
C
          DO I = 1, NCACTV
            ICON = ACT2CN(I)
            IF ( SFEZ(ICON) .NE. ZERO ) THEN
              IKT = NFREE0 - NFRSLK + I
              U(IKT) = SFEZ(ICON)
C
              NZUVEC = NZUVEC + 1
              IRWUVC(NZUVEC) = IKT
              UVEC(NZUVEC) = U(IKT)
C
            ENDIF
          enddo
C
C
C         COMPUTE THE V VECTOR
C
          DO I= 1, K
            JSHUR = ABS(INSHUR(I))
            IF ( INSHUR(I) .LT. 0 .AND. JSHUR .LE. MINEQL ) THEN
C             SLACK VAR JSHUR WAS FIXED AT THE I-TH SCHUR-COMPLEMENT
C             UPDATE.
C
              V(I) = SFEZ(MEQUAL+JSHUR)
            ENDIF
          enddo
C
        ELSE
C         ( JVAR > MINEQL + 1 ) 
C         THE ORIGINAL PROBLEM VARIABLE JVAR IS TO BE FREED.
C
C         ...ENTER THE PORTION OF U CORRESPONDING TO THE CONSTRAINT
C            SEGMENT OF THE KT MATRIX. ALSO, ENTER THE PORTION
C            OF THE V VECTOR CORRESPONDING TO INEQUALITY CONSTRAINTS
C            MADER ACTIVE BY SCUR-COMPLEMENT UPDATES.
C
C         JREAL <-- INDEX OF JVAR AS A VARIABLE IN THE ORIGINAL PROBLEM.
C                   JREAL IS THE INDEX USED FOR JVAR IN AMAT AND HMAT!
C       
          JREAL = JVAR - ( MINEQL + 1 )
C
          DO INONZ = JSTRA(JREAL), JSTRA(JREAL+1) - 1
            IIROW = IROWA(INONZ)
            ISLK = IIROW - MEQUAL
            IF ( CON2KT(IIROW) .GT. 0 ) THEN
C
C             CONSTRAINT IIROW WAS IN THE KKT MATRIX AT X0, FILL IN
C             ITS U VECTOR ENTRY.
C
              IKT = NFREE0 - NFRSLK + CON2KT(IIROW)
C
              U(IKT) = AMAT(INONZ)
C
              NZUVEC = NZUVEC + 1
              IRWUVC(NZUVEC) = IKT
              UVEC(NZUVEC) = U(IKT)
            ELSEIF ( ISLK .GT. 0 .AND. IFREKT(ISLK) .LT. 0 ) THEN
              KK = ABS(IFREKT(ISLK))
C
C             CONSTRAINT NUMBER MEQUAL+ISLK WAS MADE ACTIVE BI FIXING 
C             ISLK AT SCHUR-COMPLEMENT UPDATE KK. PUT THE APPROPRIATE
C             JACOBIAN ENTRY IN THE V VECTOR FOR JREAL.
C
              V(KK) = AMAT(INONZ)
            ENDIF
C
          enddo
C
C
C         ...ENTER THE CONTRIBUTIONS OF COL JVAR OF THE FREE PORTION
C            OF THE HESSIAN TO U AND V. ALSO, SIGMA <-- HESSIAN DIAGONAL
C            FOR JVAR.
C
C         SEARCH ALL THE COLS OF THE HESSIAN WITH INDEX < JREAL TO
C         FIND ENTRIES IN ROW JREAL. (THESE ARE ENTRIES IN COL JREAL
C         ABOVE THE HESSIAN DIAGONAL.) THEN ACCOUNT FOR HESSIAN ENTRIES
C         ON OR BELOW THE HESSIAN DIAGONAL.
C
          DO JJCOL = 1, JREAL
            inonloop: DO INONZ = JSTRH(JJCOL), JSTRH(JJCOL+1) - 1
              IIROW = IROWH(INONZ)
C
              IF ( JJCOL .LT. JREAL ) THEN
C 
                IF ( IIROW .EQ. JREAL ) THEN
C
C                 PUT HESSIAN ENTRY (JREAL,JJCOL) == (JJCOL,REAL) INTO
C                 U OR V, DEPENDING ON THE VALUES ON IFREE0 AND IFREKT.
C
                  IOFFST = MINEQL + 1
                  IDEX = IFREE0(IOFFST+JJCOL)
                  IDXSHR = IFREKT(IOFFST+JJCOL)
C
                  IF ( IDEX .GT. 0 ) THEN           
                    IF ( IDXSHR .GT. 0 ) THEN
                      IER = 4
                      IF ( IPC .GT. 0 ) THEN
                        WRITE(IPU,*) ' UVGETN: ERROR, TRYING TO',
     $                    ' FREE A FREE VARIABLE'
                      ENDIF 
                      GO TO 999
                    ENDIF
C
                    U(IDEX-NFRSLK) = HMAT(INONZ)
C
                    NZUVEC = NZUVEC + 1
                    IRWUVC(NZUVEC) = IDEX - NFRSLK
                    UVEC(NZUVEC) = U(IDEX-NFRSLK)
C
                  ELSEIF ( IDXSHR .GT. 0 ) THEN
                    V(IDXSHR) = HMAT(INONZ)
                  ENDIF
C
                ELSEIF ( IIROW .GT. JREAL ) THEN
C                 NO MORE ENTRIES FOR ROW JREAL IN COL JJCOL.
C                 THIS ASSUMES HESSIAN ENTRIES ARE ARE STORED
C                 WITH ASCENDING ROW INDICES WITHIN EACH COLUMN.
C
                  exit inonloop
                ENDIF
              ELSE
C
C               PROCESSING COL JREAL, ALL FREE HESSIAN ENTRIES BELONG
C               EITHER IN U, V OR SIGMA.
C
                IF ( IIROW .EQ. JREAL ) THEN
C                 SIGMA <-- HESSIAN DIAGONAL JREAL
                  SIGMA = HMAT(INONZ)
                ELSE
C
C                 ACCOUNT FOR THE FREE PART OF HESSIAN COL JREAL,
C                 BELOW THE DIAGONAL.(IIROW HAD BETTER BE > JREAL HERE.)
C
                  IOFFST = MINEQL + 1
                  IDEX = IFREE0(IOFFST+IIROW)
                  IDXSHR = IFREKT(IOFFST+IIROW)
C
                  IF ( IDEX .GT. 0 ) THEN           
C
                    IF ( IDXSHR .GT. 0 ) THEN
                      IER = 5
                      IF ( IPC .GT. 0 ) THEN
                        WRITE(IPU,*) ' UVGETN: FREE VARIABLE',
     $                               ' INDEXING ERROR'
                      ENDIF 
                      GO TO 999
                    ENDIF
C
                    U(IDEX-NFRSLK) = HMAT(INONZ)
C
                    NZUVEC = NZUVEC + 1
                    IRWUVC(NZUVEC) = IDEX-NFRSLK
                    UVEC(NZUVEC) = U(IDEX-NFRSLK)
C
                  ELSEIF ( IDXSHR .GT. 0 ) THEN
                    V(IDXSHR) = HMAT(INONZ)
                  ENDIF
C
                ENDIF
C
              ENDIF
C
            enddo inonloop
C
          enddo     
C
        ENDIF
C
      ENDIF
C
C
 999  CONTINUE
C
      RETURN
      END
