      SUBROUTINE UVGET ( IUPTYP, JVAR, MINEQL, MEQUAL, NROW,
     .                  NDIM, IPC, IPU, NEQNS, NX, K, KMAX,
     .                  NFREE0, IFREE0, IFREKT, SFEZ,
     .                  MXNZAM, JSTRA, IROWA, AMAT,
     .                  MXNZHM, JSTRH, IROWH, HMAT, U, NZUVEC,
     .                  UVEC, IRWUVC, V, CHI, SIGMA, IER   )
C
C ======================================================================
C     UVGET ===>ivget    J.T. BETTS
C ======================================================================
C
C
C     ROUTINE WHICH GIVES VALUES TO SIGMA AND THE VECTORS
C     u AND v FOR A SCHUR-COMPLEMENT UPDATE WITH TYPE INDICATED
C     BY THE VALUE OF IUPTYP. IN ADDITION A SPARSE STORAGE VERSION
C     OF u, WITH NZUVEC NONZEROS, IS STORED IN UVEC AND IRWUVC.
C
C     NOTE: JVAR IS THE INDEX OF THE VARIABLE (AMONG ALL VARS) TO BE
C           FIXED OR FREED. 
C
C                     { -2, FIX JVAR AT UPPER BOUND  
C           IUPTYP =  { -1, FIX JVAR AT LOWER BOUND
C                     {  1, FREE JVAR
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
C
      INTEGER IFREE0(NX), IFREKT(NX)
      INTEGER JSTRA(NDIM+1), IROWA(MXNZAM)
      INTEGER JSTRH(NDIM+1), IROWH(MXNZHM)
      INTEGER IRWUVC(NEQNS)
      INTEGER IDEX, IDXSHR
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
     .                 ( JVAR .LE. MINEQL + 1 + NDIM ) ) ) THEN
C
        IER = 1
        IF ( IPC .GT. 0 ) THEN
          WRITE(IPU,*) ' UVGET: ERROR, JVAR VALUE NOT IN PROPER',
     .                 ' RANGE, JVAR=', JVAR
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
C       JVAR IS TO BE FIXED AT A BOUND.
C
        II = IFREE0(JVAR)
C
        IF ( .NOT. ( (  II .GE. 1 ) .AND. ( II .LE. NFREE0) 
     .               .AND. (IFREKT(JVAR) .EQ. 0 )          )  ) THEN
C
C         JVAR CANNOT BE FIXED AT A BOUND BY AN UPDATE BECAUSE
C         EITHER IT WAS NOT FREE AT THE PREVIOUS KT FACTORIZATION
C         OR IT WAS FREE AND LATER UPDATED. (IN THE LATER CASE,
C         IT IS EITHER NOW FIXED OR IT SHOULD HAVE BEEN DOWNDATED.)
C
          IER = 2
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' UVGET: FREE VAR INDICATOR ERROR WHEN',
     .                   ' FIXING VAR NO.=', JVAR
          ENDIF
C
          GO TO 999
        ENDIF
C
C       U <-- (E)II AND V REMAINS ZERO.
        U(II) = ONE
C
        NZUVEC = NZUVEC + 1
        IRWUVC(NZUVEC) = II
        UVEC(NZUVEC) = ONE
C
      ELSE
C       JVAR IS TO BE FREED.
C
        IF ( .NOT. ( ( IFREE0(JVAR) .LT. 0 ) .AND. 
     .                      ( IFREKT(JVAR) .EQ. 0 ) )  ) THEN
C
          IER = 3
          IF ( IPC .GT. 0 ) THEN
            WRITE(IPU,*) ' UVGET: FREE VAR INDICATOR ERROR WHEN',
     .                   ' FREEING VAR NO.=', JVAR
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
C         THE SLACK VARS CONSTITUTE THE FIRST MINEQL VARIABLES.)
C         SINCE SLACK VARS HAVE NO CONTRIBUTIONS TO THE HESSIAN, THE
C         ONLY NONZERO CONTRIBUTIONS TO U OR V IS THE UNIT VECTOR
C         SEGMENT OF U CORRESPONDING TO THE "SLACK IDENTITY MATRIX"
C         APPENDED TO THE FRONT OF THE INEQUALITY CONSTRAINTS.
C
          U( NFREE0 + MEQUAL + JVAR ) = ONE
C
          NZUVEC = NZUVEC + 1
          IRWUVC(NZUVEC) = NFREE0 + MEQUAL + JVAR
          UVEC(NZUVEC) = ONE
C
        ELSEIF ( JVAR .EQ. MINEQL + 1 ) THEN
C         THE FEASIBILITY (ARTIFICIAL) VARIABLE XI IS TO BE FREED.
C         SINCE 0.5*CHI*(XI**2) APPEARS IN THE OBJECTIVE, XI
C         CONTRIBUTES ONLY
C         CHI TO THE APPROPRIATE DIAGONAL IN THE HESSIAN PORTION OF
C         THE KT MATRIX (THUS, ALL OF V IS ZERO).
C         THE ONLY NONZERO PORTION OF U IS THE SCALED
C         RESIDUAL VECTOR SFEZ WHICH IS ENTERED IN THE "CONSTRAINT"
C         SEGMENT OF U.
C
          SIGMA = CHI
C
          DO I = 1, NROW
            IF ( SFEZ(I) .NE. ZERO ) THEN
              U( NFREE0 + I ) = SFEZ(I)
C
              NZUVEC = NZUVEC + 1
              IRWUVC(NZUVEC) = NFREE0 + I
              UVEC(NZUVEC) = U( NFREE0 + I )
C
            ENDIF
          enddo
C
        ELSE
C         ( JVAR > MINEQL + 1 ) 
C         THE ORIGINAL PROBLEM VARIABLE JVAR IS TO BE FREED.
C
C         ...ENTER THE PORTION OF U CORRESPONDING TO THE CONSTRAINT
C            SEGMENT OF THE KT MATRIX.
C
C         JREAL <-- INDEX OF JVAR AS A VARIABLE IN THE ORIGINAL PROBLEM.
C                   JREAL IS THE INDEX USED FOR JVAR IN AMAT AND HMAT!
C       
          JREAL = JVAR - ( MINEQL + 1 )
C
          IOFFST = NFREE0
          DO INONZ = JSTRA(JREAL), JSTRA(JREAL+1) - 1
            IIROW = IROWA(INONZ)
            U( IOFFST + IIROW ) = AMAT(INONZ)
C
            NZUVEC = NZUVEC + 1
            IRWUVC(NZUVEC) = IOFFST + IIROW
            UVEC(NZUVEC) = U( IOFFST + IIROW )
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
                        WRITE(IPU,*) ' UVGET: ERROR, TRYING TO',
     .                    ' FREE A FREE VARIABLE'
                      ENDIF 
                      GO TO 999
                    ENDIF
C
                    U(IDEX) = HMAT(INONZ)
C
                    NZUVEC = NZUVEC + 1
                    IRWUVC(NZUVEC) = IDEX
                    UVEC(NZUVEC) = U(IDEX)
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
c
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
C                 BELOW THE DIAGONAL.
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
                        WRITE(IPU,*) ' UVGET: FREE VARIABLE',
     .                               ' INDEXING ERROR'
                      ENDIF 
                      GO TO 999
                    ENDIF
C
                    U(IDEX ) = HMAT(INONZ)
C
                    NZUVEC = NZUVEC + 1
                    IRWUVC(NZUVEC) = IDEX
                    UVEC(NZUVEC) = U(IDEX)
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
