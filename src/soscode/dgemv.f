      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER(LEN=1)   TRANS
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGEMV  PERFORMS ONE OF THE MATRIX-VECTOR OPERATIONS
*
*     Y := ALPHA*A*X + BETA*Y,   OR   Y := ALPHA*A'*X + BETA*Y,
*
*  WHERE ALPHA AND BETA ARE SCALARS, X AND Y ARE VECTORS AND A IS AN
*  M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  TRANS  - CHARACTER(LEN=1).
*           ON ENTRY, TRANS SPECIFIES THE OPERATION TO BE PERFORMED AS
*           FOLLOWS:
*
*              TRANS = 'N' OR 'N'   Y := ALPHA*A*X + BETA*Y.
*
*              TRANS = 'T' OR 'T'   Y := ALPHA*A'*X + BETA*Y.
*
*              TRANS = 'C' OR 'C'   Y := ALPHA*A'*X + BETA*Y.
*
*           UNCHANGED ON EXIT.
*
*  M      - INTEGER.
*           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
*           M MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
*           CONTAIN THE MATRIX OF COEFFICIENTS.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCX ) ) WHEN TRANS = 'N' OR 'N'
*           AND AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCX ) ) OTHERWISE.
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE
*           VECTOR X.
*           UNCHANGED ON EXIT.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  BETA   - DOUBLE PRECISION.
*           ON ENTRY, BETA SPECIFIES THE SCALAR BETA. WHEN BETA IS
*           SUPPLIED AS ZERO THEN Y NEED NOT BE SET ON INPUT.
*           UNCHANGED ON EXIT.
*
*  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCY ) ) WHEN TRANS = 'N' OR 'N'
*           AND AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCY ) ) OTHERWISE.
*           BEFORE ENTRY WITH BETA NON-ZERO, THE INCREMENTED ARRAY Y
*           MUST CONTAIN THE VECTOR Y. ON EXIT, Y IS OVERWRITTEN BY THE
*           UPDATED VECTOR Y.
*
*  INCY   - INTEGER.
*           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           Y. INCY MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 2 BLAS ROUTINE.
*
*  -- WRITTEN ON 22-OCTOBER-1986.
*     JACK DONGARRA, ARGONNE NATIONAL LAB.
*     JEREMY DU CROZ, NAG CENTRAL OFFICE.
*     SVEN HAMMARLING, NAG CENTRAL OFFICE.
*     RICHARD HANSON, SANDIA NATIONAL LABS.
*
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSEIF( M.LT.0 )THEN
         INFO = 2
      ELSEIF( N.LT.0 )THEN
         INFO = 3
      ELSEIF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSEIF( INCX.EQ.0 )THEN
         INFO = 8
      ELSEIF( INCY.EQ.0 )THEN
         INFO = 11
      ENDIF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      ENDIF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     SET  LENX  AND  LENY, THE LENGTHS OF THE VECTORS X AND Y, AND SET
*     UP THE START POINTS IN  X  AND  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      ENDIF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      ENDIF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      ENDIF
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
*
*     FIRST FORM  Y := BETA*Y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO I = 1, LENY
                  Y( I ) = ZERO
               ENDDO
            ELSE
               DO I = 1, LENY
                  Y( I ) = BETA*Y( I )
               ENDDO
            ENDIF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
               ENDDO
            ELSE
               DO I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        FORM  Y := ALPHA*A*X + Y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
                  ENDDO
               ENDIF
               JX = JX + INCX
            ENDDO
         ELSE
            DO J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
                  ENDDO
               ENDIF
               JX = JX + INCX
            ENDDO
         ENDIF
      ELSE
*
*        FORM  Y := ALPHA*A'*X + Y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO  J = 1, N
               TEMP = ZERO
               DO I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
               ENDDO
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
            ENDDO
         ELSE
            DO J = 1, N
               TEMP = ZERO
               IX   = KX
               DO I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
               ENDDO
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
*
*     END OF DGEMV .
*
      END
