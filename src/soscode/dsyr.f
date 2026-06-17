      SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER(LEN=1)   UPLO
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYR   PERFORMS THE SYMMETRIC RANK 1 OPERATION
*
*     A := ALPHA*X*X' + A,
*
*  WHERE ALPHA IS A REAL SCALAR, X IS AN N ELEMENT VECTOR AND A IS AN
*  N BY N SYMMETRIC MATRIX.
*
*  PARAMETERS
*  ==========
*
*  UPLO   - CHARACTER(LEN=1).
*           ON ENTRY, UPLO SPECIFIES WHETHER THE UPPER OR LOWER
*           TRIANGULAR PART OF THE ARRAY A IS TO BE REFERENCED AS
*           FOLLOWS:
*
*              UPLO = 'U' OR 'U'   ONLY THE UPPER TRIANGULAR PART OF A
*                                  IS TO BE REFERENCED.
*
*              UPLO = 'L' OR 'L'   ONLY THE LOWER TRIANGULAR PART OF A
*                                  IS TO BE REFERENCED.
*
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE ORDER OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCX ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE N
*           ELEMENT VECTOR X.
*           UNCHANGED ON EXIT.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY WITH  UPLO = 'U' OR 'U', THE LEADING N BY N
*           UPPER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE UPPER
*           TRIANGULAR PART OF THE SYMMETRIC MATRIX AND THE STRICTLY
*           LOWER TRIANGULAR PART OF A IS NOT REFERENCED. ON EXIT, THE
*           UPPER TRIANGULAR PART OF THE ARRAY A IS OVERWRITTEN BY THE
*           UPPER TRIANGULAR PART OF THE UPDATED MATRIX.
*           BEFORE ENTRY WITH UPLO = 'L' OR 'L', THE LEADING N BY N
*           LOWER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE LOWER
*           TRIANGULAR PART OF THE SYMMETRIC MATRIX AND THE STRICTLY
*           UPPER TRIANGULAR PART OF A IS NOT REFERENCED. ON EXIT, THE
*           LOWER TRIANGULAR PART OF THE ARRAY A IS OVERWRITTEN BY THE
*           LOWER TRIANGULAR PART OF THE UPDATED MATRIX.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, N ).
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
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
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
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSEIF( N.LT.0 )THEN
         INFO = 2
      ELSEIF( INCX.EQ.0 )THEN
         INFO = 5
      ELSEIF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      ENDIF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR  ', INFO )
         RETURN
      ENDIF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     SET THE START POINT IN X IF THE INCREMENT IS NOT UNITY.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSEIF( INCX.NE.1 )THEN
         KX = 1
      ENDIF
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH THE TRIANGULAR PART
*     OF A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        FORM  A  WHEN A IS STORED IN UPPER TRIANGLE.
*
         IF( INCX.EQ.1 )THEN
            DO J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
                  ENDDO
               ENDIF
            ENDDO
         ELSE
            JX = KX
            DO J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
                  ENDDO
               ENDIF
               JX = JX + INCX
            ENDDO
         ENDIF
      ELSE
*
*        FORM  A  WHEN A IS STORED IN LOWER TRIANGLE.
*
         IF( INCX.EQ.1 )THEN
            DO J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
                  ENDDO
               ENDIF
            ENDDO
         ELSE
            JX = KX
            DO J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
                  ENDDO
               ENDIF
               JX = JX + INCX
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
*
*     END OF DSYR  .
*
      END
