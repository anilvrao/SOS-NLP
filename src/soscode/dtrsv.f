      SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, LDA, N
      CHARACTER(LEN=1)   DIAG, TRANS, UPLO
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DTRSV  SOLVES ONE OF THE SYSTEMS OF EQUATIONS
*
*     A*X = B,   OR   A'*X = B,
*
*  WHERE B AND X ARE N ELEMENT VECTORS AND A IS AN N BY N UNIT, OR
*  NON-UNIT, UPPER OR LOWER TRIANGULAR MATRIX.
*
*  NO TEST FOR SINGULARITY OR NEAR-SINGULARITY IS INCLUDED IN THIS
*  ROUTINE. SUCH TESTS MUST BE PERFORMED BEFORE CALLING THIS ROUTINE.
*
*  PARAMETERS
*  ==========
*
*  UPLO   - CHARACTER(LEN=1).
*           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX IS AN UPPER OR
*           LOWER TRIANGULAR MATRIX AS FOLLOWS:
*
*              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
*
*              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
*
*           UNCHANGED ON EXIT.
*
*  TRANS  - CHARACTER(LEN=1).
*           ON ENTRY, TRANS SPECIFIES THE EQUATIONS TO BE SOLVED AS
*           FOLLOWS:
*
*              TRANS = 'N' OR 'N'   A*X = B.
*
*              TRANS = 'T' OR 'T'   A'*X = B.
*
*              TRANS = 'C' OR 'C'   A'*X = B.
*
*           UNCHANGED ON EXIT.
*
*  DIAG   - CHARACTER(LEN=1).
*           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT
*           TRIANGULAR AS FOLLOWS:
*
*              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
*
*              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
*                                  TRIANGULAR.
*
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE ORDER OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY WITH  UPLO = 'U' OR 'U', THE LEADING N BY N
*           UPPER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE UPPER
*           TRIANGULAR MATRIX AND THE STRICTLY LOWER TRIANGULAR PART OF
*           A IS NOT REFERENCED.
*           BEFORE ENTRY WITH UPLO = 'L' OR 'L', THE LEADING N BY N
*           LOWER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE LOWER
*           TRIANGULAR MATRIX AND THE STRICTLY UPPER TRIANGULAR PART OF
*           A IS NOT REFERENCED.
*           NOTE THAT WHEN  DIAG = 'U' OR 'U', THE DIAGONAL ELEMENTS OF
*           A ARE NOT REFERENCED EITHER, BUT ARE ASSUMED TO BE UNITY.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, N ).
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCX ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE N
*           ELEMENT RIGHT-HAND SIDE VECTOR B. ON EXIT, X IS OVERWRITTEN
*           WITH THE SOLUTION VECTOR X.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
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
      LOGICAL            NOUNIT
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
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSEIF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSEIF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSEIF( N.LT.0 )THEN
         INFO = 4
      ELSEIF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSEIF( INCX.EQ.0 )THEN
         INFO = 8
      ENDIF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSV ', INFO )
         RETURN
      ENDIF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     SET UP THE START POINT IN X IF THE INCREMENT IS NOT UNITY. THIS
*     WILL BE  ( N - 1 )*INCX  TOO SMALL FOR DESCENDING LOOPS.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSEIF( INCX.NE.1 )THEN
         KX = 1
      ENDIF
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        FORM  X := INV( A )*X.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
                     ENDDO
                  ENDIF
                  JX = JX - INCX
               ENDDO
            ENDIF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               JX = KX
               DO J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
                     ENDDO
                  ENDIF
                  JX = JX + INCX
               ENDDO
            ENDIF
         ENDIF
      ELSE
*
*        FORM  X := INV( A' )*X.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO J = 1, N
                  TEMP = X( J )
                  DO I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
                  ENDDO
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
               ENDDO
            ELSE
               JX = KX
               DO J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
                  ENDDO
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
               ENDDO
            ENDIF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO J = N, 1, -1
                  TEMP = X( J )
                  DO I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
                   ENDDO
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
               ENDDO
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
                  ENDDO
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
               ENDDO
            ENDIF
         ENDIF
      ENDIF
*
      RETURN
*
*     END OF DTRSV .
*
      END
