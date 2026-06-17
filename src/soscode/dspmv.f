      SUBROUTINE DSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, N
      CHARACTER(LEN=1)   UPLO
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   AP( * ), X( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSPMV  PERFORMS THE MATRIX-VECTOR OPERATION
*
*     Y := ALPHA*A*X + BETA*Y,
*
*  WHERE ALPHA AND BETA ARE SCALARS, X AND Y ARE N ELEMENT VECTORS AND
*  A IS AN N BY N SYMMETRIC MATRIX, SUPPLIED IN PACKED FORM.
*
*  PARAMETERS
*  ==========
*
*  UPLO   - CHARACTER(LEN=1).
*           ON ENTRY, UPLO SPECIFIES WHETHER THE UPPER OR LOWER
*           TRIANGULAR PART OF THE MATRIX A IS SUPPLIED IN THE PACKED
*           ARRAY AP AS FOLLOWS:
*
*              UPLO = 'U' OR 'U'   THE UPPER TRIANGULAR PART OF A IS
*                                  SUPPLIED IN AP.
*
*              UPLO = 'L' OR 'L'   THE LOWER TRIANGULAR PART OF A IS
*                                  SUPPLIED IN AP.
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
*  AP     - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( ( N*( N + 1 ) )/2 ).
*           BEFORE ENTRY WITH UPLO = 'U' OR 'U', THE ARRAY AP MUST
*           CONTAIN THE UPPER TRIANGULAR PART OF THE SYMMETRIC MATRIX
*           PACKED SEQUENTIALLY, COLUMN BY COLUMN, SO THAT AP( 1 )
*           CONTAINS A( 1, 1 ), AP( 2 ) AND AP( 3 ) CONTAIN A( 1, 2 )
*           AND A( 2, 2 ) RESPECTIVELY, AND SO ON.
*           BEFORE ENTRY WITH UPLO = 'L' OR 'L', THE ARRAY AP MUST
*           CONTAIN THE LOWER TRIANGULAR PART OF THE SYMMETRIC MATRIX
*           PACKED SEQUENTIALLY, COLUMN BY COLUMN, SO THAT AP( 1 )
*           CONTAINS A( 1, 1 ), AP( 2 ) AND AP( 3 ) CONTAIN A( 2, 1 )
*           AND A( 3, 1 ) RESPECTIVELY, AND SO ON.
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
*  BETA   - DOUBLE PRECISION.
*           ON ENTRY, BETA SPECIFIES THE SCALAR BETA. WHEN BETA IS
*           SUPPLIED AS ZERO THEN Y NEED NOT BE SET ON INPUT.
*           UNCHANGED ON EXIT.
*
*  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCY ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY Y MUST CONTAIN THE N
*           ELEMENT VECTOR Y. ON EXIT, Y IS OVERWRITTEN BY THE UPDATED
*           VECTOR Y.
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
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
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
         INFO = 6
      ELSEIF( INCY.EQ.0 )THEN
         INFO = 9
      ENDIF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSPMV ', INFO )
         RETURN
      ENDIF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     SET UP THE START POINTS IN  X  AND  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF THE ARRAY AP
*     ARE ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH AP.
*
*     FIRST FORM  Y := BETA*Y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO I = 1, N
                  Y( I ) = ZERO
               ENDDO
            ELSE
               DO I = 1, N
                  Y( I ) = BETA*Y( I )
               ENDDO
            ENDIF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
               ENDDO
            ELSE
               DO I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        FORM  Y  WHEN AP CONTAINS THE UPPER TRIANGLE.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               K     = KK
               DO I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + AP( K )*X( I )
                  K      = K      + 1
               ENDDO
               Y( J ) = Y( J ) + TEMP1*AP( KK + J - 1 ) + ALPHA*TEMP2
               KK     = KK     + J
            ENDDO
         ELSE
            JX = KX
            JY = KY
            DO J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO K = KK, KK + J - 2
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + AP( K )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
               ENDDO
               Y( JY ) = Y( JY ) + TEMP1*AP( KK + J - 1 ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + J
            ENDDO
         ENDIF
      ELSE
*
*        FORM  Y  WHEN AP CONTAINS THE LOWER TRIANGLE.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*AP( KK )
               K      = KK           + 1
               DO I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + AP( K )*X( I )
                  K      = K      + 1
               ENDDO
               Y( J ) = Y( J ) + ALPHA*TEMP2
               KK     = KK     + ( N - J + 1 )
            ENDDO
         ELSE
            JX = KX
            JY = KY
            DO J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*AP( KK )
               IX      = JX
               IY      = JY
               DO K = KK + 1, KK + N - J
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + AP( K )*X( IX )
               ENDDO
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + ( N - J + 1 )
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
*
*     END OF DSPMV .
*
      END
