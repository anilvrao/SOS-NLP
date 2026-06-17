      SUBROUTINE LAMERR( NX, EQLMTL, IFREKT, ALMSHR,
     .                   ALMBDA, ALMERR                  )
C
C ======================================================================
C     LAMERR===>lamerr   J.T. BETTS
C ======================================================================
C
C     ROUTINE TO COMPUTE WHETHER OR NOT THERE IS TOO LARGE A
C     DISCREPENCY BETWEEN THE LAGRANGE MULTIPLIERS
C     COMPUTED AS PART OF THE KT SYSTEM AND THOSE COMPUTED
C                                              T
C     DIRECTLY FROM ALMBDA = (G(X+P))FX - (A)FX PI
C
      LOGICAL ALMERR
C
      INTEGER NX, I
      INTEGER IFREKT(NX)
C
      DOUBLE PRECISION EQLMTL, ALMBDA(NX), ALMSHR(NX)
C
C                                              T
C     DIRECTLY FROM ALMBDA = (G(X+P))FX - (A)FX PI
C
      ALMERR = .FALSE.
C
      DO I = 1, NX
        IF ( IFREKT(I) .LT. 0 ) THEN
C
C         A VARIABLE FIXED BY AN UPDATE. THUS, ALMSHR
C         HAS BEEN COMPUTED AS PART OF THE KT MATRIX SOLUTION.
C
          IF ( ABS( ALMBDA(I) - ALMSHR(I) ) .GT. EQLMTL ) THEN
C
C           TOO LARGE A DISCREPENCY BETWEEN THE LAGRANGE MULTIPLIERS
C           COMPUTED AS PART OF THE KT SYSTEM AND THOSE COMPUTED
C                                                    T
C           DIRECTLY FROM ALMBDA = (G(X+P))FX - (A)FX PI
C
            ALMERR = .TRUE.
C
          ENDIF
C
        ENDIF
      enddo
C
C
      RETURN
      END
