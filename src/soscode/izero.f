      SUBROUTINE  IZERO(N,SX,INCX)
 
C
C     FILLS A VECTOR, X, WITH ZEROES.
C
C     ------------------------------------------------------------------
C
      INTEGER              SX(*)
      INTEGER              I,INCX,IX,N
 
C     ------------------------------------------------------------------
 
      IF(N.LE.0)RETURN
 
      IX = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        SX(IX) = 0
        IX = IX + INCX
   10 CONTINUE
 
C     ------------------------------------------------------------------
 
      RETURN
      END
