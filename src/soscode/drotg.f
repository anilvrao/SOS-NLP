      SUBROUTINE DROTG(SA,SB,SC,SS)
C
C     DESIGNED BY C.L.LAWSON, JPL, 1977 SEPT 08
C
C
C     CONSTRUCT THE GIVENS TRANSFORMATION
C
C         ( SC  SS )
C     G = (        ) ,    SC**2 + SS**2 = 1 ,
C         (-SS  SC )
C
C     WHICH ZEROS THE SECOND ENTRY OF THE 2-VECTOR  (SA,SB)**T .
C
C     THE QUANTITY R = (+/-)SQRT(SA**2 + SB**2) OVERWRITES SA IN
C     STORAGE.  THE VALUE OF SB IS OVERWRITTEN BY A VALUE Z WHICH
C     ALLOWS SC AND SS TO BE RECOVERED BY THE FOLLOWING ALGORITHM@
C           IF Z=1  SET  SC=0.  AND  SS=1.
C           IF ABS(Z) .LT. 1  SET  SC=SQRT(1-Z**2)  AND  SS=Z
C           IF ABS(Z) .GT. 1  SET  SC=1/Z  AND  SS=SQRT(1-SC**2)
C
C     NORMALLY, THE SUBPROGRAM SROT(N,SX,INCX,SY,INCY,SC,SS) WILL
C     NEXT BE CALLED TO APPLY THE TRANSFORMATION TO A 2 BY N MATRIX.
C
C ------------------------------------------------------------------
C
      DOUBLE PRECISION SA, SB, SC, SS, R, V, U
C
      IF (ABS(SA) .GT. ABS(SB)) THEN
C
C *** HERE ABS(SA) .GT. ABS(SB) ***
C
        U = SA + SA
        V = SB / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF SA
C
        R = SQRT(.25D0 + V**2) * U
C
C     NOTE THAT SC IS POSITIVE
C
        SC = SA / R
        SS = V * (SC + SC)
        SB = SS
        SA = R
        RETURN
      ENDIF
C
C *** HERE ABS(SA) .LE. ABS(SB) ***
C
      IF (SB .EQ. 0.) THEN
C
C *** HERE SA = SB = 0. ***
C
        SC = 1.
        SS = 0.
        RETURN
      ENDIF
      U = SB + SB
      V = SA / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF SB
C     (R IS IMMEDIATELY STORED IN SA)
C
      SA = SQRT(.25D0 + V**2) * U
C
C     NOTE THAT SS IS POSITIVE
C
      SS = SB / SA
      SC = V * (SS + SS)
      IF (SC .NE. 0.D0) THEN
        SB = 1. / SC
      ELSE
        SB = 1.
      ENDIF
C
      RETURN
      END
