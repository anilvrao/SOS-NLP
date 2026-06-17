C
C
C
      SUBROUTINE INRTCK( K, KMAX, IPU, IPC, SHURC, KPVT,
     $                    ITPSHR, INRTER )
C
C ======================================================================
C     INRTCK===>inrtck   J.T. BETTS
C ======================================================================
C
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D0 )
C
      INTEGER I, K, IPC, IPU, INLOOP, INRTER
      INTEGER KP, KMAX, KPVT(KMAX), ITPSHR(KMAX)
      INTEGER NPOSHR, NNGSHR, NZRSHR
      INTEGER NPOSQP, NNEGQP
C
      DOUBLE PRECISION SHURC(KMAX,KMAX), DIAG
      DOUBLE PRECISION A, B, C, DISCR, ALAM1, ALAM2
C
C
      INRTER = 0
C
C
C     DETERMINE THE INERTIA OF THE SCHUR-COMPLEMENT MATRIX.
C                                      T          T
C     SINCE U1*P1*...*UN*PN*SHURC*PN*UN *...*P1*U1  = D,
C     WHERE THE U'S ARE BLOCK "ELEMENTARY ELIMINATORS"
C     AND THE P'S ARE ELEMENTARY PERMUTATION MATRICES, SHURC HAS
C     THE SAME INERTIA AS THE BLOCK DIAGONAL MATRIX D (WITH
C     1 BY 1 AND 2 BY 2 BLOCKS).
C
      NPOSHR = 0
      NNGSHR = 0
      NZRSHR = 0
      INLOOP = 1
 100  CONTINUE
      IF ( INLOOP .GT. K ) GO TO 110
C
        KP = KPVT(INLOOP)
C
        IF ( KP .GT. 0 ) THEN
C         1 BY 1 PIVOT
C
          DIAG = SHURC(INLOOP,INLOOP)
C
          INLOOP = INLOOP + 1
C
          IF ( DIAG .GT. ZERO ) THEN
            NPOSHR = NPOSHR + 1
          ELSEIF ( DIAG .LT. ZERO ) THEN
            NNGSHR = NNGSHR + 1
          ELSE
            NZRSHR = NZRSHR + 1
          ENDIF
C
        ELSE
C
C         2 BY 2 PIVOT
C
          A = SHURC(INLOOP,INLOOP)
          B = SHURC(INLOOP+1,INLOOP+1)
          C = SHURC(INLOOP,INLOOP+1)
C
C
          DISCR = (A-B)**2 + 4.0D0*C**2
C
          ALAM1 = 0.5D0* ( A + B + SQRT(DISCR) )
          ALAM2 = 0.5D0* ( A + B - SQRT(DISCR) )
C
          IF ( ALAM1 .GT. ZERO ) THEN
            NPOSHR = NPOSHR + 1
          ELSEIF ( ALAM1 .LT. ZERO ) THEN
            NNGSHR = NNGSHR + 1
          ELSE
            NZRSHR = NZRSHR + 1
          ENDIF
C
          IF ( ALAM2 .GT. ZERO ) THEN
            NPOSHR = NPOSHR + 1
          ELSEIF ( ALAM2 .LT. ZERO ) THEN
            NNGSHR = NNGSHR + 1
          ELSE
            NZRSHR = NZRSHR + 1
          ENDIF
C
          INLOOP = INLOOP + 2
        ENDIF
C
      GOTO 100
 110  CONTINUE
C
C
C     DETERMINE THE INERTIA CHANGE COUNTS REQUIRED BY
C     ACTIVE SET CHANGES IN THE QP. (POSITIVE FOR DROPS
C     AND NEGATIVE FOR ADDS.)
C
      NPOSQP = 0
      NNEGQP = 0
      DO I = 1, K
C
        IF ( ITPSHR(I) .GT. 0 ) THEN
          NPOSQP = NPOSQP + 1
        ELSEIF ( ITPSHR(I) .LT. 0 ) THEN
          NNEGQP = NNEGQP + 1
        ENDIF
C
      enddo
C
C
C     DETERMINE WHETHER THE INERTIA OF THE SCHUR-COMPLEMENT
C     MATRIX IS THAT REQUIRED BY THE ACTIVE SET CHANGES IN THE
C     QP. 
C
      IF ( NPOSHR .NE. NPOSQP .OR. NNGSHR .NE. NNEGQP ) THEN
        INRTER = 1
      ENDIF
C
      IF ( IPC .GT. 100 ) THEN
        CALL IVPRIN( K, KPVT, 'KPVT', 4, IPU )
C
        CALL IVPRIN( K, ITPSHR, 'ITPSHR', 8, IPU )
      ENDIF
C
      IF ( IPC .GT. 0 ) THEN
        WRITE(IPU,900) NPOSHR, NNGSHR, NZRSHR, NPOSQP, NNEGQP
 900    FORMAT( 2X, 'INERTIA OF SCHUR COMPLEMENT = ( ', I5, ',', I5,
     $          ',', I5, ' )' , /, 2X,
     $          'INERTIA REQUIRED BY QP = ( ',  I5, ',', I5,
     $          ',   0 )'  )
C
      ENDIF
C
C
      RETURN
      END
