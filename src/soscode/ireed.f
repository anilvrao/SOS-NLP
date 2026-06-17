

      
      SUBROUTINE IREED (IPRU,IARRAY,LENG)
C
C ======================================================================
C     IREED ===>ireed    J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IARRAY(LENG)
C
      PARAMETER (MAXREC=100000)
C
      NRECS = LENG/MAXREC
C
      LASTRC = LENG - NRECS*MAXREC
C
      DO I = 1,NRECS
        ISTRT = 1 + (I-1)*MAXREC
        IEND = ISTRT + MAXREC - 1
        READ(IPRU,'(i25)') (IARRAY(II),II=ISTRT,IEND)
      ENDDO
C
      ISTRT = 1 + MAX(0,(NRECS-1))*MAXREC
      IEND = ISTRT + LASTRC - 1
      READ(IPRU,'(i25)') (IARRAY(II),II=ISTRT,IEND)
C
      RETURN
      END
