

      SUBROUTINE RITE  (IPRU,ARRAY,LENG)
C
C ======================================================================
C     RITE  ===>rite     J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION ARRAY(LENG)
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
        WRITE(IPRU,'(d25.16)') (ARRAY(II),II=ISTRT,IEND)
      ENDDO
C
      ISTRT = 1 + MAX(0,(NRECS-1))*MAXREC
      IEND = ISTRT + LASTRC - 1
      WRITE(IPRU,'(d25.16)') (ARRAY(II),II=ISTRT,IEND)
C
      RETURN
      END
