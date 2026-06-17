

      SUBROUTINE VECRYT(MAXNUM,N,AVEC,IWRK,TITLE,IPU)
C
C ======================================================================
C     VECRYT===>vecryt   J.T. BETTS
C ======================================================================
C     MODIFICATIONS:   13-APR-99 DSN, ADJUSTIBLE ARRAY DIMENSION FIX
C                                     FOR VAX-VMS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         INPUT:
C
C           MAXNUM MAX. NUMBER OF ELEMENTS PRINTED 
C           N      DIMENSION OF AVEC
C           AVEC   OUTPUT VECTOR
C           IWRK   INTEGER WORK ARRAY (N)
C           IPU    OUTPUT UNIT NO.
C
C        THIS ROUTINE WRITES THE N VECTOR AVEC IN COMPRESSED
C        FORMAT ON OUTPUT UNIT IPU.  
C
      DIMENSION  AVEC(*),IWRK(*)
      CHARACTER(LEN=*)  TITLE
C
      IF(N.EQ.0) RETURN
C
C         SORT AVEC ARRAY TO CONSTRUCT ORDER ARRAY AND THEN
C         RESTORE THE ORIGINAL ORDER OF AVEC.
C
      IOPT = 1
      ISTBLE = 0
      CALL HDSRTN(AVEC,N,IOPT,ISTBLE,IWRK,IERSRT)
      IF(IERSRT.NE.0) RETURN
C
      MAXPRN = MIN(N,MAXNUM)
      WRITE(IPU,1002) TITLE,MAXPRN,MAXPRN
C
C         WRITE OUT THE MAXPRN LARGEST AND SMALLEST ELEMENTS
C
      DO I = 1,MAXPRN
        WRITE(IPU,1001) IWRK(I),AVEC(IWRK(I))
     $                 ,IWRK(N-I+1),AVEC(IWRK(N-I+1))
      ENDDO
C
      RETURN
 1001 FORMAT(T3,'*',T6,2(5X,'(',I6,':',1PG16.8,')'),T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,A,T106,'*',/T3,'*',T106,'*'/
     $    T3,'*',T11,'Smallest ',I2,19X,'Largest ',I2,T106,'*'
     $    /T3,'*',T106,'*')
      END
