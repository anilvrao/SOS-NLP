
      SUBROUTINE MAGRYT(MAXNUM,N,AVEC,WORK,IWRK,TITLE,IPU)
C
C ======================================================================
C     MAGRYT===>magryt   J.T. BETTS
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
C           AVEC   OUTPUT VECTOR (N)
C           WORK   REAL WORK ARRAY (N)
C           IWRK   INTEGER WORK ARRAY (N)
C           IPU    OUTPUT UNIT NO.
C
C        THIS ROUTINE WRITES THE MAGNITUDES OF THE N VECTOR 
C        AVEC IN COMPRESSED FORMAT ON OUTPUT UNIT IPU.  
C
      DIMENSION  AVEC(*),WORK(*),IWRK(*)
      CHARACTER(LEN=*)  TITLE
      PARAMETER (ZERO=0.0D0)
C
      IF(N.EQ.0) RETURN
C
C         LOAD ABSOLUTE VALUE OF AVEC INTO WORK
C
      DO I = 1,N
        WORK(I) = ABS(AVEC(I))
      ENDDO
C
C         SORT WORK ARRAY TO CONSTRUCT ORDER ARRAY 
C
      IOPT = 0
      ISTBLE = 0
      CALL HDSRTN(WORK,N,IOPT,ISTBLE,IWRK,IERSRT)
C
C         LOCATE SMALLEST NONZERO ELEMENT
C
      ISMALL = 0
      DO I = 1,N
        IF(WORK(I).NE.ZERO) THEN
          ISMALL = I
          EXIT
        ENDIF
      ENDDO
C
      IF(IERSRT.NE.0.OR.ISMALL.EQ.0) RETURN
C
      NONZER = MIN(N-ISMALL+1,MAXNUM)
      WRITE(IPU,1002) TITLE,NONZER,NONZER
C
C         WRITE OUT THE NONZER LARGEST AND SMALLEST ELEMENTS
C
      DO I = 1,NONZER
        WRITE(IPU,1001) IWRK(ISMALL+I-1),ABS(AVEC(IWRK(ISMALL+I-1)))
     $               ,IWRK(N-I+1),ABS(AVEC(IWRK(N-I+1)))
      ENDDO
C
      RETURN
 1001 FORMAT(T3,'*',T6,2(5X,'(',I6,':',1PG16.8,')'),T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,A,T106,'*',/T3,'*',T106,'*'/
     $    T3,'*',T11,'Smallest ',I2,19X,'Largest ',I2,T106,'*'
     $    /T3,'*',T106,'*')
      END
