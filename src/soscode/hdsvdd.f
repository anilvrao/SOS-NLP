      SUBROUTINE HDSVDD(A,NDIMA,N,M,JOB,NDIMV,WORK,NWORK,NVALER,SIGMA,U,
     +                  V,IER)
C
C*********************************************************************
C
C PURPOSE   HDSVDD COMPUTES THE SINGULAR VALUE DECOMPOSITION (SVD) OF
C           A WHERE A IS AN N BY M REAL GENERAL MATRIX.
C
C METHOD    HDSVDD USES LINPACK SUBROUTINE DSVDC
C
C REMARKS   HDSVDD ALWAYS COMPUTES THE MIN(N,M) SINGULAR VALUES.  THE
C           USER CONTROLS WHETHER THE LEFT AND/OR RIGHT SINGULAR
C           VECTORS ARE COMPUTED.  DEPENDING ON THE END USE OF THE
C           SVD,THE USER CAN SELECT COMPUTING NONE OF THE,
C           THE FIRST MIN(N,M),OR ALL OF THE N LEFT SINGULAR
C           VECTORS.  THE USER CAN ALSO SELECT COMPUTING NONE
C           OF OR ALL OF THE M RIGHT SINGULAR VECTORS.
C           HDSVDS REQUIRES ONLY THE FIRST MIN(N,M) LEFT SINGULAR
C           VECTORS AND ALL OF THE M RIGHT SINGULAR VECTORS TO
C           COMPUTE THE LEAST SQUARES SOLUTION.
C
C
C USAGE     DOUBLE PRECISION A(NDIMA,M),SIGMA(MIN(N,M)),HOLD(NHOLD)
C           DOUBLE PRECISION U(NDIMA,L),V(NDIMV,M)
C           CALL HDSVDD (A,NDIMA,N,M,JOB,NDIMV,NHOLD,HOLD,SIGMA,U,V,IER)
C
C INPUT     A         TWO-DIMENSIONAL ARRAY WHICH STORES THE N BY M
C                     MATRIX.
C           NDIMA     DIMENSIONAL CONSTANT FOR A AND U,NDIMA .GE. N
C
C           N         NUMBER OF ROWS IN A,1 .LE. N
C
C           M         NUMBER OF COLUMNS IN A,1 .LE. M
C
C           JOB       SINGULAR VECTOR COMPUTATION CODE
C                     JOB=00  COMPUTE NO LEFT AND NO RIGHT SINGULAR
C                     VECTORS
C
C                     JOB=01  COMPUTE NO LEFT AND M RIGHT SINGULAR
C                     VECTORS
C
C                     JOB=10  COMPUTE ALL N OF THE LEFT SINGULAR VECTORS
C                     AND NO RIGHT SINGULAR VECTORS
C
C                     JOB=11  COMPUTE ALL N OF THE LEFT SINGULAR VECTORS
C                     AND M RIGHT SINGULAR VECTORS
C
C                     JOB=20  COMPUTE THE FIRST MIN(N,M) LEFT SINGULAR
C                     VECTORS AND NO RIGHT SINGULAR VECTORS
C
C                     JOB=21  COMPUTE THE FIRST MIN(N,M) LEFT SINGULAR
C                     VECTORS AND M RIGHT SINGULAR VECTORS
C
C                     WHEN JOB = 10 OR 11 THE ARRAY U MUST HAVE N
C                     COLUMNS
C                     (I.E.,L = N).  WHEN JOB .GE. 20 THE ARRAY U
C                     MUST HAVE MIN(N,M) COLUMNS (I.E.,L = MIN(N,M)).
C                     WHEN JOB .LE. 01 THE ARRAY U IS UNUSED
C
C           NDIMV     DIMENSIONAL CONSTANT FOR ARRAY V,NDIMV .GE. M
C
C WORKING   WORK      WORK VECTOR OF LENGTH NWORK
C
C           NWORK     THE LENGTH OF THE VECTOR WORK WHICH MUST BE AT
C                     LEAST N + MAX(N,M)
C
C  OUTPUT   A         IF IER = -1,-2,-3,-4,-5 OR -6,THEN
C                     A IS UNCHANGED. OTHERWISE A HAS BEEN DESTROYED.
C
C           SIGMA     IF  IER = 0 THEN SIGMA CONTAINS,IN DESCENDING
C                     ORDER,THE MIN(N,M) SINGULAR VALUES OF A
C
C           U         IF IER = 0 AND JOB = 10 OR 11 THEN THE COLUMNS
C                     OF U CONTAIN THE N LEFT SINGULAR VECTORS.  IF
C                     IER = 0 AND JOB .GE. 20 THEN THE COLUMNS OF
C                     U CONTAIN THE FIRST MIN(N,M) LEFT SINGULAR
C                     VECTORS.  OTHERWISE U IS UNUSED.
C
C           V         IF IER = 0 AND JOB = 01,11,OR 21 THEN THE
C                     COLUMNS OF V CONTAIN THE M RIGHT SINGULAR
C                     VECTORS.  OTHERWISE V IS UNUSED.
C
C           NVALER    IF IER=-9, THEN NVALER GIVES THE NUMBER OF
C                     SINGULAR VALUES THAT HDSVDD FAILED TO
C                     COMPUTE.
C
C           IER       SUCCESS/ERROR CODE. RESULTS HAVE NOT BEEN
C                     COMPUTED FOR IER < 0  POSSIBLE RETURN VALUES
C                     ARE
C
C                     IER =  0,NORMAL RETURN
C                         = -1  N .LE. 0
C                         = -2  M .LE. 0
C                         = -3  NDIMA < N
C                         = -4  NDIMV < M
C                         = -5  JOB IS NOT ONE OF 00,01,10,11,
C                         = -8  FAILURE IN COMPUTING THE SINGULAR VALUE
C                               DECOMPOSITION. THE ENTRIES IN SIGMA
C                               WILL BE SET TO HDMCON(1)
C                         = -9  NOT ALL SIGMAS COMPUTED
C
C     SUBPROGAM HISTORY
C     -----------------
C
C
C*********************************************************************
C
C        LOCAL VARIABLE
C
C     .. SCALAR ARGUMENTS ..
      INTEGER           IER, JOB, M, N, NDIMA, NDIMV, NVALER, NWORK
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION  A(NDIMA,*), SIGMA(*), U(NDIMA,*), V(NDIMV,*),
     +                  WORK(*)
C     ..
C     .. LOCAL SCALARS ..
      CHARACTER(LEN=8)  NAM
      INTEGER           I, MAXNM, NEED
C     ..
C     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION  HDMCON
      EXTERNAL          HDMCON
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL          DSVDC, HHERR
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC         MAX, MIN
C     ..
C     .. DATA STATEMENTS ..
      DATA              NAM/'HDSVDD'/
C     ..
C
C
C------------------- TEST INITIAL PARAMETERS ------------------
C
 
      MAXNM = MAX(N,M)
      NEED = N + MAXNM
      IER = 0
      NVALER = 0
      IF (N.LE.0) IER = -1
      IF (M.LE.0) IER = -2
      IF (NDIMA.LT.N) IER = -3
      IF (NDIMV.LT.M) IER = -4
      IF (JOB.NE.0 .AND. JOB.NE.1 .AND. JOB.NE.10 .AND. JOB.NE.11 .AND.
     +    JOB.NE.20 .AND. JOB.NE.21) IER = -5
      IF (NWORK.LT.NEED) THEN
         IER = -6
         CALL HHERR(2,NAM,IER,NEED)
         SIGMA(1) = HDMCON(1)
 
      ELSEIF (IER.LT.0) THEN
         CALL HHERR(1,NAM,IER,1)
         SIGMA(1) = HDMCON(1)
 
      ELSE
 
         CALL DSVDC(A,NDIMA,N,M,SIGMA,WORK(N+1),U,NDIMA,V,NDIMV,WORK,
     +              JOB,IER)
 
         IF (IER.GT.0)THEN
            IF (IER.EQ.MIN(M,N)) THEN
               IER = -8
 
            ELSE
               NVALER = IER
               IER = -9
               DO I = 1, NVALER
                  SIGMA(I) = HDMCON(1)
               ENDDO
            ENDIF
 
            CALL HHERR(3,NAM,IER,1)
         ENDIF
 
      ENDIF
 
      RETURN
      END
