
      SUBROUTINE MNMSPR(IOPT,M,N,A,IROW,JCOLST,X)
C
C ======================================================================
C     MNMSPR===>mnmspr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  COMPUTE SPARSE MATRIX ROW OR COLUMN NORMS
C
C         INPUT:
C
C            IOPT   INTEGER OPTION CODE
C                   = 1  ROW NORMS (EUCLIDEAN OR 2 NORM)
C                   = 2  COLUMN NORMS (EUCLIDEAN OR 2 NORM)
C                   = -1 ROW NORMS (MAX ABS. VALUE OR INFINITY NORM)
C                   = -2 COLUMN NORMS (MAX ABS. VALUE OR INFINITY NORM)
C            M      NUMBER OF ROWS IN A
C            N      NUMBER OF COLUMNS IN A
C            A      M X N MATRIX STORED AS A VECTOR OF LENGTH NZ 
C            IROW   INTEGER ARRAY OF LENGTH NZ, CONTAINING ROW INDEX
C                   OF NONZERO ELEMENT IN A 
C            JCOLST INTEGER COLUMN START ARRAY OF LENGTH N+1
C                   NOTE NZ = JCOLST(N+1) - 1
C
C         OUTPUT:
C
C            X      VECTOR CONTAINING ROW NORMS (IOPT=1)
C                   VECTOR CONTAINING COLUMN NORMS (IOPT=2)
C
C     *******************************************************
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DIMENSION A(*),X(*),IROW(*),JCOLST(*)
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
C     *******************************************************
C
      IF(IOPT.EQ.1) THEN
C
C         --- ROW NORMS (2-NORM)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
          X(1:M) = ZERO
          NZ = JCOLST(N+1) - 1
          DO K = 1,NZ
            IRK = IROW(K) 
            X(IRK) = X(IRK) + A(K)**2
          ENDDO
C
          DO K = 1,M
            X(K) = SQRT(X(K))
          ENDDO
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,M
            RNRM = ZERO
            DO J=0,N-1
              RNRM = RNRM + A(K+J*MAXROW)**2
            ENDDO
            X(K) = SQRT(RNRM)
          ENDDO
C
        ENDIF
C
      ELSEIF(IOPT.EQ.2) THEN
C
C         --- COLUMN NORMS (2-NORM)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
C
          DO K = 1,N
            CNRM = ZERO
            DO J=JCOLST(K),JCOLST(K+1)-1
              CNRM = CNRM + A(J)**2
            ENDDO
            X(K) = SQRT(CNRM)
          ENDDO
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,N
            J0 = (K-1)*MAXROW
            CNRM = ZERO
            DO J=J0+1,J0+M
              CNRM = CNRM + A(J)**2
            ENDDO
            X(K) = SQRT(CNRM)
          ENDDO
C
        ENDIF
C
      ELSEIF(IOPT.EQ.-1) THEN
C
C         --- ROW NORMS (INFINITY NORM)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
          X(1:M) = ZERO
          NZ = JCOLST(N+1) - 1
          DO K = 1,NZ
            IRK = IROW(K) 
            X(IRK) = MAX(X(IRK),ABS(A(K)))
          ENDDO
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,M
            RNRM = ZERO
            DO J=0,N-1
              RNRM = MAX(RNRM,ABS(A(K+J*MAXROW)))
            ENDDO
            X(K) = RNRM
          ENDDO
C
        ENDIF
C
      ELSEIF(IOPT.EQ.-2) THEN
C
C         --- COLUMN NORMS (INFINITY NORM)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
C
          DO K = 1,N
            CNRM = ZERO
            DO J=JCOLST(K),JCOLST(K+1)-1
              CNRM = MAX(CNRM,ABS(A(J)))
            ENDDO
            X(K) = CNRM
          ENDDO
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,N
            J0 = (K-1)*MAXROW
            CNRM = ZERO
            DO J=J0+1,J0+M
              CNRM = MAX(CNRM,ABS(A(J)))
            ENDDO
            X(K) = CNRM
          ENDDO
C
        ENDIF
C
      ENDIF
C
      RETURN
      END
