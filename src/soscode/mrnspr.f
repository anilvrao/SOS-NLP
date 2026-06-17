
      SUBROUTINE MRNSPR(IOPT,M,N,A,IROW,JCOLST,X)
C
C
C ======================================================================
C     MRNSPR===>mrnspr   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  COMPUTE SPARSE MATRIX ROW OR COLUMN NORMS
C
C         INPUT:
C
C            IOPT   INTEGER OPTION CODE
C                   = 1  ROW NORMS (SUM OF ABSOLUTE VALUES OF ELEMENTS 
C                        IN EACH ROW)
C                   = 2  COLUMN NORMS (SUM OF ABSOLUTE VALUES OF ELEMENTS 
C                        IN EACH COLUMN)
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
C         --- ROW NORMS (SUM OF ABSOLUTE VALUES)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
          x(1:m) = zero
          NZ = JCOLST(N+1) - 1
          DO K = 1,NZ
            IRK = IROW(K) 
            X(IRK) = X(IRK) + ABS(A(K))
          enddo
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,M
            XK = ZERO
            DO I=0,N-1
              KK = K + I*MAXROW
              XK = XK + ABS(A(KK))
            ENDDO
            X(K) = XK
          enddo
C
        ENDIF
C
      ELSEIF(IOPT.EQ.2) THEN
C
C         --- COLUMN NORMS (SUM OF ABSOLUTE VALUES)
C
        IF(QPOPTN.EQ.'SPARSE') THEN
C
          DO K = 1,N
            XK = ZERO
            DO I=JCOLST(K),JCOLST(K+1)-1
              XK = XK + ABS(A(I))
            ENDDO
            X(K) = XK
          enddo
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,N
            XK = ZERO
            I0 = 1 + (K-1)*MAXROW
            DO I=0,M-1
              XK = XK + ABS(A(I0+I))
            ENDDO
            X(K) = XK
          enddo
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
          enddo
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,M
            X(K) = DAMAX(N,A(K),MAXROW)
          enddo
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
            XK = ZERO
            DO I=JCOLST(K),JCOLST(K+1)-1
              XK = MAX(XK,ABS(A(I)))
            ENDDO
            X(K) = XK
          enddo
C
        ELSE
C
          MAXROW = JCOLST(1)
          DO K = 1,N
            XK = ZERO
            I0 = 1 + (K-1)*MAXROW
            DO I=0,M-1
              XK = XK + ABS(A(I0+I))
            ENDDO
            X(K) = XK
          enddo
C
        ENDIF
C
      ENDIF
C
      RETURN
      END
