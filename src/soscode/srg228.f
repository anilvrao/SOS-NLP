

      SUBROUTINE SRG228(IOPT,M,N,NANZ,AMAT,IROW,JCOL,ITITLE,IOUNIT)
C
C ======================================================================
C     SRG228===>srg228   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C
C        THIS ROUTINE WRITES THE M X N MATRIX AMAT, SPECIFIED BY
C        ROW AND COL, ON OUTPUT UNIT IOUNIT.  ITITLE
C        CONTAINS THE TITLE OF THE PRINT.
C
      DIMENSION  AMAT(*),IROW(*),JCOL(*),JC(5),A(5)
      CHARACTER(LEN=20)  ITITLE
C
      IF(M.LE.0.OR.N.LE.0) RETURN
C
C             WRITE TITLE.
C
      WRITE(IOUNIT,1001)
      WRITE(IOUNIT,1002) ITITLE
C
      IF (IOPT.EQ.1) THEN
C
C             MATRIX EXPRESSED BY IROW,JCOL.
C             FOR EACH ROW, PRINT NONZERO COLUMNS.
C
        DO IR=1,M
C
C             ROW INDEX PRINT.
C
          WRITE(IOUNIT,1003)  IR
C
C             WRITE NONZERO COLUMNS.
C
          NC = 0
          DO K=1,NANZ
            IF (IROW(K).EQ.IR) THEN
              NC = NC + 1
              A(NC) = AMAT(K)
              JC(NC) = JCOL(K)
              IF (NC.EQ.5)  THEN
                CALL RITLYN(IOUNIT,JC,A,5)
                NC = 0
              ENDIF
            ENDIF
          ENDDO
          IF (NC.GT.0) CALL RITLYN(IOUNIT,JC,A,NC)
        ENDDO
C
      ELSEIF (IOPT.EQ.2) THEN
C
C             MATRIX EXPRESSED BY JCOL,IRSTRT(IROW).
C             FOR EACH ROW, PRINT NONZERO COLUMNS.
C
        DO IR=1,M
C
C             ROW INDEX PRINT.
C
          WRITE(IOUNIT,1003)  IR
C
C             WRITE NONZERO COLS.
C
          NC = 0
          DO K=IROW(IR),IROW(IR+1)-1
            NC = NC + 1
            A(NC) = AMAT(K)
            JC(NC) = JCOL(K)
            IF (NC.EQ.5)  THEN
              CALL RITLYN(IOUNIT,JC,A,5)
              NC = 0
            ENDIF
          ENDDO
          IF (NC.GT.0) CALL RITLYN(IOUNIT,JC,A,NC)
        ENDDO
      ELSEIF (IOPT.EQ.3) THEN
C
C             MATRIX EXPRESSED BY IROW,JCSTRT(JCOL).
C             FOR EACH ROW, PRINT NONZERO COLUMNS.
C             ASSUMES ROW INDICES SORTED IN EACH COLUMN.
C
        DO IR=1,M
C
C             ROW INDEX PRINT.
C
          WRITE(IOUNIT,1003)  IR
C
C             WRITE NONZERO COLUMNS.
C
          NC = 0
          DO J=1,N
            DO K=JCOL(J),JCOL(J+1)-1
              IF (IROW(K).EQ.IR)  THEN
                NC = NC + 1
                A(NC) = AMAT(K)
                JC(NC) = J
                IF (NC.EQ.5)  THEN
                  CALL RITLYN(IOUNIT,JC,A,5)
                  NC = 0
                ENDIF
              ELSEIF (IROW(K).GT.IR) THEN
                EXIT
              ENDIF
            ENDDO
          ENDDO
          IF (NC.GT.0) CALL RITLYN(IOUNIT,JC,A,NC)
        ENDDO
      ELSEIF (IOPT.EQ.4) THEN
C
C           DENSE RECTANGULAR
C
        MAXN = M
        CALL RYTMAT(AMAT,MAXN,M,N,IOUNIT)
C
      ELSEIF (IOPT.EQ.5) THEN
C
C           DENSE SYMMETRIC LOWER TRIANGULAR
C
        CALL SYMMAT(AMAT,N,IOUNIT)
C
      ENDIF
C
      WRITE(IOUNIT,1001)
C
      RETURN
 1001 FORMAT(T3,'*',T106,'*'/2X,104('*'))
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T11,A20,T106,'*'/T3,'*',T106,'*')
 1003 FORMAT(T3,'*',T6,'ROW ',I4,T20,40('- '),T106,'*')
      END
