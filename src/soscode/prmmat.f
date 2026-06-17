

      SUBROUTINE PRMMAT(IOPT,NZ,A,IROW,JCOL,IPRM,KPERM,
     $    NPERM,IWORK,IER)
C
C ======================================================================
C     PRMMAT===>prmmat   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  APPLY ROW OR COLUMN PERMUTATION TO A SPARSE 
C                   MATRIX
C
C         INPUT:
C
C            IOPT   INTEGER OPTION CODE
C                   = 1  APPLY ROW PERMUTATION SPECIFIED BY 
C                        KPERM AND SORT INTO ROW ORDER. IF
C                   KPERM = (3,1,2),  THE RESULT WILL BE
C                   `ROW AFTER'   `ROW BEFORE'
C                   3              1
C                   1              2
C                   2              3
C                   = 2  APPLY COLUMN PERMUTATION SPECIFIED BY 
C                        KPERM AND SORT INTO COLUMN ORDER
C                   = -1 ROW PERMUTATION ON A VECTOR 
C                   = -2 COLUMN PERMUTATION ON A VECTOR
C            NZ     NUMBER OF NONZERO ELEMENTS IN A
C            A      M X N MATRIX STORED AS A VECTOR (NZ)
C            IROW   INTEGER ARRAY, CONTAINING ROW INDEX
C                   OF NONZERO ELEMENT IN A (NZ)
C            JCOL   INTEGER ARRAY, CONTAINING COLUMN
C                   INDEX OF NONZERO ELEMENT IN A (NZ)
C            IPRM   INTEGER PERMUTATION RECORD (NZ)
C            KPERM  INTEGER PERMUTATION ARRAY DEFINING ORDER
C                   OF RESULT.  (NPERM)
C            NPERM  NUMBER OF PERMUTATIONS WHEN IOPT = 1 KPERM SPECIFIES
C                   THE DESIRED ROW PERMUTATION AND MUST BE
C                   OF LENGTH (M).  WHEN IOPT = 2 KPERM SPECIFIES
C                   THE DESIRED COLUMN PERMUTATION AND MUST BE
C                   OF LENGTH (N)
C            IWORK  INTEGER WORK ARRAY OF LENGTH NZ
C
C         OUTPUT:
C
C            A      VECTOR SORTED TO CORRESPOND TO ROW 
C                   PERMUTATION WHEN IOPT = 1 OR COLUMN
C                   PERMUTATION WHEN IOPT = 2.
C
C            IROW   ROW INDICES REORDERED BY COLUMN PERMUTATION
C                   WHEN IOPT = 2, OR IN ASCENDING ORDER WHEN 
C                   A ROW PERMUTATION IS PERFORMED (IOPT = 1) 
C
C            JCOL   COLUMN INDICES IN ASCENDING ORDER WHEN A
C                   A COLUMN PERMUTATION IS PERFORMED (IOPT=2)
C                   REORDERED WHEN A ROW PERMUTATION IS MADE.
C
C            IPRM   MASTER PERMUTATION RECORD, I.E. 
C                   `AFTER'      `BEFORE'
C                   A(K)     =   A(IPRM(K))
C
C            IER    ERROR RETURN FLAG
C                   0         NORMAL RETURN
C                   .LT. 0    PROCESSING ERROR
C
C     *******************************************************
C
      DIMENSION A(NZ),IROW(NZ),JCOL(NZ),IPRM(NZ),KPERM(NPERM)
     $    ,IWORK(NZ)
C
C     *******************************************************
C
      IER = 0
C
C         CONSTRUCT THE INDEX ARRAY WHICH WILL PUT THE PERMUTATION 
C         ARRAY KPERM INTO ASCENDING ORDER -- RESTORE KPERM TO 
C         ORIGINAL ORDER 
C
      CALL HJSRTN(KPERM,NPERM,1,0,IWORK,IER)
      IF(IER.LT.0) RETURN
C
      IF(ABS(IOPT).EQ.1) THEN
C
C         ***ROW PERMUTATION***
C
C         APPLY ROW PERMUTATION TO ROW INDICES
C
        DO I = 1,NZ
          IRWI = IROW(I)
          IF(IRWI.LE.0.OR.IRWI.GT.NZ) THEN
            IER = -1
            RETURN
          ELSE
            IROW(I) = IWORK(IRWI)
          ENDIF
        ENDDO
C
C         SORT THE ROW INDEX ARRAY INTO ASCENDING ORDER AND SAVE
C         THE ORDER INDEX IN IWORK
C
        CALL HJSRTN(IROW,NZ,0,0,IWORK,IER)
        IF(IER.LT.0) RETURN
C
C         REORDER THE OTHER THREE ARRAYS INTO THE SAME ORDER 
C
        IF(IOPT.GT.0) CALL HJPRMX(JCOL,NZ,IWORK,IER)
        IF(IER.LT.0) RETURN
        CALL HJPRMX(IPRM,NZ,IWORK,IER)
        IF(IER.LT.0) RETURN
        CALL HDPRMX(A,NZ,IWORK,IER)
        IF(IER.LT.0) RETURN
C
      ELSE
C
C         ***COLUMN PERMUTATION***
C
C         APPLY COLUMN PERMUTATION TO COLUMN INDICES
C
        DO I = 1,NZ
          JCLI = JCOL(I)
          IF(JCLI.LE.0.OR.JCLI.GT.NZ) THEN
            IER = -1
            RETURN
          ELSE
            JCOL(I) = IWORK(JCLI)
          ENDIF
        ENDDO
C
C         SORT THE COLUMN INDEX ARRAY INTO ASCENDING ORDER AND SAVE
C         THE ORDER INDEX IN IWORK
C
        CALL HJSRTN(JCOL,NZ,0,0,IWORK,IER)
        IF(IER.LT.0) RETURN
C
C         REORDER THE OTHER THREE ARRAYS INTO THE SAME ORDER 
C
        IF(IOPT.GT.0) CALL HJPRMX(IROW,NZ,IWORK,IER)
        IF(IER.LT.0) RETURN
        CALL HJPRMX(IPRM,NZ,IWORK,IER)
        IF(IER.LT.0) RETURN
        CALL HDPRMX(A,NZ,IWORK,IER)
        IF(IER.LT.0) RETURN
      ENDIF
C
      RETURN
      END 
