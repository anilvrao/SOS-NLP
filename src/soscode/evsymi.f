


      SUBROUTINE EVSYMI(SYMBOL,VALSYM,NSYM,SYMI,INTSYM,FOUND,IER)
C
C ======================================================================
C     EVSYMI===>evsymi   J.T. BETTS
C ======================================================================
C
C
C     PURPOSE:  Find SYMBOL in the list SYMI and load the value in
C               VALSYM string into INTSYM.  WARNING: The length
C               of the strings SYMBOL and SYMI should be the same.
C
C     INPUT:
C
C        SYMBOL    CHARACTER STRING CONTAINING SYMBOL.
C        VALSYM    CHARACTER STRING CONTAINING VALUE.
C        SYMI      LIST OF VALID NAMES.
C
C     OUTPUT:
C
C        INTSYM    ARRAY OF VALUES TO BE UPDATED.
C        FOUND     LOGICAL SPECIFYING WHETHER SYMBOL WAS FOUND IN SYMI.
C        IER       ERROR RETURN FLAG.
C
C     ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(LEN=*)  SYMBOL,VALSYM
      CHARACTER(LEN=6)  SYMI(NSYM)
      CHARACTER(LEN=1)  TYPE
      DIMENSION  INTSYM(NSYM)
      LOGICAL  FOUND
C
C     ------------------------------------------------------------------
C
      IER = 0
      FOUND = .FALSE.
C
C             COMPARE THE INPUT SYMBOL AGAINST THE LIST OF POSSIBLE
C             CANDIDATES.
C
      DO ISYM = 1,NSYM
        IF(SYMBOL.EQ.SYMI(ISYM)) THEN
          FOUND = .TRUE.
C
C             NOW FIND THE NUMERIC VALUE IN VALSYM.
C             AFTER THE CALL THE NUMERIC VALUE IS FOUND
C             IN VALSYM(LCNRL:LCNRR), AND ITS TYPE IS 
C             GIVEN BY THE VARIABLE TYPE.
C
          CALL HHFCNR(VALSYM,' ',LCNRL,LCNRR,TYPE)
C
C             CONVERT THE CHARACTER VARIABLE TO ITS NUMERIC
C             REPRESENTATION AS INDICATED BY THE TYPE FLAG.
C
          IF (TYPE.EQ.'I')
     &      THEN
              CALL HHCNDP(VALSYM(LCNRL:LCNRR),RVAL,IERCNV)
              IF (IERCNV.EQ.0) THEN
                IVAL = INT(RVAL)
                RTST = DBLE(IVAL)
                IF (RTST .NE. RVAL)  IERCNV = -1
              ENDIF
              IF (IERCNV.EQ.0)
     &          THEN
                  INTSYM(ISYM) = IVAL
                ELSE
                  IER = -1
              ENDIF
            ELSE
              IER = -2
          ENDIF
          RETURN
        ENDIF
      ENDDO
C
      RETURN
      END
